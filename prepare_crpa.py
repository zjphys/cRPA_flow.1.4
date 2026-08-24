#!/usr/bin/env python3
"""Prepare a cRPA calculation from a completed workflow Wannier stage."""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


MARKER = ".generated-by-poscar-workflow"
COMPLETION_RE = re.compile(
    r"General timing and accounting informations? for this job",
    re.IGNORECASE,
)
NBANDS_RE = re.compile(r"\bNBANDS\s*=\s*(\d+)", re.IGNORECASE)
INCAR_ASSIGNMENT_RE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$"
)
TEMPLATE_TOKEN_RE = re.compile(r"\{\{([A-Z0-9_]+)\}\}")
TARGET_STATE_TOKEN_RE = re.compile(r"^(\d+)(?:-(\d+))?$")

REQUIRED_WANNIER_FILES = (
    "POSCAR",
    "POTCAR",
    "KPOINTS",
    "CHGCAR",
    "WAVECAR",
    "WANPROJ",
)

DEFAULT_INCAR_CRPA_TEMPLATE = """\
SYSTEM = {{SYSTEM}} cRPA
ENCUT = {{ENCUT}}
ISMEAR = 0; SIGMA = 0.1
EDIFF = {{EDIFF}}
NBANDS = {{NBANDS}}
NBANDSGW = {{NBANDSGW}}

# SOC can be enabled in an overridden INCAR_CRPA_TEMPLATE when required.
# LSORBIT = .TRUE.
ISYM = -1

LMAXMIX = 6

ALGO = CRPA
NTARGET_STATES = {{TARGET_STATES}}
ENCUTGW = {{ENCUTGW}}

KPAR = {{CRPA_KPAR}}
"""


class CrpaPreparationError(RuntimeError):
    """Raised when the cRPA stage cannot be prepared safely."""


@dataclass(frozen=True)
class CrpaSettings:
    """Resolved scientific settings written to the cRPA INCAR."""

    nbands: int
    nbandsgw: int
    encut: float
    encutgw: float
    num_wann: int
    ispin: int
    target_states: tuple[int, ...]
    kpar: int
    copied_waveder: bool


def require_file(path: Path, description: str) -> None:
    """Require a non-empty regular file."""

    if not path.is_file() or path.stat().st_size == 0:
        raise CrpaPreparationError(f"{description} is missing or empty: {path}")


def _split_incar_comment(line: str) -> str:
    positions = [
        position
        for marker in ("#", "!")
        if (position := line.find(marker)) >= 0
    ]
    if positions:
        return line[: min(positions)]
    return line


def read_optional_incar_assignment(text: str, tag: str) -> str | None:
    """Return the last effective scalar assignment, or None when absent."""

    wanted = tag.upper()
    values: list[str] = []
    for line in text.splitlines():
        code = _split_incar_comment(line)
        for field in code.split(";"):
            match = INCAR_ASSIGNMENT_RE.match(field)
            if match and match.group(1).upper() == wanted:
                value = match.group(2).strip()
                if value:
                    values.append(value)
    return values[-1] if values else None


def read_incar_assignment(text: str, tag: str) -> str:
    """Return the last effective scalar assignment for an INCAR tag."""

    value = read_optional_incar_assignment(text, tag)
    if value is None:
        raise CrpaPreparationError(f"04_wann/INCAR does not define {tag}")
    return value


def read_effective_ispin(text: str, source: str) -> int:
    """Read ISPIN, applying VASP's default of 1 when it is omitted."""

    value = read_optional_incar_assignment(text, "ISPIN")
    if value is None:
        return 1
    try:
        ispin = int(value)
    except ValueError as exc:
        raise CrpaPreparationError(f"{source} has invalid ISPIN value: {value!r}") from exc
    if ispin not in (1, 2):
        raise CrpaPreparationError(
            f"{source} has unsupported ISPIN value {ispin}; expected 1 or 2"
        )
    return ispin


def read_wanproj_header(path: Path) -> tuple[int, int, int, int]:
    """Read and validate WANPROJ's ISPIN, NKPTS, NB_TOT, and NW fields."""

    require_file(path, "Wannier restart/input WANPROJ")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        raise CrpaPreparationError(f"WANPROJ header is incomplete: {path}")
    fields = lines[1].split()
    if len(fields) != 4:
        raise CrpaPreparationError(
            f"WANPROJ header in {path} must contain ISPIN NKPTS NB_TOT NW"
        )
    try:
        ispin, nkpoints, nbands, num_wann = map(int, fields)
    except ValueError as exc:
        raise CrpaPreparationError(
            f"WANPROJ header in {path} contains non-integer dimensions"
        ) from exc
    if ispin not in (1, 2):
        raise CrpaPreparationError(
            f"WANPROJ header in {path} has unsupported ISPIN {ispin}"
        )
    if nkpoints <= 0 or nbands <= 0 or num_wann <= 0:
        raise CrpaPreparationError(
            f"WANPROJ header in {path} contains non-positive dimensions"
        )
    return ispin, nkpoints, nbands, num_wann


def read_positive_int_assignment(text: str, tag: str) -> int:
    """Read a positive integer INCAR assignment."""

    value = read_incar_assignment(text, tag)
    try:
        result = int(value)
    except ValueError as exc:
        raise CrpaPreparationError(
            f"04_wann/INCAR has invalid {tag} value: {value!r}"
        ) from exc
    if result <= 0:
        raise CrpaPreparationError(
            f"04_wann/INCAR has non-positive {tag} value: {result}"
        )
    return result


def read_positive_float_assignment(text: str, tag: str) -> float:
    """Read a positive finite floating-point INCAR assignment."""

    value = read_incar_assignment(text, tag)
    try:
        result = float(value.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise CrpaPreparationError(
            f"04_wann/INCAR has invalid {tag} value: {value!r}"
        ) from exc
    if not math.isfinite(result) or result <= 0:
        raise CrpaPreparationError(
            f"04_wann/INCAR has non-positive or non-finite {tag} value: {value!r}"
        )
    return result


def read_effective_nbands(path: Path) -> int:
    """Read the final effective NBANDS value recorded by VASP."""

    require_file(path, "completed Wannier OUTCAR")
    text = path.read_text(encoding="utf-8", errors="replace")
    if not COMPLETION_RE.search(text):
        raise CrpaPreparationError(
            f"04_wann is not complete; VASP completion marker is absent from {path}"
        )
    values = [int(match.group(1)) for match in NBANDS_RE.finditer(text)]
    if not values:
        raise CrpaPreparationError(
            f"could not find an effective 'NBANDS =' value in {path}"
        )
    if values[-1] <= 0:
        raise CrpaPreparationError(
            f"04_wann/OUTCAR contains invalid NBANDS: {values[-1]}"
        )
    return values[-1]


def validate_target_states(
    target_states: Sequence[int], num_wann: int
) -> tuple[int, ...]:
    """Validate target Wannier-state indices without changing their order."""

    if not target_states:
        raise CrpaPreparationError("at least one --target-states value is required")
    result = tuple(target_states)
    if any(value <= 0 for value in result):
        raise CrpaPreparationError("--target-states values must be positive integers")
    if len(set(result)) != len(result):
        raise CrpaPreparationError("--target-states values must be unique")
    invalid = [value for value in result if value > num_wann]
    if invalid:
        raise CrpaPreparationError(
            "target state indices exceed NUM_WANN "
            f"({num_wann}): {' '.join(str(value) for value in invalid)}"
        )
    return result


def expand_target_state_tokens(tokens: Sequence[str]) -> tuple[int, ...]:
    """Expand inclusive target-state ranges while preserving input order."""

    expanded: list[int] = []
    for token in tokens:
        match = TARGET_STATE_TOKEN_RE.fullmatch(token)
        if not match:
            raise CrpaPreparationError(
                f"invalid target-state value {token!r}; "
                "use a positive integer or inclusive range such as 1-10"
            )
        start = int(match.group(1))
        end_text = match.group(2)
        if end_text is None:
            expanded.append(start)
            continue
        end = int(end_text)
        if end < start:
            raise CrpaPreparationError(
                f"target-state range {token!r} is descending; "
                "the range start must not exceed its end"
            )
        expanded.extend(range(start, end + 1))
    return tuple(expanded)


def _system_name(poscar: Path) -> str:
    require_file(poscar, "Wannier POSCAR")
    lines = poscar.read_text(encoding="utf-8", errors="replace").splitlines()
    title = lines[0].strip() if lines else ""
    return title or "VASP calculation"


def _format_number(value: float) -> str:
    if value == 0:
        return "0"
    return f"{value:.12g}"


def render_template(template: str, values: Mapping[str, str]) -> str:
    """Render the supported cRPA template placeholders."""

    unknown = sorted(set(TEMPLATE_TOKEN_RE.findall(template)) - set(values))
    if unknown:
        raise CrpaPreparationError(
            "unrecognized INCAR_CRPA_TEMPLATE token(s): "
            + ", ".join(f"{{{{{name}}}}}" for name in unknown)
        )

    rendered = TEMPLATE_TOKEN_RE.sub(lambda match: values[match.group(1)], template)
    if not rendered.endswith("\n"):
        rendered += "\n"
    return rendered


def _install_stage(staged: Path, stage: Path, temporary_root: Path) -> None:
    """Install a fully prepared stage and restore the previous one on failure."""

    if not stage.exists():
        staged.replace(stage)
        return

    backup = temporary_root / "previous-05_crpa"
    stage.replace(backup)
    try:
        staged.replace(stage)
    except Exception:
        backup.replace(stage)
        raise
    shutil.rmtree(backup)


def prepare_crpa(
    root: Path,
    target_states: Sequence[int] | None = None,
    nbandsgw: int | None = None,
    encutgw: float | None = None,
    *,
    kpar: int = 1,
    template: str = DEFAULT_INCAR_CRPA_TEMPLATE,
    force: bool = False,
) -> tuple[Path, CrpaSettings]:
    """Create a complete workflow-owned 05_crpa directory."""

    root = root.expanduser().resolve()
    wannier = root / "04_wann"
    if not (wannier / MARKER).is_file():
        raise CrpaPreparationError(
            f"{wannier} is not a workflow-owned Wannier stage"
        )

    for name in REQUIRED_WANNIER_FILES:
        require_file(wannier / name, f"Wannier restart/input {name}")

    outcar = wannier / "OUTCAR"
    nbands = read_effective_nbands(outcar)
    incar_path = wannier / "INCAR"
    require_file(incar_path, "Wannier INCAR")
    incar_text = incar_path.read_text(encoding="utf-8", errors="replace")
    num_wann = read_positive_int_assignment(incar_text, "NUM_WANN")
    ispin = read_effective_ispin(incar_text, "04_wann/INCAR")
    wanproj_ispin, _, _, wanproj_num_wann = read_wanproj_header(
        wannier / "WANPROJ"
    )
    if wanproj_ispin != ispin:
        raise CrpaPreparationError(
            f"04_wann/INCAR ISPIN ({ispin}) does not match WANPROJ ISPIN "
            f"({wanproj_ispin})"
        )
    if wanproj_num_wann != num_wann:
        raise CrpaPreparationError(
            f"04_wann/INCAR NUM_WANN ({num_wann}) does not match WANPROJ NW "
            f"({wanproj_num_wann})"
        )
    encut = read_positive_float_assignment(incar_text, "ENCUT")
    ediff = read_positive_float_assignment(incar_text, "EDIFF")
    targets = (
        tuple(range(1, num_wann + 1))
        if target_states is None
        else validate_target_states(target_states, num_wann)
    )

    resolved_nbandsgw = nbands if nbandsgw is None else nbandsgw
    if resolved_nbandsgw <= 0:
        raise CrpaPreparationError("--nbandsgw must be positive")
    if resolved_nbandsgw > nbands:
        raise CrpaPreparationError(
            f"NBANDSGW ({resolved_nbandsgw}) cannot exceed NBANDS ({nbands})"
        )

    resolved_encutgw = (2.0 * encut / 3.0) if encutgw is None else encutgw
    if not math.isfinite(resolved_encutgw) or resolved_encutgw <= 0:
        raise CrpaPreparationError("--encutgw must be a positive finite number")
    if kpar <= 0:
        raise CrpaPreparationError("CRPA_KPAR must be a positive integer")

    stage = root / "05_crpa"
    if stage.exists():
        if not (stage / MARKER).is_file():
            raise CrpaPreparationError(
                f"{stage} exists and was not generated by this workflow; "
                "refusing to overwrite it"
            )
        if not force:
            raise CrpaPreparationError(
                f"{stage} already exists; use --force to replace it"
            )

    copied_waveder = (wannier / "WAVEDER").is_file() and (
        wannier / "WAVEDER"
    ).stat().st_size > 0
    values = {
        "SYSTEM": _system_name(wannier / "POSCAR"),
        "ENCUT": _format_number(encut),
        "EDIFF": _format_number(ediff),
        "NBANDS": str(nbands),
        "NBANDSGW": str(resolved_nbandsgw),
        "TARGET_STATES": " ".join(str(value) for value in targets),
        "ENCUTGW": _format_number(resolved_encutgw),
        "CRPA_KPAR": str(kpar),
    }
    rendered_incar = render_template(template, values)
    rendered_ispin = read_effective_ispin(rendered_incar, "rendered cRPA INCAR")
    if rendered_ispin != ispin:
        raise CrpaPreparationError(
            f"rendered cRPA INCAR ISPIN ({rendered_ispin}) does not match "
            f"04_wann ISPIN ({ispin}); custom templates must preserve ISPIN"
        )

    temporary_root = Path(tempfile.mkdtemp(prefix=".prepare-crpa.", dir=root))
    staged = temporary_root / "05_crpa"
    staged.mkdir()
    try:
        for name in REQUIRED_WANNIER_FILES:
            shutil.copy2(wannier / name, staged / name)
        if copied_waveder:
            shutil.copy2(wannier / "WAVEDER", staged / "WAVEDER")
        (staged / "INCAR").write_text(rendered_incar, encoding="utf-8")
        (staged / MARKER).touch()
        _install_stage(staged, stage, temporary_root)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    shutil.rmtree(temporary_root, ignore_errors=True)

    settings = CrpaSettings(
        nbands=nbands,
        nbandsgw=resolved_nbandsgw,
        encut=encut,
        encutgw=resolved_encutgw,
        num_wann=num_wann,
        ispin=ispin,
        target_states=targets,
        kpar=kpar,
        copied_waveder=copied_waveder,
    )
    return stage, settings


def _positive_int(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare workflow stage 05_crpa from completed 04_wann outputs."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="workflow root (default: directory containing this script)",
    )
    parser.add_argument(
        "--target-states",
        nargs="+",
        metavar="I_OR_RANGE",
        help=(
            "Wannier-state indices or inclusive ranges to exclude "
            "(default: every state from 1 through NUM_WANN)"
        ),
    )
    parser.add_argument(
        "--nbandsgw",
        type=_positive_int,
        help="QP band count (default: effective Wannier NBANDS)",
    )
    parser.add_argument(
        "--encutgw",
        type=float,
        help="response-function cutoff in eV (default: 2/3 of ENCUT)",
    )
    parser.add_argument(
        "--kpar",
        type=_positive_int,
        default=os.environ.get("CRPA_KPAR", "1"),
        help="cRPA KPAR value (default: CRPA_KPAR or 1)",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    template = os.environ.get(
        "INCAR_CRPA_TEMPLATE", DEFAULT_INCAR_CRPA_TEMPLATE
    )
    try:
        target_states = (
            None
            if args.target_states is None
            else expand_target_state_tokens(args.target_states)
        )
        stage, settings = prepare_crpa(
            args.root,
            target_states,
            args.nbandsgw,
            args.encutgw,
            kpar=args.kpar,
            template=template,
            force=args.force,
        )
    except (CrpaPreparationError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Prepared {stage}")
    print(f"NBANDS = {settings.nbands}")
    print(f"NBANDSGW = {settings.nbandsgw}")
    print(f"ENCUTGW = {_format_number(settings.encutgw)}")
    print(f"ISPIN = {settings.ispin}")
    print(
        "NTARGET_STATES = "
        + " ".join(str(value) for value in settings.target_states)
    )
    if not settings.copied_waveder:
        print(
            "WARNING: 04_wann/WAVEDER is missing; preparing cRPA without "
            "precomputed long-wave derivatives.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

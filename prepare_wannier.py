#!/usr/bin/env python3
"""Prepare a post-SCF VASP/Wannier90 calculation from workflow outputs."""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import rank_wannier_bands as ranker


MARKER = ".generated-by-poscar-workflow"
CONTROLLED_INCAR_TAGS = frozenset(
    (
        "SYSTEM",
        "KSPACING",
        "KGAMMA",
        "NBANDS",
        "ICHARG",
        "ISYM",
        "LWRITE_WANPROJ",
        "LWANNIER90_RUN",
        "NUM_WANN",
        "WANNIER90_WIN",
    )
)
INCAR_ASSIGNMENT_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=")
NBANDS_RE = re.compile(r"\bNBANDS\s*=\s*(\d+)", re.IGNORECASE)
KPATH_BEGIN_RE = re.compile(r"^\s*begin\s+kpoint_path\b", re.IGNORECASE)
KPATH_END_RE = re.compile(r"^\s*end\s+kpoint_path\b", re.IGNORECASE)
SHELL_MULTIPLICITIES = {"s": 1, "p": 3, "d": 5, "f": 7}


class WannierPreparationError(RuntimeError):
    """Raised when the Wannier stage cannot be prepared safely."""


@dataclass(frozen=True)
class WannierWindows:
    dis_froz_min: float
    dis_froz_max: float
    dis_win_min: float
    dis_win_max: float


@dataclass(frozen=True)
class EnergyGuard:
    energy: float
    band_index: int
    spin: str
    kpoint_index: int


def require_file(path: Path, description: str) -> None:
    """Require a non-empty regular file."""

    if not path.is_file() or path.stat().st_size == 0:
        raise WannierPreparationError(f"{description} is missing or empty: {path}")


def read_effective_nbands(path: Path) -> int:
    """Return the last effective NBANDS value recorded in a VASP OUTCAR."""

    require_file(path, "SCF OUTCAR")
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise WannierPreparationError(f"could not read SCF OUTCAR {path}: {exc}") from exc
    values = [int(match.group(1)) for match in NBANDS_RE.finditer(text)]
    if not values:
        raise WannierPreparationError(
            f"could not find an effective 'NBANDS =' value in SCF OUTCAR: {path}"
        )
    if values[-1] <= 0:
        raise WannierPreparationError(f"SCF OUTCAR contains invalid NBANDS: {values[-1]}")
    return values[-1]


def _is_integer_token(value: str) -> bool:
    try:
        return str(int(value)) == value.lstrip("+")
    except ValueError:
        return False


def read_poscar_atom_counts(path: Path) -> dict[str, int]:
    """Return case-insensitive VASP 5-style species counts from a POSCAR."""

    require_file(path, "SCF POSCAR")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise WannierPreparationError(f"could not read SCF POSCAR {path}: {exc}") from exc
    if len(lines) < 7:
        raise WannierPreparationError(
            f"SCF POSCAR must contain element symbols on line 6 and atom counts "
            f"on line 7: {path}"
        )

    elements = lines[5].split()
    count_fields = lines[6].split()
    if not elements or all(_is_integer_token(token) for token in elements):
        raise WannierPreparationError(
            f"VASP 4-style SCF POSCAR detected; add element symbols on line 6: {path}"
        )
    if len(elements) != len(count_fields):
        raise WannierPreparationError(
            f"SCF POSCAR has {len(elements)} element symbols on line 6 but "
            f"{len(count_fields)} atom counts on line 7: {path}"
        )

    counts: dict[str, int] = {}
    for element, field in zip(elements, count_fields):
        key = element.casefold()
        if not element.strip() or key in counts:
            raise WannierPreparationError(
                f"SCF POSCAR contains an empty or duplicate element symbol "
                f"{element!r}: {path}"
            )
        if not _is_integer_token(field) or int(field) <= 0:
            raise WannierPreparationError(
                f"SCF POSCAR atom count for {element!r} must be a positive "
                f"integer, found {field!r}: {path}"
            )
        counts[key] = int(field)
    return counts


def infer_num_wann(
    poscar: Path,
    elements: Sequence[str],
    orbitals: Sequence[str],
) -> int:
    """Infer NUM_WANN from POSCAR counts and paired orbital multiplicities."""

    paired_projections(elements, orbitals)
    atom_counts = read_poscar_atom_counts(poscar)

    total = 0
    for element, orbital in zip(elements, orbitals):
        key = element.casefold()
        if key not in atom_counts:
            raise WannierPreparationError(
                f"requested element {element!r} is absent from SCF POSCAR {poscar}"
            )
        orbital_key = ranker.normalized_name(orbital)
        multiplicity = SHELL_MULTIPLICITIES.get(orbital_key)
        if multiplicity is None:
            raise WannierPreparationError(
                f"unsupported orbital {orbital!r}; SCF LORBIT=10 ranking accepts "
                "only aggregate shells s, p, d, and f"
            )
        total += atom_counts[key] * multiplicity

    if total <= 0:
        raise WannierPreparationError("inferred NUM_WANN must be positive")
    return total


def _spin_energy_guard(
    eigenval: ranker.EigenvalData,
    spin: str,
    band_indices: Sequence[int],
    *,
    highest: bool,
) -> EnergyGuard | None:
    """Return one spin channel's energy extreme for a group of bands."""

    guard: EnergyGuard | None = None
    for band_index in band_indices:
        for kpoint_index, energy in enumerate(
            eigenval.energies[spin][band_index], start=1
        ):
            candidate = EnergyGuard(energy, band_index, spin, kpoint_index)
            if guard is None or (energy > guard.energy if highest else energy < guard.energy):
                guard = candidate
    return guard


def _format_guard(label: str, guard: EnergyGuard | None) -> str:
    if guard is None:
        return f"{label}=none"
    return (
        f"{label}={_format_energy(guard.energy)} eV "
        f"(band {guard.band_index}, spin {guard.spin}, "
        f"k-point {guard.kpoint_index})"
    )


def _largest_contiguous_interval(
    ranked_band_indices: Sequence[int],
) -> tuple[int, ...]:
    """Choose the longest consecutive run, preferring the best-ranked run."""

    rank_positions = {
        band_index: position
        for position, band_index in enumerate(ranked_band_indices)
    }
    ordered = sorted(rank_positions)
    runs: list[list[int]] = []
    for band_index in ordered:
        if not runs or band_index != runs[-1][-1] + 1:
            runs.append([band_index])
        else:
            runs[-1].append(band_index)
    return tuple(
        min(
            runs,
            key=lambda run: (
                -len(run),
                min(rank_positions[band_index] for band_index in run),
                run[0],
            ),
        )
    )


def calculate_windows(
    eigenval: ranker.EigenvalData,
    target_band_indices: Mapping[str, Sequence[int]],
    num_wann: int,
    frozen_margin: float = 0.1,
) -> WannierWindows:
    """Calculate one window from rank-ordered selections in each spin channel."""

    if not math.isfinite(frozen_margin) or frozen_margin < 0:
        raise WannierPreparationError("--frozen-margin must be finite and nonnegative")
    if num_wann <= 0:
        raise WannierPreparationError("NUM_WANN must be positive")

    expected_spins = set(eigenval.spin_channels)
    actual_spins = set(target_band_indices)
    if actual_spins != expected_spins:
        raise WannierPreparationError(
            "SCF PROCAR target spin channels do not match EIGENVAL: "
            f"targets={sorted(actual_spins)}, EIGENVAL={sorted(expected_spins)}"
        )

    ranked_targets_by_spin: dict[str, tuple[int, ...]] = {}
    window_bands_by_spin: dict[str, tuple[int, ...]] = {}
    target_minima: list[float] = []
    target_maxima: list[float] = []
    lower_guards: list[EnergyGuard] = []
    upper_guards: list[EnergyGuard] = []
    for spin in eigenval.spin_channels:
        ranked_target_bands = tuple(target_band_indices[spin])
        target_bands = tuple(sorted(set(ranked_target_bands)))
        if len(target_bands) != num_wann:
            raise WannierPreparationError(
                f"selected spin {spin} target bands {list(target_bands)} contain "
                f"{len(target_bands)} unique indices, but NUM_WANN is {num_wann}"
            )
        if target_bands[0] < 1 or target_bands[-1] > eigenval.nbands:
            raise WannierPreparationError(
                f"selected spin {spin} target bands {list(target_bands)} are "
                f"outside EIGENVAL bands 1..{eigenval.nbands}"
            )
        window_bands = _largest_contiguous_interval(ranked_target_bands)
        ranked_targets_by_spin[spin] = ranked_target_bands
        window_bands_by_spin[spin] = window_bands

        target_energies = [
            energy
            for band_index in window_bands
            for energy in eigenval.energies[spin][band_index]
        ]
        if not target_energies or not all(
            math.isfinite(value) for value in target_energies
        ):
            raise WannierPreparationError(
                f"spin {spin} target energy range is empty or non-finite"
            )
        target_minima.append(min(target_energies))
        target_maxima.append(max(target_energies))

        lower_guard = _spin_energy_guard(
            eigenval, spin, range(1, window_bands[0]), highest=True
        )
        upper_guard = _spin_energy_guard(
            eigenval,
            spin,
            range(window_bands[-1] + 1, eigenval.nbands + 1),
            highest=False,
        )
        if lower_guard is not None:
            lower_guards.append(lower_guard)
        if upper_guard is not None:
            upper_guards.append(upper_guard)

    target_min = min(target_minima)
    target_max = max(target_maxima)
    lower_guard = max(lower_guards, key=lambda item: item.energy, default=None)
    upper_guard = min(upper_guards, key=lambda item: item.energy, default=None)
    candidate_min = max(
        target_min,
        lower_guard.energy if lower_guard is not None else -math.inf,
    )
    candidate_max = min(
        target_max,
        upper_guard.energy if upper_guard is not None else math.inf,
    )
    frozen_min = candidate_min + frozen_margin
    frozen_max = candidate_max - frozen_margin
    context = (
        "ranked target bands="
        + ", ".join(
            f"{spin}:{list(bands)}"
            for spin, bands in ranked_targets_by_spin.items()
        )
        + "; window bands="
        + ", ".join(
            f"{spin}:{list(bands)}" for spin, bands in window_bands_by_spin.items()
        )
        + f"; window target union=[{_format_energy(target_min)}, "
        f"{_format_energy(target_max)}] eV; "
        f"{_format_guard('lower guard', lower_guard)}; "
        f"{_format_guard('upper guard', upper_guard)}; "
        f"candidate=[{_format_energy(candidate_min)}, "
        f"{_format_energy(candidate_max)}] eV; "
        f"margin={_format_energy(frozen_margin)} eV"
    )
    if frozen_min >= frozen_max:
        raise WannierPreparationError(
            "the isolated frozen window is empty or inverted after shrinking; "
            + context
        )

    for spin in eigenval.spin_channels:
        target_set = frozenset(window_bands_by_spin[spin])
        for kpoint_offset in range(eigenval.nkpoints):
            inside = [
                band_index
                for band_index in range(1, eigenval.nbands + 1)
                if frozen_min
                <= eigenval.energies[spin][band_index][kpoint_offset]
                <= frozen_max
            ]
            non_targets = [band for band in inside if band not in target_set]
            if len(inside) > num_wann:
                raise WannierPreparationError(
                    f"non-target bands enter the frozen window ({non_targets}) and cause it to "
                    f"contain {len(inside)} states at spin {spin}, "
                    f"k-point {kpoint_offset + 1}, exceeding NUM_WANN={num_wann}; "
                    f"bands inside={inside}; {context}"
                )
            if non_targets:
                raise WannierPreparationError(
                    "non-target bands enter the inclusive frozen window at "
                    f"spin {spin}, k-point {kpoint_offset + 1}: {non_targets}; "
                    + context
                )

    return WannierWindows(
        dis_froz_min=frozen_min,
        dis_froz_max=frozen_max,
        dis_win_min=target_min - 5.0,
        dis_win_max=target_max + 5.0,
    )


def paired_projections(
    elements: Sequence[str], orbitals: Sequence[str]
) -> tuple[str, ...]:
    """Return positionally paired Wannier90 projection specifications."""

    if not elements or not orbitals:
        raise WannierPreparationError(
            "at least one element and one orbital are required"
        )
    if len(elements) != len(orbitals):
        raise WannierPreparationError(
            "--elements and --orbitals must contain the same number of values "
            "for positional pairing"
        )
    projections: list[str] = []
    for element, orbital in zip(elements, orbitals):
        if (
            not element.strip()
            or any(character.isspace() or character in ':"' for character in element)
        ):
            raise WannierPreparationError(
                f"invalid projection element {element!r}; use a single element token"
            )
        if (
            not orbital.strip()
            or any(character.isspace() or character in ':"' for character in orbital)
        ):
            raise WannierPreparationError(
                f"invalid projection orbital {orbital!r}; use a single orbital token"
            )
        if ranker.normalized_name(orbital) not in SHELL_MULTIPLICITIES:
            raise WannierPreparationError(
                f"unsupported projection orbital {orbital!r}; use an aggregate "
                "LORBIT=10 shell: s, p, d, or f"
            )
        projections.append(f"{element}:{orbital}")
    return tuple(projections)


def extract_kpoint_path(text: str, source: str = "KPATH.wannier90") -> str:
    """Extract exactly one complete Wannier90 kpoint_path block."""

    lines = text.splitlines()
    starts = [index for index, line in enumerate(lines) if KPATH_BEGIN_RE.match(line)]
    ends = [index for index, line in enumerate(lines) if KPATH_END_RE.match(line)]
    if len(starts) != 1 or len(ends) != 1 or starts[0] >= ends[0]:
        raise WannierPreparationError(
            f"{source} must contain exactly one ordered "
            "'begin kpoint_path'/'end kpoint_path' block"
        )
    return "\n".join(lines[starts[0] : ends[0] + 1])


def _split_incar_comment(line: str) -> tuple[str, str]:
    comment_positions = [
        position for marker in ("#", "!") if (position := line.find(marker)) >= 0
    ]
    if not comment_positions:
        return line, ""
    position = min(comment_positions)
    return line[:position], line[position:]


def _without_controlled_assignments(line: str) -> str | None:
    code, comment = _split_incar_comment(line)
    if not code.strip():
        return line
    kept: list[str] = []
    for field in code.split(";"):
        match = INCAR_ASSIGNMENT_RE.match(field)
        if match and match.group(1).upper() in CONTROLLED_INCAR_TAGS:
            continue
        if field.strip():
            kept.append(field.strip())
    if kept:
        rebuilt = "; ".join(kept)
        if comment:
            rebuilt += " " + comment.lstrip()
        return rebuilt
    if comment:
        return comment
    return None


def _strip_controlled_incar_settings(text: str) -> list[str]:
    lines = text.splitlines()
    result: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        code, _ = _split_incar_comment(line)
        match = INCAR_ASSIGNMENT_RE.match(code)
        if match and match.group(1).upper() == "WANNIER90_WIN":
            quote_count = code.count('"')
            index += 1
            while quote_count % 2 == 1 and index < len(lines):
                quote_count += lines[index].count('"')
                index += 1
            continue
        filtered = _without_controlled_assignments(line)
        if filtered is not None:
            result.append(filtered)
        index += 1
    while result and not result[-1].strip():
        result.pop()
    return result


def _format_energy(value: float) -> str:
    if value == 0:
        return "0"
    return f"{value:.10g}"


def rewrite_incar(
    scf_text: str,
    system: str,
    nbands: int,
    num_wann: int,
    projections: Sequence[str],
    windows: WannierWindows,
    kpoint_path: str,
) -> str:
    """Clone SCF settings and append the controlled Wannier settings."""

    if nbands <= 0 or num_wann <= 0:
        raise WannierPreparationError("NBANDS and NUM_WANN must be positive")
    if num_wann > nbands:
        raise WannierPreparationError(
            f"NUM_WANN ({num_wann}) cannot exceed Wannier NBANDS ({nbands})"
        )
    if not system.strip() or "\n" in system or "\r" in system:
        raise WannierPreparationError("the Wannier SYSTEM name must be one line")
    if not projections:
        raise WannierPreparationError("at least one Wannier projection is required")
    kpoint_path = extract_kpoint_path(kpoint_path, "k-point path")

    inherited = _strip_controlled_incar_settings(scf_text)
    controlled = [
        f"SYSTEM = {system.strip()} Wannier",
        "ICHARG = 11",
        "ISYM = -1",
        f"NBANDS = {nbands}",
        "LWRITE_WANPROJ = .TRUE.",
        "LWANNIER90_RUN = .TRUE.",
        f"NUM_WANN = {num_wann}",
        'WANNIER90_WIN = "',
        "write_hr = true",
        "num_iter = 0",
        "dis_num_iter = 1000",
        f"dis_win_min = {_format_energy(windows.dis_win_min)}",
        f"dis_win_max = {_format_energy(windows.dis_win_max)}",
        f"dis_froz_min = {_format_energy(windows.dis_froz_min)}",
        f"dis_froz_max = {_format_energy(windows.dis_froz_max)}",
        "begin projections",
        *projections,
        "end projections",
        "bands_plot = true",
        kpoint_path,
        "bands_num_points = 20",
        '"',
    ]
    if inherited:
        return "\n".join((*inherited, "", *controlled)) + "\n"
    return "\n".join(controlled) + "\n"


def _vaspkit_argv(command: str) -> list[str]:
    argv = shlex.split(command)
    if not argv:
        raise WannierPreparationError("the VASPKIT command is empty")
    return argv


def run_vaspkit(
    command: str,
    arguments: Sequence[str],
    directory: Path,
    description: str,
    input_text: str | None = None,
) -> None:
    """Run one VASPKIT operation and provide a concise failure."""

    try:
        completed = subprocess.run(
            [*_vaspkit_argv(command), *arguments],
            cwd=directory,
            input=input_text,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise WannierPreparationError(
            f"could not run VASPKIT for {description}: {exc}"
        ) from exc
    if completed.returncode != 0:
        raise WannierPreparationError(
            f"VASPKIT failed during {description} with exit status "
            f"{completed.returncode}"
        )


def _system_name(poscar: Path) -> str:
    require_file(poscar, "SCF POSCAR")
    try:
        title = poscar.read_text(encoding="utf-8", errors="replace").splitlines()[0]
    except (OSError, IndexError) as exc:
        raise WannierPreparationError(
            f"could not read the system name from {poscar}: {exc}"
        ) from exc
    return title.strip() or "VASP calculation"


def _prepare_vaspkit_outputs(
    generation_directory: Path,
    scf_poscar: Path,
    vaspkit: str,
    kpr: float,
) -> tuple[str, Path]:
    shutil.copy2(scf_poscar, generation_directory / "POSCAR")
    run_vaspkit(
        vaspkit,
        (),
        generation_directory,
        "3D Wannier90 k-path generation (task 304, option 3)",
        "304\n3\n",
    )
    kpath_file = generation_directory / "KPATH.wannier90"
    require_file(kpath_file, "VASPKIT task 304 output KPATH.wannier90")
    kpoint_path = extract_kpoint_path(
        kpath_file.read_text(encoding="utf-8", errors="replace"),
        str(kpath_file),
    )

    run_vaspkit(
        vaspkit,
        (),
        generation_directory,
        "Gamma-centered KPOINTS generation (task 102)",
        f"102\n2\n{kpr:.10g}\n",
    )
    kpoints = generation_directory / "KPOINTS"
    require_file(kpoints, "VASPKIT task 102 output KPOINTS")
    return kpoint_path, kpoints


def prepare_wannier(
    root: Path,
    elements: Sequence[str],
    orbitals: Sequence[str],
    num_bands: int | None,
    kpr: float,
    vaspkit: str,
    force: bool = False,
    frozen_margin: float = 0.1,
) -> tuple[Path, ranker.EnergyFrontier, int, int]:
    """Create a complete workflow-owned 04_wann directory."""

    root = root.expanduser().resolve()
    if num_bands is not None and num_bands <= 0:
        raise WannierPreparationError("--num-bands must be positive")
    if not math.isfinite(kpr) or kpr <= 0:
        raise WannierPreparationError("--kpr must be a positive finite number")
    if not math.isfinite(frozen_margin) or frozen_margin < 0:
        raise WannierPreparationError("--frozen-margin must be finite and nonnegative")
    projections = paired_projections(elements, orbitals)

    scf_directory = root / "01_scf"
    dos_directory = root / "02_dos"
    required_scf = {
        "INCAR": "SCF INCAR",
        "POSCAR": "SCF POSCAR",
        "POTCAR": "SCF POTCAR",
        "CHGCAR": "SCF CHGCAR",
        "OUTCAR": "SCF OUTCAR",
        "PROCAR": "SCF PROCAR",
    }
    for name, description in required_scf.items():
        require_file(scf_directory / name, description)
    require_file(dos_directory / "EIGENVAL", "DOS EIGENVAL")
    if num_bands is None:
        num_bands = infer_num_wann(scf_directory / "POSCAR", elements, orbitals)
    procar = ranker.parse_procar(scf_directory / "PROCAR")
    ranked_by_spin = ranker.rank_procar(
        procar, scf_directory / "POSCAR", elements, orbitals
    )
    selected_by_spin = {
        spin: ranker.select_top_bands(ranking, num_bands)
        for spin, ranking in ranked_by_spin.items()
    }
    eigenval = ranker.parse_eigenval(dos_directory / "EIGENVAL")
    ranker.validate_procar_eigenval(procar, eigenval)
    results = ranker.add_spin_energy_statistics(selected_by_spin, eigenval)
    frontier = ranker.calculate_total_frontier(results)
    windows = calculate_windows(
        eigenval,
        {
            spin: [item.band_index for item in selected]
            for spin, selected in selected_by_spin.items()
        },
        num_bands,
        frozen_margin,
    )

    scf_nbands = read_effective_nbands(scf_directory / "OUTCAR")
    wannier_nbands = 2 * scf_nbands
    if num_bands > wannier_nbands:
        raise WannierPreparationError(
            f"NUM_WANN ({num_bands}) exceeds doubled SCF NBANDS "
            f"({wannier_nbands})"
        )

    stage = root / "04_wann"
    if stage.exists():
        if not (stage / MARKER).is_file():
            raise WannierPreparationError(
                f"{stage} exists and was not generated by this workflow; "
                "refusing to overwrite it"
            )
        if not force:
            raise WannierPreparationError(
                f"{stage} already exists; use --force to replace it"
            )

    temporary_root = Path(
        tempfile.mkdtemp(prefix=".prepare-wannier.", dir=root)
    )
    staged = temporary_root / "04_wann"
    generation = temporary_root / "vaspkit"
    staged.mkdir()
    generation.mkdir()
    try:
        kpoint_path, generated_kpoints = _prepare_vaspkit_outputs(
            generation,
            scf_directory / "POSCAR",
            vaspkit,
            kpr,
        )
        scf_incar = (scf_directory / "INCAR").read_text(
            encoding="utf-8", errors="replace"
        )
        incar = rewrite_incar(
            scf_incar,
            _system_name(scf_directory / "POSCAR"),
            wannier_nbands,
            num_bands,
            projections,
            windows,
            kpoint_path,
        )
        for name in ("POSCAR", "POTCAR", "CHGCAR"):
            shutil.copy2(scf_directory / name, staged / name)
        shutil.copy2(generated_kpoints, staged / "KPOINTS")
        (staged / "INCAR").write_text(incar, encoding="utf-8")
        ranker.write_csv(staged / "wannier_band_ranking.csv", results)
        (staged / MARKER).touch()

        if stage.exists():
            shutil.rmtree(stage)
        staged.replace(stage)
    except Exception:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    shutil.rmtree(temporary_root, ignore_errors=True)
    print(
        "Frozen window (eV): "
        f"[{_format_energy(windows.dis_froz_min)}, "
        f"{_format_energy(windows.dis_froz_max)}]"
    )
    return stage, frontier, wannier_nbands, num_bands


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare workflow stage 04_wann from completed SCF and DOS outputs."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="workflow root (default: directory containing this script)",
    )
    parser.add_argument("--elements", nargs="+", required=True)
    parser.add_argument(
        "--orbitals",
        nargs="+",
        required=True,
        help="positionally paired aggregate shells; accepted values: s, p, d, f",
    )
    parser.add_argument(
        "--num-bands",
        type=int,
        help=(
            "number of top-ranked bands and NUM_WANN (default: infer from "
            "01_scf/POSCAR atom counts and paired orbital multiplicities)"
        ),
    )
    parser.add_argument(
        "--kpr",
        type=float,
        default=0.04,
        help="VASPKIT reciprocal resolution for the Gamma mesh (default: 0.04)",
    )
    parser.add_argument(
        "--frozen-margin",
        type=float,
        default=0.1,
        help="inward frozen-window margin in eV (default: 0.1)",
    )
    parser.add_argument(
        "--vaspkit",
        default=os.environ.get("VASPKIT_BIN", "vaspkit"),
        help="VASPKIT command (default: VASPKIT_BIN or vaspkit)",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        stage, frontier, wannier_nbands, num_wann = prepare_wannier(
            args.root,
            args.elements,
            args.orbitals,
            args.num_bands,
            args.kpr,
            args.vaspkit,
            args.force,
            args.frozen_margin,
        )
    except (WannierPreparationError, ranker.PBandError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Prepared {stage}")
    print(f"NBANDS = {wannier_nbands}")
    print(f"NUM_WANN = {num_wann}")
    print(
        "Total energy frontier (eV): "
        f"[{frontier.min_value:.10g}, {frontier.max_value:.10g}]"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

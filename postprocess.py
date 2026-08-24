#!/usr/bin/env python3
"""Plot element- or orbital-projected bands and DOS from a completed workflow.

VASPKIT task 213 is run in ``03_band`` and task 113 in ``02_dos``.  Matching
``PBAND_<element>.dat`` and ``PDOS_<element>.dat`` files are combined into one
band/DOS figure.  Spin-polarized ``*_UP.dat`` and ``*_DW.dat`` pairs are also
recognized.  By default colors represent elements; ``--orbital-element`` uses
the named orbital columns for one selected element instead.  ``--wannier-bands``
overlays interpolated bands from a completed ``04_wann`` calculation.
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot element- or orbital-projected bands and DOS after workflow.sh."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="workflow directory (default: directory containing this script)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="projected_band_dos",
        help="output path without an extension (default: ROOT/projected_band_dos)",
    )
    parser.add_argument(
        "--format",
        dest="formats",
        action="append",
        choices=("png", "pdf", "svg"),
        help="output format; repeat for several formats (default: png and pdf)",
    )
    parser.add_argument("--emin", type=float, default=-5.0, help="minimum energy (eV)")
    parser.add_argument("--emax", type=float, default=5.0, help="maximum energy (eV)")
    parser.add_argument(
        "--dos-max",
        type=float,
        default=None,
        help="maximum DOS axis magnitude (default: determine from visible data)",
    )
    parser.add_argument(
        "--marker-scale",
        type=float,
        default=35.0,
        help="projected-band marker area scale (default: 35)",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        help="plot only these element symbols (default: all matching elements)",
    )
    parser.add_argument(
        "--orbital-element",
        metavar="ELEMENT",
        help=(
            "plot orbital projections for one element instead of element "
            "projections; its total DOS is always included"
        ),
    )
    parser.add_argument(
        "--orbitals",
        nargs="+",
        help=(
            "orbital columns to plot with --orbital-element, for example "
            "'s p d f' or 'dxy dz2 dx2-y2' (default: all common orbital columns)"
        ),
    )
    parser.add_argument("--title", default=None, help="optional figure title")
    parser.add_argument("--dpi", type=int, default=600, help="PNG resolution")
    parser.add_argument(
        "--reuse-data",
        action="store_true",
        help="do not run VASPKIT; reuse existing PBAND/PDOS files",
    )
    parser.add_argument(
        "--wannier-bands",
        action="store_true",
        help="overlay Fermi-aligned interpolated bands from ROOT/04_wann",
    )
    parser.add_argument(
        "--vaspkit",
        default=os.environ.get("VASPKIT_BIN", "vaspkit"),
        help="VASPKIT command (default: VASPKIT_BIN or vaspkit)",
    )
    args = parser.parse_args()
    if args.emin >= args.emax:
        parser.error("--emin must be smaller than --emax")
    if args.dos_max is not None and args.dos_max <= 0:
        parser.error("--dos-max must be positive")
    if args.marker_scale <= 0:
        parser.error("--marker-scale must be positive")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.orbitals and not args.orbital_element:
        parser.error("--orbitals requires --orbital-element")
    if args.orbital_element and args.elements:
        parser.error("--orbital-element cannot be combined with --elements")
    return args


def require_file(path: Path, description: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{description} is missing or empty: {path}")


def run_vaspkit(command: str, task: str, directory: Path) -> None:
    for name in ("INCAR", "POSCAR", "DOSCAR"):
        require_file(directory / name, f"VASPKIT input {name}")
    if task == "213":
        for name in ("EIGENVAL", "KPOINTS", "PROCAR"):
            require_file(directory / name, f"VASPKIT projected-band input {name}")

    argv = shlex.split(command)
    if not argv:
        raise RuntimeError("the VASPKIT command is empty")
    if shutil.which(argv[0]) is None and not Path(argv[0]).is_file():
        raise RuntimeError(
            f"VASPKIT executable '{argv[0]}' was not found; load VASPKIT, set "
            "VASPKIT_BIN, or use --reuse-data"
        )

    print(f"Running VASPKIT task {task} in {directory}")
    result = subprocess.run(
        [*argv, "-task", task],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-20:])
        raise RuntimeError(
            f"VASPKIT task {task} failed with exit code {result.returncode}:\n{tail}"
        )


def load_numeric_table(path: Path, minimum_columns: int) -> np.ndarray:
    require_file(path, "data file")
    try:
        data = np.loadtxt(path, comments=("#", "@", "&"), ndmin=2)
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"could not read numeric data from {path}: {exc}") from exc
    if data.shape[1] < minimum_columns:
        raise RuntimeError(
            f"{path} has {data.shape[1]} columns; expected at least {minimum_columns}"
        )
    if not np.all(np.isfinite(data)):
        raise RuntimeError(f"{path} contains non-finite values")
    return data


def read_effective_ispin(path: Path) -> int:
    """Read the last effective ISPIN assignment, using VASP's default of 1."""

    require_file(path, "Wannier INCAR")
    assignments: list[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        code = re.split(r"[#!]", raw_line, maxsplit=1)[0]
        for field in code.split(";"):
            match = re.match(r"^\s*ISPIN\s*=\s*(\S+)\s*$", field, re.IGNORECASE)
            if match:
                assignments.append(match.group(1))
    if not assignments:
        return 1
    try:
        ispin = int(assignments[-1])
    except ValueError as exc:
        raise RuntimeError(
            f"{path} has invalid ISPIN value: {assignments[-1]!r}"
        ) from exc
    if ispin not in (1, 2):
        raise RuntimeError(f"{path} has unsupported ISPIN={ispin}; expected 1 or 2")
    return ispin


def read_last_fermi_energy(path: Path) -> float:
    """Read the final VASP E-fermi value from a completed OUTCAR."""

    require_file(path, "Wannier OUTCAR")
    matches = re.findall(
        r"E-fermi\s*:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
        path.read_text(encoding="utf-8", errors="replace"),
        flags=re.IGNORECASE,
    )
    if not matches:
        raise RuntimeError(f"no E-fermi value could be read from {path}")
    fermi = float(matches[-1])
    if not np.isfinite(fermi):
        raise RuntimeError(f"the final E-fermi value in {path} is not finite")
    return fermi


def load_wannier_band_blocks(path: Path) -> list[np.ndarray]:
    """Parse blank-separated Wannier90 bands as independent x/energy arrays."""

    require_file(path, "Wannier90 band data")
    blocks: list[np.ndarray] = []
    rows: list[tuple[float, float]] = []

    def finish_block() -> None:
        if rows:
            blocks.append(np.asarray(rows, dtype=float))
            rows.clear()

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1
    ):
        stripped = raw_line.strip()
        if not stripped:
            finish_block()
            continue
        if stripped.startswith(("#", "@", "&")):
            continue
        fields = stripped.split()
        if len(fields) < 2:
            raise RuntimeError(
                f"{path}:{line_number} has {len(fields)} column(s); expected at least 2"
            )
        try:
            x_value, energy = float(fields[0]), float(fields[1])
        except ValueError as exc:
            raise RuntimeError(
                f"{path}:{line_number} contains non-numeric k-distance or energy"
            ) from exc
        if not np.isfinite(x_value) or not np.isfinite(energy):
            raise RuntimeError(f"{path}:{line_number} contains non-finite values")
        rows.append((x_value, energy))
    finish_block()

    if not blocks:
        raise RuntimeError(f"{path} contains no readable Wannier90 band blocks")
    reference_x = blocks[0][:, 0]
    if len(reference_x) < 2 or float(np.ptp(reference_x)) <= 0.0:
        raise RuntimeError(f"{path} uses an empty or zero-length Wannier k-path")
    for index, block in enumerate(blocks[1:], start=2):
        if block.shape[0] != reference_x.shape[0] or not np.allclose(
            block[:, 0], reference_x
        ):
            raise RuntimeError(
                f"band block {index} in {path} does not use the same k-grid as band 1"
            )
    return blocks


def load_wannier_bands(wannier_directory: Path) -> tuple[dict[str, list[np.ndarray]], float]:
    """Load the ISPIN-appropriate Wannier90 band files and final Fermi energy."""

    ispin = read_effective_ispin(wannier_directory / "INCAR")
    if ispin == 1:
        paths = {"plain": wannier_directory / "wannier90_band.dat"}
    else:
        paths = {
            "up": wannier_directory / "wannier90.1_band.dat",
            "down": wannier_directory / "wannier90.2_band.dat",
        }
    bands = {spin: load_wannier_band_blocks(path) for spin, path in paths.items()}
    reference_x = next(iter(bands.values()))[0][:, 0]
    for spin, spin_bands in list(bands.items())[1:]:
        x_values = spin_bands[0][:, 0]
        if x_values.shape != reference_x.shape or not np.allclose(x_values, reference_x):
            raise RuntimeError(
                f"Wannier90 {spin} bands do not use the same k-grid as spin up"
            )
    fermi = read_last_fermi_energy(wannier_directory / "OUTCAR")
    return bands, fermi


def align_wannier_bands(
    bands: dict[str, list[np.ndarray]],
    fermi: float,
    target_min: float,
    target_max: float,
) -> dict[str, list[np.ndarray]]:
    """Shift energies by EF and map the Wannier k-distance to a target domain."""

    reference = next(iter(bands.values()))[0][:, 0]
    source_min = float(np.min(reference))
    source_max = float(np.max(reference))
    source_span = source_max - source_min
    if source_span <= 0.0:
        raise RuntimeError("Wannier90 bands use a zero-length k-path")
    scale = (target_max - target_min) / source_span
    return {
        spin: [
            np.column_stack(
                (
                    target_min + (band[:, 0] - source_min) * scale,
                    band[:, 1] - fermi,
                )
            )
            for band in spin_bands
        ]
        for spin, spin_bands in bands.items()
    }


def read_column_names(path: Path, column_count: int) -> list[str]:
    candidates: list[list[str]] = []
    inferred_total_candidates: list[list[str]] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = raw_line.lstrip()
        if not stripped.startswith(("#", "@", "&")):
            continue
        fields = stripped[1:].split()
        if len(fields) == column_count:
            candidates.append(fields)
        elif len(fields) == column_count - 1 and fields:
            inferred_total_candidates.append([*fields, "tot"])
    if candidates:
        return candidates[-1]
    return inferred_total_candidates[-1] if inferred_total_candidates else []


def normalized_name(name: str) -> str:
    return "".join(character for character in name.upper() if character.isalnum())


def is_total_column(name: str) -> bool:
    normalized = normalized_name(name)
    return "TOT" in normalized or "TOTAL" in normalized


def total_band_weight(path: Path, data: np.ndarray) -> np.ndarray:
    names = read_column_names(path, data.shape[1])
    total_columns = [
        index
        for index, name in enumerate(names)
        if index >= 2 and is_total_column(name)
    ]
    if total_columns:
        weight = data[:, total_columns[-1]]
    else:
        weight = np.sum(data[:, 2:], axis=1)
    return np.clip(weight, 0.0, None)


def projected_dos(path: Path, data: np.ndarray) -> tuple[np.ndarray, np.ndarray | None]:
    names = read_column_names(path, data.shape[1])
    total_columns = [
        index
        for index, name in enumerate(names)
        if index >= 1 and is_total_column(name)
    ]
    if not total_columns:
        return np.abs(np.sum(data[:, 1:], axis=1)), None

    up_column = next(
        (
            index
            for index in total_columns
            if "UP" in normalized_name(names[index])
        ),
        None,
    )
    down_column = next(
        (
            index
            for index in total_columns
            if "DOWN" in normalized_name(names[index])
            or "DW" in normalized_name(names[index])
        ),
        None,
    )
    if up_column is not None and down_column is not None:
        return np.abs(data[:, up_column]), np.abs(data[:, down_column])
    return np.abs(data[:, total_columns[-1]]), None


def orbital_column_names(path: Path, data: np.ndarray, first_column: int) -> dict[str, tuple[str, int]]:
    names = read_column_names(path, data.shape[1])
    if not names:
        raise RuntimeError(
            f"{path} has no readable column header; orbital-resolved plotting "
            "requires named columns"
        )
    columns: dict[str, tuple[str, int]] = {}
    for index, name in enumerate(names[first_column:], start=first_column):
        if is_total_column(name):
            continue
        key = normalized_name(name)
        if key:
            columns[key] = (name, index)
    if not columns:
        raise RuntimeError(f"{path} has no named orbital columns")
    return columns


def split_spin_name(name: str) -> tuple[str, str | None]:
    normalized = normalized_name(name)
    for marker, spin in (
        ("SPINDOWN", "down"),
        ("DOWN", "down"),
        ("SPINDW", "down"),
        ("DW", "down"),
        ("SPINUP", "up"),
        ("UP", "up"),
    ):
        if normalized.endswith(marker) and len(normalized) > len(marker):
            return normalized[: -len(marker)], spin
        if normalized.startswith(marker) and len(normalized) > len(marker):
            return normalized[len(marker) :], spin
    return normalized, None


def dos_orbital_columns(
    path: Path, data: np.ndarray
) -> dict[str, tuple[str, int | None, int | None, int | None]]:
    raw_columns = orbital_column_names(path, data, first_column=1)
    grouped: dict[str, dict[str, tuple[str, int]]] = {}
    for display_name, index in raw_columns.values():
        orbital, spin = split_spin_name(display_name)
        grouped.setdefault(orbital, {})[spin or "plain"] = (display_name, index)

    result: dict[str, tuple[str, int | None, int | None, int | None]] = {}
    for orbital, columns in grouped.items():
        plain = columns.get("plain")
        up = columns.get("up")
        down = columns.get("down")
        display_name = (plain or up or down)[0]
        result[orbital] = (
            display_name,
            plain[1] if plain else None,
            up[1] if up else None,
            down[1] if down else None,
        )
    return result


def select_orbitals(
    band_columns: dict[str, tuple[str, int]],
    dos_columns: dict[str, tuple[str, int | None, int | None, int | None]],
    requested: list[str] | None,
) -> list[str]:
    matching = [orbital for orbital in band_columns if orbital in dos_columns]
    if requested:
        requested_keys = [normalized_name(name) for name in requested]
        missing = [
            name for name, key in zip(requested, requested_keys) if key not in matching
        ]
        if missing:
            raise RuntimeError(
                "no matching PBAND/PDOS orbital column(s): " + ", ".join(missing)
            )
        matching = requested_keys
    else:
        shell_order = ("S", "P", "D", "F")
        matching = [
            *(orbital for orbital in shell_order if orbital in matching),
            *(orbital for orbital in matching if orbital not in shell_order),
        ]
    if not matching:
        raise RuntimeError("no matching orbital columns were found in PBAND and PDOS")
    return matching


def orbital_dos(
    data: np.ndarray,
    columns: tuple[str, int | None, int | None, int | None],
) -> tuple[np.ndarray, np.ndarray | None]:
    _, plain_column, up_column, down_column = columns
    if up_column is not None and down_column is not None:
        return np.abs(data[:, up_column]), np.abs(data[:, down_column])
    if plain_column is not None:
        return np.abs(data[:, plain_column]), None
    if up_column is not None:
        return np.abs(data[:, up_column]), None
    if down_column is not None:
        return np.abs(data[:, down_column]), None
    raise RuntimeError("internal error: orbital has no DOS column")


GREEK_KLABELS = {
    "ALPHA": r"\alpha",
    "BETA": r"\beta",
    "GAMMA": r"\Gamma",
    "DELTA": r"\Delta",
    "EPSILON": r"\epsilon",
    "ZETA": r"\zeta",
    "ETA": r"\eta",
    "THETA": r"\Theta",
    "IOTA": r"\iota",
    "KAPPA": r"\kappa",
    "LAMBDA": r"\Lambda",
    "MU": r"\mu",
    "NU": r"\nu",
    "XI": r"\Xi",
    "OMICRON": r"O",
    "PI": r"\Pi",
    "RHO": r"\rho",
    "SIGMA": r"\Sigma",
    "TAU": r"\tau",
    "UPSILON": r"\Upsilon",
    "PHI": r"\Phi",
    "CHI": r"\chi",
    "PSI": r"\Psi",
    "OMEGA": r"\Omega",
}


def klabel_expression(label: str) -> str:
    label = label.strip().strip("$")
    upper = label.upper()
    if upper in {"Γ", r"\GAMMA"}:
        return r"\Gamma"
    if upper in GREEK_KLABELS:
        return GREEK_KLABELS[upper]
    if label.startswith("\\"):
        return label

    indexed = re.fullmatch(r"([A-Za-z]+)_?(\d+)", label)
    if indexed:
        base, index = indexed.groups()
        base_expression = GREEK_KLABELS.get(base.upper(), base)
        return rf"{base_expression}_{{{index}}}"

    primed = re.fullmatch(r"([A-Za-z]+)('+)", label)
    if primed:
        base, primes = primed.groups()
        base_expression = GREEK_KLABELS.get(base.upper(), base)
        return base_expression + "^{" + r"\prime" * len(primes) + "}"
    return label.replace(" ", r"\,")


def latex_klabel(label: str) -> str:
    expressions = [klabel_expression(part) for part in label.split("|")]
    return "$" + r"\mid ".join(expressions) + "$"


def read_klabels(path: Path) -> tuple[list[float], list[str]]:
    require_file(path, "VASPKIT KLABELS file")
    coordinates: list[float] = []
    labels: list[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw_line.split()
        if len(fields) < 2:
            continue
        try:
            coordinate = float(fields[-1])
        except ValueError:
            continue
        label = " ".join(fields[:-1])
        if label.upper().startswith("UNDEFINED"):
            label = ""
        if coordinates and np.isclose(coordinate, coordinates[-1]):
            if label and label != labels[-1]:
                labels[-1] = f"{labels[-1]}|{label}" if labels[-1] else label
            continue
        coordinates.append(coordinate)
        labels.append(label)
    if not coordinates:
        raise RuntimeError(f"no k-point labels could be read from {path}")
    return coordinates, [latex_klabel(label) if label else "" for label in labels]


def discover_files(directory: Path, prefix: str) -> dict[str, dict[str, Path]]:
    files: dict[str, dict[str, Path]] = {}
    for path in sorted(directory.glob(f"{prefix}_*.dat")):
        element = path.stem[len(prefix) + 1 :]
        spin = "plain"
        spin_match = re.search(r"_(UP|DW|DOWN)$", element, flags=re.IGNORECASE)
        if spin_match:
            spin = "up" if spin_match.group(1).upper() == "UP" else "down"
            element = element[: spin_match.start()]
        if element.upper() in {"SUM", "USER"}:
            continue
        key = element.lower()
        if spin in files.setdefault(key, {}):
            raise RuntimeError(
                f"multiple {spin} {prefix} files were found for element {element}"
            )
        files[key][spin] = path
    return files


def element_label(files: dict[str, Path], prefix: str) -> str:
    path = next(iter(files.values()))
    label = path.stem.removeprefix(f"{prefix}_")
    return re.sub(r"_(?:UP|DW|DOWN)$", "", label, flags=re.IGNORECASE)


def validate_spin_files(files: dict[str, Path], description: str) -> None:
    channels = set(files)
    if channels == {"plain"} or channels == {"up", "down"}:
        return
    if "plain" in channels:
        raise RuntimeError(
            f"{description} mixes a non-spin file with spin-resolved files"
        )
    missing = "*_DW.dat" if "up" in channels else "*_UP.dat"
    raise RuntimeError(f"{description} is missing its matching {missing} file")


def select_elements(
    band_files: dict[str, dict[str, Path]],
    dos_files: dict[str, dict[str, Path]],
    requested: list[str] | None,
) -> list[str]:
    matching = [element for element in band_files if element in dos_files]
    if requested:
        lookup = {element.lower(): element for element in requested}
        missing = [symbol for symbol in requested if symbol.lower() not in matching]
        if missing:
            raise RuntimeError(
                "no matching PBAND/PDOS files for element(s): " + ", ".join(missing)
            )
        matching = [element for element in matching if element in lookup]
    if not matching:
        raise RuntimeError(
            "no matching PBAND/PDOS element files were found"
        )
    expected_channels: set[str] | None = None
    for element in matching:
        validate_spin_files(band_files[element], f"PBAND data for {element}")
        validate_spin_files(dos_files[element], f"PDOS data for {element}")
        if set(band_files[element]) != set(dos_files[element]):
            raise RuntimeError(
                f"PBAND and PDOS spin channels do not match for element {element}"
            )
        channels = set(band_files[element])
        if expected_channels is None:
            expected_channels = channels
        elif channels != expected_channels:
            raise RuntimeError("selected elements do not use the same spin channels")
    return matching


def band_slices(path: Path, x: np.ndarray) -> list[slice]:
    text = path.read_text(encoding="utf-8", errors="replace")
    layout = re.search(
        r"NKPTS\s*&\s*NBANDS\s*:\s*(\d+)\s+(\d+)",
        text,
        flags=re.IGNORECASE,
    )
    if layout:
        number_kpoints, number_bands = map(int, layout.groups())
        expected_rows = number_kpoints * number_bands
        if len(x) != expected_rows:
            raise RuntimeError(
                f"{path} declares {number_kpoints} k-points and {number_bands} "
                f"bands ({expected_rows} rows), but contains {len(x)} numeric rows"
            )
        return [
            slice(index * number_kpoints, (index + 1) * number_kpoints)
            for index in range(number_bands)
        ]

    number_bands = len(
        re.findall(r"^\s*#\s*Band-Index\s*:", text, flags=re.IGNORECASE | re.MULTILINE)
    )
    if number_bands:
        number_kpoints, remainder = divmod(len(x), number_bands)
        if remainder:
            raise RuntimeError(
                f"{path} contains {len(x)} rows that cannot be divided into "
                f"{number_bands} Band-Index blocks"
            )
        return [
            slice(index * number_kpoints, (index + 1) * number_kpoints)
            for index in range(number_bands)
        ]

    # Compatibility fallback for headerless tables: a decrease in cumulative
    # k-path distance marks the start of the next band.
    starts = [0, *(np.flatnonzero(np.diff(x) < 0.0) + 1).tolist()]
    stops = [*starts[1:], len(x)]
    return [slice(start, stop) for start, stop in zip(starts, stops)]


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 9,
            "mathtext.fontset": "stix",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.65,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "figure.dpi": 150,
        }
    )


def plot_projected_band_dos(
    band_files: dict[str, dict[str, Path]],
    dos_files: dict[str, dict[str, Path]],
    elements: list[str],
    orbital_element: str | None,
    requested_orbitals: list[str] | None,
    labels_path: Path,
    output_stem: Path,
    formats: list[str],
    emin: float,
    emax: float,
    dos_max: float | None,
    marker_scale: float,
    title: str | None,
    dpi: int,
    wannier_bands: dict[str, list[np.ndarray]] | None = None,
    wannier_fermi: float | None = None,
) -> list[Path]:
    band_data = {
        element: {
            spin: load_numeric_table(path, minimum_columns=3)
            for spin, path in band_files[element].items()
        }
        for element in elements
    }
    dos_data = {
        element: {
            spin: load_numeric_table(path, minimum_columns=2)
            for spin, path in dos_files[element].items()
        }
        for element in elements
    }
    tick_positions, tick_labels = read_klabels(labels_path)

    spin_channels = (
        ("plain",)
        if "plain" in band_data[elements[0]]
        else ("up", "down")
    )
    if wannier_bands is not None:
        if set(wannier_bands) != set(spin_channels):
            raise RuntimeError(
                "Wannier90 and projected-band spin channels do not match: "
                f"Wannier90={sorted(wannier_bands)}, "
                f"PBAND={sorted(spin_channels)}"
            )
        if wannier_fermi is None or not np.isfinite(wannier_fermi):
            raise RuntimeError("a finite Wannier Fermi energy is required")
    for spin in spin_channels:
        reference = band_data[elements[0]][spin]
        for element in elements:
            data = band_data[element][spin]
            if data.shape[0] != reference.shape[0] or not np.allclose(
                data[:, :2], reference[:, :2]
            ):
                raise RuntimeError(
                    f"{band_files[element][spin]} does not use the same {spin} band grid"
                )
    x_reference = band_data[elements[0]][spin_channels[0]][:, 0]
    for spin in spin_channels[1:]:
        x = band_data[elements[0]][spin][:, 0]
        if x.shape != x_reference.shape or not np.allclose(x, x_reference):
            raise RuntimeError("spin-up and spin-down PBAND files use different k grids")

    series: list[str]
    display_labels: dict[str, str]
    band_weights: dict[str, dict[str, np.ndarray]]
    dos_series: list[str]
    dos_labels: dict[str, str]
    total_dos_key: str | None = None
    traces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray | None]] = {}
    any_spin = set(spin_channels) == {"up", "down"}
    visible_max = 0.0
    if orbital_element is not None:
        element = elements[0]
        band_columns = {
            spin: orbital_column_names(
                band_files[element][spin], band_data[element][spin], first_column=2
            )
            for spin in spin_channels
        }
        dos_columns = {
            spin: dos_orbital_columns(dos_files[element][spin], dos_data[element][spin])
            for spin in spin_channels
        }
        first_spin = spin_channels[0]
        series = select_orbitals(
            band_columns[first_spin], dos_columns[first_spin], requested_orbitals
        )
        missing_by_spin = {
            spin: [
                orbital
                for orbital in series
                if orbital not in band_columns[spin] or orbital not in dos_columns[spin]
            ]
            for spin in spin_channels
        }
        if any(missing_by_spin.values()):
            details = "; ".join(
                f"{spin}: {', '.join(missing)}"
                for spin, missing in missing_by_spin.items()
                if missing
            )
            raise RuntimeError(f"orbital columns differ between spin channels ({details})")
        display_element = element_label(band_files[element], "PBAND")
        display_labels = {
            orbital: f"{display_element}-{band_columns[first_spin][orbital][0]}"
            for orbital in series
        }
        dos_labels = dict(display_labels)
        band_weights = {
            orbital: {
                spin: np.clip(
                    band_data[element][spin][:, band_columns[spin][orbital][1]],
                    0.0,
                    None,
                )
                for spin in spin_channels
            }
            for orbital in series
        }

        sorted_dos = {}
        for spin in spin_channels:
            data = dos_data[element][spin]
            sorted_dos[spin] = data[np.argsort(data[:, 0])]
        dos_energy = sorted_dos[first_spin][:, 0]
        for spin in spin_channels[1:]:
            spin_energy = sorted_dos[spin][:, 0]
            if spin_energy.shape != dos_energy.shape or not np.allclose(
                spin_energy, dos_energy
            ):
                raise RuntimeError("spin-up and spin-down PDOS files use different energy grids")
        visible = (dos_energy >= emin) & (dos_energy <= emax)
        if not np.any(visible):
            raise RuntimeError(
                f"PDOS data for {display_element} has no points in {emin:g}..{emax:g} eV"
            )
        total_dos_key = "__TOTAL_DOS__"
        if any_spin:
            total_up = projected_dos(
                dos_files[element]["up"], sorted_dos["up"]
            )[0]
            total_down = projected_dos(
                dos_files[element]["down"], sorted_dos["down"]
            )[0]
        else:
            total_up, total_down = projected_dos(
                dos_files[element][first_spin], sorted_dos[first_spin]
            )
        visible_max = max(visible_max, float(np.max(total_up[visible])))
        if total_down is not None:
            any_spin = True
            visible_max = max(visible_max, float(np.max(total_down[visible])))
        traces[total_dos_key] = (dos_energy, total_up, total_down)
        dos_labels[total_dos_key] = f"{display_element}-total"
        dos_series = [total_dos_key, *series]
        for orbital in series:
            if any_spin:
                dos_up = orbital_dos(
                    sorted_dos["up"], dos_columns["up"][orbital]
                )[0]
                dos_down = orbital_dos(
                    sorted_dos["down"], dos_columns["down"][orbital]
                )[0]
            else:
                dos_up, dos_down = orbital_dos(
                    sorted_dos[first_spin], dos_columns[first_spin][orbital]
                )
            visible_max = max(visible_max, float(np.max(dos_up[visible])))
            if dos_down is not None:
                any_spin = True
                visible_max = max(visible_max, float(np.max(dos_down[visible])))
            traces[orbital] = (dos_energy, dos_up, dos_down)
    else:
        series = elements
        display_labels = {
            element: element_label(band_files[element], "PBAND")
            for element in elements
        }
        dos_labels = {
            element: f"{display_labels[element]}-total" for element in elements
        }
        dos_series = series
        band_weights = {
            element: {
                spin: total_band_weight(
                    band_files[element][spin], band_data[element][spin]
                )
                for spin in spin_channels
            }
            for element in elements
        }
        for element, element_data in dos_data.items():
            sorted_dos = {
                spin: data[np.argsort(data[:, 0])]
                for spin, data in element_data.items()
            }
            first_spin = spin_channels[0]
            dos_energy = sorted_dos[first_spin][:, 0]
            for spin in spin_channels[1:]:
                spin_energy = sorted_dos[spin][:, 0]
                if spin_energy.shape != dos_energy.shape or not np.allclose(
                    spin_energy, dos_energy
                ):
                    raise RuntimeError(
                        f"spin-up and spin-down PDOS files for {display_labels[element]} "
                        "use different energy grids"
                    )
            if any_spin:
                dos_up = projected_dos(
                    dos_files[element]["up"], sorted_dos["up"]
                )[0]
                dos_down = projected_dos(
                    dos_files[element]["down"], sorted_dos["down"]
                )[0]
            else:
                dos_up, dos_down = projected_dos(
                    dos_files[element][first_spin], sorted_dos[first_spin]
                )
            visible = (dos_energy >= emin) & (dos_energy <= emax)
            if not np.any(visible):
                raise RuntimeError(
                    f"PDOS data for {display_labels[element]} has no points in "
                    f"{emin:g}..{emax:g} eV"
                )
            visible_max = max(visible_max, float(np.max(dos_up[visible])))
            if dos_down is not None:
                any_spin = True
                visible_max = max(visible_max, float(np.max(dos_down[visible])))
            traces[element] = (dos_energy, dos_up, dos_down)
    if dos_max is None:
        dos_max = 1.05 * visible_max if visible_max > 0 else 1.0

    configure_style()
    colors = plt.get_cmap("tab10").colors
    fig = plt.figure(figsize=(5.3, 3.8))
    grid = GridSpec(1, 2, width_ratios=(2.25, 1.0), wspace=0.08)
    ax_band = fig.add_subplot(grid[0, 0])
    ax_dos = fig.add_subplot(grid[0, 1], sharey=ax_band)

    for spin in spin_channels:
        x = band_data[elements[0]][spin][:, 0]
        energy = band_data[elements[0]][spin][:, 1]
        linestyle = "--" if spin == "down" else "-"
        for section in band_slices(band_files[elements[0]][spin], x):
            ax_band.plot(
                x[section],
                energy[section],
                color="#303030",
                linewidth=0.45,
                linestyle=linestyle,
                zorder=1,
            )
        for index, key in enumerate(series):
            ax_band.scatter(
                x,
                energy,
                s=band_weights[key][spin] * marker_scale,
                color=colors[index % len(colors)],
                alpha=0.55,
                edgecolors="none",
                rasterized=True,
                label=display_labels[key] if spin == spin_channels[0] else None,
                zorder=2 + index,
            )
    wannier_handles: list[Line2D] = []
    if wannier_bands is not None:
        dft_min = float(np.min(x_reference))
        dft_max = float(np.max(x_reference))
        aligned_bands = align_wannier_bands(
            wannier_bands, wannier_fermi, dft_min, dft_max
        )
        wannier_color = "#b2182b"
        for spin in spin_channels:
            linestyle = "--" if spin == "down" else "-"
            for band in aligned_bands[spin]:
                ax_band.plot(
                    band[:, 0],
                    band[:, 1],
                    color=wannier_color,
                    linewidth=0.9,
                    linestyle=linestyle,
                    zorder=20,
                )
            label = (
                "Wannier"
                if spin == "plain"
                else f"Wannier spin {'up' if spin == 'up' else 'down'}"
            )
            wannier_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=wannier_color,
                    linewidth=0.9,
                    linestyle=linestyle,
                    label=label,
                )
            )
    for position in tick_positions:
        ax_band.axvline(position, color="#c8c8c8", linewidth=0.55, zorder=0)
    ax_band.axhline(0.0, color="#555555", linewidth=0.65, linestyle="--")
    ax_band.set_xlim(float(np.min(x_reference)), float(np.max(x_reference)))
    ax_band.set_ylim(emin, emax)
    ax_band.set_xticks(tick_positions, tick_labels)
    ax_band.set_ylabel(r"$E-E_{\mathrm{F}}$ (eV)")
    band_handles, band_labels = ax_band.get_legend_handles_labels()
    ax_band.legend(
        handles=[*band_handles, *wannier_handles],
        labels=[*band_labels, *(handle.get_label() for handle in wannier_handles)],
        loc="best",
        fontsize=7,
        frameon=False,
        markerscale=0.8,
    )

    for key in dos_series:
        dos_energy, dos_up, dos_down = traces[key]
        if key == total_dos_key:
            color = "#202020"
            linewidth = 1.15
        else:
            color_index = series.index(key)
            color = colors[color_index % len(colors)]
            linewidth = 0.9
        ax_dos.plot(
            dos_up,
            dos_energy,
            color=color,
            linewidth=linewidth,
            label=dos_labels[key],
        )
        if dos_down is not None:
            ax_dos.plot(
                -dos_down,
                dos_energy,
                color=color,
                linewidth=linewidth,
                linestyle="--",
            )
    ax_dos.axhline(0.0, color="#555555", linewidth=0.65, linestyle="--")
    if any_spin:
        ax_dos.axvline(0.0, color="#777777", linewidth=0.5)
        spin_legend = ax_dos.legend(
            handles=[
                Line2D([0], [0], color="#555555", label="Spin up"),
                Line2D([0], [0], color="#555555", linestyle="--", label="Spin down"),
            ],
            loc="lower right",
            fontsize=7,
            frameon=False,
        )
        ax_dos.add_artist(spin_legend)
        ax_dos.set_xlim(-dos_max, dos_max)
    else:
        ax_dos.set_xlim(0.0, dos_max)
    ax_dos.legend(loc="upper right", fontsize=7, frameon=False)
    ax_dos.set_xlabel("DOS (states/eV)")
    ax_dos.tick_params(labelleft=False)

    for axis in (ax_band, ax_dos):
        axis.tick_params(which="major", length=3.0, width=0.6)
        axis.tick_params(which="minor", length=1.6, width=0.45)
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax_dos.xaxis.set_minor_locator(AutoMinorLocator(2))
    if title:
        fig.suptitle(title, y=0.985, fontsize=10)
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.14, top=0.92)
    else:
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.14, top=0.97)

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for extension in formats:
        path = output_stem.with_suffix(f".{extension}")
        save_options = {"bbox_inches": "tight", "pad_inches": 0.03}
        if extension == "png":
            save_options["dpi"] = dpi
        fig.savefig(path, **save_options)
        written.append(path)
    plt.close(fig)
    return written


def main() -> int:
    args = parse_args()
    root = args.root.expanduser().resolve()
    band_dir = root / "03_band"
    dos_dir = root / "02_dos"

    try:
        if not args.reuse_data:
            run_vaspkit(args.vaspkit, "213", band_dir)
            run_vaspkit(args.vaspkit, "113", dos_dir)

        labels_file = band_dir / "KLABELS"
        require_file(labels_file, "VASPKIT k-point labels")
        band_files = discover_files(band_dir, "PBAND")
        dos_files = discover_files(dos_dir, "PDOS")
        requested_elements = (
            [args.orbital_element] if args.orbital_element is not None else args.elements
        )
        elements = select_elements(band_files, dos_files, requested_elements)
        wannier_bands = None
        wannier_fermi = None
        if args.wannier_bands:
            wannier_bands, wannier_fermi = load_wannier_bands(root / "04_wann")

        output = Path(args.output).expanduser()
        if not output.is_absolute():
            output = root / output
        formats = list(dict.fromkeys(args.formats or ["png", "pdf"]))
        written = plot_projected_band_dos(
            band_files,
            dos_files,
            elements,
            args.orbital_element,
            args.orbitals,
            labels_file,
            output,
            formats,
            args.emin,
            args.emax,
            args.dos_max,
            args.marker_scale,
            args.title,
            args.dpi,
            wannier_bands,
            wannier_fermi,
        )
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        "Elements: "
        + ", ".join(
            element_label(band_files[element], "PBAND") for element in elements
        )
    )
    if args.orbital_element:
        print("Projection: orbital-resolved")
    if args.wannier_bands:
        print(f"Wannier bands: E-fermi = {wannier_fermi:g} eV")
    print("Saved " + ", ".join(str(path) for path in written))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

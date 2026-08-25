#!/usr/bin/env python3
"""Rank full-zone SCF PROCAR projections and report DOS-grid band ranges."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, TextIO


PROCAR_LAYOUT_RE = re.compile(
    r"#\s*of\s+k-points\s*:\s*(\d+)\s+"
    r"#\s*of\s+bands\s*:\s*(\d+)\s+"
    r"#\s*of\s+ions\s*:\s*(\d+)",
    re.IGNORECASE,
)
PROCAR_FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
PROCAR_COORDINATE_SEPARATOR = r"(?:\s+|(?=[+-]))"
KPOINT_RE = re.compile(
    r"^\s*k-point\s+(\d+)\s*:\s*"
    rf"({PROCAR_FLOAT_PATTERN}){PROCAR_COORDINATE_SEPARATOR}"
    rf"({PROCAR_FLOAT_PATTERN}){PROCAR_COORDINATE_SEPARATOR}"
    rf"({PROCAR_FLOAT_PATTERN})\s+"
    rf"weight\s*=\s*({PROCAR_FLOAT_PATTERN})\s*$",
    re.IGNORECASE,
)
PROCAR_BAND_RE = re.compile(
    r"^\s*band\s+(\d+)\s+#\s*energy\s+([-+0-9.eEdD]+)\s+"
    r"#\s*occ\.\s*([-+0-9.eEdD]+)\s*$",
    re.IGNORECASE,
)
SPIN_COMPONENT_RE = re.compile(r"^\s*spin\s+component\s+(\d+)\s*$", re.IGNORECASE)
SHELL_NAMES = frozenset(("s", "p", "d", "f"))


class PBandError(RuntimeError):
    """Raised when projection or EIGENVAL input is inconsistent or malformed."""


@dataclass(frozen=True)
class ProcarBand:
    band_index: int
    energy: float
    occupation: float
    ion_projections: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class ProcarKPoint:
    kpoint_index: int
    coordinates: tuple[float, float, float]
    weight: float
    bands: tuple[ProcarBand, ...]


@dataclass(frozen=True)
class ProcarData:
    path: Path
    columns: tuple[str, ...]
    nkpoints: int
    nbands: int
    nions: int
    spin_channels: tuple[str, ...]
    kpoints: dict[str, tuple[ProcarKPoint, ...]]


@dataclass(frozen=True)
class RankedBand:
    rank: int
    band_index: int
    bz_weighted_projection: float
    n_kpoints: int
    kpoint_weight_sum: float


@dataclass(frozen=True)
class EigenvalData:
    path: Path
    nkpoints: int
    nbands: int
    spin_channels: tuple[str, ...]
    energies: dict[str, dict[int, tuple[float, ...]]]


@dataclass(frozen=True)
class BandResult:
    rank: int
    band_index: int
    spin: str
    bz_weighted_projection: float
    projection_n_kpoints: int
    kpoint_weight_sum: float
    dos_n_kpoints: int
    min_value: float
    max_value: float
    bandwidth: float


@dataclass(frozen=True)
class EnergyFrontier:
    min_value: float
    max_value: float
    width: float
    number_of_bands: int


def normalized_name(name: str) -> str:
    """Return a case- and punctuation-insensitive orbital name."""

    return "".join(character for character in name.lower() if character.isalnum())


def _float_field(value: str, path: Path, line_number: int, name: str) -> float:
    try:
        parsed = float(value.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise PBandError(
            f"{path}:{line_number}: invalid {name} value {value!r}"
        ) from exc
    if not math.isfinite(parsed):
        raise PBandError(f"{path}:{line_number}: non-finite {name} value")
    return parsed


def _next_content(lines: Sequence[str], cursor: int) -> int:
    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    return cursor


def parse_procar(path: Path) -> ProcarData:
    """Parse shell-resolved non-spin or collinear-spin VASP PROCAR data."""

    if not path.is_file():
        raise PBandError(f"SCF PROCAR file does not exist: {path}")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PBandError(f"could not read {path}: {exc}") from exc

    layout_index = next(
        (index for index, line in enumerate(lines) if PROCAR_LAYOUT_RE.search(line)),
        None,
    )
    if layout_index is None:
        raise PBandError(f"{path}: no PROCAR k-point/band/ion layout header found")
    layout = PROCAR_LAYOUT_RE.search(lines[layout_index])
    assert layout is not None
    nkpoints, nbands, nions = map(int, layout.groups())
    if min(nkpoints, nbands, nions) <= 0:
        raise PBandError(f"{path}:{layout_index + 1}: layout counts must be positive")

    cursor = layout_index + 1
    all_kpoints: list[ProcarKPoint] = []
    columns: tuple[str, ...] | None = None
    while True:
        cursor = _next_content(lines, cursor)
        if cursor >= len(lines):
            break
        spin_match = SPIN_COMPONENT_RE.match(lines[cursor])
        if spin_match is not None:
            if len(all_kpoints) % nkpoints:
                raise PBandError(
                    f"{path}:{cursor + 1}: spin component header interrupts a k grid"
                )
            expected_component = len(all_kpoints) // nkpoints + 1
            if int(spin_match.group(1)) != expected_component:
                raise PBandError(
                    f"{path}:{cursor + 1}: expected spin component "
                    f"{expected_component}, found {spin_match.group(1)}"
                )
            cursor = _next_content(lines, cursor + 1)
            if cursor >= len(lines):
                raise PBandError(f"{path}: spin component has no k-point blocks")
        kpoint_match = KPOINT_RE.match(lines[cursor])
        if kpoint_match is None:
            raise PBandError(
                f"{path}:{cursor + 1}: expected a k-point block; "
                "noncollinear or lm/phase-resolved PROCAR is unsupported"
            )
        kpoint_index = int(kpoint_match.group(1))
        expected_kpoint = len(all_kpoints) % nkpoints + 1
        if kpoint_index != expected_kpoint:
            raise PBandError(
                f"{path}:{cursor + 1}: expected k-point {expected_kpoint}, "
                f"found {kpoint_index}"
            )
        coordinates = tuple(
            _float_field(value, path, cursor + 1, "k-point coordinate")
            for value in kpoint_match.groups()[1:4]
        )
        weight = _float_field(kpoint_match.group(5), path, cursor + 1, "k-point weight")
        if weight < 0:
            raise PBandError(f"{path}:{cursor + 1}: k-point weight must be nonnegative")
        cursor += 1

        bands: list[ProcarBand] = []
        for expected_band in range(1, nbands + 1):
            cursor = _next_content(lines, cursor)
            if cursor >= len(lines):
                raise PBandError(f"{path}: unexpected end while reading band {expected_band}")
            band_match = PROCAR_BAND_RE.match(lines[cursor])
            if band_match is None or int(band_match.group(1)) != expected_band:
                found = band_match.group(1) if band_match else lines[cursor].strip()
                raise PBandError(
                    f"{path}:{cursor + 1}: expected band {expected_band}, found {found!r}"
                )
            energy = _float_field(band_match.group(2), path, cursor + 1, "band energy")
            occupation = _float_field(
                band_match.group(3), path, cursor + 1, "band occupation"
            )
            cursor = _next_content(lines, cursor + 1)
            if cursor >= len(lines):
                raise PBandError(f"{path}: missing projection header for band {expected_band}")
            header = tuple(lines[cursor].split())
            if len(header) < 3 or normalized_name(header[0]) != "ion":
                raise PBandError(f"{path}:{cursor + 1}: expected 'ion ... tot' header")
            found_columns = tuple(normalized_name(name) for name in header[1:])
            if found_columns[-1] != "tot" or any(not name for name in found_columns):
                raise PBandError(f"{path}:{cursor + 1}: invalid PROCAR projection header")
            shell_columns = tuple(name for name in found_columns[:-1])
            if any(name not in SHELL_NAMES for name in shell_columns):
                raise PBandError(
                    f"{path}:{cursor + 1}: only LORBIT=10 shell columns s/p/d/f "
                    f"are supported; found {', '.join(shell_columns)}"
                )
            if len(set(found_columns)) != len(found_columns):
                raise PBandError(f"{path}:{cursor + 1}: duplicate projection column")
            if columns is None:
                columns = found_columns
            elif columns != found_columns:
                raise PBandError(f"{path}:{cursor + 1}: projection columns changed")
            cursor += 1

            ion_rows: list[tuple[float, ...]] = []
            for expected_ion in range(1, nions + 1):
                cursor = _next_content(lines, cursor)
                if cursor >= len(lines):
                    raise PBandError(f"{path}: missing ion {expected_ion} projection row")
                fields = lines[cursor].split()
                if not fields or fields[0] != str(expected_ion):
                    raise PBandError(
                        f"{path}:{cursor + 1}: expected ion {expected_ion} projection row"
                    )
                if len(fields) != len(found_columns) + 1:
                    raise PBandError(
                        f"{path}:{cursor + 1}: ion row has {len(fields) - 1} values; "
                        f"expected {len(found_columns)}"
                    )
                ion_rows.append(
                    tuple(
                        _float_field(value, path, cursor + 1, "projection")
                        for value in fields[1:]
                    )
                )
                cursor += 1

            cursor = _next_content(lines, cursor)
            total_fields = lines[cursor].split() if cursor < len(lines) else []
            has_total_row = bool(
                total_fields and normalized_name(total_fields[0]) == "tot"
            )
            # VASP omits the redundant sum-over-ions row for a single ion.
            if has_total_row:
                if len(total_fields) != len(found_columns) + 1:
                    raise PBandError(
                        f"{path}:{cursor + 1}: malformed total projection row"
                    )
                for value in total_fields[1:]:
                    _float_field(value, path, cursor + 1, "total projection")
                cursor += 1
            elif nions != 1:
                if cursor >= len(lines):
                    raise PBandError(f"{path}: missing total projection row")
                raise PBandError(f"{path}:{cursor + 1}: expected total projection row")
            bands.append(
                ProcarBand(expected_band, energy, occupation, tuple(ion_rows))
            )

        all_kpoints.append(
            ProcarKPoint(kpoint_index, coordinates, weight, tuple(bands))
        )

    if len(all_kpoints) not in (nkpoints, 2 * nkpoints):
        raise PBandError(
            f"{path}: expected {nkpoints} non-spin or {2 * nkpoints} collinear-spin "
            f"k-point blocks, found {len(all_kpoints)}"
        )
    spin_channels = ("none",) if len(all_kpoints) == nkpoints else ("up", "down")
    kpoints = {
        spin: tuple(all_kpoints[index * nkpoints : (index + 1) * nkpoints])
        for index, spin in enumerate(spin_channels)
    }
    assert columns is not None
    return ProcarData(
        path, columns, nkpoints, nbands, nions, spin_channels, kpoints
    )


def _next_nonempty_line(
    lines: Sequence[str], start: int, path: Path, description: str
) -> tuple[int, str]:
    index = start
    while index < len(lines) and not lines[index].strip():
        index += 1
    if index >= len(lines):
        raise PBandError(f"{path}: unexpected end of file while reading {description}")
    return index, lines[index]


def _integer_field(value: str, path: Path, line_number: int, name: str) -> int:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise PBandError(
            f"{path}:{line_number}: invalid {name} value {value!r}"
        ) from exc
    if not math.isfinite(parsed) or not parsed.is_integer():
        raise PBandError(
            f"{path}:{line_number}: {name} must be an integer, found {value!r}"
        )
    return int(parsed)


def parse_eigenval(path: Path) -> EigenvalData:
    """Parse band energies from a VASP EIGENVAL file."""

    if not path.is_file():
        raise PBandError(f"EIGENVAL file does not exist: {path}")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PBandError(f"could not read {path}: {exc}") from exc
    if len(lines) < 6:
        raise PBandError(f"{path}: incomplete EIGENVAL header")

    layout_fields = lines[5].split()
    if len(layout_fields) < 3:
        raise PBandError(
            f"{path}:6: expected NELECT, NKPTS, and NBANDS in EIGENVAL header"
        )
    nkpoints = _integer_field(layout_fields[1], path, 6, "NKPTS")
    nbands = _integer_field(layout_fields[2], path, 6, "NBANDS")
    if nkpoints <= 0 or nbands <= 0:
        raise PBandError(f"{path}: NKPTS and NBANDS must be positive")

    cursor = 6
    spin_channels: tuple[str, ...] | None = None
    collected: dict[str, dict[int, list[float]]] = {}
    for kpoint_number in range(1, nkpoints + 1):
        cursor, kpoint_line = _next_nonempty_line(
            lines, cursor, path, f"k-point {kpoint_number}"
        )
        try:
            kpoint_values = tuple(float(field) for field in kpoint_line.split())
        except ValueError as exc:
            raise PBandError(
                f"{path}:{cursor + 1}: invalid k-point row"
            ) from exc
        if len(kpoint_values) != 4:
            raise PBandError(
                f"{path}:{cursor + 1}: k-point row has {len(kpoint_values)} "
                "columns; expected 4"
            )
        if not all(math.isfinite(value) for value in kpoint_values):
            raise PBandError(f"{path}:{cursor + 1}: non-finite k-point value")
        cursor += 1

        for expected_band in range(1, nbands + 1):
            cursor, band_line = _next_nonempty_line(
                lines,
                cursor,
                path,
                f"Band-Index {expected_band} at k-point {kpoint_number}",
            )
            fields = band_line.split()
            if len(fields) == 3:
                row_channels = ("none",)
                energy_positions = (1,)
            elif len(fields) == 5:
                row_channels = ("up", "down")
                energy_positions = (1, 2)
            else:
                raise PBandError(
                    f"{path}:{cursor + 1}: band row has {len(fields)} columns; "
                    "expected 3 for non-spin or 5 for collinear spin"
                )
            if spin_channels is None:
                spin_channels = row_channels
                collected = {
                    spin: {band: [] for band in range(1, nbands + 1)}
                    for spin in spin_channels
                }
            elif row_channels != spin_channels:
                raise PBandError(
                    f"{path}:{cursor + 1}: inconsistent spin layout in band rows"
                )

            band_index = _integer_field(
                fields[0], path, cursor + 1, "band index"
            )
            if band_index != expected_band:
                raise PBandError(
                    f"{path}:{cursor + 1}: found band index {band_index}; "
                    f"expected {expected_band} at k-point {kpoint_number}"
                )
            try:
                numeric_fields = tuple(float(field) for field in fields[1:])
            except ValueError as exc:
                raise PBandError(
                    f"{path}:{cursor + 1}: invalid numeric band row"
                ) from exc
            if not all(math.isfinite(value) for value in numeric_fields):
                raise PBandError(f"{path}:{cursor + 1}: non-finite band value")
            for spin, position in zip(spin_channels, energy_positions):
                collected[spin][band_index].append(float(fields[position]))
            cursor += 1

    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    if cursor != len(lines):
        raise PBandError(
            f"{path}:{cursor + 1}: unexpected data after declared EIGENVAL blocks"
        )
    if spin_channels is None:
        raise PBandError(f"{path}: EIGENVAL contains no band rows")

    return EigenvalData(
        path=path,
        nkpoints=nkpoints,
        nbands=nbands,
        spin_channels=spin_channels,
        energies={
            spin: {
                band: tuple(values) for band, values in band_energies.items()
            }
            for spin, band_energies in collected.items()
        },
    )


def read_poscar_ion_map(path: Path) -> tuple[dict[str, tuple[int, ...]], int]:
    """Map VASP 5 POSCAR element labels to zero-based PROCAR ion indices."""

    if not path.is_file():
        raise PBandError(f"SCF POSCAR file does not exist: {path}")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PBandError(f"could not read {path}: {exc}") from exc
    if len(lines) < 7:
        raise PBandError(f"{path}: POSCAR must contain symbols and counts on lines 6-7")
    symbols = lines[5].split()
    count_fields = lines[6].split()
    if not symbols or all(field.lstrip("+").isdigit() for field in symbols):
        raise PBandError(f"{path}: VASP 4 POSCAR is unsupported; add element symbols")
    if len(symbols) != len(count_fields):
        raise PBandError(f"{path}: POSCAR symbol/count lengths differ")
    mapping: dict[str, tuple[int, ...]] = {}
    cursor = 0
    for symbol, field in zip(symbols, count_fields):
        key = symbol.casefold()
        if key in mapping:
            raise PBandError(f"{path}: duplicate element symbol {symbol!r}")
        try:
            count = int(field)
        except ValueError as exc:
            raise PBandError(f"{path}: invalid atom count {field!r}") from exc
        if str(count) != field.lstrip("+") or count <= 0:
            raise PBandError(f"{path}: atom count for {symbol!r} must be positive")
        mapping[key] = tuple(range(cursor, cursor + count))
        cursor += count
    return mapping, cursor


def resolve_shells(data: ProcarData, requested: Sequence[str]) -> tuple[int, ...]:
    """Resolve unique aggregate shell requests to PROCAR numeric columns."""

    if not requested:
        raise PBandError("at least one orbital shell is required")
    available = {name: index for index, name in enumerate(data.columns)}
    selected: list[int] = []
    seen: set[str] = set()
    for raw_name in requested:
        key = normalized_name(raw_name)
        if key not in SHELL_NAMES:
            raise PBandError(
                f"unsupported orbital {raw_name!r}; SCF LORBIT=10 ranking accepts "
                "only aggregate shells s, p, d, and f"
            )
        if key in seen:
            raise PBandError(f"orbital shell {raw_name!r} was requested more than once")
        seen.add(key)
        if key not in available:
            shells = ", ".join(name for name in data.columns if name != "tot")
            raise PBandError(
                f"{data.path}: orbital shell {raw_name!r} is absent; "
                f"available shells: {shells}"
            )
        selected.append(available[key])
    return tuple(selected)


def resolve_elements(
    data: ProcarData, poscar: Path, requested: Sequence[str]
) -> tuple[int, ...]:
    """Resolve unique element requests to zero-based PROCAR ion indices."""

    if not requested:
        raise PBandError("at least one element is required")
    mapping, nions = read_poscar_ion_map(poscar)
    if nions != data.nions:
        raise PBandError(
            f"{data.path}: PROCAR contains {data.nions} ions but {poscar} contains {nions}"
        )
    selected: list[int] = []
    seen: set[str] = set()
    for element in requested:
        key = element.casefold()
        if key in seen:
            raise PBandError(f"element {element!r} was requested more than once")
        seen.add(key)
        if key not in mapping:
            available = ", ".join(mapping)
            raise PBandError(
                f"requested element {element!r} is absent from {poscar}; "
                f"available elements: {available}"
            )
        selected.extend(mapping[key])
    return tuple(selected)


def rank_procar(
    data: ProcarData,
    poscar: Path,
    elements: Sequence[str],
    orbitals: Sequence[str],
) -> dict[str, list[RankedBand]]:
    """Rank every spin channel by normalized Brillouin-zone projection weight."""

    ion_indices = resolve_elements(data, poscar, elements)
    column_indices = resolve_shells(data, orbitals)
    rankings: dict[str, list[RankedBand]] = {}
    for spin in data.spin_channels:
        kpoints = data.kpoints[spin]
        weight_sum = math.fsum(kpoint.weight for kpoint in kpoints)
        if not math.isfinite(weight_sum) or weight_sum <= 0:
            raise PBandError(
                f"{data.path}: spin {spin} k-point weights have nonpositive total"
            )
        scores: dict[int, float] = {}
        for band_index in range(1, data.nbands + 1):
            weighted_terms: list[float] = []
            for kpoint in kpoints:
                band = kpoint.bands[band_index - 1]
                projection = math.fsum(
                    band.ion_projections[ion][column]
                    for ion in ion_indices
                    for column in column_indices
                )
                weighted_terms.append(kpoint.weight * projection)
            scores[band_index] = math.fsum(weighted_terms) / weight_sum
        ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
        rankings[spin] = [
            RankedBand(rank, band_index, score, data.nkpoints, weight_sum)
            for rank, (band_index, score) in enumerate(ordered, start=1)
        ]
    return rankings


def add_spin_energy_statistics(
    selected_by_spin: dict[str, Sequence[RankedBand]], eigenval: EigenvalData
) -> list[BandResult]:
    """Join independently ranked SCF spin channels to DOS-grid energies."""

    if set(selected_by_spin) != set(eigenval.spin_channels):
        raise PBandError(
            "SCF PROCAR spin channels do not match DOS EIGENVAL: "
            f"PROCAR={sorted(selected_by_spin)}, "
            f"EIGENVAL={sorted(eigenval.spin_channels)}"
        )
    results: list[BandResult] = []
    for spin in eigenval.spin_channels:
        for item in selected_by_spin[spin]:
            energies = eigenval.energies[spin].get(item.band_index)
            if energies is None or len(energies) != eigenval.nkpoints:
                raise PBandError(
                    f"{eigenval.path}: ranked band {item.band_index}, spin {spin}, "
                    f"is absent or incomplete"
                )
            minimum = min(energies)
            maximum = max(energies)
            results.append(
                BandResult(
                    item.rank,
                    item.band_index,
                    spin,
                    item.bz_weighted_projection,
                    item.n_kpoints,
                    item.kpoint_weight_sum,
                    eigenval.nkpoints,
                    minimum,
                    maximum,
                    maximum - minimum,
                )
            )
    return results


def validate_procar_eigenval(data: ProcarData, eigenval: EigenvalData) -> None:
    """Require SCF projections and DOS energies to use the same bands/spins."""

    if data.nbands != eigenval.nbands:
        raise PBandError(
            f"SCF PROCAR has {data.nbands} bands but DOS EIGENVAL has "
            f"{eigenval.nbands}"
        )
    if data.spin_channels != eigenval.spin_channels:
        raise PBandError(
            "SCF PROCAR spin channels do not match DOS EIGENVAL: "
            f"PROCAR={list(data.spin_channels)}, "
            f"EIGENVAL={list(eigenval.spin_channels)}"
        )


def select_top_bands(
    ranked: Sequence[RankedBand], number_of_bands: int | None
) -> list[RankedBand]:
    """Return all ranked bands or only the requested highest-ranked subset."""

    if number_of_bands is None:
        return list(ranked)
    if number_of_bands <= 0:
        raise PBandError("number of bands must be positive")
    if number_of_bands > len(ranked):
        raise PBandError(
            f"requested {number_of_bands} bands, but only {len(ranked)} "
            "ranked bands are available"
        )
    return list(ranked[:number_of_bands])


def calculate_total_frontier(results: Sequence[BandResult]) -> EnergyFrontier:
    """Calculate the combined energy range of all selected bands and spins."""

    if not results:
        raise PBandError("cannot calculate an energy frontier without band results")
    min_value = min(item.min_value for item in results)
    max_value = max(item.max_value for item in results)
    return EnergyFrontier(
        min_value=min_value,
        max_value=max_value,
        width=max_value - min_value,
        number_of_bands=len({item.band_index for item in results}),
    )


def write_table(
    results: Sequence[BandResult],
    elements: Sequence[str],
    orbitals: Sequence[str],
    stream: TextIO,
) -> None:
    """Write a human-readable ranking table."""

    stream.write(f"Elements: {', '.join(elements)}\n")
    stream.write(f"Orbitals: {', '.join(orbitals)}\n")
    headers = (
        "Rank",
        "Band",
        "Spin",
        "BZ-weighted projection",
        "SCF NKPTS",
        "K-weight sum",
        "DOS NKPTS",
        "Min value (eV)",
        "Max value (eV)",
        "Bandwidth (eV)",
    )
    rows = [
        (
            str(item.rank),
            str(item.band_index),
            item.spin,
            f"{item.bz_weighted_projection:.10g}",
            str(item.projection_n_kpoints),
            f"{item.kpoint_weight_sum:.10g}",
            str(item.dos_n_kpoints),
            f"{item.min_value:.10g}",
            f"{item.max_value:.10g}",
            f"{item.bandwidth:.10g}",
        )
        for item in results
    ]
    widths = [
        max(len(headers[index]), *(len(row[index]) for row in rows))
        for index in range(len(headers))
    ]
    stream.write(
        "  ".join(header.rjust(width) for header, width in zip(headers, widths))
        + "\n"
    )
    stream.write("  ".join("-" * width for width in widths) + "\n")
    for row in rows:
        stream.write(
            "  ".join(value.rjust(width) for value, width in zip(row, widths))
            + "\n"
        )
    frontier = calculate_total_frontier(results)
    band_label = "band" if frontier.number_of_bands == 1 else "bands"
    stream.write("\n")
    stream.write(
        f"Total energy frontier ({frontier.number_of_bands} selected {band_label}, "
        "all spins):\n"
    )
    stream.write(f"  Min value (eV): {frontier.min_value:.10g}\n")
    stream.write(f"  Max value (eV): {frontier.max_value:.10g}\n")
    stream.write(f"  Frontier width (eV): {frontier.width:.10g}\n")


def write_csv(path: Path, results: Sequence[BandResult]) -> None:
    """Write the ranking in machine-readable form."""

    frontier = calculate_total_frontier(results)
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                (
                    "rank",
                    "band_index",
                    "spin",
                    "bz_weighted_projection",
                    "scf_n_kpoints",
                    "kpoint_weight_sum",
                    "dos_n_kpoints",
                    "min_value",
                    "max_value",
                    "bandwidth",
                    "total_frontier_min_value",
                    "total_frontier_max_value",
                    "total_frontier_width",
                )
            )
            for item in results:
                writer.writerow(
                    (
                        item.rank,
                        item.band_index,
                        item.spin,
                        f"{item.bz_weighted_projection:.17g}",
                        item.projection_n_kpoints,
                        f"{item.kpoint_weight_sum:.17g}",
                        item.dos_n_kpoints,
                        f"{item.min_value:.17g}",
                        f"{item.max_value:.17g}",
                        f"{item.bandwidth:.17g}",
                        f"{frontier.min_value:.17g}",
                        f"{frontier.max_value:.17g}",
                        f"{frontier.width:.17g}",
                    )
                )
    except OSError as exc:
        raise PBandError(f"could not write CSV file {path}: {exc}") from exc


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank bands by k-weighted full-zone SCF PROCAR shell projections "
            "and report DOS-stage EIGENVAL extrema and bandwidths."
        )
    )
    parser.add_argument(
        "--scf-directory",
        type=Path,
        default=Path("01_scf"),
        help="directory containing POSCAR and PROCAR (default: 01_scf)",
    )
    parser.add_argument(
        "--elements",
        nargs="+",
        required=True,
        help="POSCAR elements whose ion projections are combined",
    )
    parser.add_argument(
        "--orbitals",
        nargs="+",
        required=True,
        help="aggregate shell columns to sum; accepted values are s, p, d, and f",
    )
    parser.add_argument(
        "--num-bands",
        type=int,
        default=None,
        help=(
            "analyze only the N bands with the largest weighted projection "
            "(default: all bands)"
        ),
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="optional path for a machine-readable CSV copy",
    )
    parser.add_argument(
        "--dos-directory",
        type=Path,
        default=None,
        help=(
            "directory containing EIGENVAL (default: sibling 02_dos of the "
            "resolved SCF directory)"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        scf_directory = args.scf_directory.resolve()
        dos_directory = (
            args.dos_directory.resolve()
            if args.dos_directory is not None
            else scf_directory.parent / "02_dos"
        )
        procar = parse_procar(scf_directory / "PROCAR")
        ranked_by_spin = rank_procar(
            procar, scf_directory / "POSCAR", args.elements, args.orbitals
        )
        selected_by_spin = {
            spin: select_top_bands(ranked, args.num_bands)
            for spin, ranked in ranked_by_spin.items()
        }
        eigenval = parse_eigenval(dos_directory / "EIGENVAL")
        validate_procar_eigenval(procar, eigenval)
        results = add_spin_energy_statistics(selected_by_spin, eigenval)
        write_table(results, args.elements, args.orbitals, sys.stdout)
        if args.csv is not None:
            write_csv(args.csv, results)
    except PBandError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

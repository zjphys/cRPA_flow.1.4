from __future__ import annotations

import contextlib
import csv
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

VERSION_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VERSION_DIR))
import rank_wannier_bands as ranker

POSCAR = """Test material
1.0
1 0 0
0 1 0
0 0 1
Fe O
1 1
Direct
0 0 0
0.5 0.5 0.5
"""


def procar_text(channels) -> str:
    """Build shell-only PROCAR: channel -> kpoint -> (weight, bands -> ions)."""
    nkpoints = len(channels[0])
    nbands = len(channels[0][0][1])
    nions = len(channels[0][0][1][0])
    lines = [
        "PROCAR lm decomposed",
        f" # of k-points: {nkpoints} # of bands: {nbands} # of ions: {nions}",
    ]
    for channel in channels:
        for kpoint_index, (weight, bands) in enumerate(channel, start=1):
            lines.extend(("", f" k-point {kpoint_index} : 0 0 0 weight = {weight}"))
            for band_index, ions in enumerate(bands, start=1):
                lines.extend((
                    "", f" band {band_index} # energy {band_index - 2}.0 # occ. 1.0",
                    "", " ion s p d tot",
                ))
                totals = [0.0, 0.0, 0.0]
                for ion_index, values in enumerate(ions, start=1):
                    totals = [left + right for left, right in zip(totals, values)]
                    lines.append(
                        f" {ion_index} " + " ".join(str(value) for value in (*values, sum(values)))
                    )
                lines.append(" tot " + " ".join(str(value) for value in (*totals, sum(totals))))
    return "\n".join(lines) + "\n"


def eigenval_text(kpoint_bands) -> str:
    nkpoints = len(kpoint_bands)
    nbands = len(kpoint_bands[0])
    lines = ["header"] * 5 + [f" 2 {nkpoints} {nbands}"]
    for kpoint_index, bands in enumerate(kpoint_bands):
        lines.extend(("", f" {kpoint_index}.0 0 0 {1.0 / nkpoints}"))
        for band_index, energy in enumerate(bands, start=1):
            if isinstance(energy, tuple):
                lines.append(f" {band_index} {energy[0]} {energy[1]} 1 1")
            else:
                lines.append(f" {band_index} {energy} 1")
    return "\n".join(lines) + "\n"


class RankWannierBandsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.scf = self.root / "01_scf"
        self.dos = self.root / "02_dos"
        self.scf.mkdir()
        self.dos.mkdir()
        (self.scf / "POSCAR").write_text(POSCAR, encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_procar(self, channels) -> Path:
        path = self.scf / "PROCAR"
        path.write_text(procar_text(channels), encoding="utf-8")
        return path

    def test_weighted_full_zone_ranking_and_normalization(self) -> None:
        path = self.write_procar([[
            (1.8, [[(0, 0, 1), (0, 0, 0)], [(0, 0, 0), (0, 0, 0)]]),
            (0.2, [[(0, 0, 0), (0, 0, 0)], [(0, 0, 1), (0, 0, 0)]]),
        ]])
        data = ranker.parse_procar(path)
        ranked = ranker.rank_procar(data, self.scf / "POSCAR", ["fe"], ["D"])["none"]
        self.assertEqual([item.band_index for item in ranked], [1, 2])
        self.assertAlmostEqual(ranked[0].bz_weighted_projection, 0.9)
        self.assertAlmostEqual(ranked[1].bz_weighted_projection, 0.1)
        self.assertAlmostEqual(ranked[0].kpoint_weight_sum, 2.0)

    def test_sums_elements_and_shells_as_cross_product(self) -> None:
        path = self.write_procar([[(1, [
            [(0.1, 0.2, 0.3), (0.4, 0.5, 0.6)],
            [(0, 0.1, 0), (0, 0.1, 0)],
        ])]])
        ranked = ranker.rank_procar(
            ranker.parse_procar(path), self.scf / "POSCAR", ["Fe", "O"], ["p", "d"]
        )["none"]
        self.assertEqual(ranked[0].band_index, 1)
        self.assertAlmostEqual(ranked[0].bz_weighted_projection, 1.6)

    def test_collinear_spins_are_ranked_independently(self) -> None:
        path = self.write_procar([
            [(1, [[(0, 0, 0.8), (0, 0, 0)], [(0, 0, 0.1), (0, 0, 0)]])],
            [(1, [[(0, 0, 0.1), (0, 0, 0)], [(0, 0, 0.9), (0, 0, 0)]])],
        ])
        data = ranker.parse_procar(path)
        rankings = ranker.rank_procar(data, self.scf / "POSCAR", ["Fe"], ["d"])
        self.assertEqual(data.spin_channels, ("up", "down"))
        self.assertEqual(rankings["up"][0].band_index, 1)
        self.assertEqual(rankings["down"][0].band_index, 2)

    def test_accepts_explicit_vasp_spin_component_headers(self) -> None:
        path = self.write_procar([
            [(1, [[(0, 0, 0.8), (0, 0, 0)]])],
            [(1, [[(0, 0, 0.2), (0, 0, 0)]])],
        ])
        text = path.read_text(encoding="utf-8")
        first = text.index(" k-point 1")
        second = text.index(" k-point 1", first + 1)
        text = text[:second] + " spin component 2\n" + text[second:]
        text = text[:first] + " spin component 1\n" + text[first:]
        path.write_text(text, encoding="utf-8")
        self.assertEqual(ranker.parse_procar(path).spin_channels, ("up", "down"))

    def test_accepts_adjacent_signed_kpoint_coordinates(self) -> None:
        path = self.write_procar([[(1, [[(0, 0, 1), (0, 0, 0)]])]])
        text = path.read_text(encoding="utf-8").replace(
            "k-point 1 : 0 0 0 weight",
            "k-point 1 : 0.38461538-0.00000000 0.00000000 weight",
        )
        path.write_text(text, encoding="utf-8")

        data = ranker.parse_procar(path)

        self.assertEqual(
            data.kpoints["none"][0].coordinates,
            (0.38461538, -0.0, 0.0),
        )

    def test_accepts_single_ion_procar_without_total_rows(self) -> None:
        path = self.write_procar([[
            (0.25, [[(0.9, 0, 0)], [(0, 0.8, 0)]]),
            (0.75, [[(0.7, 0, 0)], [(0, 0.6, 0)]]),
        ]])
        text = "\n".join(
            line for line in path.read_text(encoding="utf-8").splitlines()
            if not line.lstrip().startswith("tot ")
        )
        path.write_text(text + "\n", encoding="utf-8")

        data = ranker.parse_procar(path)

        self.assertEqual(data.nions, 1)
        self.assertEqual(data.nkpoints, 2)
        self.assertEqual(data.nbands, 2)
        self.assertAlmostEqual(
            data.kpoints["none"][1].bands[1].ion_projections[0][1], 0.6
        )

    def test_still_requires_total_rows_for_multiple_ions(self) -> None:
        path = self.write_procar([[(1, [[(0, 0, 1), (0, 0, 0)]])]])
        text = "\n".join(
            line for line in path.read_text(encoding="utf-8").splitlines()
            if not line.lstrip().startswith("tot ")
        )
        path.write_text(text + "\n", encoding="utf-8")

        with self.assertRaisesRegex(ranker.PBandError, "total projection row"):
            ranker.parse_procar(path)

    def test_ties_use_ascending_band_index(self) -> None:
        path = self.write_procar([[(1, [
            [(0, 0, 0.5), (0, 0, 0)], [(0, 0, 0.5), (0, 0, 0)],
        ])]])
        ranked = ranker.rank_procar(
            ranker.parse_procar(path), self.scf / "POSCAR", ["Fe"], ["d"]
        )["none"]
        self.assertEqual([item.band_index for item in ranked], [1, 2])

    def test_rejects_components_duplicates_and_missing_requests(self) -> None:
        data = ranker.parse_procar(self.write_procar([[(1, [[(0, 0, 1), (0, 0, 0)]])]]))
        for orbitals, message in ((["dxy"], "only aggregate"), (["d", "D"], "more than once"), (["f"], "absent")):
            with self.subTest(orbitals=orbitals), self.assertRaisesRegex(ranker.PBandError, message):
                ranker.rank_procar(data, self.scf / "POSCAR", ["Fe"], orbitals)
        with self.assertRaisesRegex(ranker.PBandError, "absent from"):
            ranker.rank_procar(data, self.scf / "POSCAR", ["Mn"], ["d"])

    def test_rejects_negative_and_zero_total_weights(self) -> None:
        text = procar_text([[(1, [[(0, 0, 1), (0, 0, 0)]])]])
        path = self.scf / "PROCAR"
        path.write_text(text.replace("weight = 1", "weight = -1"), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "nonnegative"):
            ranker.parse_procar(path)
        path.write_text(text.replace("weight = 1", "weight = 0"), encoding="utf-8")
        data = ranker.parse_procar(path)
        with self.assertRaisesRegex(ranker.PBandError, "nonpositive total"):
            ranker.rank_procar(data, self.scf / "POSCAR", ["Fe"], ["d"])
        path.write_text(text.replace("weight = 1", "weight = NaN"), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "expected a k-point block"):
            ranker.parse_procar(path)

    def test_rejects_lm_extra_tables_and_bad_layout(self) -> None:
        text = procar_text([[(1, [[(0, 0, 1), (0, 0, 0)]])]])
        path = self.scf / "PROCAR"
        path.write_text(text.replace("ion s p d tot", "ion s py pz px dxy tot"), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "only LORBIT=10"):
            ranker.parse_procar(path)
        path.write_text(text + "ion s p d tot\n", encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "noncollinear"):
            ranker.parse_procar(path)
        path.write_text("not a procar\n", encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "layout header"):
            ranker.parse_procar(path)
        path.write_text("\n".join(text.splitlines()[:-2]) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "missing ion|unexpected end"):
            ranker.parse_procar(path)

    def test_rejects_poscar_procar_ion_mismatch(self) -> None:
        data = ranker.parse_procar(self.write_procar([[(1, [[(0, 0, 1), (0, 0, 0)]])]]))
        (self.scf / "POSCAR").write_text(POSCAR.replace("1 1", "1 2"), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "contains 3"):
            ranker.rank_procar(data, self.scf / "POSCAR", ["Fe"], ["d"])

    def test_validates_dos_band_layout(self) -> None:
        data = ranker.parse_procar(self.write_procar([[(1, [[(0, 0, 1), (0, 0, 0)]])]]))
        (self.dos / "EIGENVAL").write_text(eigenval_text([[0, 1]]), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "1 bands.*2"):
            ranker.validate_procar_eigenval(data, ranker.parse_eigenval(self.dos / "EIGENVAL"))
        (self.dos / "EIGENVAL").write_text(eigenval_text([[(0.0, 0.1)]]), encoding="utf-8")
        with self.assertRaisesRegex(ranker.PBandError, "spin channels"):
            ranker.validate_procar_eigenval(data, ranker.parse_eigenval(self.dos / "EIGENVAL"))

    def test_cli_table_and_csv_report_weighted_fields(self) -> None:
        self.write_procar([[(1, [
            [(0, 0, 0.2), (0, 0, 0)], [(0, 0, 0.8), (0, 0, 0)],
        ])]])
        (self.dos / "EIGENVAL").write_text(eigenval_text([[-1, 0.5], [-0.5, 1.5]]), encoding="utf-8")
        csv_path = self.root / "ranking.csv"
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            result = ranker.main([
                "--scf-directory", str(self.scf), "--elements", "Fe",
                "--orbitals", "d", "--num-bands", "1", "--csv", str(csv_path),
            ])
        self.assertEqual(result, 0)
        self.assertIn("BZ-weighted projection", stdout.getvalue())
        with csv_path.open(encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["band_index"], "2")
        self.assertAlmostEqual(float(row["bz_weighted_projection"]), 0.8)
        self.assertEqual(row["scf_n_kpoints"], "1")
        self.assertEqual(row["kpoint_weight_sum"], "1")

    def test_cli_defaults_to_scf_and_dos_siblings(self) -> None:
        self.write_procar([[(1, [[(0, 0, 0.5), (0, 0, 0)]])]])
        (self.dos / "EIGENVAL").write_text(eigenval_text([[-2.5]]), encoding="utf-8")
        original = Path.cwd()
        stdout = io.StringIO()
        try:
            os.chdir(self.root)
            with contextlib.redirect_stdout(stdout):
                result = ranker.main(["--elements", "Fe", "--orbitals", "d"])
        finally:
            os.chdir(original)
        self.assertEqual(result, 0)
        self.assertIn("-2.5", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()

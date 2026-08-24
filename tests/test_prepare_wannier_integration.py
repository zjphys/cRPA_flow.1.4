#!/usr/bin/env python3

from __future__ import annotations

import shlex
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

SOURCE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE_DIR))

import prepare_wannier as preparer
from tests.test_rank_wannier_bands import procar_text


POSCAR = """Integration material
1.0
1 0 0
0 1 0
0 0 1
Mn Sb
1 1
Direct
0 0 0
0.5 0.5 0.5
"""

EIGENVAL = """header
header
header
header
header
  4  2  2

  0.0 0.0 0.0 0.5
  1 -1.0 1.0
  2  0.5 0.0

  0.5 0.0 0.0 0.5
  1 -0.5 1.0
  2  1.5 0.0
"""


class PrepareWannierIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.scf = self.root / "01_scf"
        self.dos = self.root / "02_dos"
        self.band = self.root / "03_band"
        self.scf.mkdir()
        self.dos.mkdir()
        self.band.mkdir()

        (self.scf / "INCAR").write_text(
            """SYSTEM = Integration SCF
ENCUT = 520
PREC = Accurate
EDIFF = 1e-6
KSPACING = 0.22
KGAMMA = .TRUE.
NBANDS = 36
""",
            encoding="utf-8",
        )
        (self.scf / "POSCAR").write_text(POSCAR, encoding="utf-8")
        (self.scf / "POTCAR").write_text("POTCAR data\n", encoding="utf-8")
        (self.scf / "CHGCAR").write_text("CHGCAR data\n", encoding="utf-8")
        (self.scf / "OUTCAR").write_text(
            "first run NBANDS = 36\nfinal setup NBANDS = 40\n",
            encoding="utf-8",
        )
        (self.scf / "PROCAR").write_text(
            procar_text([[
                (0.75, [
                    [(0.1, 0.2, 0.8), (0.1, 0.2, 0.1)],
                    [(0.1, 0.7, 0.3), (0.1, 0.7, 0.1)],
                ]),
                (0.25, [
                    [(0.1, 0.2, 0.8), (0.1, 0.2, 0.1)],
                    [(0.1, 0.7, 0.3), (0.1, 0.7, 0.1)],
                ]),
            ]]),
            encoding="utf-8",
        )
        (self.dos / "EIGENVAL").write_text(EIGENVAL, encoding="utf-8")
        for name in ("INCAR", "POSCAR", "DOSCAR", "EIGENVAL", "KPOINTS", "PROCAR"):
            (self.band / name).write_text(f"{name} input\n", encoding="utf-8")

        self.mock = self.root / "mock_vaspkit.py"
        self.mock.write_text(
            textwrap.dedent(
                """\
                from pathlib import Path
                import sys

                cwd = Path.cwd()
                log = Path(__file__).with_name("vaspkit.log")
                inputs = sys.stdin.read().splitlines() if not sys.argv[1:] else []
                with log.open("a", encoding="utf-8") as handle:
                    handle.write(" ".join(sys.argv[1:]) + "\\n")
                    if inputs:
                        handle.write("stdin: " + " ".join(inputs) + "\\n")

                if inputs == ["304", "3"]:
                    (cwd / "KPATH.wannier90").write_text(
                        "begin kpoint_path\\n"
                        "G 0 0 0 X 0.5 0 0\\n"
                        "end kpoint_path\\n",
                        encoding="utf-8",
                    )
                elif not sys.argv[1:]:
                    if len(inputs) != 3 or inputs[:2] != ["102", "2"]:
                        raise SystemExit(3)
                    (cwd / "KPOINTS").write_text(
                        f"Mock Gamma KPR {inputs[2]}\\n0\\nGamma\\n6 6 6\\n0 0 0\\n",
                        encoding="utf-8",
                    )
                else:
                    raise SystemExit(2)
                """
            ),
            encoding="utf-8",
        )
        self.vaspkit = f"{shlex.quote(sys.executable)} {shlex.quote(str(self.mock))}"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_prepares_complete_stage_and_force_replaces_it(self) -> None:
        stage, frontier, nbands, num_wann = preparer.prepare_wannier(
            self.root,
            ("Mn", "Sb"),
            ("d", "p"),
            2,
            0.05,
            self.vaspkit,
        )

        self.assertEqual(stage, self.root / "04_wann")
        self.assertEqual(nbands, 80)
        self.assertEqual(num_wann, 2)
        self.assertAlmostEqual(frontier.min_value, -1.0)
        self.assertAlmostEqual(frontier.max_value, 1.5)
        for name in (
            "INCAR",
            "POSCAR",
            "POTCAR",
            "CHGCAR",
            "KPOINTS",
            "wannier_band_ranking.csv",
            preparer.MARKER,
        ):
            self.assertTrue((stage / name).is_file(), name)
        self.assertEqual(
            (stage / "POSCAR").read_text(encoding="utf-8"),
            (self.scf / "POSCAR").read_text(encoding="utf-8"),
        )

        incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("ENCUT = 520", incar)
        self.assertIn("PREC = Accurate", incar)
        self.assertIn("EDIFF = 1e-6", incar)
        self.assertNotIn("KSPACING", incar)
        self.assertNotIn("KGAMMA", incar)
        self.assertIn("NBANDS = 80", incar)
        self.assertIn("NUM_WANN = 2", incar)
        self.assertIn("Mn:d\nSb:p", incar)
        self.assertIn("dis_froz_min = -0.9", incar)
        self.assertIn("dis_froz_max = 1.4", incar)
        self.assertIn("dis_win_min = -6", incar)
        self.assertIn("dis_win_max = 6.5", incar)
        self.assertIn("begin kpoint_path", incar)
        self.assertIn("Mock Gamma KPR 0.05", (stage / "KPOINTS").read_text())

        log = (self.root / "vaspkit.log").read_text(encoding="utf-8")
        self.assertNotIn("-task 213", log)
        self.assertIn("stdin: 304 3", log)

        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "already exists"
        ):
            preparer.prepare_wannier(
                self.root,
                ("Mn", "Sb"),
                ("d", "p"),
                2,
                0.06,
                self.vaspkit,
            )

        preparer.prepare_wannier(
            self.root,
            ("Mn", "Sb"),
            ("d", "p"),
            2,
            0.06,
            self.vaspkit,
            force=True,
            frozen_margin=0.1,
        )
        self.assertIn(
            "Mock Gamma KPR 0.06",
            (stage / "KPOINTS").read_text(encoding="utf-8"),
        )
        regenerated_incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("dis_froz_min = -0.9", regenerated_incar)
        self.assertIn("dis_froz_max = 1.4", regenerated_incar)

    def test_noncontiguous_ranking_uses_largest_run_without_reducing_selection(
        self,
    ) -> None:
        eigenval = """header
header
header
header
header
  4  2  5

  0.0 0.0 0.0 0.5
  1 -3.0 1.0
  2 -1.0 1.0
  3 0.5 0.0
  4 2.0 0.0
  5 4.0 0.0

  0.5 0.0 0.0 0.5
  1 -2.5 1.0
  2 -0.5 1.0
  3 1.0 0.0
  4 2.5 0.0
  5 4.5 0.0
"""
        (self.dos / "EIGENVAL").write_text(eigenval, encoding="utf-8")
        bands = [
            [(0, 0, weight), (0, 0, 0)]
            for weight in (0.1, 0.9, 0.8, 0.2, 0.7)
        ]
        (self.scf / "PROCAR").write_text(
            procar_text([[(0.5, bands), (0.5, bands)]]), encoding="utf-8"
        )

        stage, frontier, _, num_wann = preparer.prepare_wannier(
            self.root,
            ("Mn",),
            ("d",),
            3,
            0.05,
            self.vaspkit,
        )

        self.assertEqual(3, num_wann)
        self.assertAlmostEqual(-1.0, frontier.min_value)
        self.assertAlmostEqual(4.5, frontier.max_value)
        incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("NUM_WANN = 3", incar)
        self.assertIn("dis_froz_min = -0.9", incar)
        self.assertIn("dis_froz_max = 0.9", incar)
        self.assertIn("dis_win_min = -6", incar)
        self.assertIn("dis_win_max = 6", incar)
        ranking = (stage / "wannier_band_ranking.csv").read_text(encoding="utf-8")
        self.assertIn("1,2,none,", ranking)
        self.assertIn("2,3,none,", ranking)
        self.assertIn("3,5,none,", ranking)

    def test_prepares_from_spin_resolved_scf_procar(self) -> None:
        spin_eigenval = """header
header
header
header
header
  4  2  2

  0.0 0.0 0.0 0.5
  1 -1.0 -0.8 1.0 1.0
  2 0.5 0.7 0.0 0.0

  0.5 0.0 0.0 0.5
  1 -0.5 -0.3 1.0 1.0
  2 1.5 1.7 0.0 0.0
"""
        (self.dos / "EIGENVAL").write_text(spin_eigenval, encoding="utf-8")
        with (self.scf / "INCAR").open("a", encoding="utf-8") as handle:
            handle.write("ISPIN = 2\n")

        up_bands = [
            [(0, 0.1, 0.8), (0, 0.1, 0.8)],
            [(0, 0.1, 0.1), (0, 0.1, 0.1)],
        ]
        down_bands = [
            [(0, 0.1, 0.1), (0, 0.1, 0.1)],
            [(0, 0.8, 0.1), (0, 0.8, 0.1)],
        ]
        (self.scf / "PROCAR").write_text(
            procar_text([
                [(0.5, up_bands), (0.5, up_bands)],
                [(0.5, down_bands), (0.5, down_bands)],
            ]),
            encoding="utf-8",
        )

        stage, frontier, _, num_wann = preparer.prepare_wannier(
            self.root,
            ("Mn", "Sb"),
            ("d", "p"),
            1,
            0.05,
            self.vaspkit,
        )

        self.assertAlmostEqual(-1.0, frontier.min_value)
        self.assertAlmostEqual(1.7, frontier.max_value)
        self.assertEqual(num_wann, 1)
        incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("ISPIN = 2", incar)
        self.assertIn("NUM_WANN = 1", incar)
        self.assertIn("dis_froz_min = -0.2", incar)
        self.assertIn("dis_froz_max = 0.4", incar)
        self.assertIn("dis_win_min = -6", incar)
        self.assertIn("dis_win_max = 6.7", incar)
        ranking = (stage / "wannier_band_ranking.csv").read_text(encoding="utf-8")
        self.assertIn("1,1,up,", ranking)
        self.assertIn("1,2,down,", ranking)
        self.assertNotIn("1,2,up,", ranking)
        self.assertNotIn("1,1,down,", ranking)
        log = (self.root / "vaspkit.log").read_text(encoding="utf-8")
        self.assertNotIn("-task 213", log)

    def test_infers_num_wann_when_num_bands_is_omitted(self) -> None:
        stage, _, _, num_wann = preparer.prepare_wannier(
            self.root,
            ("mn",),
            ("s",),
            None,
            0.05,
            self.vaspkit,
        )

        self.assertEqual(num_wann, 1)
        self.assertIn(
            "NUM_WANN = 1",
            (stage / "INCAR").read_text(encoding="utf-8"),
        )

    def test_requires_scf_procar(self) -> None:
        (self.scf / "PROCAR").unlink()
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "SCF PROCAR.*missing or empty"
        ):
            preparer.prepare_wannier(
                self.root,
                ("Mn", "Sb"),
                ("d", "p"),
                2,
                0.04,
                self.vaspkit,
            )

    def test_force_refuses_unowned_stage(self) -> None:
        stage = self.root / "04_wann"
        stage.mkdir()
        manual_input = stage / "manual-input"
        manual_input.write_text("keep me\n", encoding="utf-8")

        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "not generated by this workflow"
        ):
            preparer.prepare_wannier(
                self.root,
                ("Mn", "Sb"),
                ("d", "p"),
                2,
                0.04,
                self.vaspkit,
                force=True,
            )
        self.assertEqual(manual_input.read_text(encoding="utf-8"), "keep me\n")


if __name__ == "__main__":
    unittest.main()

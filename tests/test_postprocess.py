from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


VERSION_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VERSION_DIR))

import postprocess


class PostprocessSpinTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.band_dir = self.root / "03_band"
        self.dos_dir = self.root / "02_dos"
        self.band_dir.mkdir()
        self.dos_dir.mkdir()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_spin_data(self, spin: str, energy_shift: float) -> None:
        (self.band_dir / f"PBAND_Fe_{spin}.dat").write_text(
            "# K-Path Energy s tot\n"
            "# NKPTS & NBANDS: 2 1\n"
            "# Band-Index: 1\n"
            f"0.0 {-1.0 + energy_shift} 0.2 0.2\n"
            f"1.0 {1.0 + energy_shift} 0.4 0.4\n",
            encoding="utf-8",
        )
        (self.dos_dir / f"PDOS_Fe_{spin}.dat").write_text(
            "# Energy s tot\n"
            "-2.0 0.1 0.1\n"
            "0.0 0.8 0.8\n"
            "2.0 0.2 0.2\n",
            encoding="utf-8",
        )

    def test_discovers_up_and_down_as_one_element(self) -> None:
        self.write_spin_data("UP", 0.0)
        self.write_spin_data("DW", 0.3)

        files = postprocess.discover_files(self.band_dir, "PBAND")

        self.assertEqual({"fe"}, set(files))
        self.assertEqual({"up", "down"}, set(files["fe"]))

    def test_spin_resolved_files_generate_a_plot(self) -> None:
        self.write_spin_data("UP", 0.0)
        self.write_spin_data("DW", 0.3)
        (self.band_dir / "KLABELS").write_text(
            "GAMMA 0.0\nX 1.0\n", encoding="utf-8"
        )

        result = subprocess.run(
            [
                sys.executable,
                str(VERSION_DIR / "postprocess.py"),
                "--root",
                str(self.root),
                "--reuse-data",
                "--format",
                "png",
                "--dpi",
                "50",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

        self.assertEqual(0, result.returncode, result.stdout)
        self.assertTrue((self.root / "projected_band_dos.png").is_file())
        self.assertIn("Elements: Fe", result.stdout)

    def test_spin_resolved_orbital_plot(self) -> None:
        self.write_spin_data("UP", 0.0)
        self.write_spin_data("DW", 0.3)
        (self.band_dir / "KLABELS").write_text(
            "GAMMA 0.0\nX 1.0\n", encoding="utf-8"
        )

        result = subprocess.run(
            [
                sys.executable,
                str(VERSION_DIR / "postprocess.py"),
                "--root",
                str(self.root),
                "--reuse-data",
                "--format",
                "png",
                "--dpi",
                "50",
                "--orbital-element",
                "Fe",
                "--orbitals",
                "s",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

        self.assertEqual(0, result.returncode, result.stdout)
        self.assertTrue((self.root / "projected_band_dos.png").is_file())
        self.assertIn("Projection: orbital-resolved", result.stdout)

    def test_incomplete_spin_pair_is_rejected(self) -> None:
        self.write_spin_data("UP", 0.0)
        band_files = postprocess.discover_files(self.band_dir, "PBAND")
        dos_files = postprocess.discover_files(self.dos_dir, "PDOS")

        with self.assertRaisesRegex(RuntimeError, r"missing.*_DW\.dat"):
            postprocess.select_elements(band_files, dos_files, None)

    def test_kpoint_labels_are_converted_to_latex(self) -> None:
        labels_path = self.band_dir / "KLABELS"
        labels_path.write_text(
            "GAMMA 0.0\n"
            "X 1.0\n"
            "M_1 2.0\n"
            "DELTA 3.0\n"
            "X 4.0\n"
            "M 4.0\n",
            encoding="utf-8",
        )

        coordinates, labels = postprocess.read_klabels(labels_path)

        self.assertEqual([0.0, 1.0, 2.0, 3.0, 4.0], coordinates)
        self.assertEqual(
            [r"$\Gamma$", "$X$", "$M_{1}$", r"$\Delta$", r"$X\mid M$"],
            labels,
        )


class PostprocessWannierTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.band_dir = self.root / "03_band"
        self.dos_dir = self.root / "02_dos"
        self.wannier_dir = self.root / "04_wann"
        self.band_dir.mkdir()
        self.dos_dir.mkdir()
        self.wannier_dir.mkdir()
        (self.band_dir / "KLABELS").write_text(
            "GAMMA 0.0\nX 1.0\n", encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def wannier_text(energy_shift: float = 0.0) -> str:
        return (
            f"0.0 {4.0 + energy_shift}\n"
            f"2.0 {6.0 + energy_shift}\n\n"
            f"0.0 {5.0 + energy_shift}\n"
            f"2.0 {7.0 + energy_shift}\n"
        )

    def write_plain_projected_data(self) -> None:
        (self.band_dir / "PBAND_Fe.dat").write_text(
            "# K-Path Energy s tot\n"
            "# NKPTS & NBANDS: 2 1\n"
            "# Band-Index: 1\n"
            "0.0 -1.0 0.2 0.2\n"
            "1.0 1.0 0.4 0.4\n",
            encoding="utf-8",
        )
        (self.dos_dir / "PDOS_Fe.dat").write_text(
            "# Energy s tot\n-2.0 0.1 0.1\n0.0 0.8 0.8\n2.0 0.2 0.2\n",
            encoding="utf-8",
        )

    def write_spin_projected_data(self) -> None:
        for suffix, shift in (("UP", 0.0), ("DW", 0.2)):
            (self.band_dir / f"PBAND_Fe_{suffix}.dat").write_text(
                "# K-Path Energy s tot\n"
                "# NKPTS & NBANDS: 2 1\n"
                "# Band-Index: 1\n"
                f"0.0 {-1.0 + shift} 0.2 0.2\n"
                f"1.0 {1.0 + shift} 0.4 0.4\n",
                encoding="utf-8",
            )
            (self.dos_dir / f"PDOS_Fe_{suffix}.dat").write_text(
                "# Energy s tot\n-2.0 0.1 0.1\n0.0 0.8 0.8\n2.0 0.2 0.2\n",
                encoding="utf-8",
            )

    def run_postprocess(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERSION_DIR / "postprocess.py"),
                "--root",
                str(self.root),
                "--reuse-data",
                "--wannier-bands",
                "--format",
                "png",
                "--dpi",
                "50",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def test_effective_ispin_uses_default_and_last_assignment(self) -> None:
        incar = self.wannier_dir / "INCAR"
        incar.write_text("ENCUT = 400\n", encoding="utf-8")
        self.assertEqual(1, postprocess.read_effective_ispin(incar))
        incar.write_text(
            "ISPIN = 1 ; ENCUT = 400\nISPIN = 2 ! effective\n",
            encoding="utf-8",
        )
        self.assertEqual(2, postprocess.read_effective_ispin(incar))

    def test_band_blocks_are_preserved_and_aligned(self) -> None:
        path = self.wannier_dir / "wannier90_band.dat"
        path.write_text(self.wannier_text(), encoding="utf-8")
        blocks = postprocess.load_wannier_band_blocks(path)
        aligned = postprocess.align_wannier_bands(
            {"plain": blocks}, fermi=5.0, target_min=0.0, target_max=1.0
        )

        self.assertEqual(2, len(blocks))
        np.testing.assert_allclose(aligned["plain"][0][:, 0], [0.0, 1.0])
        np.testing.assert_allclose(aligned["plain"][0][:, 1], [-1.0, 1.0])

    def test_non_spin_wannier_overlay_generates_plot(self) -> None:
        self.write_plain_projected_data()
        (self.wannier_dir / "INCAR").write_text("ENCUT = 400\n", encoding="utf-8")
        (self.wannier_dir / "OUTCAR").write_text(
            "E-fermi : 4.5\nE-fermi : 5.0\n", encoding="utf-8"
        )
        (self.wannier_dir / "wannier90_band.dat").write_text(
            self.wannier_text(), encoding="utf-8"
        )

        result = self.run_postprocess()

        self.assertEqual(0, result.returncode, result.stdout)
        self.assertTrue((self.root / "projected_band_dos.png").is_file())
        self.assertIn("Wannier bands: E-fermi = 5 eV", result.stdout)

    def test_spin_wannier_overlay_generates_plot(self) -> None:
        self.write_spin_projected_data()
        (self.wannier_dir / "INCAR").write_text("ISPIN = 2\n", encoding="utf-8")
        (self.wannier_dir / "OUTCAR").write_text(
            "E-fermi : 5.0\n", encoding="utf-8"
        )
        (self.wannier_dir / "wannier90.1_band.dat").write_text(
            self.wannier_text(), encoding="utf-8"
        )
        (self.wannier_dir / "wannier90.2_band.dat").write_text(
            self.wannier_text(0.2), encoding="utf-8"
        )

        result = self.run_postprocess()

        self.assertEqual(0, result.returncode, result.stdout)
        self.assertTrue((self.root / "projected_band_dos.png").is_file())

    def test_missing_spin_file_is_reported(self) -> None:
        self.write_spin_projected_data()
        (self.wannier_dir / "INCAR").write_text("ISPIN = 2\n", encoding="utf-8")
        (self.wannier_dir / "OUTCAR").write_text(
            "E-fermi : 5.0\n", encoding="utf-8"
        )
        (self.wannier_dir / "wannier90.1_band.dat").write_text(
            self.wannier_text(), encoding="utf-8"
        )

        result = self.run_postprocess()

        self.assertNotEqual(0, result.returncode)
        self.assertIn("wannier90.2_band.dat", result.stdout)
        self.assertIn("missing or empty", result.stdout)

    def test_malformed_and_inconsistent_band_blocks_are_rejected(self) -> None:
        path = self.wannier_dir / "wannier90_band.dat"
        path.write_text("0.0 1.0\n1.0 nope\n", encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError, "non-numeric"):
            postprocess.load_wannier_band_blocks(path)

        path.write_text(
            "0.0 1.0\n1.0 2.0\n\n0.0 3.0\n2.0 4.0\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(RuntimeError, "same k-grid"):
            postprocess.load_wannier_band_blocks(path)

    def test_missing_fermi_and_unsupported_ispin_are_rejected(self) -> None:
        (self.wannier_dir / "INCAR").write_text("ISPIN = 3\n", encoding="utf-8")
        with self.assertRaisesRegex(RuntimeError, "unsupported ISPIN"):
            postprocess.load_wannier_bands(self.wannier_dir)

        (self.wannier_dir / "INCAR").write_text("ISPIN = 1\n", encoding="utf-8")
        (self.wannier_dir / "wannier90_band.dat").write_text(
            self.wannier_text(), encoding="utf-8"
        )
        (self.wannier_dir / "OUTCAR").write_text(
            "no Fermi value here\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(RuntimeError, "no E-fermi"):
            postprocess.load_wannier_bands(self.wannier_dir)


if __name__ == "__main__":
    unittest.main()

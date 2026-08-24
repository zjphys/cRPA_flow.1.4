#!/usr/bin/env python3

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

SOURCE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE_DIR))

import prepare_wannier as preparer
import rank_wannier_bands as ranker


KPOINT_PATH = """header
begin kpoint_path
G 0 0 0 X 0.5 0 0
end kpoint_path
footer
"""


def eigenval_data(
    energies: dict[str, dict[int, tuple[float, ...]]]
) -> ranker.EigenvalData:
    first_spin = next(iter(energies.values()))
    return ranker.EigenvalData(
        path=Path("EIGENVAL"),
        nkpoints=len(next(iter(first_spin.values()))),
        nbands=len(first_spin),
        spin_channels=tuple(energies),
        energies=energies,
    )


class PrepareWannierUnitTests(unittest.TestCase):
    @staticmethod
    def write_poscar(path: Path, elements: str, counts: str) -> None:
        path.write_text(
            "Test material\n1.0\n1 0 0\n0 1 0\n0 0 1\n"
            f"{elements}\n{counts}\nDirect\n0 0 0\n",
            encoding="utf-8",
        )

    def test_reads_last_effective_nbands_from_outcar(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "OUTCAR"
            path.write_text(
                " dimension x,y,z NGX = 1\n NBANDS= 40\n"
                " restart\n NBANDS =     48\n",
                encoding="utf-8",
            )
            self.assertEqual(preparer.read_effective_nbands(path), 48)

    def test_rejects_outcar_without_nbands(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "OUTCAR"
            path.write_text("no matching setting\n", encoding="utf-8")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "could not find"
            ):
                preparer.read_effective_nbands(path)

    def test_pairs_elements_and_orbitals_by_position(self) -> None:
        self.assertEqual(
            preparer.paired_projections(["Mn", "Sb"], ["d", "p"]),
            ("Mn:d", "Sb:p"),
        )
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "same number"
        ):
            preparer.paired_projections(["Mn", "Sb"], ["d"])

    def test_reads_vasp5_poscar_atom_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "POSCAR"
            self.write_poscar(path, "U O", "2 4")
            self.assertEqual(preparer.read_poscar_atom_counts(path), {"u": 2, "o": 4})

    def test_rejects_malformed_and_vasp4_poscars(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "POSCAR"
            self.write_poscar(path, "2 4", "Direct")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "VASP 4-style"
            ):
                preparer.read_poscar_atom_counts(path)

            self.write_poscar(path, "U O", "2")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "element symbols"
            ):
                preparer.read_poscar_atom_counts(path)

    def test_infers_shell_num_wann_from_paired_elements(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "POSCAR"
            self.write_poscar(path, "U Mn Sb", "2 3 4")
            self.assertEqual(
                preparer.infer_num_wann(path, ("u",), ("F",)),
                14,
            )
            self.assertEqual(
                preparer.infer_num_wann(path, ("Mn", "Sb"), ("d", "p")),
                27,
            )

    def test_rejects_component_orbital(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            poscar = directory / "POSCAR"
            self.write_poscar(poscar, "U", "3")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "aggregate.*shell"
            ):
                preparer.infer_num_wann(poscar, ("U",), ("fxyz",))

    def test_inference_rejects_unknown_non_shell_orbital(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            poscar = directory / "POSCAR"
            self.write_poscar(poscar, "U", "1")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "aggregate.*shell"
            ):
                preparer.infer_num_wann(poscar, ("U",), ("custom",))

    def test_inference_rejects_element_absent_from_poscar(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "POSCAR"
            self.write_poscar(path, "O", "2")
            with self.assertRaisesRegex(
                preparer.WannierPreparationError, "absent from SCF POSCAR"
            ):
                preparer.infer_num_wann(path, ("U",), ("f",))

    def test_num_bands_cli_option_is_optional(self) -> None:
        args = preparer.parse_args(["--elements", "U", "--orbitals", "f"])
        self.assertIsNone(args.num_bands)

    def test_calculates_isolated_windows_from_neighboring_bands(self) -> None:
        eigenval = eigenval_data(
            {
                "none": {
                    1: (-2.0, -1.0),
                    2: (-0.5, -1.5),
                    3: (2.0, 1.0),
                    4: (1.5, 3.0),
                }
            }
        )
        windows = preparer.calculate_windows(eigenval, {"none": (2, 3)}, 2)
        self.assertAlmostEqual(windows.dis_froz_min, -0.9)
        self.assertAlmostEqual(windows.dis_froz_max, 1.4)
        self.assertAlmostEqual(windows.dis_win_min, -6.5)
        self.assertAlmostEqual(windows.dis_win_max, 7.0)

    def test_handles_target_block_at_eigenval_edges(self) -> None:
        eigenval = eigenval_data(
            {"none": {1: (-1.0, -0.5), 2: (0.5, 1.0), 3: (2.0, 2.5)}}
        )
        lower = preparer.calculate_windows(
            eigenval,
            {"none": (1, 2)},
            2,
            0.1,
        )
        self.assertAlmostEqual(lower.dis_froz_min, -0.9)
        self.assertAlmostEqual(lower.dis_froz_max, 0.9)

        upper = preparer.calculate_windows(
            eigenval,
            {"none": (2, 3)},
            2,
            0.1,
        )
        self.assertAlmostEqual(upper.dis_froz_min, 0.6)
        self.assertAlmostEqual(upper.dis_froz_max, 2.4)

    def test_uses_largest_contiguous_target_interval(self) -> None:
        eigenval = eigenval_data(
            {
                "none": {
                    1: (-4.0, -3.5),
                    2: (-2.0, -1.0),
                    3: (0.0, 1.0),
                    4: (3.0, 3.5),
                    5: (5.0, 5.5),
                }
            }
        )
        windows = preparer.calculate_windows(
            eigenval,
            {"none": (2, 3, 5)},
            3,
        )

        self.assertAlmostEqual(-1.9, windows.dis_froz_min)
        self.assertAlmostEqual(0.9, windows.dis_froz_max)
        self.assertAlmostEqual(-7.0, windows.dis_win_min)
        self.assertAlmostEqual(6.0, windows.dis_win_max)

    def test_equal_length_intervals_prefer_highest_ranked_band(self) -> None:
        eigenval = eigenval_data(
            {
                "none": {
                    1: (-6.0,),
                    2: (-4.0,),
                    3: (-3.0,),
                    4: (-1.0,),
                    5: (1.0,),
                    6: (2.0,),
                    7: (4.0,),
                }
            }
        )
        windows = preparer.calculate_windows(
            eigenval,
            {"none": (5, 2, 3, 6)},
            4,
        )

        self.assertAlmostEqual(1.1, windows.dis_froz_min)
        self.assertAlmostEqual(1.9, windows.dis_froz_max)
        self.assertAlmostEqual(-4.0, windows.dis_win_min)
        self.assertAlmostEqual(7.0, windows.dis_win_max)

    def test_selected_band_outside_window_interval_is_not_a_window_target(self) -> None:
        eigenval = eigenval_data(
            {
                "none": {
                    1: (-2.0,),
                    2: (0.0,),
                    3: (1.0,),
                    4: (2.0,),
                    5: (1.0,),
                }
            }
        )
        with self.assertRaises(preparer.WannierPreparationError) as caught:
            preparer.calculate_windows(
                eigenval,
                {"none": (2, 3, 5)},
                3,
                0.0,
            )
        self.assertIn("non-target bands enter", str(caught.exception))
        self.assertIn("[5]", str(caught.exception))
        self.assertIn("window bands=none:[2, 3]", str(caught.exception))

    def test_rejects_empty_isolated_window(self) -> None:
        eigenval = eigenval_data({"none": {1: (1.0, 1.4)}})
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "empty or inverted"
        ):
            preparer.calculate_windows(eigenval, {"none": (1,)}, 1, 0.2)

    def test_rejects_invalid_frozen_margin(self) -> None:
        eigenval = eigenval_data({"none": {1: (-1.0, 1.0)}})
        for margin in (-0.1, math.inf, math.nan):
            with self.subTest(margin=margin), self.assertRaisesRegex(
                preparer.WannierPreparationError, "finite and nonnegative"
            ):
                preparer.calculate_windows(
                    eigenval, {"none": (1,)}, 1, margin
                )

    def test_zero_margin_rejects_non_target_on_inclusive_boundary(self) -> None:
        eigenval = eigenval_data(
            {"none": {1: (-2.0, -1.0), 2: (-1.0, 1.0)}}
        )
        with self.assertRaises(preparer.WannierPreparationError) as caught:
            preparer.calculate_windows(eigenval, {"none": (2,)}, 1, 0.0)
        self.assertIn("non-target bands enter", str(caught.exception))
        self.assertIn("exceeding NUM_WANN=1", str(caught.exception))

    def test_zero_margin_is_accepted_without_neighbor_invasion(self) -> None:
        eigenval = eigenval_data({"none": {1: (-1.0, 1.0)}})
        windows = preparer.calculate_windows(
            eigenval, {"none": (1,)}, 1, 0.0
        )
        self.assertEqual(windows.dis_froz_min, -1.0)
        self.assertEqual(windows.dis_froz_max, 1.0)

    def test_parses_frozen_margin_cli_option(self) -> None:
        args = preparer.parse_args(
            [
                "--elements",
                "Mn",
                "--orbitals",
                "d",
                "--num-bands",
                "5",
                "--frozen-margin",
                "0.15",
            ]
        )
        self.assertAlmostEqual(args.frozen_margin, 0.15)

    def test_validates_spin_channels_independently(self) -> None:
        eigenval = eigenval_data(
            {
                "up": {1: (-2.0, -2.0), 2: (-1.0, 1.0)},
                "down": {1: (-0.9, -0.9), 2: (-1.0, 1.0)},
            }
        )
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "spin down"
        ):
            preparer.calculate_windows(
                eigenval, {"up": (2,), "down": (2,)}, 1, 0.0
            )

    def test_intersects_independent_spin_target_and_guard_ranges(self) -> None:
        eigenval = eigenval_data(
            {
                "up": {
                    1: (-1.2, -0.9),
                    2: (-1.5, -0.4),
                    3: (0.2, 0.8),
                    4: (1.3, 1.6),
                    5: (2.5, 2.8),
                },
                "down": {
                    1: (-2.0, -1.8),
                    2: (-1.3, -1.0),
                    3: (-0.8, 0.1),
                    4: (1.0, 2.0),
                    5: (1.5, 2.4),
                },
            }
        )

        windows = preparer.calculate_windows(
            eigenval,
            {"up": (2, 3), "down": (3, 4)},
            2,
            0.1,
        )

        self.assertAlmostEqual(-0.8, windows.dis_froz_min)
        self.assertAlmostEqual(1.2, windows.dis_froz_max)
        self.assertAlmostEqual(-6.5, windows.dis_win_min)
        self.assertAlmostEqual(7.0, windows.dis_win_max)

    def test_rejects_mismatched_procar_and_eigenval_spin_channels(self) -> None:
        eigenval = eigenval_data(
            {"up": {1: (-1.0,)}, "down": {1: (-0.5,)}}
        )

        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "do not match EIGENVAL"
        ):
            preparer.calculate_windows(eigenval, {"up": (1,)}, 1)

    def test_selects_largest_contiguous_interval_independently_per_spin(self) -> None:
        eigenval = eigenval_data(
            {
                "up": {
                    1: (-3.0,),
                    2: (-1.0,),
                    3: (1.0,),
                    4: (3.0,),
                    5: (4.0,),
                },
                "down": {
                    1: (-1.0,),
                    2: (1.0,),
                    3: (3.0,),
                    4: (4.0,),
                    5: (5.0,),
                },
            }
        )

        windows = preparer.calculate_windows(
            eigenval,
            {"up": (2, 3, 5), "down": (4, 1, 2)},
            3,
        )

        self.assertAlmostEqual(-0.9, windows.dis_froz_min)
        self.assertAlmostEqual(0.9, windows.dis_froz_max)
        self.assertAlmostEqual(-6.0, windows.dis_win_min)
        self.assertAlmostEqual(6.0, windows.dis_win_max)

    def test_rejects_empty_cross_spin_guard_intersection(self) -> None:
        eigenval = eigenval_data(
            {
                "up": {1: (1.0,), 2: (0.0,)},
                "down": {1: (0.0,), 2: (-1.0,)},
            }
        )

        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "empty or inverted"
        ):
            preparer.calculate_windows(
                eigenval, {"up": (2,), "down": (1,)}, 1, 0.0
            )

    def test_extracts_exactly_one_kpoint_path(self) -> None:
        self.assertEqual(
            preparer.extract_kpoint_path(KPOINT_PATH),
            "begin kpoint_path\nG 0 0 0 X 0.5 0 0\nend kpoint_path",
        )
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "exactly one"
        ):
            preparer.extract_kpoint_path("begin kpoint_path\nG 0 0 0\n")
        with self.assertRaisesRegex(
            preparer.WannierPreparationError, "exactly one"
        ):
            preparer.extract_kpoint_path(KPOINT_PATH + KPOINT_PATH)

    def test_rewrites_scf_incar_and_preserves_unrelated_settings(self) -> None:
        source = """SYSTEM = Example SCF
ENCUT = 500
PREC = Accurate
ISMEAR = 0; SIGMA = 0.05
KSPACING = 0.22
KGAMMA = .TRUE.
NBANDS = 40
# keep this comment
"""
        windows = preparer.WannierWindows(-1.3, 1.8, -6.5, 7.0)
        result = preparer.rewrite_incar(
            source,
            "Example",
            96,
            2,
            ("Mn:d", "Sb:p"),
            windows,
            KPOINT_PATH,
        )

        self.assertIn("ENCUT = 500", result)
        self.assertIn("PREC = Accurate", result)
        self.assertIn("ISMEAR = 0; SIGMA = 0.05", result)
        self.assertIn("# keep this comment", result)
        self.assertNotIn("KSPACING", result)
        self.assertNotIn("KGAMMA", result)
        self.assertEqual(result.count("NBANDS ="), 1)
        self.assertIn("NBANDS = 96", result)
        self.assertIn("NUM_WANN = 2", result)
        self.assertIn("ICHARG = 11", result)
        self.assertIn("ISYM = -1", result)
        self.assertIn("Mn:d\nSb:p", result)
        self.assertIn("dis_froz_min = -1.3", result)
        self.assertIn("dis_froz_max = 1.8", result)
        self.assertIn("dis_win_min = -6.5", result)
        self.assertIn("dis_win_max = 7", result)
        self.assertIn("begin kpoint_path", result)
        self.assertIn("bands_num_points = 20", result)

    def test_removes_old_multiline_wannier_block(self) -> None:
        source = """ENCUT = 400
WANNIER90_WIN = "
old = true
"
PREC = high
"""
        result = preparer.rewrite_incar(
            source,
            "Example",
            20,
            1,
            ("Mn:d",),
            preparer.WannierWindows(-1.0, 1.0, -6.0, 6.0),
            KPOINT_PATH,
        )
        self.assertNotIn("old = true", result)
        self.assertIn("ENCUT = 400", result)
        self.assertIn("PREC = high", result)
        self.assertEqual(result.count("WANNIER90_WIN"), 1)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

SOURCE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE_DIR))

import prepare_crpa as preparer


POSCAR = """cRPA integration material
1.0
1 0 0
0 1 0
0 0 1
Mn Se
1 2
Direct
0 0 0
0.25 0.25 0.25
0.75 0.75 0.75
"""


class PrepareCrpaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.wannier = self.root / "04_wann"
        self.wannier.mkdir()
        (self.wannier / preparer.MARKER).touch()
        (self.wannier / "POSCAR").write_text(POSCAR, encoding="utf-8")
        (self.wannier / "INCAR").write_text(
            """SYSTEM = cRPA integration material Wannier
ENCUT = 500
EDIFF = 1D-06
NBANDS = 72
NUM_WANN = 10
WANNIER90_WIN = "
begin projections
Mn:d
end projections
"
""",
            encoding="utf-8",
        )
        (self.wannier / "OUTCAR").write_text(
            """initial setup NBANDS = 72
parallel setup NBANDS = 80
General timing and accounting informations for this job:
""",
            encoding="utf-8",
        )
        for name in preparer.REQUIRED_WANNIER_FILES:
            if name == "POSCAR":
                continue
            content = (
                "# ISPIN NKPTS NB_TOT NW\n1 1 80 10\n"
                if name == "WANPROJ"
                else f"{name} restart data\n"
            )
            (self.wannier / name).write_text(content, encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_prepares_defaults_and_preserves_target_order(self) -> None:
        stage, settings = preparer.prepare_crpa(
            self.root, [3, 1, 2], kpar=8
        )

        self.assertEqual(stage, self.root / "05_crpa")
        self.assertEqual(settings.nbands, 80)
        self.assertEqual(settings.nbandsgw, 80)
        self.assertTrue(math.isclose(settings.encutgw, 1000.0 / 3.0))
        self.assertEqual(settings.target_states, (3, 1, 2))
        self.assertEqual(settings.num_wann, 10)
        self.assertEqual(settings.ispin, 1)
        self.assertEqual(settings.kpar, 8)
        self.assertFalse(settings.copied_waveder)

        incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("SYSTEM = cRPA integration material cRPA", incar)
        self.assertIn("NBANDS = 80", incar)
        self.assertIn("NBANDSGW = 80", incar)
        self.assertNotIn("ISPIN", incar)
        self.assertIn("NTARGET_STATES = 3 1 2", incar)
        self.assertIn("ENCUTGW = 333.333333333", incar)
        self.assertIn("KPAR = 8", incar)
        self.assertTrue((stage / preparer.MARKER).is_file())
        self.assertTrue((stage / "WANPROJ").is_file())
        self.assertFalse((stage / "OUTCAR").exists())
        self.assertFalse((stage / "WAVEDER").exists())

    def test_defaults_target_states_to_all_wannier_states(self) -> None:
        stage, settings = preparer.prepare_crpa(self.root)

        self.assertEqual(settings.target_states, tuple(range(1, 11)))
        self.assertIn(
            "NTARGET_STATES = 1 2 3 4 5 6 7 8 9 10",
            (stage / "INCAR").read_text(encoding="utf-8"),
        )

    def test_target_states_cli_option_is_optional(self) -> None:
        args = preparer.parse_args([])
        self.assertIsNone(args.target_states)

    def test_accepts_overrides_and_copies_waveder(self) -> None:
        (self.wannier / "WAVEDER").write_text(
            "long-wave derivatives\n", encoding="utf-8"
        )

        stage, settings = preparer.prepare_crpa(
            self.root,
            [1, 10],
            nbandsgw=64,
            encutgw=280.5,
            kpar=4,
        )

        self.assertEqual(settings.nbandsgw, 64)
        self.assertEqual(settings.encutgw, 280.5)
        self.assertTrue(settings.copied_waveder)
        self.assertEqual(
            (stage / "WAVEDER").read_text(encoding="utf-8"),
            "long-wave derivatives\n",
        )
        incar = (stage / "INCAR").read_text(encoding="utf-8")
        self.assertIn("NBANDSGW = 64", incar)
        self.assertIn("ENCUTGW = 280.5", incar)

    def test_rejects_invalid_target_states(self) -> None:
        cases = (
            ([], "at least one"),
            ([0, 1], "positive"),
            ([1, 1], "unique"),
            ([1, 11], "exceed NUM_WANN"),
        )
        for values, message in cases:
            with self.subTest(values=values):
                with self.assertRaisesRegex(
                    preparer.CrpaPreparationError, message
                ):
                    preparer.prepare_crpa(self.root, values)

    def test_expands_target_state_ranges_and_mixed_values(self) -> None:
        self.assertEqual(
            preparer.expand_target_state_tokens(("1-3", "5", "7-8")),
            (1, 2, 3, 5, 7, 8),
        )

    def test_rejects_malformed_or_descending_target_ranges(self) -> None:
        for token, message in (
            ("1:3", "invalid target-state"),
            ("3-1", "descending"),
            ("1-", "invalid target-state"),
        ):
            with self.subTest(token=token):
                with self.assertRaisesRegex(
                    preparer.CrpaPreparationError, message
                ):
                    preparer.expand_target_state_tokens((token,))

    def test_rejects_invalid_overrides(self) -> None:
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "cannot exceed NBANDS"
        ):
            preparer.prepare_crpa(self.root, [1], nbandsgw=81)
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "encutgw"
        ):
            preparer.prepare_crpa(self.root, [1], encutgw=float("nan"))
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "CRPA_KPAR"
        ):
            preparer.prepare_crpa(self.root, [1], kpar=0)

    def test_requires_completed_workflow_owned_wannier(self) -> None:
        (self.wannier / "OUTCAR").write_text(
            "NBANDS = 80\nunfinished\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "not complete"
        ):
            preparer.prepare_crpa(self.root, [1])

        (self.wannier / "OUTCAR").write_text(
            "NBANDS = 80\n"
            "General timing and accounting informations for this job:\n",
            encoding="utf-8",
        )
        (self.wannier / preparer.MARKER).unlink()
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "not a workflow-owned"
        ):
            preparer.prepare_crpa(self.root, [1])

    def test_requires_every_restart_input(self) -> None:
        (self.wannier / "WANPROJ").unlink()
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "WANPROJ"
        ):
            preparer.prepare_crpa(self.root, [1])

    def test_preserves_spin_polarized_wannier_mode(self) -> None:
        with (self.wannier / "INCAR").open("a", encoding="utf-8") as handle:
            handle.write("ISPIN = 2\n")
        (self.wannier / "WANPROJ").write_text(
            "# ISPIN NKPTS NB_TOT NW\n2 1 80 10\n", encoding="utf-8"
        )

        stage, settings = preparer.prepare_crpa(
            self.root,
            [1, 10],
            template=preparer.DEFAULT_INCAR_CRPA_TEMPLATE + "ISPIN = 2\n",
        )

        self.assertEqual(settings.ispin, 2)
        self.assertIn(
            "ISPIN = 2", (stage / "INCAR").read_text(encoding="utf-8")
        )
        self.assertEqual(settings.target_states, (1, 10))

    def test_rejects_wannier_incar_wanproj_spin_mismatch(self) -> None:
        with (self.wannier / "INCAR").open("a", encoding="utf-8") as handle:
            handle.write("ISPIN = 2\n")

        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "does not match WANPROJ ISPIN"
        ):
            preparer.prepare_crpa(self.root, [1])

    def test_rejects_malformed_or_inconsistent_wanproj_header(self) -> None:
        for header, message in (
            ("not a WANPROJ header\n", "header is incomplete"),
            ("# header\n1 1 bad 10\n", "non-integer"),
            ("# header\n1 1 80 9\n", "does not match WANPROJ NW"),
        ):
            with self.subTest(header=header):
                (self.wannier / "WANPROJ").write_text(header, encoding="utf-8")
                with self.assertRaisesRegex(preparer.CrpaPreparationError, message):
                    preparer.prepare_crpa(self.root, [1])

    def test_spin_custom_template_must_preserve_ispin(self) -> None:
        with (self.wannier / "INCAR").open("a", encoding="utf-8") as handle:
            handle.write("ISPIN = 2\n")
        (self.wannier / "WANPROJ").write_text(
            "# ISPIN NKPTS NB_TOT NW\n2 1 80 10\n", encoding="utf-8"
        )
        for template in (
            "ALGO = CRPA\nNTARGET_STATES = {{TARGET_STATES}}\n",
            "ISPIN = 1\nALGO = CRPA\nNTARGET_STATES = {{TARGET_STATES}}\n",
        ):
            with self.subTest(template=template), self.assertRaisesRegex(
                preparer.CrpaPreparationError, "must preserve ISPIN"
            ):
                preparer.prepare_crpa(self.root, [1], template=template)

    def test_refuses_unowned_stage(self) -> None:
        stage = self.root / "05_crpa"
        stage.mkdir()
        (stage / "keep").write_text("user data\n", encoding="utf-8")

        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "refusing to overwrite"
        ):
            preparer.prepare_crpa(self.root, [1], force=True)
        self.assertEqual(
            (stage / "keep").read_text(encoding="utf-8"), "user data\n"
        )

    def test_force_replaces_only_owned_stage(self) -> None:
        stage, _ = preparer.prepare_crpa(self.root, [1])
        (stage / "stale-output").write_text("old\n", encoding="utf-8")

        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "already exists"
        ):
            preparer.prepare_crpa(self.root, [2])

        replaced, settings = preparer.prepare_crpa(
            self.root, [2], force=True
        )
        self.assertEqual(settings.target_states, (2,))
        self.assertFalse((replaced / "stale-output").exists())
        self.assertIn(
            "NTARGET_STATES = 2",
            (replaced / "INCAR").read_text(encoding="utf-8"),
        )

    def test_rejects_unknown_template_tokens_before_writing_stage(self) -> None:
        with self.assertRaisesRegex(
            preparer.CrpaPreparationError, "UNKNOWN"
        ):
            preparer.prepare_crpa(
                self.root,
                [1],
                template="ALGO = CRPA\nVALUE = {{UNKNOWN}}\n",
            )
        self.assertFalse((self.root / "05_crpa").exists())

    def test_renders_an_overridden_template(self) -> None:
        template = """\
# site-specific cRPA
SYSTEM = {{SYSTEM}}
ALGO = CRPA
NBANDS = {{NBANDS}}
NTARGET_STATES = {{TARGET_STATES}}
KPAR = {{CRPA_KPAR}}
"""
        stage, _ = preparer.prepare_crpa(
            self.root, [4, 2], kpar=6, template=template
        )
        self.assertEqual(
            (stage / "INCAR").read_text(encoding="utf-8"),
            """# site-specific cRPA
SYSTEM = cRPA integration material
ALGO = CRPA
NBANDS = 80
NTARGET_STATES = 4 2
KPAR = 6
""",
        )


if __name__ == "__main__":
    unittest.main()

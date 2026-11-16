"""Integration coverage for the independent platform apps (no robot required)."""
import ast
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PlatformApps(unittest.TestCase):
    def test_generate_both_platforms(self):
        for folder, tag, forbidden in [
            ("DNABOT-OT2", "OT2", ("flex_", "robot_type", "DECK_LAYOUT")),
            (
                "DNABOT-Flex",
                "FLEX",
                (
                    "p20_single_gen2",
                    "p300_multi_gen2",
                    "magneticModuleV1",
                    "magneticModuleV2",
                    "robot_type",
                    "DECK_LAYOUT",
                ),
            ),
        ]:
            with self.subTest(platform=tag), tempfile.TemporaryDirectory() as output:
                app = ROOT / folder
                for path in app.glob("*.py"):
                    ast.parse(path.read_text(encoding="utf-8"))
                if tag == "OT2":
                    settings_text = (app / "default_settings_2_3_OT2.yaml").read_text()
                    self.assertIn("magneticModuleV1", settings_text)
                    self.assertIn("magneticModuleV2", settings_text)
                    self.assertNotIn("magnetic module gen1", settings_text)
                for path in (app / "template_opentrons_scripts").glob("*.py"):
                    text = path.read_text(encoding="utf-8")
                    ast.parse(text)
                    for token in forbidden:
                        self.assertNotIn(token, text, str(path))
                    if tag == "OT2" and "def get_low_volume_push_out" in text:
                        self.assertNotIn("return 7", text, str(path))
                        self.assertIn("return 2", text, str(path))
                    if tag == "FLEX" and "def get_low_volume_push_out" in text:
                        self.assertIn(
                            "return 7 if float(volume) < 5 else 2", text, str(path)
                        )
                    if path.name.startswith("1_"):
                        self.assertIn(
                            "aspirate_height = get_premix_aspirate_height(volume)", text
                        )
                        self.assertIn("dispense_height = 2", text)
                        self.assertNotIn("get_premix_dispense_height", text)
                        self.assertNotIn("final_dispense_height", text)
                        self.assertIn("Point(x=x_offset, y=0, z=0)", text)
                result = subprocess.run(
                    [
                        sys.executable,
                        str(app / f"dnabot_app_2_3_{tag}.py"),
                        "nogui",
                        "--construct_path",
                        str(app / "examples/constructs_temp.csv"),
                        "--source_paths",
                        str(app / "examples/BIOLEGIO_BASIC_RBS_EXT_SET.csv"),
                        str(app / "examples/part_plate_2_230419.csv"),
                        "--output_dir",
                        output,
                    ],
                    cwd=output,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                scripts = sorted(Path(output).glob("*.py"))
                self.assertEqual(len(scripts), 4)
                for script in scripts:
                    tree = ast.parse(script.read_text(encoding="utf-8"))
                    assignments = {
                        n.targets[0].id: n.value
                        for n in tree.body
                        if isinstance(n, ast.Assign)
                        and isinstance(n.targets[0], ast.Name)
                    }
                    self.assertEqual(
                        ast.literal_eval(assignments["metadata"])["apiLevel"], "2.21"
                    )
                    if tag == "OT2":
                        magnet_model = ast.literal_eval(assignments["__HARDWARE"])[
                            "mag_deck"
                        ]["id"]
                        self.assertIn(
                            magnet_model, ("magneticModuleV1", "magneticModuleV2")
                        )
                    if tag == "FLEX":
                        self.assertEqual(
                            ast.literal_eval(assignments["requirements"])["robotType"],
                            "Flex",
                        )
                    else:
                        self.assertNotIn("requirements", assignments)
                    if script.name.startswith("4_"):
                        groups = ast.literal_eval(assignments["spotting_tuples"])
                        self.assertEqual(groups[0][1][:4], ("A1", "B1", "C1", "A2"))
                        self.assertTrue(all(len(group[0]) <= 12 for group in groups))
                        self.assertTrue(
                            all(volume == 40 for group in groups for volume in group[2])
                        )
                self.assertTrue(
                    list(Path(output, "metainformation").glob("*clip_run_info.csv"))
                )
                deck_report = Path(output, "metainformation", "deck.md").read_text()
                self.assertEqual(deck_report.count("## "), 4)
                purification_report = deck_report.split("## 2_", 1)[1].split(
                    "## 3_", 1
                )[0]
                if tag == "OT2":
                    self.assertNotIn("thermocycler", purification_report)

    def test_12_well_mapping_repeats_for_each_plate(self):
        for folder, tag in [("DNABOT-OT2", "OT2"), ("DNABOT-Flex", "FLEX")]:
            spec = importlib.util.spec_from_file_location(
                "plates_" + tag, ROOT / folder / f"mplates_{tag}.py"
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            expected = [f"{row}{col}" for col in range(1, 5) for row in "ABC"]
            self.assertEqual(
                [module.final_12wellplate(i) for i in range(1, 97)], expected * 8
            )


if __name__ == "__main__":
    unittest.main()

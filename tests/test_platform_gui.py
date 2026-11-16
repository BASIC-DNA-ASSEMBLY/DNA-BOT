"""Headless GUI smoke tests; requires the apps' PySide6 dependency."""
import importlib.util
import os
from pathlib import Path
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
try:
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication
    import yaml
except ImportError:
    QApplication = None

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipIf(QApplication is None, "PySide6 is not installed")
class PlatformGui(unittest.TestCase):
    def test_load_and_save_settings(self):
        app = QApplication.instance() or QApplication([])
        for folder, tag in [("DNABOT-OT2", "OT2"), ("DNABOT-Flex", "FLEX")]:
            with self.subTest(platform=tag):
                root = ROOT / folder
                settings = yaml.safe_load(
                    (root / f"default_settings_2_3_{tag}.yaml").read_text()
                )
                settings.update(etoh_well="A11", soc_column=1)
                spec = importlib.util.spec_from_file_location(
                    "gui_" + tag, root / f"dnabot_gui_2_3_{tag}.py"
                )
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                errors = []
                alternate_source_plate = settings["labwares"]["clip_source_plate"][
                    "options"
                ][1]

                def save():
                    try:
                        window = next(
                            w
                            for w in app.topLevelWidgets()
                            if isinstance(w, module.GUI)
                        )
                        self.assertFalse(hasattr(window, "robot_type"))
                        self.assertEqual(window.single_pipette.count(), 1)
                        self.assertEqual(window.multi_pipette.count(), 1)
                        expected_magnets = (
                            ["Gen1", "Gen2"]
                            if tag == "OT2"
                            else settings["hardware"]["mag_deck"]["options"]
                        )
                        actual_magnets = [
                            window.magdeck.itemText(index)
                            for index in range(window.magdeck.count())
                        ]
                        self.assertEqual(actual_magnets, expected_magnets)
                        if tag == "OT2":
                            self.assertEqual(window.magdeck.currentText(), "Gen2")
                            window.magdeck.setCurrentText("Gen1")
                        else:
                            self.assertEqual(
                                window.magdeck.currentText(),
                                settings["hardware"]["mag_deck"]["id"],
                            )
                        for key, configuration in settings["labwares"].items():
                            expected = [
                                str(value) for value in configuration["options"]
                            ]
                            actual = [
                                window.labware_widgets[key].itemText(index)
                                for index in range(window.labware_widgets[key].count())
                            ]
                            self.assertEqual(actual, expected, key)
                        window.labware_widgets["clip_source_plate"].setCurrentText(
                            alternate_source_plate
                        )
                        window.construct_path = str(
                            root / "examples/constructs_temp.csv"
                        )
                        window.source_paths = [
                            str(root / "examples/part_plate_2_230419.csv")
                        ]
                        window.generate()
                    except BaseException as exc:
                        errors.append(exc)
                        app.quit()

                QTimer.singleShot(0, save)
                gui = module.GUI(None, settings)
                if errors:
                    raise errors[0]
                self.assertFalse(gui.quit_status)
                self.assertNotIn("robot_type", settings["hardware"])
                self.assertEqual(
                    settings["labwares"]["clip_source_plate"]["id"],
                    alternate_source_plate,
                )
                if tag == "OT2":
                    self.assertEqual(
                        settings["hardware"]["mag_deck"]["id"],
                        "magneticModuleV1",
                    )
                gui.deleteLater()
                app.processEvents()


if __name__ == "__main__":
    unittest.main()

# DNABOT Flex 2.3

Independent Flex app, GUI, settings, deck helpers, custom labware and four protocol templates. No files from the parent DNABOT installation are required. The existing shared app is retained separately for reference.

The Flex app has **its own virtual environment** in this folder (`DNABOT-Flex/.venv`), separate from the OT-2 environment at the repo root. Flex uses `opentrons==10.1.0`; the OT-2 app stays on `opentrons==9.0.0`, because 9.1 and newer reject OT-2 protocols. Keeping them apart lets the Flex side update Opentrons (and later Python) without affecting the OT-2.

Set it up once with Python 3.12, from inside `DNABOT-Flex`:

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe dnabot_app_2_3_FLEX.py
```

`.venv` is git-ignored. In VS Code, open the `DNABOT-Flex` folder itself (File → Open Folder) and choose **Python: Select Interpreter** → `.venv`; VS Code remembers this for the folder.

Generate the bundled example without opening the GUI:

```powershell
python dnabot_app_2_3_FLEX.py nogui --construct_path examples/constructs_temp.csv --source_paths examples/BIOLEGIO_BASIC_RBS_EXT_SET.csv examples/part_plate_2_230419.csv --output_dir examples/generated
```

The app writes four numbered protocols and a `metainformation` folder containing reagent calculations, source locations, assembly mappings, well information and `deck.md`. The deck report includes candidate tiprack locations; the protocol determines the occupied racks. `--help` and `nogui --help` describe settings overrides. The GUI uses ethanol well A11 and SOC column 1; CLI options can override these.

Hardware: Flex 1-channel 50 µL, Flex 8-channel 1000 µL, separate mounts, Gen2 thermocycler at B1, trash at A3, and magnetic block at D1 for purification. Clip source plates use D2, C2 and D1, with the tube rack at C1. Assembly takes place on the thermocycler. Copy the required definitions from `labware` into your Opentrons App custom labware library before analysing protocols.

All templates declare API 2.21 and `requirements = {"robotType": "Flex"}`, preserving the established API baseline. API upgrades can now be made independently in this app's four templates; keep metadata and filenames consistent. `requirements.txt` pins `opentrons==10.1.0` (Protocol API up to 2.31) for simulation and validation. Update it here when the Flex robots are updated.

Validation:

```powershell
python validate.py --generate-only
python validate.py
```

GUI load/save, CLI generation, Python syntax and generated variable checks passed. All four bundled examples passed `validate.py` simulation with Python 3.12 / Opentrons 10.1.0. Analyse all four generated protocols in the Flex App before laboratory use. Deck geometry, liquid handling and agar depth still need verification on your robot.

Changes carried through during the split:

- Removed robot selection and all alternate-platform hardware paths; standardised Flex instrument load names.
- Fixed broken clip-template indentation, missing helper arguments/imports, duplicate module loading and template sample globals that overwrote generated settings.
- Clip premixing aspirates progressively higher as the working volume increases, drawing upper liquid down through the tube. Every high-speed dispense remains 2 mm above the bottom with alternating lateral offsets to resuspend DNA condensate.
- Low-volume mode is explicitly switched back before larger mixing strokes.
- Purification uses 200 µL tips for sample handling and 1000 µL tips for washes. A magnetic block cannot disengage: the protocol prompts for removal of the empty mixing plate from C1 and manual movement of the magnetic plate to C1 for elution, then back to the block for separation. No gripper is required. Plate moves use `move_labware(..., use_gripper=False)` so the software tracks locations.
- Corrected repeating 12-well target coordinates and use of the 12-well spotting tuples. Each group of up to 12 constructs prompts for a fresh agar plate. Spotting uses 20 µL aliquots and the inherited 13 mm agar depth.
- The inherited 150 µL thermocycler outgrowth setting is rejected by API 2.21, whose block-volume argument is limited to 100 µL. Step 4 explicitly pauses for external incubation at 37 C for 60 minutes after adding 100 µL SOC, then asks for the plate to be returned. Heat shock remains on the thermocycler. Review this workflow change before laboratory use.
- The eight-channel SOC transfer requires medium in all eight wells of the selected deepwell column, and competent cells in all eight wells of each occupied transformation column.

Limits: 96 constructs, 48 unique clips, three source plates, and four assembly tipracks. Runs exceeding the assembly deck's tip capacity are rejected. Clip/transformation protocols prompt for tiprack refills when required.

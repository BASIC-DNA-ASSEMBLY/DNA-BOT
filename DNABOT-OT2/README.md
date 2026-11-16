# DNABOT OT2 2.3

Independent OT-2 app, GUI, settings, deck helpers, custom labware and four protocol templates. No files from the parent DNABOT installation are required. The existing shared app is retained separately for reference.

Use Python 3.10 and install the app dependencies:

```powershell
python -m pip install -r requirements.txt
python dnabot_app_2_3_OT2.py
```

Generate the bundled example without opening the GUI:

```powershell
python dnabot_app_2_3_OT2.py nogui --construct_path examples/constructs_temp.csv --source_paths examples/BIOLEGIO_BASIC_RBS_EXT_SET.csv examples/part_plate_2_230419.csv --output_dir examples/generated
```

The app writes four numbered protocols and a `metainformation` folder containing reagent calculations, source locations, assembly mappings, well information and `deck.md`. The deck report includes candidate tiprack locations; the protocol determines the occupied racks. `--help` and `nogui --help` describe settings overrides. The GUI uses ethanol well A11 and SOC column 1; CLI options can override these.

Hardware: P20 single Gen2, P300 multi Gen2, separate mounts, thermocycler Gen1/Gen2, and magnetic module Gen1/Gen2. Clip source plates use slots 2, 5 and 1; the tube rack uses slot 4. Assembly takes place on the thermocycler. Copy the required definitions from `labware` into your Opentrons App custom labware library before analysing protocols.

All templates declare API 2.21, preserving the established baseline. API upgrades can now be made independently in this app's four templates; keep metadata and filenames consistent. The generator does not depend on the Opentrons Python package. Install a compatible Opentrons simulator separately if needed.

Validation:

```powershell
python validate.py --generate-only
python validate.py
```

All four bundled examples passed simulation with Opentrons Python 8.2.0 / API 2.21. GUI load/save and CLI generation also passed automated checks. This is software validation; deck geometry, liquid handling and agar depth still need verification on your robot.

Changes carried through during the split:

- Removed robot selection and all alternate-platform hardware paths.
- Fixed broken clip-template indentation, missing helper arguments/imports, duplicate module loading and template sample globals that overwrote generated settings.
- Clip premixing aspirates progressively higher as the working volume increases, drawing upper liquid down through the tube. Every high-speed dispense remains 2 mm above the bottom with alternating lateral offsets to resuspend DNA condensate.
- Corrected repeating 12-well target coordinates and use of the 12-well spotting tuples. Each group of up to 12 constructs prompts for a fresh agar plate. Spotting uses 20 µL aliquots and the inherited 13 mm agar depth.
- The inherited 150 µL thermocycler outgrowth setting is rejected by API 2.21, whose block-volume argument is limited to 100 µL. Step 4 explicitly pauses for external incubation at 37 C for 60 minutes after adding 100 µL SOC, then asks for the plate to be returned. Heat shock remains on the thermocycler. Review this workflow change before laboratory use.
- The eight-channel SOC transfer requires medium in all eight wells of the selected deepwell column, and competent cells in all eight wells of each occupied transformation column.

Limits: 96 constructs, 48 unique clips, three source plates, and four assembly tipracks. Runs exceeding the assembly deck's tip capacity are rejected. Clip/transformation protocols prompt for tiprack refills when required.

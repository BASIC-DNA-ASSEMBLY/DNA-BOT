# DNA-BOT

**Feb 2022 - DNA-BOT has now been updated to produce scripts that run on the Opentrons OT2 in APIv2.** It also incorporates options for different hardware configurations by producing multiple scripts for each of the 4 steps. The Thermocycler module can now be optionally used in the step 1 clip reactions as well as the step 4 transformations.
In addition, different labware inputs are available through the Graphical User Interface (GUI) or Command Line Interface (CLI).

This work was originally published in [OUP Synthetic Biology](https://academic.oup.com/synbio/article/5/1/ysaa010/5869449)

## Getting Started

Users looking to implement the DNA-BOT workflow are encouraged to consult the [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf). If you are looking to contribute to this project, please raise an issue or pull request. Otherwise, feel free to reach out to [geoffbaldwin](mailto:g.baldwin@imperial.ac.uk).

dnabot can be run in 2 modes:

- with a graphical interface: see `Running the dnabot app` section in [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf). dnabot was developed using Python v3.7. Refer to [requirements.txt](requirements.txt).
- without a graphical interface: you need to specify all settings through the command line, you can see the instructions in the sections below.

## Installation from conda

```bash
conda create --name <myenv>
conda activate <myenv>
conda install -c conda-forge -c brsynth dnabot
```

`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.

## Usage

### Using the GUI

```bash
conda activate <myenv>
python -m dnabot.dnabot_app --help
python -m dnabot.dnabot_app
```

### Using the CLI only

```bash
conda activate <myenv>
python -m dnabot.dnabot_app nogui --help
python -m dnabot.dnabot_app nogui \
    --construct_path /path/to/constructs.csv \
    --source_paths /path/to/linker_parts_coord.csv /path/to/user_parts_coord.csv \
    --output_dir /path/to/output/dir
```

## Command line arguments

### GUI mode

```bash
usage: dnabot_app.py [-h] [--default_settings_file DEFAULT_SETTINGS_FILE] {nogui} ...

DNA assembly using BASIC on OpenTrons.

positional arguments:
  {nogui}               Optional, switch to define settings from the terminal instead of the graphical interface. Type "python dnabot_app.py
                        nogui -h" for more info.

optional arguments:
  -h, --help            show this help message and exit
  --default_settings_file DEFAULT_SETTINGS_FILE
                        Optional, file providing labware IDs and parameter to be used. Default:
                        /Users/tduigou/code/dnabot/dnabot/default_settings.yaml.
```

### No GUI mode

```bash
usage: dnabot_app.py nogui [-h] --construct_path CONSTRUCT_PATH --source_paths SOURCE_PATHS [SOURCE_PATHS ...]
                           [--etoh_well ETOH_WELL] [--soc_column SOC_COLUMN] [--output_dir OUTPUT_DIR]
                           [--template_dir TEMPLATE_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --construct_path CONSTRUCT_PATH
                        File listing constructs to be implemented.
  --source_paths SOURCE_PATHS [SOURCE_PATHS ...]
                        File(s) listing parts to be used in constructs.
  --etoh_well ETOH_WELL
                        Coordinates of the well plate providing ethanol for the purification step. Default: A11
  --soc_column SOC_COLUMN
                        Coordinate of the column plate providing SOC media for the transformation step. Default: 1
  --output_dir OUTPUT_DIR
                        Output directory. Default: same directory than the one containing the 'construct_path' file
  --template_dir TEMPLATE_DIR
                        Template directory. Default: 'template_ot2_scripts' located next to the present script.
```

## Change default values

Use the `--default_settings_file` argument to set different default values. This option is 
available either using the GUI or the CLI interface.

```bash
conda activate <myenv>
python -m dnabot.dnabot_app --default_settings_file /path/to/custom/default_settings.yaml.
```

The default settings file should follow the structure below (yaml file). The 
labware IDs to be used can be updated with the `labwares` section, while the 
parameters for the seperation step are listed in the `parameters` section.

```yaml
labwares:

  # Pipettes #############################################
  # Opentrons P20 Single-Channel Electronic Pipette
  p20_single:
    id: p20_single_gen2

  # Opentrons P300 8-Channel Electronic Pipette
  p300_multi:
    id: p300_multi_gen2

  # Modules ###############################################
  # Opentrons magnetic module (step: purification)
  mag_deck:
    id: magdeck
    # id: magnetic module gen2  # BRS

  # Tip racks #############################################
  # Opentrons 20μL tips rack
  96_tiprack_20ul:
    id: opentrons_96_tiprack_20ul
    # id: tipone_3dprinted_96_tiprack_20ul  # BRS

  # Opentrons 300μL tips rack
  96_tiprack_300ul:
    id: opentrons_96_tiprack_300ul
    # id: tipone_yellow_3dprinted_96_tiprack_300ul  # BRS

  # Plates ################################################
  # Opentrons 4-in-1 tubes rack for 1.5 ml eppendorf tubes (steps: clip, assembly, transformation)
  24_tuberack_1500ul:
    id: e14151500starlab_24_tuberack_1500ul
    # id: opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap  # BRS

  # Clip reaction source plate (steps: clip)
  clip_source_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: green_96_wellplate_200ul_pcr  # BRS (transparent)

  # Clip reaction plate (steps: clip, purif, assembly)
  clip_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: black_96_wellplate_200ul_pcr  # BRS (black)

  # Mix plate (step: purification)
  mix_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: black_96_wellplate_200ul_pcr  # BRS (black)
  
  # Final assembly plate (steps: assembly, transformation)
  final_assembly_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: black_96_wellplate_200ul_pcr  # BRS (black)
  
  # Transformation plate with thermocycler (step: transformation)
  transfo_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: black_96_wellplate_200ul_pcr  # BRS (black)

  # Transformation plate without thermocycler (step: transformation)
  transfo_plate_wo_thermo:
    id: 4ti0960rig_96_wellplate_200ul
    # id: green_96_wellplate_200ul_pcr  # BRS (transparent)

  # Agar plate (transformation step)
  agar_plate:
    id: 4ti0960rig_96_wellplate_200ul
    # id: thermoomnitrayfor96spots_96_wellplate_50ul  # BRS

  # Reservoir plate 21 mL 12 channels (step: purification)
  12_reservoir_21000ul:
    id: 4ti0131_12_reservoir_21000ul
    # id: citadel_12_wellplate_22000ul  # BRS

  # 96 deep well plate 2 mL wells (steps; purification, transormation)
  96_deepwellplate_2ml:
    id: 4ti0136_96_wellplate_2200ul
    # id: transparent_96_wellplate_2ml_deep  # BRS
  
  # Corning 12 Well Plate 6.9 mL Flat (step: transformation)
  12_corning_wellplate:
    id: corning_12_wellplate_6.9ml_flat
    # id: sarstedtcplatte_12_wellplate_6640ul  # BRS
  

parameters:

  # Clip reaction step ####################################
  # Keep the thermocycler lid closed at 4°C at the end of execution?
  # 1 for yes, 0 for no
  clip_keep_thermo_lid_closed:
    value: 0

  # Purification step #####################################
  # Magnetic module height (mm) - purification step
  purif_magdeck_height: 
    value: 20
    # value: 10.8  # BRS

  # Washing time (min) - purification step
  purif_wash_time:
    value: 0.5

  # Bead ratio - purification step
  purif_bead_ratio:
    value: 1.8

  # Incubation time (min) - purification step
  purif_incubation_time:
    value: 5

  # Settling time (min) - purification step
  purif_settling_time:
    value: 2
    # value: 6  # BRS

  # Drying time (min) - purification step
  purif_drying_time:
    value: 5
    # value: 15  # BRS

  # Elution time (min) - purif step
  purif_elution_time:
    value: 2
    # value: 5  # BRS
  
  # Transformation step ###################################
  # Incubation temperature
  transfo_incubation_temp:
    value: 4
    # value: 8  # BRS

  # Incubation time (min)
  transfo_incubation_time:
    value: 20
    # value: 30  # BRS

```

## For developers

### Development installation

After a git clone:

```bash
conda env create -f environment.yaml -n <dev_env>
conda develop -n <dev_env> .
```

You may be prompted to install `conda-build` in your base environment (`conda install conda-build`).
The default conda environment name will be `dnabot-dev` if not specified by `-n <dev-env>`.

Test your installation with:

```bash
conda activate <dev_env>
python -m dnabot.dnabot_app nogui --help
python -m pytest tests
```

To uninstall:

```bash
conda deactivate
conda env remove -n <dev_env>
```

### Tests

You need to install *pytest* if it's not done yet (`conda install pytest`).

```bash
conda install pytest
python -m pytest tests
```

## Authors

The update of DNABOT to APIv2 and improvements to the front end involved the work of several people:

- Thomas Duigou - [tduigou](https://github.com/tduigou)
- Geoff Baldwin - [geoffbaldwin](https://github.com/geoffbaldwin)
- Gizem Buldum  - [gizembuldum](https://github.com/gizembuldum)

Initial work revising the template scripts to run in Opentrons APIv2 was done by a team of Masters students from the 
MRes in Systems and Synthetic Biology at Imperial College London, thanks to:

- Xin Luo
- Ruihan Bai
- Zhenhua Wu
- Lianne Wu
- Ting An Lee
- Xiangming Xu

The original code for DNABOT produced OT2 scripts that ran in APIv1 and was authored by:

- **Matthew C Haines** - [hainesm6](https://github.com/hainesm6)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

- Marko Storch for all his help with [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf) and DNA-BOT manuscript.
- Geoff Baldwin for all his help with the DNA-BOT manuscript.

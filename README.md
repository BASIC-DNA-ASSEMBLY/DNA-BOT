# dnabot

Material for [DNA assembly using BASIC on OpenTrons (DNA-BOT)](https://www.biorxiv.org/content/10.1101/832139v1).

## Getting Started

Users looking to implement the DNA-BOT workflow are encouraged to consult the [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf). If you are looking to contribute to this project, please raise an issue or pull request. Otherwise, feel free to reach out to [hainesm6](mailto:hainesm6@gmail.com).

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

### Change default values

Use the `--labware_settings_file` argument to set different default values. This option is 
available either using the GUI or the CLI interface.

```bash
conda activate <myenv>
python -m dnabot.dnabot_app --labware_settings_file /path/to/custom/labware_settings.yaml.
```

The labware settings file should follow the structure below (yaml file)
where the "id“ values could be updated with the labware IDs that the user
is planning to use.

```yaml
labwares:

  # Opentrons P10 Single-Channel Electronic Pipette
  p10_single:
    id: p20_single_gen2

  # Opentrons P300 8-Channel Electronic Pipette
  p300_multi:
    id: p300_multi_gen2

  # Opentrons magnetic module
  mag_deck:
    id: magdeck

  # Opentrons 4-in-1 tubes rack for 1.5 ml eppendorf tubes
  24_tuberack_1500ul:
    id: e14151500starlab_24_tuberack_1500ul

  # Opentrons 10μL tips rack
  96_tiprack_10ul:
    id: opentrons_96_tiprack_20ul

  # Opentrons 300μL tips rack
  96_tiprack_300ul:
    id: opentrons_96_tiprack_300ul

  # 96 well rigid PCR plate (clip and transformation steps)
  96_wellplate_200ul_pcr_step_14:
    id: 4ti0960rig_96_wellplate_200ul

  # 96 well rigid PCR plate (purification step)
  96_wellplate_200ul_pcr_step_2:
    id: 4ti0960rig_96_wellplate_200ul

  # 96 well rigid PCR plate (assembly step)
  96_wellplate_200ul_pcr_step_3:
    id: 4ti0960rig_96_wellplate_200ul

  # Reservoir plate 21 mL 12 channels
  12_reservoir_21000ul:
    id: 4ti0131_12_reservoir_21000ul

  # 96 deep well plate 2 mL wells
  96_deepwellplate_2ml:
    id: 4ti0136_96_wellplate_2200ul
```


## Command line arguments

```
 optional arguments:
   -h, --help            show this help message and exit
   --construct_path CONSTRUCT_PATH
                         Construct CSV file.
   --source_paths SOURCE_PATHS [SOURCE_PATHS ...]
                         Source CSV files.
   --etoh_well ETOH_WELL
                         Well coordinate for Ethanol. Default: A11
   --soc_column SOC_COLUMN
                         Column coordinate for SOC. Default: 1
   --output_dir OUTPUT_DIR
                         Output directory. Default: same directory than the one
                         containing the "construct_path" file
   --template_dir TEMPLATE_DIR
                         Template directory. Default: "template_ot2_scripts"
                         located next to the present script.
```

## For developers

### Development installation

After a git clone:

```bash
cd <repository>
conda env create -f environment.yaml -n <dev_env>
conda develop -n <dev_env> .
```

You may be prompted to install `conda-build` in your base environment (`conda install conda-build`).
The default conda environment name will be `dnabot-dev` if not specified by `-n <dev-env>`.

Test your installation with:

```bash
conda activate <dev_env>
python -m dnabot.dnabot_app nogui --help
```

To uninstall:

```bash
conda deactivate
conda env remove -n <dev_env>
```

### Tests

You need to install *pytest* if it's not done yet (`conda install pytest`).

```bash
cd <repository>
conda install pytest
python -m pytest tests
```

## Authors

* **Matthew C Haines** - [hainesm6](https://github.com/hainesm6)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Marko Storch for all his help with [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf) and DNA-BOT manuscript.
* Geoff Baldwin for all his help with the DNA-BOT manuscript.
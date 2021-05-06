# dnabot

Material for [DNA assembly using BASIC on OpenTrons (DNA-BOT)](https://www.biorxiv.org/content/10.1101/832139v1).

## Getting Started

Users looking to implement the DNA-BOT workflow are encouraged to consult the [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf). If you are looking to contribute to this project, please raise an issue or pull request. Otherwise, feel free to reach out to [hainesm6](mailto:hainesm6@gmail.com).

dnabot can be run in 2 modes:
- with a graphical interface: see `Running the dnabot app` section in [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf). dnabot was developed using Python v3.7. Refer to [requirements.txt](requirements.txt).
- without a graphical interface: you need to specify all settings through the command line, you can see the instructions in the sections below.

## Installation

```bash
conda create --name <myenv>
conda activate <myenv>
conda install -c conda-forge dnabot
```

`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.

## Usage

After a git clone:

```bash
conda activate <myenv>
python dnabot/dnabot_app.py nogui --help
python dnabot/dnabot_app.py nogui \
    --construct_path path/to/constructs.csv \
    --source_paths /path/to/linker_parts_coord.csv /path/to/user_parts_coord.csv \
    --output_dir path/to/output/dir
```

## Command line arguments

```bash
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
conda activate <dev_env>
conda develop -n <dev_env> .
```

You may be prompted to install *conda-build* in your base environment (`conda install conda-build`).
The default conda environment name will be `dev_dnabot` if not specified by `-n <dev_env>`.

Test your installation with:

```bash
conda activate <dev_env>
python dnabot/dnabot_app.py nogui --help
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
python -m pytest
```

## Authors

* **Matthew C Haines** - [hainesm6](https://github.com/hainesm6)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Marko Storch for all his help with [instructions](docs/DNA_BOT_instructions_v1.0.0.pdf) and DNA-BOT manuscript.
* Geoff Baldwin for all his help with the DNA-BOT manuscript.
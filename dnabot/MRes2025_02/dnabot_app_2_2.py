str# -*- coding: utf-8 -*-
"""
Originally created on Thu Apr 11 14:26:07 2019 by @author: mh2210

@authors: geoffbaldwin, mh2210, gizembuldum, tduigou
"""
#line below started to throw a syntax error - no idea why
#from __future__ import annotations  # Enable the "hint" feature for objects

import os
import sys
import csv
import argparse
import pandas as pd
import numpy as np
import json
import tkinter as tk
import yaml
from pathlib import Path

#add dnabot module to syspath
abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, abs_path)

import mplates
import slots_2_1 as slots 

# Constant str
TEMPLATE_DIR_NAME = 'template_opentrons_scripts'

# CLIP_TEMP_FNAME_1 = '1_UNI_clip_template_APIv2_21.py'
# #CLIP_TEMP_FNAME_2 = '1_OT-2_clip_template_APIv2_21.py'
# #CLIP_TEMP_FNAME_4 = 'clip_template_Thermocycler_Gen2_APIv2_19.py'

# MAGBEAD_TEMP_FNAME_1 = '2_Flex_purification_template_APIv2_21.py'
# #MAGBEAD_TEMP_FNAME_2 = '2_OT-2_purification_template_APIv2_21.py'

# F_ASSEMBLY_TEMP_FNAME_1 = '3_Flex_assembly_template_APIv2_21.py'
# #F_ASSEMBLY_TEMP_FNAME_2 = '3_OT-2_assembly_template_APIv2_21.py'
# #F_ASSEMBLY_TEMP_FNAME_3 = 'assembly_template_Thermocycler_Gen1_APIv2.8.py'
# #F_ASSEMBLY_TEMP_FNAME_4 = 'assembly_template_Thermocycler_Gen2_APIv2.8.py'

# TRANSFORMATION_TEMP_FNAME_1 = '4_Flex_transformation_template_12wellplate_APIv2_21.py'
# TRANSFORMATION_TEMP_FNAME_2 = '4_OT-2_transformation_template_12wellplate_APIv2_21.py'
# #TRANSFORMATION_TEMP_FNAME_3 = 'transformation_template_Thermocycler_Gen1_APIv2.8.py'
# #TRANSFORMATION_TEMP_FNAME_4 = 'transformation_template_Thermocycler_Gen2_APIv2.8.py'

# CLIP_FNAME_1 = '1_Flex_clip_APIv2_21.py'
# CLIP_FNAME_2 = '1_OT-2_clip_APIv2_21.py'

# MAGBEAD_FNAME_1 = '2_Flex_purification_APIv2_21.py'
# MAGBEAD_FNAME_2 = '2_OT-2_purification_APIv2_21.py'

# F_ASSEMBLY_FNAME_1 = '3_Flex_assembly_APIv2_21.py'
# F_ASSEMBLY_FNAME_2 = '3_OT-2_assembly_APIv2_21.py'

# TRANSFORMATION_FNAME_1 = '4_Flex_transformation_12wellplate_APIv2_21.py'
# TRANSFORMATION_FNAME_2 = '4_OT-2_transformation_12wellplate_APIv2_21.py'
# #TRANSFORMATION_FNAME_3 = '4_transformation_ot2_Thermocycler_APIv2.8.py'
# #TRANSFORMATION_FNAME_4 = '4_transformation_ot2_Thermocycler_12wellplate_APIv2.8.py'

CLIPS_INFO_FNAME = 'clip_run_info.csv'
FINAL_ASSEMBLIES_INFO_FNAME = 'final_assembly_run_info.csv'
WELL_OUTPUT_FNAME = 'wells.txt'
DECK_OUTPUT_FNAME = "deck.md"

# Constant floats/ints
CLIP_DEAD_VOL = 60
CLIP_VOL = 30
T4_BUFF_VOL = 3
BSAI_VOL = 1
T4_LIG_VOL = 0.5
CLIP_MAST_WATER = 15.5
#ngDNA_PART_PER_CLIP = 200
MIN_VOL = 1
MAX_CONSTRUCTS = 96
MAX_CLIPS = 48
FINAL_ASSEMBLIES_PER_CLIP = 15
#DEFAULT_PART_VOL = 1
MAX_SOURCE_PLATES = 6
MAX_FINAL_ASSEMBLY_TIPRACKS = 7

# Constant dicts for 96 and 12 well plate formats
SPOTTING_VOLS_DICT = {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5}
SPOTTING_VOLS_DICT_12 = {2: 40, 3: 40, 4: 40, 5: 40, 6: 40, 7: 40}

# Constant lists


# Settings
DEFAULT_SETTINGS_FILE = Path(__file__).resolve().parent / 'default_settings_2_1_final.yaml'


def __cli():
    """Command line interface.

    :returns: CLI arguments
    :rtype: <argparse.Namespace>
    """
    desc = "DNA assembly using BASIC on OpenTrons."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--default_settings_file',
                        help='Optional, file providing Hardware IDs, labware IDs and parameter to be used. '
                             'Default: ' + str(DEFAULT_SETTINGS_FILE) +'.',
                        default= DEFAULT_SETTINGS_FILE)
    # Specific options for collecting settings from command line
    subparsers = parser.add_subparsers(help='Optional, switch to define settings from the terminal '
                                            'instead of the graphical interface. '
                                            'Type "python dnabot_app_2_2.py nogui -h" for more info.')
    parser_nogui = subparsers.add_parser('nogui')
    parser_nogui.add_argument('--robot_type',
                              help="Robot type to target when generating protocols.",
                              choices=['OT-2', 'Flex'],
                              required=True)
    parser_nogui.add_argument('--construct_path',
                              help='File listing constructs to be implemented.',
                              required=True)
    parser_nogui.add_argument('--source_paths',
                              help='File(s) listing parts to be used in constructs.',
                              nargs='+', required=True)
    parser_nogui.add_argument('--etoh_well',
                              help="Coordinates of the well plate providing ethanol "
                                   "for the purification step. Default: A11",
                              default='A11', type=str)
    parser_nogui.add_argument('--soc_column',
                              help="Coordinate of the column plate providing SOC media "
                                   "for the transformation step. Default: 1",
                              default=1, type=int)
    parser_nogui.add_argument('--output_dir',
                              help="Output directory. Default: same directory than the "
                                   "one containing the 'construct_path' file",
                              default=None, type=str or None)
    parser_nogui.add_argument('--template_dir',
                              help="Template directory. Default: 'template_opentrons_scripts' "
                                   "located next to the present script.",
                              default=None, type=str or None)
    parser_nogui.add_argument('--premix_linkers',
                              help="Override clip linker premixing. Choices: Yes or No. Defaults to the YAML value.",
                              choices=['Yes', 'No'],
                              default=None)
    parser_nogui.add_argument('--premix_parts',
                              help="Override clip part premixing. Choices: Yes or No. Defaults to the YAML value.",
                              choices=['Yes', 'No'],
                              default=None)
    parser_nogui.add_argument('--linkers_volume',
                              help="Override linker premix volume in uL. Defaults to the YAML value.",
                              default=None, type=float)
    parser_nogui.add_argument('--parts_volume',
                              help="Override part premix volume in uL. Defaults to the YAML value.",
                              default=None, type=float)
    parser_nogui.add_argument('--thermo_temp',
                              help="Override clip thermocycler hold temperature in C. Defaults to the YAML value.",
                              default=None, type=float)
    # Makes life easier to decide if we should switch to GUI or not
    parser.set_defaults(nogui=False)
    parser_nogui.set_defaults(nogui=True)
    return parser.parse_args()


def __get_settings_from_file(file_path: str) -> None:
    with open(file_path) as ifh:
        return yaml.safe_load(ifh)


def __info_from_gui(user_settings: dict) -> dict:
    """Pop GUI to collect user inputs

    Parameters
    ----------
    settings : dict
        default hardware, labware and parameter settings

    Returns
    -------
    dict
        actual settings to be used
    """
    user_settings = {
        'hardware': user_settings['hardware'],
        'labwares': user_settings['labwares'],
        'parameters': user_settings['parameters'],
        'construct_path': None,
        'sources_paths': None,
        'etoh_well': None,
        'soc_column': None,
        'robot_type': None
    }

    # Obtain user input
    print("Requesting user input, if not visible checked minimized windows.")
    # Collect info
    import dnabot_gui_2_2 as gui
    root = tk.Tk() 
    gui_inst = gui.GUI(root, user_settings)
    root.destroy()
    # User asked to quit?
    if gui_inst.quit_status:
        sys.exit("User specified 'QUIT' during app.")
    # Collect info
    user_settings = gui_inst.user_settings

    return user_settings


def main():
    
    # Parse args if any
    args = __cli()

    # Settings on labwares
    user_settings = __get_settings_from_file(args.default_settings_file)

    if args.nogui:
        robot_type = args.robot_type
        user_settings['hardware']['robot_type']['id'] = robot_type
        etoh_well = args.etoh_well
        soc_column = args.soc_column
        hardware_settings = user_settings['hardware']
        labware_settings = user_settings['labwares']
        parameter_settings = user_settings['parameters']
        if args.premix_linkers is not None:
            parameter_settings['premix_linkers']['value'] = args.premix_linkers
        if args.premix_parts is not None:
            parameter_settings['premix_parts']['value'] = args.premix_parts
        if args.linkers_volume is not None:
            parameter_settings['linkers_volume']['value'] = args.linkers_volume
        if args.parts_volume is not None:
            parameter_settings['parts_volume']['value'] = args.parts_volume
        if args.thermo_temp is not None:
            parameter_settings['thermo_temp']['value'] = args.thermo_temp
        construct_path = os.path.abspath(args.construct_path)
        sources_paths = [os.path.abspath(path) for path in args.source_paths]
        template_dir = os.path.abspath(args.template_dir) if args.template_dir else None
        if args.output_dir is not None:
            output_dir = args.output_dir
        else:
            output_dir = os.path.dirname(construct_path)
    else:
        user_inputs = __info_from_gui(user_settings)
        robot_type = user_inputs['robot_type']
        etoh_well = user_inputs['etoh_well']
        soc_column = user_inputs['soc_column']
        hardware_settings = user_inputs['hardware']  #changed from user_settings to user_inputs so that changed hardware is registered as default
        labware_settings = user_inputs['labwares']  # update labwares IDs
        parameter_settings = user_inputs['parameters']  # update parameters
        construct_path = user_inputs['construct_path']
        if construct_path is None:
            print('No file provided for constructs. Exit.')
            sys.exit()
        sources_paths = user_inputs['sources_paths']
        if sources_paths is None or len(sources_paths) == 0:
            print('No files prodived for sources. Exit.')
            sys.exit()
        output_dir = os.path.dirname(construct_path)
        template_dir = None

    # Select deck positions based on robot type
    # print("DEBUG: robot_type =", robot_type)
    # print("DEBUG: single_pipette =", hardware_settings['single_pipette']['id'])





    # if robot_type=='OT-2':
    #     if hardware_settings['single_pipette']['id'] == 'p20_single_gen2':
    #         pass
    #     else:
    #         raise ValueError("Invalid pipette for robot type. OT-2 requires 'p20_single_gen2'.")
        

    



    if robot_type=='OT-2':
        SOURCE_DECK_POS = ['2', '5', '1']
    elif robot_type=='Flex':
        SOURCE_DECK_POS = ["D2", "C2", "D1"]
    else:
        raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")

    # Check pipette compatibility with robot type
    if robot_type=='OT-2':
        if hardware_settings['single_pipette']['id']==('p20_single_gen2'):
            pass
        else:
            raise ValueError("Invalid pipette for robot type. OT-2 requires 'p20_single_gen2'.")

    if robot_type=='Flex':
        if hardware_settings['single_pipette']['id'] in ('Flex_1channel_50', 'flex_1channel_50'):
            pass
        else:
            raise ValueError("Invalid single pipette for robot type. Flex requires 'Flex_1channel_50'.")

    if robot_type=='OT-2':
        if hardware_settings['multi_pipette']['id']==('p300_multi_gen2'):
            pass
        else:
            raise ValueError("Invalid multi pipette for robot type. OT-2 requires 'p300_multi_gen2'.")
    if robot_type=='Flex':
        if hardware_settings['multi_pipette']['id'] in ('Flex_8channel_1000', 'flex_8channel_1000'):
            pass
        else:
            raise ValueError("Invalid multi pipette for robot type. Flex requires 'Flex_8channel_1000'.")
 # Check Mag plate compatibility with robot type
    if robot_type=='Flex':
        if hardware_settings['mag_deck']['id']==('magneticBlockV1'):
            pass
        else:
            raise ValueError("Invalid purification magnet for robot type. Flex requires 'magneticBlockV1'.")

    if robot_type=='OT-2':
        if hardware_settings['mag_deck']['id'] in ('magnetic module gen1', 'magnetic module gen2', 'magneticModuleV1', 'magneticModuleV2'):
            pass
        else:
            raise ValueError("Invalid purification magnet for robot type. OT-2 requires magnetic module Gen1 or Gen2.")

    # Args checking
    if len(sources_paths) > len(SOURCE_DECK_POS):
        raise ValueError('Number of source plates exceeds deck positions.')

    # Path to template directory
    if template_dir is not None:
        # Just to comment this case: only way to fall here is that the variable
        # has been set throught the command line arguments, nothing to do.
        template_dir_path = template_dir
        pass
    elif __name__ == '__main__':
        # Alternatively, try to automatically deduce the path relatively to the main script path
        script_path = os.path.abspath(__file__)
        template_dir_path = os.path.abspath(os.path.join(script_path, '..', TEMPLATE_DIR_NAME))
    else:
        # Fallback
        generator_dir = os.getcwd()
        template_dir_path = os.path.abspath(os.path.join(generator_dir, TEMPLATE_DIR_NAME))

    # Dealing with output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)

    # Prefix name
    construct_base = os.path.basename(construct_path)
    construct_base = os.path.splitext(construct_base)[0]
    print('User input successfully collected.')

    # Process input csv files
    print('Processing input csv files...')
    constructs_list = generate_constructs_list(construct_path)
    clips_df = generate_clips_df(constructs_list)
    sources_dict = generate_sources_dict(sources_paths, SOURCE_DECK_POS)

    # Calculate OT-2 script variables
    print('Calculating OT-2 variables...')
    clips_dict = generate_clips_dict(
        clips_df,
        sources_dict
        )
    magbead_sample_number = clips_df['number'].sum()
    final_assembly_dict = generate_final_assembly_dict(
        constructs_list,
        clips_df
        )
    final_assembly_tipracks = calculate_final_assembly_tipracks(
        final_assembly_dict
        )
    spotting_tuples = generate_spotting_tuples(
        constructs_list,
        SPOTTING_VOLS_DICT
        )
    spotting_tuples_12 = generate_spotting_tuples_12(
        constructs_list,
        SPOTTING_VOLS_DICT_12
        )
# Write OT-2 scripts
    print('Writing files...')

    if robot_type == 'Flex':
        generate_opentrons_script(
            '1_Flex_clip_APIv2_21.py',
            os.path.join(template_dir_path, '1_UNI_clip_template_APIv2_21_final.py'),            clips_dict=clips_dict,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '2_Flex_purification_APIv2_21.py',
            os.path.join(template_dir_path, '2_UNI_purification_template_APIv2_21_final.py'),            sample_number=magbead_sample_number,
            ethanol_well=etoh_well,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '3_Flex_assembly_APIv2_21.py',
            os.path.join(template_dir_path, '3_UNI_assembly_template_APIv2_21_final.py'),            final_assembly_dict=final_assembly_dict,
            tiprack_num=final_assembly_tipracks,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '4_Flex_transformation_12wellplate_APIv2_21.py',
            os.path.join(template_dir_path, '4_UNI_transformation_template_12wellplate_APIv2_21.py'),
            spotting_tuples=spotting_tuples,
            soc_well=f"A{soc_column}",
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
    
    elif robot_type == 'OT-2':
        generate_opentrons_script(
            '1_OT-2_clip_APIv2_21.py',
            os.path.join(template_dir_path, '1_UNI_clip_template_APIv2_21_final.py'),            clips_dict=clips_dict,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '2_OT-2_purification_APIv2_21.py',
            os.path.join(template_dir_path, '2_UNI_purification_template_APIv2_21_final.py'),            sample_number=magbead_sample_number,
            ethanol_well=etoh_well,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '3_OT-2_assembly_APIv2_21.py',
            os.path.join(template_dir_path, '3_UNI_assembly_template_APIv2_21_final.py'),            final_assembly_dict=final_assembly_dict,
            tiprack_num=final_assembly_tipracks,
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)
        generate_opentrons_script(
            '4_OT-2_transformation_12wellplate_APIv2_21.py',
            os.path.join(template_dir_path, '4_UNI_transformation_template_12wellplate_APIv2_21.py'),
            spotting_tuples=spotting_tuples,
            soc_well=f"A{soc_column}",
            __HARDWARE=hardware_settings,
            __LABWARES=labware_settings,
            __PARAMETERS=parameter_settings)   
    else:
        pass

    # # generate_opentrons_script(
    # #     CLIP_FNAME_1,
    # #     os.path.join(template_dir_path, CLIP_TEMP_FNAME_1),
    # #     clips_dict=clips_dict,
    # #     __HARDWARE=hardware_settings,
    # #     __LABWARES=labware_settings,
    # #     __PARAMETERS=parameter_settings)
    # # generate_opentrons_script(
    # #     CLIP_FNAME_2,
    # #     os.path.join(template_dir_path, CLIP_TEMP_FNAME_2),
    # #     clips_dict=clips_dict,
    # #     __HARDWARE=hardware_settings,
    # #     __LABWARES=labware_settings,
    # #     __PARAMETERS=parameter_settings)

       
    # generate_opentrons_script(
    #     MAGBEAD_FNAME_1,
    #     os.path.join(template_dir_path, MAGBEAD_TEMP_FNAME_1),
    #     sample_number=magbead_sample_number,
    #     ethanol_well=etoh_well,
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)
    
    # generate_opentrons_script(
    #     MAGBEAD_FNAME_2,
    #     os.path.join(template_dir_path, MAGBEAD_TEMP_FNAME_2),
    #     sample_number=magbead_sample_number,
    #     ethanol_well=etoh_well,
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)
    
    # generate_opentrons_script(
    #     F_ASSEMBLY_FNAME_1,
    #     os.path.join(template_dir_path, F_ASSEMBLY_TEMP_FNAME_1),
    #     final_assembly_dict=final_assembly_dict,
    #     tiprack_num=final_assembly_tipracks,
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)

    # generate_opentrons_script(
    #     F_ASSEMBLY_FNAME_2,
    #     os.path.join(template_dir_path, F_ASSEMBLY_TEMP_FNAME_2),
    #     final_assembly_dict=final_assembly_dict,
    #     tiprack_num=final_assembly_tipracks,
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)   
    
    # generate_opentrons_script(
    #     TRANSFORMATION_FNAME_1,
    #     os.path.join(template_dir_path, TRANSFORMATION_TEMP_FNAME_1),
    #     spotting_tuples=spotting_tuples,
    #     soc_well=f"A{soc_column}",
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)

    # generate_opentrons_script(
    #     TRANSFORMATION_FNAME_2,
    #     os.path.join(template_dir_path, TRANSFORMATION_TEMP_FNAME_2),
    #     spotting_tuples=spotting_tuples,
    #     soc_well=f"A{soc_column}",
    #     __HARDWARE=hardware_settings,
    #     __LABWARES=labware_settings,
    #     __PARAMETERS=parameter_settings)

    # Write information scripts
    metainfo_dir = Path().resolve() / "metainformation"
    metainfo_dir.mkdir(exist_ok=True)
    master_mix_df = generate_master_mix_df(clips_df['number'].sum())
    sources_paths_df = generate_sources_paths_df(sources_paths, SOURCE_DECK_POS)
    dfs_to_csv(
        metainfo_dir / f"{construct_base}_{CLIPS_INFO_FNAME}",
        index=False,
        MASTER_MIX=master_mix_df,
        SOURCE_PLATES=sources_paths_df,
        CLIP_REACTIONS=clips_df
        )
    with open(metainfo_dir / f"{construct_base}_{FINAL_ASSEMBLIES_INFO_FNAME}", "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for final_assembly_well, construct_clips in final_assembly_dict.items():
            csvwriter.writerow([final_assembly_well, construct_clips])
    with open(metainfo_dir / f"{construct_base}_{WELL_OUTPUT_FNAME}", "w") as f:
        f.write('Magbead ethanol well: {}'.format(etoh_well))
        f.write('\n')
        f.write('SOC column: {}'.format(soc_column))

    # Write deck position info
    # with open(metainfo_dir / f"{construct_base}_{DECK_OUTPUT_FNAME}", "w") as ofh:
    #     for fname in (CLIP_FNAME_1, CLIP_FNAME_2):
    #         deck = slots.get_positions_from_clip(fname)
    #         s = slots.format_deck_info(deck, section = f"Clip reaction script: {fname}")
    #         ofh.write(s)
    #     for fname in (MAGBEAD_FNAME_1, MAGBEAD_FNAME_2):
    #         deck = slots.get_positions_from_purif(fname)
    #         s = slots.format_deck_info(deck, section = f"Purification script: {fname}")
    #         ofh.write(s)
    #     for fname in (F_ASSEMBLY_FNAME_1, F_ASSEMBLY_FNAME_2):
    #         deck = slots.get_positions_from_assembly(fname)
    #         s = slots.format_deck_info(deck, section = f"Assembly script: {fname}")
    #         ofh.write(s)
    #     for fname in (TRANSFORMATION_FNAME_1, TRANSFORMATION_FNAME_2):
    #         deck = slots.get_positions_from_transfo(fname)
    #         s = slots.format_deck_info(deck, section = f"Transformation script: {fname}")
    #         ofh.write(s)
    print('BOT-2 generator successfully completed!')


def generate_constructs_list(path):
    """Generates a list of dataframes corresponding to each construct. Each
    dataframe lists components of the CLIP reactions required.

    """

    def process_construct(construct):
        """Processes an individual construct into a dataframe of CLIP reactions
        outlining prefix linkers, parts and suffix linkers.

        """

        def interogate_linker(linker):
            """Interrogates linker to determine if the suffix linker is a UTR
            linker.

            """
            if linker.startswith('U'):
                return linker.split('-')[0] + '-S'
            else:
                return linker + "-S"

        clips_info = {
            'prefixes': [],
            'parts': [],
            'suffixes': []
            }
        for i, sequence in enumerate(construct):
            if i % 2 != 0:
                clips_info['parts'].append(sequence)
                clips_info['prefixes'].append(
                    construct[i - 1] + '-P')
                if i == len(construct) - 1:
                    suffix_linker = interogate_linker(construct[0])
                    clips_info['suffixes'].append(suffix_linker)
                else:
                    suffix_linker = interogate_linker(construct[i + 1])
                    clips_info['suffixes'].append(suffix_linker)
        return pd.DataFrame.from_dict(clips_info)

    constructs_list = []
    with open(path, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        for index, construct in enumerate(csv_reader):
            if index != 0:  # Checks if row is header.
                construct = list(filter(None, construct))
                if not construct[1:]:
                    break
                else:
                    constructs_list.append(process_construct(construct[1:]))

    # Errors
    if len(constructs_list) > MAX_CONSTRUCTS:
        raise ValueError(
            'Number of constructs exceeds maximum. Reduce construct number in construct.csv.')
    else:
        return constructs_list


def generate_clips_df(constructs_list):
    """Generates a dataframe containing information about all the unique CLIP
    reactions required to synthesise the constructs in constructs_list.

    """
    merged_construct_dfs = pd.concat(constructs_list, ignore_index=True)
    unique_clips_df = merged_construct_dfs.drop_duplicates()
    unique_clips_df = unique_clips_df.reset_index(drop=True)
    clips_df = unique_clips_df.copy()

    # Error
    if len(unique_clips_df.index) > MAX_CLIPS:
        raise ValueError(
            'Number of CLIP reactions exceeds 48. Reduce number of constructs in construct.csv.')

    # Count number of each CLIP reaction
    clip_count = np.zeros(len(clips_df.index))
    for i, unique_clip in unique_clips_df.iterrows():
        for _, clip in merged_construct_dfs.iterrows():
            if unique_clip.equals(clip):
                clip_count[i] = clip_count[i] + 1
    clip_count = clip_count // FINAL_ASSEMBLIES_PER_CLIP + 1
    clips_df['number'] = [int(i) for i in clip_count.tolist()]

    # Associate well/s for each CLIP reaction
    clips_df['mag_well'] = pd.Series(['0'] * len(clips_df.index),
                                     index=clips_df.index)
    for index, number in clips_df['number'].items():
        if index == 0:
            mag_wells = []
            for x in range(number):
                mag_wells.append(mplates.final_well(x + 1 + 48))
            clips_df.at[index, 'mag_well'] = tuple(mag_wells)
        else:
            mag_wells = []
            for x in range(number):
                well_count = clips_df.loc[
                    :index - 1, 'number'].sum() + x + 1 + 48
                mag_wells.append(mplates.final_well(well_count))
            clips_df.at[index, 'mag_well'] = tuple(mag_wells)
    return clips_df


def generate_sources_dict(paths, SOURCE_DECK_POS):
    """Imports csvs files containing a series of parts/linkers with
    corresponding information into a dictionary where the key corresponds with
    part/linker and the value contains a tuple of corresponding information.
    Makes sure there are at least 5 columns (4 indexes) of information
    Deck position is added as the 5th column.

    Args:
        paths (list): list of strings each corresponding to a path for a
                      sources csv file.
    """

    sources_dict = {}
    for deck_index, path in enumerate(paths):
        with open(path, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            for index, source in enumerate(csv_reader):
                if index != 0:
                    csv_values = source[1:]  # Extract values excluding the first column
                    
                    # Ensure the list has at least 5 columns before inserting
                    while len(csv_values) < 4:
                        csv_values.append("")  # Fill missing columns with empty strings
                    
                    # Insert SOURCE_DECK_POS[deck_index] at index 4 (column 5)
                    csv_values.insert(4, SOURCE_DECK_POS[deck_index])
                    
                    sources_dict[str(source[0])] = tuple(csv_values)
    # print(sources_dict)
    return sources_dict
    
def generate_clips_dict(clips_df, sources_dict):
    """Using clips_df and sources_dict, returns a clips_dict which acts as the
    sole variable for the opentrons script "clip.ot2.py".

    """
    max_part_vol = CLIP_VOL - (T4_BUFF_VOL + BSAI_VOL + T4_LIG_VOL
                               + CLIP_MAST_WATER + 2)
    clips_dict = {'prefixes_wells': [], 'prefixes_plates': [],
                  'suffixes_wells': [], 'suffixes_plates': [],
                  'parts_wells': [], 'parts_plates': [], 'parts_vols': [],
                  'water_vols': []}

    # Generate clips_dict from args
    try:
        # Iterate through each row in the clips dataframe
        for _, clip_info in clips_df.iterrows():
            # Extract prefix linker and store corresponding well and plate information
            prefix_linker = clip_info['prefixes']
            clips_dict['prefixes_wells'].append([sources_dict[prefix_linker][0]] * clip_info['number'])
            clips_dict['prefixes_plates'].append([sources_dict[prefix_linker][4]] * clip_info['number'])
            
            # Extract suffix linker and store corresponding well and plate information
            suffix_linker = clip_info['suffixes']
            clips_dict['suffixes_wells'].append([sources_dict[suffix_linker][0]] * clip_info['number'])
            clips_dict['suffixes_plates'].append([sources_dict[suffix_linker][4]] * clip_info['number'])
            # Extract part and store corresponding well and plate information
            part = clip_info['parts']
            clips_dict['parts_wells'].append([sources_dict[part][0]] * clip_info['number'])
            clips_dict['parts_plates'].append([sources_dict[part][4]] * clip_info['number'])
            
            # Check if the concentration or part size is missing in sources_dict
            if not sources_dict[part][1] or not sources_dict[part][2]:
                sys.exit(f'Error: Missing values for part {part}. Both concentration and part_size are required.')
            else:
                # Calculate part volume using part size and concentration
                part_vol = round((float(sources_dict[part][2]) * 50) / float(sources_dict[part][1]), 2)
                
                # Ensure part_vol is within the allowed range
                if part_vol < MIN_VOL:
                    raise ValueError("Part concentration is too high, you need 50 ng per kb of plasmid and final volumes must be between 1 and 8 µL.")
                elif part_vol > max_part_vol:
                    raise ValueError("Part concentration is too low, you need 50 ng per kb of plasmid and final part volumes must be between 1 and 8 µL.")
                
                # Calculate the required water volume to reach the total volume
                water_vol = max_part_vol - part_vol +2
                
                # Store calculated part volume and water volume in clips_dict
                clips_dict['parts_vols'].append([part_vol] * clip_info['number'])
                clips_dict['water_vols'].append([water_vol] * clip_info['number'])
    except KeyError:
        # Handle missing part or linker in sources.csv
        sys.exit('likely part/linker not listed in sources.csv')

    # Flatten the lists in clips_dict to remove nested lists
    for key, value in clips_dict.items():
        clips_dict[key] = [item for sublist in value for item in sublist]

    # Return the final dictionary containing all clip information
    return clips_dict


def generate_final_assembly_dict(constructs_list, clips_df):
    """Using constructs_list and clips_df, returns a dictionary of final
    assemblies with keys defining destination plate well positions and values
    indicating which clip reaction wells are used.

    """
    final_assembly_dict = {}
    clips_count = np.zeros(len(clips_df.index))
    for construct_index, construct_df in enumerate(constructs_list):
        construct_well_list = []
        for _, clip in construct_df.iterrows():
            clip_info = clips_df[(clips_df['prefixes'] == clip['prefixes']) &
                                 (clips_df['parts'] == clip['parts']) &
                                 (clips_df['suffixes'] == clip['suffixes'])]
            clip_wells = clip_info.at[clip_info.index[0], 'mag_well']
            clip_num = int(clip_info.index[0])
            clip_well = clip_wells[int(clips_count[clip_num] //
                                       FINAL_ASSEMBLIES_PER_CLIP)]
            clips_count[clip_num] = clips_count[clip_num] + 1
            construct_well_list.append(clip_well)
        final_assembly_dict[mplates.final_well(
            construct_index + 1)] = construct_well_list
    return final_assembly_dict


def calculate_final_assembly_tipracks(final_assembly_dict):
    """Calculates the number of final assembly tipracks required ensuring
    no more than MAX_FINAL_ASSEMBLY_TIPRACKS are used.

    """
    final_assembly_lengths = []
    for values in final_assembly_dict.values():
        final_assembly_lengths.append(len(values))
    master_mix_tips = len(list(set(final_assembly_lengths)))
    total_tips = master_mix_tips + sum(final_assembly_lengths)
    final_assembly_tipracks = total_tips // 96 + (
        1 if total_tips % 96 > 0 else 0)
    if final_assembly_tipracks > MAX_FINAL_ASSEMBLY_TIPRACKS:
        raise ValueError(
            'Final assembly tiprack number exceeds number of slots. Reduce number of constructs in constructs.csv')
    else:
        return final_assembly_tipracks


def generate_spotting_tuples(constructs_list, spotting_vols_dict):
    """Using constructs_list, generates a spotting tuple
    (Refer to 'transformation_spotting_template.py') for every column of
    constructs, assuming the 1st construct is located in well A1 and wells
    increase linearly. Target wells locations are equivalent to construct well
    locations and spotting volumes are defined by spotting_vols_dict.

    Args:
        spotting_vols_dict (dict): Part number defined by keys, spottting
            volumes defined by corresponding value.

    """
    # Calculate wells and volumes
    wells = [mplates.final_well(x + 1) for x in range(len(constructs_list))]
    vols = [SPOTTING_VOLS_DICT[len(construct_df.index)]
            for construct_df in constructs_list]

    # Package spotting tuples
    spotting_tuple_num = len(constructs_list)//8 + (1
                                                    if len(constructs_list) % 8 > 0 else 0)
    spotting_tuples = []
    for x in range(spotting_tuple_num):
        if x == spotting_tuple_num - 1:
            tuple_wells = tuple(wells[8*x:])
            tuple_vols = tuple(vols[8*x:])
        else:
            tuple_wells = tuple(wells[8*x:8*x + 8])
            tuple_vols = tuple(vols[8*x:8*x + 8])
        spotting_tuples.append((tuple_wells, tuple_wells, tuple_vols))
    return spotting_tuples

# Introduced 12 well plate format for spotting
def generate_spotting_tuples_12(constructs_list, spotting_vols_dict):
    """Using constructs_list, generates a spotting tuple
    (Refer to 'transformation_spotting_template.py') for every column of
    constructs, assuming the 1st construct is located in well A1 and wells
    increase linearly. Target wells locations are equivalent to construct well
    locations and spotting volumes are defined by spotting_vols_dict.

    Args:
        spotting_vols_dict (dict): Part number defined by keys, spottting
            volumes defined by corresponding value.

    """
    # Calculate wells and volumes
    wells = [mplates.final_well(x + 1) for x in range(len(constructs_list))]
    vols = [SPOTTING_VOLS_DICT_12[len(construct_df.index)]
            for construct_df in constructs_list]

    spot_wells = [mplates.final_12wellplate(x + 1) for x in range(len(constructs_list))]

    # Package spotting tuples
    spotting_tuple_num = len(constructs_list)//12 + (1
                                                    if len(constructs_list) % 12 > 0 else 0)
    spotting_tuples_12 = []

    for x in range(spotting_tuple_num): 
        if x == spotting_tuple_num - 1:
            tuple_wells = tuple(wells[12*x:])
            tuple_vols = tuple(vols[12*x:])
            tuple_spot_wells = tuple(spot_wells[12*x:])
        else:
            tuple_wells = tuple(wells[12*x:12*x + 12])
            tuple_vols = tuple(vols[12*x:12*x + 12])
            tuple_spot_wells = tuple(spot_wells[12*x:12*x + 12])
        spotting_tuples_12.append((tuple_wells, tuple_spot_wells, tuple_vols))
    return spotting_tuples_12

#OLD write code
# def generate_opentrons_script(opentrons_script_path, template_path, **kwargs):
#     """Generates an ot2 script named 'opentrons_script_path', where kwargs are
#     written as global variables at the top of the script. For each kwarg, the
#     keyword defines the variable name while the value defines the name of the
#     variable. The remainder of template file is subsequently written below.

#     """
#     with open(opentrons_script_path, 'w') as wf:
#         with open(template_path, 'r') as rf:
#             for index, line in enumerate(rf):
#                 if line[:3] == 'def':
#                     function_start = index
#                     break
#                 else:
#                     wf.write(line)
#             for key, value in kwargs.items():
#                 wf.write('{}='.format(key))
#                 if type(value) == dict:
#                     wf.write(json.dumps(value))
#                 elif type(value) == str:
#                     wf.write("'{}'".format(value))
#                 else:
#                     wf.write(str(value))
#                 wf.write('\n')
#             wf.write('\n')
#         with open(template_path, 'r') as rf:
#             for index, line in enumerate(rf):
#                 if index >= function_start - 1:
#                     wf.write(line)

# NEW code to write before requirements={} so robot type in requirements can be specified from __hardware definition
def generate_opentrons_script(opentrons_script_path, template_path, **kwargs):
    """Generates an ot2 script named 'opentrons_script_path', where kwargs are
    written as global variables at the top of the script. For each kwarg, the
    keyword defines the variable name while the value defines the name of the
    variable. The remainder of template file is subsequently written below.

    """
    with open(opentrons_script_path, 'w') as wf:
        with open(template_path, 'r') as rf:
            lines = rf.readlines()

        hardware = kwargs.get('__HARDWARE', {})
        robot_type = hardware.get('robot_type', {}).get('id')

        # Inject generated globals after the import block so module imports are
        # resolved before the template-specific dictionaries appear.
        import_end_index = 0
        for index, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('import ') or stripped.startswith('from '):
                import_end_index = index + 1

        # Find the first line that starts with 'requirements'
        requirements_index = None
        for index, line in enumerate(lines):
            if line.strip().startswith('requirements'):
                requirements_index = index
                break

        for line in lines[:import_end_index]:
            wf.write(line)

        wf.write('\n')
        if robot_type == 'Flex':
            wf.write('requirements = {"robotType": "Flex"}\n\n')

        for key, value in kwargs.items():
            wf.write('{}='.format(key))
            if isinstance(value, dict):
                wf.write(json.dumps(value))
            elif isinstance(value, str):
                wf.write("'{}'".format(value))
            else:
                wf.write(str(value))
            wf.write('\n')
        wf.write('\n')

        if requirements_index is None:
            for line in lines[import_end_index:]:
                wf.write(line)
            return

        skip_requirements_block = False
        for index, line in enumerate(lines[import_end_index:], start=import_end_index):
            stripped = line.strip()
            if index == requirements_index:
                skip_requirements_block = True
                continue
            if skip_requirements_block:
                if stripped == '':
                    skip_requirements_block = False
                continue
            wf.write(line)

def generate_master_mix_df(clip_number):
    """Generates a dataframe detailing the components required in the clip
    reaction master mix.

    """
    COMPONENTS = {'Component': ['Promega T4 DNA Ligase buffer, 10X',
                                'Water', 'NEB BsaI-HFv2',
                                'Promega T4 DNA Ligase']}
    VOL_COLUMN = 'Volume (uL)'
    master_mix_df = pd.DataFrame.from_dict(COMPONENTS)
    master_mix_df[VOL_COLUMN] = (clip_number + CLIP_DEAD_VOL/CLIP_VOL) * \
        np.array([T4_BUFF_VOL,
                  CLIP_MAST_WATER,
                  BSAI_VOL,
                  T4_LIG_VOL])
    return master_mix_df


def generate_sources_paths_df(paths, SOURCE_DECK_POS):
    """Generates a dataframe detailing source plate information.

    Args:
        paths (list): list of strings specifying paths to source plates.
        SOURCE_DECK_POS (list): list of strings specifying candidate deck positions.

    """
    source_plates_dict = {'Deck position': [], 'Source plate': [], 'Path': []}
    for index, path in enumerate(paths):
        source_plates_dict['Deck position'].append(SOURCE_DECK_POS[index])
        source_plates_dict['Source plate'].append(os.path.basename(path))
        source_plates_dict['Path'].append(path)
    return pd.DataFrame(source_plates_dict)


def dfs_to_csv(path, index=True, **kw_dfs):
    """Generates a csv file defined by path, where kw_dfs are
    written one after another with each key acting as a title. If index=True,
    df indexes are written to the csv file.

    """
    with open(path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for key, value in kw_dfs.items():
            csvwriter.writerow([str(key)])
            value.to_csv(csvfile, index=index)
            csvwriter.writerow('')


if __name__ == '__main__':
    main()

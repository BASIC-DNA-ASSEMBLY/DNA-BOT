# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:26:07 2019

@author: mh2210

TO DO
    - add in new transformation protocol
    - ammend meta information
    - add robot comments 
    - add different file locations for different protocol versions
    - add clip input functionality
    - ammend instruction manuel
    - new GUI??
 
"""
import os
import sys

#add dnabot module to syspath
abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, abs_path)

import csv
import argparse
import pandas as pd
import numpy as np
import json
import tkinter as tk
import dnabot_gui as gui
import regex as re

# Constant str
TEMPLATE_DIR_NAME = 'template_ot2_scripts'
CLIP_TEMP_FNAME_1 = 'clip_template_APIv1.py'
CLIP_TEMP_FNAME_2 = 'clip_template_APIv2.8.py'
CLIP_TEMP_FNAME_3 = 'clip_template_Thermocycler_module_APIv2.8.py'

MAGBEAD_TEMP_FNAME_1 = 'purification_template_APIv1.py'
MAGBEAD_TEMP_FNAME_2 = 'purification_template_APIv2.8.py'

F_ASSEMBLY_TEMP_FNAME_1 = 'assembly_template_APIv1.py'
F_ASSEMBLY_TEMP_FNAME_2 = 'assembly_template_APIv2.8.py'
F_ASSEMBLY_TEMP_FNAME_3 = 'assembly_template_Thermocycler_module_APIv2.8.py'

TRANS_SPOT_TEMP_FNAME_1 = 'transformation_template_APIv1.py'
TRANS_SPOT_TEMP_FNAME_2 = 'transformation_template_APIv2.8.py'
TRANS_SPOT_TEMP_FNAME_3 = 'transformation_template_Thermocycler_module_APIv2.8.py'

CLIP_FNAME_1 = 'A_clip_ot2_APIv1'                           # removed '.py'
CLIP_FNAME_2 = 'A_clip_ot2_APIv2.8'                         # removed '.py'
CLIP_FNAME_3 = 'A_clip_ot2_Thermocycler_APIv2.8'            # removed '.py'

MAGBEAD_FNAME_1 = 'B_purification_ot2_APIv1'                # removed '.py'
MAGBEAD_FNAME_2 = 'B_purification_ot2_APIv2.8'              # removed '.py'

F_ASSEMBLY_FNAME_1 = 'C_assembly_ot2_APIv1'                 # removed '.py'
F_ASSEMBLY_FNAME_2 = 'C_assembly_ot2_APIv2.8'               # removed '.py'
F_ASSEMBLY_FNAME_3 = 'C_assembly_ot2_Thermocycler_APIv2.8'  # removed '.py'

TRANS_SPOT_FNAME_1 = 'D_transformation_ot2_APIv1.py'
TRANS_SPOT_FNAME_2 = 'D_transformation_ot2_APIv2.8.py'
TRANS_SPOT_FNAME_3 = 'D_transformation_ot2_Thermocycler_APIv2.8.py'

CLIPS_INFO_FNAME = 'clip_run_info.csv'
FINAL_ASSEMBLIES_INFO_FNAME = 'final_assembly_run_info.csv'
WELL_OUTPUT_FNAME = 'wells.txt'

V1_PATH = "\\APIv1\\"
V2_8_PATH = "\\APIv2.8\\"
V2_8_TC_PATH = "\\Thermocycler_APIv2.8\\"

# Constant floats/ints
CLIP_DEAD_VOL = 60
CLIP_VOL = 30
T4_BUFF_VOL = 3
BSAI_VOL = 1
T4_LIG_VOL = 0.5
CLIP_MAST_WATER = 15.5
PART_PER_CLIP = 200
MIN_VOL = 1
# MAX_CONSTRUCTS = 96               # not required
MAX_CLIPS_PER_PLATE = 48            # Max clips per clip plate
MAX_CLIPS_TOTAL = 96*2              # 48
MAX_ASSEMBLIES_PER_PLATE = 96
FINAL_ASSEMBLIES_PER_CLIP = 15
DEFAULT_PART_VOL = 1
MAX_SOURCE_PLATES = 6
MAX_FINAL_ASSEMBLY_TIPRACKS = 7

# Constant dicts
SPOTTING_VOLS_DICT = {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5}

# Constant lists
SOURCE_DECK_POS = ['2', '5', '8', '7', '10', '11']      # NB for thermocycler protocols, the thermocycler takes up slots 7, 8, 10, 11


def __cli():
    """Command line interface.

    :returns: CLI arguments
    :rtype: <argparse.Namespace>
    """
    desc = "DNA assembly using BASIC on OpenTrons."
    parser = argparse.ArgumentParser(description=desc)

    # Specific options for collecting settings from command line
    subparsers = parser.add_subparsers(help='Optional, to define settings from the terminal instead of the graphical '
                                            'interface. Type "python dnabot_app.py nogui -h" for more info.')
    parser_nogui = subparsers.add_parser('nogui')
    parser_nogui.add_argument('--construct_path', help='Construct CSV file.', required=True)
    parser_nogui.add_argument('--source_paths', help='Source CSV files.', nargs='+', required=True)
    parser_nogui.add_argument('--etoh_well', help='Well coordinate for Ethanol. Default: A11', default='A11', type=str)
    parser_nogui.add_argument('--soc_column', help='Column coordinate for SOC. Default: 1', default=1, type=int)
    parser_nogui.add_argument('--output_dir',
                              help='Output directory. Default: same directory than the one containing the '
                                   '"construct_path" file',
                              default=None, type=str or None)
    parser_nogui.add_argument('--template_dir',
                              help='Template directory. Default: "template_ot2_scripts" located next to the present '
                                   'script.',
                              default=None, type=str or None)
    
    parser.set_defaults(nogui=False)
    parser_nogui.set_defaults(nogui=True)
    return parser.parse_args()


def __info_from_gui():
    """Pop GUI to collect user inputs.

    :returns user_inputs: info collected
    :rtype: dict
    """
    user_inputs = {
        'construct_path': None,
        'sources_paths': None,
        'etoh_well': None,
        'soc_column': None
    }

    # Obtain user input
    print("Requesting user input, if not visible checked minimized windows.")
    root = tk.Tk()
    dnabotinst = gui.DnabotApp(root)
    root.mainloop()
    root.destroy()
    if dnabotinst.quit_status:
        sys.exit("User specified 'QUIT' during app.")
    # etoh_well and soc_column are silently collected by the gui
    user_inputs['etoh_well'] = dnabotinst.etoh_well
    user_inputs['soc_column'] = dnabotinst.soc_column
    # construct file path
    root = tk.Tk()
    user_inputs['construct_path'] = gui.UserDefinedPaths(root, 'Construct csv file').output
    root.destroy()

    # part & linker file paths
    root = tk.Tk()
    user_inputs['sources_paths'] = gui.UserDefinedPaths(root, 'Sources csv files', multiple_files=True).output
    root.destroy()

    return user_inputs


def main():
    # Settings
    # args = __cli()
    
    # if args.nogui:
    #     etoh_well = args.etoh_well
    #     soc_column = args.soc_column
    #     construct_path = args.construct_path
    #     sources_paths = args.source_paths
    #     output_dir = args.output_dir
    #     template_dir = args.template_dir
    # else:
    #     user_inputs = __info_from_gui()
    #     etoh_well = user_inputs['etoh_well']
    #     soc_column = user_inputs['soc_column']
    #     construct_path = user_inputs['construct_path']
    #     sources_paths = user_inputs['sources_paths']
    #     output_dir = os.path.dirname(construct_path)
    #     template_dir = None

    #### TEST FILES - TO BE DELETED ####
    etoh_well = 'A11'
    soc_column = 1
    construct_path = 'C:\\Users\\ljh119\\OneDrive - Imperial College London\\648_build\\DNA-BOT\\648_constructs\\multistage_builds\\stage_2\\stage2_constructs.csv'
    sources_paths = ['C:\\Users\\ljh119\\OneDrive - Imperial College London\\648_build\\DNA-BOT\\648_constructs\\multistage_builds\\stage_2\\stage2_parts.csv']
    output_dir = 'C:\\Users\\ljh119\\OneDrive - Imperial College London\\648_build\\DNA-BOT\\648_constructs\\multistage_builds\\stage_2'
    template_dir = None

    ####################################

    # Args checking
    if len(sources_paths) > len(SOURCE_DECK_POS):
        raise ValueError('Number of source plates exceeds deck positions.')

    # Path to template directory
    if template_dir is not None:
        # Just to comment this case: only way to fall here is that the variable has been set throught the command
        # line arguments, nothing to do.
        template_dir_path = template_dir
        pass
    elif __name__ == '__main__':
        # Alternatively, try to automatically deduce the path relatively to the main script path
        script_path = os.path.abspath(__file__)
        template_dir_path = os.path.abspath(os.path.join(script_path, '..', TEMPLATE_DIR_NAME))
    else:
        # Fallback - if template directory is neither given in cli nor the file is run directly, 
        # check for the template dir in the current wd
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
    constructs_list = generate_constructs_list(construct_path)      # returns a list of individual dfs for each construct in which each row is a required clip reaction
    clips_df = generate_clips_df(constructs_list)                   # takes the constructs_list and returns a df of all unique clip reactions with the number of times each one is required to be made and the location(s) of it
    sources_dict = generate_sources_dict(sources_paths)             # returns a dict with source id as keys and location as values

    # output clips csv
    output_file = os.path.join(output_dir, 'clips_df.csv')
    clips_df.to_csv(output_file, index=False)

    # calculate OT2 script variables
    print('Calculating OT-2 variables...')
    clips_dict_list = generate_clips_dict_list(clips_df, sources_dict)      # takes the clips df and the sources dict and produces a list of dictionaries of tuples of the parts, prefixes, and suffixes' wells and plate locations
    
    magbead_sample_number_total = clips_df['number'].sum()                  # the total number of clips is built in sets of >48, these are combined in pairs into up to 2 plates each of >96 clips, the total number of clips possible is 192 (see error trapping in clip_dict_list_generator)
    magbead_sample_list = [96 for i in range(magbead_sample_number_total//96)] + [magbead_sample_number_total % 96 + 1] # creates a list of values where each value represents the number of mag purifications for a single step


    clips_dict = generate_clips_dict(clips_df, sources_dict)                # one dict of everything - original version - to be used by final assembly generator functions
    final_assembly_dict_list = generate_final_assembly_dict_list(constructs_list, clips_df)
    
    ################### UNCOMMENT ###################
    # spotting_tuples = generate_spotting_tuples(constructs_list,
    #                                            SPOTTING_VOLS_DICT)
    #################################################

    print('Writing files...')
    # Write OT2 scripts
    for clip_plate in range(len(clips_dict_list)):
        sub_clip_dict = clips_dict_list[clip_plate]
                
        generate_ot2_script(CLIP_FNAME_1 + '_SCRIPT_{}.py'.format(clip_plate+1), os.path.join(
            template_dir_path, CLIP_TEMP_FNAME_1), clips_dict=sub_clip_dict)
        generate_ot2_script(CLIP_FNAME_2 + '_SCRIPT_{}.py'.format(clip_plate+1), os.path.join(
            template_dir_path, CLIP_TEMP_FNAME_2), clips_dict=sub_clip_dict)
        generate_ot2_script(CLIP_FNAME_3 + '_SCRIPT_{}.py'.format(clip_plate+1), os.path.join(
            template_dir_path, CLIP_TEMP_FNAME_3), clips_dict=sub_clip_dict)
    
        
    ######## NEED TO AMMEND META DATA ##########
        
    for i, magbead_sample_number in enumerate(magbead_sample_list):
        
        generate_ot2_script(MAGBEAD_FNAME_1 + '_SCRIPT_{}.py'.format(i+1), os.path.join(
            template_dir_path, MAGBEAD_TEMP_FNAME_1),
            sample_number=magbead_sample_number,
            ethanol_well=etoh_well)
        generate_ot2_script(MAGBEAD_FNAME_2 + '_SCRIPT_{}.py'.format(i+1), os.path.join(
            template_dir_path, MAGBEAD_TEMP_FNAME_2),
            sample_number=magbead_sample_number,
            ethanol_well=etoh_well)
    
    ############### If I want to add in the clips as a new source point I need to append the new clips plate at this point ##########################

    for assembly_plate in range(len(final_assembly_dict_list)):
        final_assembly_dict = final_assembly_dict_list[assembly_plate]
        final_assembly_tipracks = calculate_final_assembly_tipracks(final_assembly_dict)

        generate_ot2_script(F_ASSEMBLY_FNAME_1 + '_SCRIPT_{}.py'.format(assembly_plate+1), os.path.join(
            template_dir_path, F_ASSEMBLY_TEMP_FNAME_1),
            final_assembly_dict=final_assembly_dict,
            tiprack_num=final_assembly_tipracks)
        generate_ot2_script(F_ASSEMBLY_FNAME_2 + '_SCRIPT_{}.py'.format(assembly_plate+1), os.path.join(
            template_dir_path, F_ASSEMBLY_TEMP_FNAME_2),
            final_assembly_dict=final_assembly_dict,
            tiprack_num=final_assembly_tipracks)
        generate_ot2_script(F_ASSEMBLY_FNAME_3 + '_SCRIPT_{}.py'.format(assembly_plate+1), os.path.join(
            template_dir_path, F_ASSEMBLY_TEMP_FNAME_3),
            final_assembly_dict=final_assembly_dict,
            tiprack_num=final_assembly_tipracks)
    
    # generate_ot2_script(TRANS_SPOT_FNAME_1, os.path.join(
    #     template_dir_path, TRANS_SPOT_TEMP_FNAME_1),
    #     spotting_tuples=spotting_tuples,
    #     soc_well=f"A{soc_column}")
    # generate_ot2_script(TRANS_SPOT_FNAME_2, os.path.join(
    #     template_dir_path, TRANS_SPOT_TEMP_FNAME_2),
    #     spotting_tuples=spotting_tuples,
    #     soc_well=f"A{soc_column}")
    # generate_ot2_script(TRANS_SPOT_FNAME_3, os.path.join(
    #     template_dir_path, TRANS_SPOT_TEMP_FNAME_3),
    #     spotting_tuples=spotting_tuples,
    #     soc_well=f"A{soc_column}")

    # Write non-OT2 scripts
    if 'metainformation' in os.listdir():
        pass
    else:
        os.makedirs('metainformation')
    os.chdir('metainformation')
    master_mix_df = generate_master_mix_df(clips_df['number'].sum())
    sources_paths_df = generate_sources_paths_df(sources_paths, SOURCE_DECK_POS)
    dfs_to_csv(construct_base + '_' + CLIPS_INFO_FNAME, index=False,
               MASTER_MIX=master_mix_df, SOURCE_PLATES=sources_paths_df,
               CLIP_REACTIONS=clips_df)
    with open(construct_base + '_' + FINAL_ASSEMBLIES_INFO_FNAME,
              'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for final_assembly_well, construct_clips in final_assembly_dict.items():
            csvwriter.writerow([final_assembly_well, construct_clips])
    with open(construct_base + '_' + WELL_OUTPUT_FNAME, 'w') as f:
        f.write('Magbead ethanol well: {}'.format(etoh_well))
        f.write('\n')
        f.write('SOC column: {}'.format(soc_column))
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

        clips_info = {'prefixes': [], 'parts': [],
                      'suffixes': []}
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
    
    return constructs_list


def generate_clips_df(constructs_list):
    """Generates a dataframe containing information about all unique CLIP reactions required to synthesise the constructs 
    in constructs_list. Mag_well and Clip plate info identified and added for each clip. Multiple entries added when required

    """
    merged_construct_dfs = pd.concat(constructs_list, ignore_index=True)
    unique_clips_df = merged_construct_dfs.drop_duplicates()
    unique_clips_df = unique_clips_df.reset_index(drop=True)
    clips_df = unique_clips_df.copy()

    def count_unique_clips(clips_df, merged_construct_dfs):
        ''' Count number of each CLIP reaction
            Iterates through unique clip list, counts number of instances of each unique clip, adds this as a new column'''
        clip_count = np.zeros(len(clips_df.index))
        for i, unique_clip in clips_df.iterrows():
            for _, clip in merged_construct_dfs.iterrows():
                if unique_clip.equals(clip):
                    clip_count[i] += 1
        clip_count = clip_count // FINAL_ASSEMBLIES_PER_CLIP + 1
        clips_df['number'] = [int(i) for i in clip_count.tolist()]
        return clips_df
    
    clips_df = count_unique_clips(unique_clips_df, merged_construct_dfs)
    
    # Associate well/s and plate/s for each CLIP reaction
    clips_df['mag_well'] = pd.Series(['0'] * len(clips_df.index), index=clips_df.index)
    clips_df['plate'] = pd.Series(['0'] * len(clips_df.index), index=clips_df.index)
    clip_count = 0
    
    for unique_clip_count, clip_number in clips_df['number'].iteritems():   # clip number is the number of each unique clip that must be built
        mag_wells = []
        plates = []         

        for well in range(clip_count, clip_count + clip_number):    # generate mag well and plate value(s) for a given unique clip
            mag_wells.append(tip_counter(well % 96))
            plates.append(1 + well//96)                             # plus one for one indexing
            
        clips_df.at[unique_clip_count, 'mag_well'] = tuple(mag_wells)
        clips_df.at[unique_clip_count, 'plate'] = tuple(plates)
        clip_count += clip_number       # increase count of all previous clips by amount for current unique clip

    return clips_df


def generate_sources_dict(paths):
    """Imports csvs files containing a series of parts/linkers with
    corresponding information into a dictionary where the key corresponds with
    part/linker and the value contains a tuple of corresponding information.

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
                    csv_values = source[1:]
                    csv_values.append(SOURCE_DECK_POS[deck_index])
                    sources_dict[str(source[0])] = tuple(csv_values)
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
        for _, clip_info in clips_df.iterrows():
            prefix_linker = clip_info['prefixes']
            clips_dict['prefixes_wells'].append([sources_dict[prefix_linker][0]]
                                                * clip_info['number'])
            clips_dict['prefixes_plates'].append(
                [handle_2_columns(sources_dict[prefix_linker])[2]] * clip_info['number'])
            suffix_linker = clip_info['suffixes']
            clips_dict['suffixes_wells'].append([sources_dict[suffix_linker][0]]
                                                * clip_info['number'])
            clips_dict['suffixes_plates'].append(
                [handle_2_columns(sources_dict[suffix_linker])[2]] * clip_info['number'])
            part = clip_info['parts']
            clips_dict['parts_wells'].append([sources_dict[part][0]]
                                             * clip_info['number'])
            clips_dict['parts_plates'].append([handle_2_columns(sources_dict[part])[2]]
                                              * clip_info['number'])
            if not sources_dict[part][1]:                               # add default vols for part and water if not user defined
                clips_dict['parts_vols'].append([DEFAULT_PART_VOL] *
                                                clip_info['number'])
                clips_dict['water_vols'].append([max_part_vol - DEFAULT_PART_VOL]
                                                * clip_info['number'])
            else:                                                       # add bespoke vols for part and water if part has a declared conc
                part_vol = round(
                    PART_PER_CLIP / float(sources_dict[part][1]), 1)
                if part_vol < MIN_VOL:
                    part_vol = MIN_VOL
                elif part_vol > max_part_vol:
                    part_vol = max_part_vol
                water_vol = max_part_vol - part_vol
                clips_dict['parts_vols'].append(
                    [part_vol] * clip_info['number'])
                clips_dict['water_vols'].append(
                    [water_vol] * clip_info['number'])
                    
    except KeyError:
        sys.exit('likely part/linker not listed in sources.csv')
    
    for key, value in clips_dict.items():                               # unlist all sublist in clips dict to yield a single list as the value for every key
        clips_dict[key] = [item for sublist in value for item in sublist]

    return clips_dict


def generate_clips_dict_list(clips_df, sources_dict):
    '''Subsets the clips df into chunks of 48, runs the generate_clips_dict function 
    for each and returns a list of the resulting sub clips dicts'''

    def clip_df_long_format(clip_df):
        ''' Takes clips df and returns a long format df with one row per mag well 
        rather than one row per unique clip'''

        df = clip_df.copy()
        long_clip_df = pd.DataFrame(columns = df.columns)

        index = 0

        for row, clip in clip_df.iterrows():
            for clip_rep in range(clip.number):                         # copy row and replace mag well and plate lists with well and plate for current rep
                data = clip.copy()
                data.number = 1
                data.mag_well = clip.mag_well[clip_rep]
                data.plate = clip.plate[clip_rep]
                long_clip_df.loc[index] = data
                index += 1

        return long_clip_df
    
    long_clip_df = clip_df_long_format(clips_df)

    # CLIP_COUNT = len(clips_df)
    CLIP_COUNT = clips_df['number'].sum()
    CLIP_PLATE_COUNT = CLIP_COUNT // MAX_CLIPS_PER_PLATE + 1            # plus one to include final partially full plate

    # Error 
    if clips_df['number'].sum() > MAX_CLIPS_TOTAL:
        raise ValueError(
            'Number of CLIP reactions exceeds {}. Reduce number of constructs in construct.csv.'.format(MAX_CLIPS_TOTAL))

    clips_dict_list = []

    for plate in range(CLIP_PLATE_COUNT):
        subset_lower = (plate * MAX_CLIPS_PER_PLATE)                    # set upper and lower bounds for subset of clips for a given plate
        subset_upper = subset_lower + MAX_CLIPS_PER_PLATE

        if subset_upper > CLIP_COUNT:                                   # set total number number of clips as upper bound if plate incomplete
            subset_upper = CLIP_COUNT
    
        # sub_clip_df = clips_df.iloc[subset_lower:subset_upper, :]
        sub_clip_df = long_clip_df.iloc[subset_lower:subset_upper, :]
        sub_clip_dict = generate_clips_dict(sub_clip_df, sources_dict)
        clips_dict_list.append(sub_clip_dict)                           # generate and append sub_clip_dict to list - allows for multiple clip reactions
        print(sub_clip_dict, '\n\n')
    return clips_dict_list


def generate_final_assembly_dict(constructs_list, clips_df):
    """Using constructs_list and clips_df, returns a dictionary of final
    assemblies with keys defining destination plate well positions and values
    indicating which clip reaction wells are used.

    """
    # final_assembly_dict = {}                                                      # old version generates dictionary
    final_assembly_dict_keys = []
    final_assembly_dict_values = []

    clips_count = np.zeros(len(clips_df.index))

    for construct_index, construct_df in enumerate(constructs_list):                # for each construct df in constructs_list
        construct_well_list = []
        construct_plate_list = []

        for _, clip in construct_df.iterrows():                                     # for each clip in construct
            clip_info = clips_df[(clips_df['prefixes'] == clip['prefixes']) &       # find clips in clips_df that match required clip
                                 (clips_df['parts'] == clip['parts']) &
                                 (clips_df['suffixes'] == clip['suffixes'])]
            
            clip_num = int(clip_info.index[0])                                      # row index of clip in clips_df
            clip_wells = clip_info.at[clip_num, 'mag_well']                         # list of all mag_wells for clip
            clip_plates = clip_info.at[clip_num, 'plate']                           # list of all plates for clip

            chosen_well = int(clips_count[clip_num] // FINAL_ASSEMBLIES_PER_CLIP)   # next viable mag well for clip
            clip_well = clip_wells[chosen_well]
            clip_plate = clip_plates[chosen_well]
            construct_well_list.append(clip_well)
            construct_plate_list.append(clip_plate)

            clips_count[clip_num] = clips_count[clip_num] + 1

        # final_assembly_dict[tip_counter(construct_index)] = [construct_well_list, construct_plate_list]   # generate dictionary
        final_assembly_dict_keys.append(tip_counter(construct_index))
        final_assembly_dict_values.append([construct_well_list, construct_plate_list])

    # return final_assembly_dict
    return final_assembly_dict_keys, final_assembly_dict_values                     # return list of dict keys and values


def generate_final_assembly_dict_list(constructs_list, clips_df):
    ''' Runs the generate_final_assembly_dict function producing ORDERED lists of 
    keys and values for the assembly source and destination (index of list denotes 
    the pair). Subsets key and value lists into chunks of up to 96 (max assemblies 
    per plate). A sub dictionary is generated for each chunk and then function  
    returns a list of the resulting sub assembly dicts
    
    This method ensures that subsequent assembly scripts do not reuse empty clip 
    wells while also keeping the correct destination well locations '''

    ASSEMBLY_COUNT = len(constructs_list)
    ASSEMBLY_PLATE_COUNT = ASSEMBLY_COUNT // MAX_ASSEMBLIES_PER_PLATE + 1           # plus one to include final partially full plate

    # # Error ########################### NO ERROR AS NO MAX ASSEMBLIES
    # if clips_df['number'].sum() > MAX_CLIPS_TOTAL:
    #     raise ValueError(
    #         'Number of CLIP reactions exceeds {}. Reduce number of constructs in construct.csv.'.format(MAX_CLIPS_TOTAL))

    # final_assembly_dict = generate_final_assembly_dict(constructs_list, clips_df)
    final_assembly_dict_keys, final_assembly_dict_values = generate_final_assembly_dict(constructs_list, clips_df)

    assembly_dict_list = []

    for plate in range(ASSEMBLY_PLATE_COUNT):
        subset_lower = (plate * MAX_ASSEMBLIES_PER_PLATE)               # set upper and lower bounds for subset of assemblies for a given plate
        subset_upper = subset_lower + MAX_ASSEMBLIES_PER_PLATE

        if subset_upper > ASSEMBLY_COUNT:                               # set total number number of assemblies as upper bound if plate incomplete
            subset_upper = ASSEMBLY_COUNT

        # sub_assembly_df = constructs_list[subset_lower:subset_upper]
        # sub_assembly_dict = generate_final_assembly_dict(sub_assembly_df, clips_df)

        keys = [tip_counter(i) for i in list(range(subset_upper - subset_lower))]
        values = final_assembly_dict_values[subset_lower:subset_upper]
        sub_assembly_dict = {keys[i]: values[i] for i in range(len(keys))}

        assembly_dict_list.append(sub_assembly_dict)                    # generate and append sub_clip_dict to list - allows for multiple clip reactions

    return assembly_dict_list


def calculate_final_assembly_tipracks(final_assembly_dict):
    """Calculates the number of final assembly tipracks required ensuring
    no more than MAX_FINAL_ASSEMBLY_TIPRACKS are used.

    """
    final_assembly_lens = []
    # final_assembly_dict = final_assembly_dict[0]

    for values in final_assembly_dict.values():
        final_assembly_lens.append(len(values[0]))                      # create a list of number of clips per assembly 
        
    master_mix_tips = len(list(set(final_assembly_lens)))               # number of different build lengths in construct build
    total_tips = master_mix_tips + sum(final_assembly_lens) 
    final_assembly_tipracks = (total_tips - 1) // 96 + 1
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
    wells = [tip_counter(x) for x in range(len(constructs_list))]
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


def generate_ot2_script(ot2_script_path, template_path, **kwargs):
    """Generates an ot2 script named 'ot2_script_path', where kwargs are
    written as global variables at the top of the script. For each kwarg, the
    keyword defines the variable name while the value defines the name of the
    variable. The remainder of template file is subsequently written below.

    """
    with open(ot2_script_path, 'w') as wf:
        with open(template_path, 'r') as rf:
            for index, line in enumerate(rf):
                if line[:3] == 'def':
                    function_start = index
                    break
                else:
                    wf.write(line)
            for key, value in kwargs.items():
                wf.write('{}='.format(key))
                if type(value) == dict:
                    wf.write(json.dumps(value))
                elif type(value) == str:
                    wf.write("'{}'".format(value))
                else:
                    wf.write(str(value))
                wf.write('\n')
            wf.write('\n')
        with open(template_path, 'r') as rf:
            for index, line in enumerate(rf):
                if index >= function_start - 1:
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


def generate_sources_paths_df(paths, deck_positions):
    """Generates a dataframe detailing source plate information.

    Args:
        paths (list): list of strings specifying paths to source plates.
        deck_positions (list): list of strings specifying candidate deck positions.

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


def handle_2_columns(datalist):
    """This function has the intent of changing:
    ('A8', '2') => ('A8', '', '2')
    ('A8', '', '2') => ('A8', '', '2')
    [('E2', '5')] => [('E2', '', '5')]
    [('G1', '', '5')] => [('G1', '', '5')]
    with the purpose of handling 2 column csv part file inputs,
    as at times when 2 column csv files are input it creates tuples
    of length 2 instead of 3
    """
    return_list = 0
    if isinstance(datalist,list):
        datalist = datalist[0]
        return_list = 1
    if len(datalist) == 2:
        datalist = list(datalist)
        datalist.insert(1,"")
        datalist = tuple(datalist)
    if return_list:
        mylist = [""]
        mylist[0] = datalist
        return mylist
    return datalist


def counter(rows):
    def inner(n):
        """ Takes either a value or a well location and converts to other fomat """
        
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
        inv_row_dict = {v: k for k, v in row_dict.items()}

        if type(n) == int:
            row = row_dict[n % rows]
            col = 1 + n // rows
            return row + f'{col}'
            # return row + f'{col:02d}' # for if 2 sf number required (i.e. 'A01' rather than 'A1')

        elif type(n) == str:
            row, col = re.findall('\d+|\D+', n)
            col = int(col) - 1

            return col*rows + inv_row_dict[row]

    return inner
tip_counter = counter(8)


if __name__ == '__main__':
    main()
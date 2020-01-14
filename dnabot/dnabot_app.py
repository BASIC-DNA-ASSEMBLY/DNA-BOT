# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:26:07 2019

@author: mh2210
"""
import os
import csv
import pandas as pd
import numpy as np
import json
import sys
import dnabot_gui as gui
import tkinter as tk
import mplates

# Constant str
TEMPLATE_DIR_NAME = 'template_ot2_scripts'
CLIP_TEMP_FNAME = 'clip_template.py'
MAGBEAD_TEMP_FNAME = 'purification_template.py'
F_ASSEMBLY_TEMP_FNAME = 'assembly_template.py'
TRANS_SPOT_TEMP_FNAME = 'transformation_template.py'
CLIP_FNAME = '1_clip.ot2.py'
MAGBEAD_FNAME = '2_purification.ot2.py'
F_ASSEMBLY_FNAME = '3_assembly.ot2.py'
TRANS_SPOT_FNAME = '4_transformation.ot2.py'
CLIPS_INFO_FNAME = 'clip_run_info.csv'
FINAL_ASSEMBLIES_INFO_FNAME = 'final_assembly_run_info.csv'
WELL_OUTPUT_FNAME = 'wells.txt'

# Constant floats/ints
CLIP_DEAD_VOL = 60
CLIP_VOL = 30
T4_BUFF_VOL = 3
BSAI_VOL = 1
T4_LIG_VOL = 0.5
CLIP_MAST_WATER = 15.5
PART_PER_CLIP = 200
MIN_VOL = 1
MAX_CONSTRUCTS = 96
MAX_CLIPS = 48
FINAL_ASSEMBLIES_PER_CLIP = 15
DEFAULT_PART_VOL = 1
MAX_SOURCE_PLATES = 6
MAX_FINAL_ASSEMBLY_TIPRACKS = 7

# Constant dicts
SPOTTING_VOLS_DICT = {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5}

# Constant lists
SOURCE_DECK_POS = ['2', '5', '8', '7', '10', '11']


def main():
    # Parent directories
    generator_dir = os.getcwd()
    template_dir_path = os.path.join(generator_dir, TEMPLATE_DIR_NAME)

    # Obtain user input
    print("Requesting user input, if not visible checked minimized windows.")
    root = tk.Tk()
    dnabotinst = gui.DnabotApp(root)
    root.mainloop()
    root.destroy()
    if dnabotinst.quit_status:
        sys.exit("User specified 'QUIT' during app.")
    root = tk.Tk()
    construct_path = gui.UserDefinedPaths(root, 'Construct csv file')
    root.destroy()
    root = tk.Tk()
    sources_paths = gui.UserDefinedPaths(root, 'Sources csv files',
                                         multiple_files=True)
    if len(sources_paths.output) > len(SOURCE_DECK_POS):
        raise ValueError(
            'Number of source plates exceeds deck positions.')
    root.destroy()
    os.chdir(os.path.dirname(construct_path.output))
    construct_base = os.path.basename(construct_path.output)
    construct_base = os.path.splitext(construct_base)[0]
    print('User input successfully collected.')

    # Process input csv files
    print('Processing input csv files...')
    constructs_list = generate_constructs_list(construct_path.output)
    clips_df = generate_clips_df(constructs_list)
    sources_dict = generate_sources_dict(sources_paths.output)

    # calculate OT2 script variables
    print('Calculating OT-2 variables...')
    clips_dict = generate_clips_dict(clips_df, sources_dict)
    magbead_sample_number = clips_df['number'].sum()
    final_assembly_dict = generate_final_assembly_dict(constructs_list,
                                                       clips_df)
    final_assembly_tipracks = calculate_final_assembly_tipracks(
        final_assembly_dict)
    spotting_tuples = generate_spotting_tuples(constructs_list,
                                               SPOTTING_VOLS_DICT)

    print('Writing files...')
    # Write OT2 scripts
    generate_ot2_script(CLIP_FNAME, os.path.join(
        template_dir_path, CLIP_TEMP_FNAME), clips_dict=clips_dict)
    generate_ot2_script(MAGBEAD_FNAME, os.path.join(
        template_dir_path, MAGBEAD_TEMP_FNAME),
        sample_number=magbead_sample_number,
        ethanol_well=dnabotinst.etoh_well)
    generate_ot2_script(F_ASSEMBLY_FNAME, os.path.join(
        template_dir_path, F_ASSEMBLY_TEMP_FNAME),
        final_assembly_dict=final_assembly_dict,
        tiprack_num=final_assembly_tipracks)
    generate_ot2_script(TRANS_SPOT_FNAME, os.path.join(
        template_dir_path, TRANS_SPOT_TEMP_FNAME),
        spotting_tuples=spotting_tuples,
        soc_well=f"A{dnabotinst.soc_column}")

    # Write non-OT2 scripts
    if 'metainformation' in os.listdir():
        pass
    else:
        os.makedirs('metainformation')
    os.chdir('metainformation')
    master_mix_df = generate_master_mix_df(clips_df['number'].sum())
    sources_paths_df = generate_sources_paths_df(
        sources_paths.output, SOURCE_DECK_POS)
    dfs_to_csv(construct_base + '_' + CLIPS_INFO_FNAME, index=False,
               MASTER_MIX=master_mix_df, SOURCE_PLATES=sources_paths_df,
               CLIP_REACTIONS=clips_df)
    with open(construct_base + '_' + FINAL_ASSEMBLIES_INFO_FNAME,
              'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for final_assembly_well, construct_clips in final_assembly_dict.items():
            csvwriter.writerow([final_assembly_well, construct_clips])
    with open(construct_base + '_' + WELL_OUTPUT_FNAME, 'w') as f:
        f.write('Magbead ethanol well: {}'.format(dnabotinst.etoh_well))
        f.write('\n')
        f.write('SOC column: {}'.format(dnabotinst.soc_column))
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
            """Interogates linker to determine if the suffix linker is a UTR
            linker.

            """
            if len(linker) >= 4:
                if linker[:3] == 'UTR':
                    return linker[:4] + '-S'
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

    # Errors
    if len(constructs_list) > MAX_CONSTRUCTS:
        raise ValueError(
            f"Number of constructs exceeds maximum ({MAX_CONSTRUCTS}). Reduce construct number in construct.csv.")
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
            f"Number of CLIP reactions exceeds {MAX_CLIPS}. Reduce number of constructs in construct.csv.")

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
    for index, number in clips_df['number'].iteritems():
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
                [sources_dict[prefix_linker][2]] * clip_info['number'])
            suffix_linker = clip_info['suffixes']
            clips_dict['suffixes_wells'].append([sources_dict[suffix_linker][0]]
                                                * clip_info['number'])
            clips_dict['suffixes_plates'].append(
                [sources_dict[suffix_linker][2]] * clip_info['number'])
            part = clip_info['parts']
            clips_dict['parts_wells'].append([sources_dict[part][0]]
                                             * clip_info['number'])
            clips_dict['parts_plates'].append([sources_dict[part][2]]
                                              * clip_info['number'])
            if not sources_dict[part][1]:
                clips_dict['parts_vols'].append([DEFAULT_PART_VOL] *
                                                clip_info['number'])
                clips_dict['water_vols'].append([max_part_vol - DEFAULT_PART_VOL]
                                                * clip_info['number'])
            else:
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
    for key, value in clips_dict.items():
        clips_dict[key] = [item for sublist in value for item in sublist]
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
    final_assembly_lens = []
    for values in final_assembly_dict.values():
        final_assembly_lens.append(len(values))
    master_mix_tips = len(list(set(final_assembly_lens)))
    total_tips = master_mix_tips + sum(final_assembly_lens)
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


if __name__ == '__main__':
    main()

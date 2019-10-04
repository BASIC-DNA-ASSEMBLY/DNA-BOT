# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:26:07 2019

@author: mh2210
"""
#from tkinter import *
import os
import csv
import pandas as pd
import numpy as np
import json

# Constant str
CSV_PATH = r"C:\Users\mh2210\Box Sync\Work_box\1_Post_doctoral_research\Experiments\2019-04-18_linker_testing\Assembly_1" # Path to CSV files temparily given as a constant.
CONSTRUCTS_CSV_FNAME = 'constructs.csv'
SOURCES_CSV_FNAME = 'sources.csv'
CLIPS_FNAME = 'clip_reactions.xlsx'
MASTER_MIX_FNAME = 'clips_master_mix.xlsx'
TEMPLATE_DIR_NAME = 'template_ot2_scripts'
CLIP_TEMP_FNAME = 'clip_template.py'
MAGBEAD_TEMP_FNAME = 'magbead_template.py'
F_ASSEMBLY_TEMP_FNAME = 'final_assembly_template.py'
TRANS_SPOT_TEMP_FNAME = 'transformation_spotting_template.py'
CLIP_FNAME = 'clip.ot2.py'
MAGBEAD_FNAME = 'magbead.ot2.py'
F_ASSEMBLY_FNAME = 'final_assembly.ot2.py'
TRANS_SPOT_FNAME = 'transformation_spotting.ot2.py'

# Constant floats/ints
CLIP_MAST_DEAD_VOL = 60
CLIP_VOL = 30
T4_BUFF_VOL = 3
BSAI_VOL = 1
T4_LIG_VOL = 0.5
CLIP_MAST_WATER = 15.5
PART_PER_CLIP = 200
FINAL_ASSEMBLIES_PER_CLIP = 20
DEFAULT_PART_VOL = 1
MAX_SOURCE_PLATES = 6


def main():
    #    user_specified_path()
    generator_dir = os.getcwd()
    os.chdir(CSV_PATH)
    template_dir_path = os.path.join(generator_dir, TEMPLATE_DIR_NAME)
    
    # Process input csv files and write non-OT2 files
    constructs_list = generate_constructs_list(CONSTRUCTS_CSV_FNAME)
    clips_df = generate_clips_df(constructs_list)
    sources_dict = generate_sources_dict(SOURCES_CSV_FNAME)
    clips_df.to_excel(CLIPS_FNAME, index=False)
    generate_master_mix(clips_df['number'].sum(), MASTER_MIX_FNAME)
    
    
    # calculate OT2 script variables
    clips_dict = generate_clips_dict(clips_df, sources_dict)
    magbead_sample_number = clips_df['number'].sum()
    
    # Write OT2 scripts
    generate_ot2_script(CLIP_FNAME, os.path.join(
            template_dir_path, CLIP_TEMP_FNAME), clips_dict=clips_dict)
    generate_ot2_script(MAGBEAD_FNAME, os.path.join(
            template_dir_path, MAGBEAD_TEMP_FNAME),
    sample_number=magbead_sample_number)
    

def user_specifications():
    pass


def generate_constructs_list(path):
    """Generates a list of dataframes corresponding to each construct. Each 
    dataframe lists components of the CLIP reactions required."""

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

        clips_dict = {'prefixes': [], 'parts': [],
                      'suffixes': []}
        for i, sequence in enumerate(construct):
            if i % 2 != 0:
                clips_dict['parts'].append(sequence)
                clips_dict['prefixes'].append(
                    construct[i - 1] + '-P')
                if i == len(construct) - 1:
                    suffix_linker = interogate_linker(construct[0])
                    clips_dict['suffixes'].append(suffix_linker)
                else:
                    suffix_linker = interogate_linker(construct[i + 1])
                    clips_dict['suffixes'].append(suffix_linker)           
        return pd.DataFrame.from_dict(clips_dict)

    constructs_list = []
    with open(path, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        for construct in csv_reader:
            constructs_list.append(process_construct(construct))     
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
    if len(unique_clips_df.index) > 48:
        raise ValueError(
            'Number of CLIP reactions exceeds 48. Reduce number of constructs in construct.csv.')

    # Count number of each CLIP reaction
    clip_count = np.zeros(len(clips_df.index))
    for i, unique_clip in unique_clips_df.iterrows():
        for j, clip in merged_construct_dfs.iterrows():
            if unique_clip.equals(clip):
                clip_count[i] = clip_count[i] + 1            
    clip_count = clip_count//FINAL_ASSEMBLIES_PER_CLIP + 1
    clips_df['number'] = [int(i) for i in clip_count.tolist()]

    # Associate well/s for each CLIP reaction
    clips_df['mag_well'] = pd.Series(['0'] * len(clips_df.index),
                                              index=clips_df.index)
    for index, number in clips_df['number'].iteritems():
        if index == 0:
            mag_wells = []
            for x in range(number):
                mag_wells.append(final_well(x + 1 + 48))
            clips_df.at[index, 'mag_well'] = mag_wells
        else:
            mag_wells = []
            for x in range(number):
                well_count = clips_df.loc[
                    :index - 1, 'number'].sum() + x + 1 + 48
                mag_wells.append(final_well(well_count))
            clips_df.at[index, 'mag_well'] = mag_wells         
    return clips_df


def generate_sources_dict(path):
    """Imports a csv file containing a series of parts/linkers with 
    corresponding information into a dictionary where the key corresponds with
    part/linker and the value contains a list of corresponding information.
    """
    sources_dict = {}
    with open(path, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        for index, source in enumerate(csv_reader):
            if index != 0:
                sources_dict[str(source[0])] = source[1:]           
    return sources_dict


def generate_master_mix(clip_number, master_mix_fname):
    """Generates a .xlsx file detailing the components required in the clip 
    reaction master mix.
    """
    COMPONENTS = {'Component': ['Promega T4 DNA Ligase buffer, 10X',
                  'Water', 'NEB BsaI-HFv2',
                  'Promega T4 DNA Ligase']}
    VOL_COLUMN = 'Volume (uL)'
    master_mix_df = pd.DataFrame.from_dict(COMPONENTS)
    master_mix_df[VOL_COLUMN] = clip_number * np.array([T4_BUFF_VOL,
                 CLIP_MAST_WATER,
                 BSAI_VOL,
                 T4_LIG_VOL])
    master_mix_df.to_excel(master_mix_fname, index = False)


def generate_clips_dict(clips_df, sources_dict):
    """Using clips_df and sources_dict, returns a clips_dict which acts as the 
    sole variable for the opentrons script "clip.ot2.py".
    """
    max_part_vol = CLIP_VOL - (T4_BUFF_VOL + BSAI_VOL + T4_LIG_VOL
                               + CLIP_MAST_WATER - 2)
    clips_dict = {'prefixes_wells': [], 'prefixes_plates': [],
                 'suffixes_wells': [], 'suffixes_plates': [],
                 'parts_wells': [], 'parts_plates': [], 'parts_vols': [],
                 'water_vols': []}
    for index, clip_info in clips_df.iterrows():
        prefix_linker = clip_info['prefixes']
        clips_dict['prefixes_wells'].append(sources_dict[prefix_linker][0])
        clips_dict['prefixes_plates'].append(int(sources_dict[prefix_linker][1]))
        
        suffix_linker = clip_info['suffixes']
        clips_dict['suffixes_wells'].append(sources_dict[suffix_linker][0])
        clips_dict['suffixes_plates'].append(int(sources_dict[suffix_linker][1]))
        
        part = clip_info['parts']
        clips_dict['parts_wells'].append(sources_dict[part][0])
        clips_dict['parts_plates'].append(int(sources_dict[part][1]))
        if not sources_dict[part][2]:
            clips_dict['parts_vols'].append(DEFAULT_PART_VOL)
            clips_dict['water_vols'].append(max_part_vol - DEFAULT_PART_VOL)
        else:
            part_vol = (round(PART_PER_CLIP/float(sources_dict[part][2]), 1)
            if round(PART_PER_CLIP/float(sources_dict[part][2]), 1)
            <= max_part_vol else max_part_vol)
            water_vol = max_part_vol - part_vol
            clips_dict['parts_vols'].append(part_vol)
            clips_dict['water_vols'].append(water_vol)
            
    # Error
    source_plate_number = max(clips_dict['prefixes_plates'] + 
                              clips_dict['suffixes_plates'] + 
                              clips_dict['parts_plates'])
    if source_plate_number > MAX_SOURCE_PLATES:
        raise ValueError('Number of source plates is {}, maximum allowed is {}'.format(source_plate_number, MAX_SOURCE_PLATES))
    else:
        return clips_dict


def generate_ot2_script(ot2_script_fname, template_path, **kwargs):
    """Generates an ot2 script named 'ot2_script_fname', where kwargs are 
    written as global variables at the top of the script. For each kwarg, the 
    keyword defines the variable name while the value defines the name of the 
    variable. The remainder of template file is subsequently written below.        
    """
    with open(ot2_script_fname, 'w') as wf:
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
                else:
                    wf.write(str(value))
            wf.write('\n')
        with open(template_path, 'r') as rf:   
             for index, line in enumerate(rf):
                if index >= function_start - 1:
                    wf.write(line)


def final_well(sample_number):
    """Determines well containing the final sample from sample number.
    """
    letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    final_well_column = sample_number // 8 + \
        (1 if sample_number % 8 > 0 else 0)
    final_well_row = letter[sample_number - (final_well_column - 1) * 8 - 1]
    return final_well_row + str(final_well_column)


if __name__ == '__main__':
    main()
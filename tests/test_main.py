# -*- coding: utf-8 -*-

import pytest
import subprocess
from pathlib import Path
from os.path import join


in_file_construct = Path(__file__).resolve().parent / 'inputs' / 'constructs.csv'
in_file_linker = Path(__file__).resolve().parent / 'inputs' / 'linker_parts_coords.csv'
in_file_user = Path(__file__).resolve().parent / 'inputs' / 'user_parts_coords.csv'
in_default_settings = Path(__file__).resolve().parent / "inputs" / "default_settings.yaml"

files_to_tests = [
    '1_clip_ot2_APIv1.py',
    '1_clip_ot2_APIv2.8.py',
    '1_clip_ot2_Thermocycler_APIv2.8.py',
    '2_purification_ot2_APIv1.py',
    '2_purification_ot2_APIv2.8.py',
    '3_assembly_ot2_APIv1.py',
    '3_assembly_ot2_APIv2.8.py',
    '3_assembly_ot2_Thermocycler_APIv2.8.py',
    '4_transformation_ot2_APIv1.py',
    '4_transformation_ot2_APIv2.8.py',
    '4_transformation_ot2_Thermocycler_APIv2.8.py',
    '4_transformation_ot2_Thermocycler_12wellplate_APIv2.8.py',
    join('metainformation','constructs_wells.txt'),
    join('metainformation','constructs_final_assembly_run_info.csv'),
    join('metainformation','constructs_deck.md')
    ]


def test_output(tmpdir):

    import difflib

    _ = subprocess.run([
        'python','-m', 'dnabot.dnabot_app',
        '--default_settings_file', in_default_settings,
        'nogui',
        '--construct_path', in_file_construct,
        '--source_paths', in_file_linker, in_file_user,
        '--output_dir', tmpdir,
        ],stdout=subprocess.PIPE)

    for fname in files_to_tests:
        out_file = tmpdir / fname
        with open(out_file) as out_fh:
            out_str = out_fh.readlines()
        ref_file = Path(__file__).resolve().parent / 'outputs' / fname
        with open(ref_file) as ref_fh:
            ref_str = ref_fh.readlines()
    
        try:
            assert out_str == ref_str
        except AssertionError as e:
            print(f'Difference in {fname}...')
            for diff in difflib.context_diff(out_str, ref_str):
                print(diff)
            raise e
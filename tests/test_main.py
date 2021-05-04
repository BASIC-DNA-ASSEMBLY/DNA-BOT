# -*- coding: utf-8 -*-

import os
import pytest
import subprocess

path = os. getcwd()

tool = os.path.join(path, 'dnabot', 'dnabot_app.py')
in_file_construct = os.path.join(path, 'tests', 'inputs', 'constructs.csv')
in_file_linker = os.path.join(path, 'tests', 'inputs', 'linker_parts_coords.csv')
in_file_user = os.path.join(path, 'tests', 'inputs', 'user_parts_coords.csv')

files_to_tests=['1_clip.ot2.py', '2_purification.ot2.py', '3_assembly.ot2.py', '4_transformation.ot2.py', os.path.join('metainformation','constructs_wells.txt'), os.path.join('metainformation','constructs_final_assembly_run_info.csv'), os.path.join('metainformation','constructs_clip_run_info.csv')]


def test_output(tmpdir):

    #compare the contents of output files actual vs expected
    #python dnabot/dnabot_app.py nogui --construct_path inputs/constructs.csv --source_paths inputs/linker_parts_coords.csv /home/kenza/Bureau/conda_packaging/DNA-BOT/inputs/user_parts_coords.csv --output_dir bot_py36_out
    
    result = subprocess.run(['python', tool, 'nogui', '--construct_path', in_file_construct, '--source_paths', in_file_linker, in_file_user, '--output_dir', tmpdir], stdout=subprocess.PIPE)

    for i in files_to_tests:
        out_file = tmpdir.join(i).strpath
        with open(out_file) as f:
            actual = f.read()

        fname = os.path.join('tests', 'outputs', i)
        with open(fname) as f:
            expected = f.read()

        print('actual:', actual)
        print('expected:', expected)
        assert actual == expected
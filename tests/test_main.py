# -*- coding: utf-8 -*-

#import os
import pytest
import subprocess
from pathlib import Path
from os.path import dirname, abspath, join


tool = join(dirname(dirname(abspath(__file__))), 'dnabot', 'dnabot_app.py')
in_file_construct = Path(__file__).resolve().parent / 'inputs' / 'constructs.csv'
in_file_linker = Path(__file__).resolve().parent / 'inputs' / 'linker_parts_coords.csv'
in_file_user = Path(__file__).resolve().parent / 'inputs' / 'user_parts_coords.csv'


files_to_tests=['1_clip.ot2.py', '2_purification.ot2.py', '3_assembly.ot2.py', '4_transformation.ot2.py', join('metainformation','constructs_wells.txt'), join('metainformation','constructs_final_assembly_run_info.csv')]


def test_output(tmpdir):

    #compare the contents of output files actual vs expected
    #python DNA-BOT/dnabot/dnabot_app.py nogui --construct_path DNA-BOT/tests/inputs/constructs.csv --source_paths DNA-BOT/tests/inputs/linker_parts_coords.csv DNA-BOT/tests/inputs/user_parts_coords.csv --output_dir tmpdir
    
    result = subprocess.run(['python', tool, 'nogui', '--construct_path', in_file_construct, '--source_paths', in_file_linker, in_file_user, '--output_dir', tmpdir], stdout=subprocess.PIPE)

    for i in files_to_tests:
        out_file = tmpdir.join(i).strpath
        with open(out_file) as f:
            actual = f.read()

        fname = Path(__file__).resolve().parent / 'outputs' / i
        with open(fname) as f:
            expected = f.read()

        print('actual:', actual)
        print('expected:', expected)
        assert actual == expected
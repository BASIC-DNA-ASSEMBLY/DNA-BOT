from opentrons.simulate import simulate, format_runlog
import sys
import os

# Add the path to the test scripts directory
test_scripts_path = os.path.join(os.path.dirname(__file__), '..', 'test_scripts', 'clip_optimisation_test', 'dnabot_run_20250730_150731', 'OT2_Scripts')
sys.path.append(test_scripts_path)

# Import the run function from the C1 assembly script
from C1_assembly_ot2_Thermocycler_APIv2_8 import run

# read the file
protocol_file = open(os.path.join(test_scripts_path, 'C1_assembly_ot2_Thermocycler_APIv2.8.py'))
# simulate() the protocol, keeping the runlog
runlog, _bundle = simulate(protocol_file, os.path.join(test_scripts_path, 'C1_assembly_ot2_Thermocycler_APIv2.8.py'))
# print the runlog
print('\n', format_runlog(runlog), '\n', sep = '')

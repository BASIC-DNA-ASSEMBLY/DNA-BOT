from opentrons import labware, instruments, modules, robot
import numpy as np


final_assembly_dict={"A1": ["A7", "H7", "A8", "C8", "G8"], "B1": ["A7", "H7", "A8", "C8", "A9"], "C1": ["A7", "H7", "A8", "C9", "G8"], "D1": ["A7", "H7", "A8", "C9", "A9"], "E1": ["A7", "G9", "C8", "H9", "B10"], "F1": ["A7", "G9", "C8", "D10", "B10"], "G1": ["A7", "G9", "C9", "H9", "B10"], "H1": ["A7", "G9", "C9", "D10", "B10"], "A2": ["A7", "H7", "F10", "C8", "G8"], "B2": ["A7", "H7", "F10", "C8", "A9"], "C2": ["A7", "H7", "F10", "C9", "G8"], "D2": ["A7", "H7", "F10", "C9", "A9"], "E2": ["A7", "G9", "C8", "H9", "H10"], "F2": ["A7", "G9", "C8", "D10", "H10"], "G2": ["B7", "G9", "C9", "H9", "H10"], "H2": ["B7", "G9", "C9", "D10", "H10"], "A3": ["B7", "H7", "B11", "C8", "G8"], "B3": ["B7", "H7", "B11", "C8", "A9"], "C3": ["B7", "H7", "B11", "C9", "G8"], "D3": ["B7", "H7", "B11", "C9", "A9"], "E3": ["B7", "G9", "C8", "H9", "C11"], "F3": ["B7", "G9", "C8", "D10", "C11"], "G3": ["B7", "G9", "C9", "H9", "C11"], "H3": ["B7", "G9", "C9", "D10", "C11"], "A4": ["B7", "D11", "A8", "C8", "G8"], "B4": ["B7", "D11", "A8", "C8", "A9"], "C4": ["B7", "D11", "A8", "C9", "G8"], "D4": ["B7", "D11", "A8", "C9", "A9"], "E4": ["C7", "E11", "D8", "H9", "B10"], "F4": ["C7", "E11", "D8", "D10", "B10"], "G4": ["C7", "E11", "D9", "H9", "B10"], "H4": ["C7", "E11", "D9", "D10", "B10"], "A5": ["C7", "D11", "F10", "D8", "G8"], "B5": ["C7", "D11", "F10", "D8", "A9"], "C5": ["C7", "D11", "F10", "D9", "G8"], "D5": ["C7", "D11", "F10", "D9", "A9"], "E5": ["C7", "E11", "D8", "H9", "H10"], "F5": ["C7", "E11", "D8", "D10", "H10"], "G5": ["C7", "E11", "D9", "H9", "H10"], "H5": ["C7", "E11", "D9", "D10", "H10"], "A6": ["C7", "D11", "B11", "D8", "G8"], "B6": ["C7", "D11", "B11", "D8", "A9"], "C6": ["D7", "D11", "B11", "D9", "G8"], "D6": ["D7", "D11", "B11", "D9", "A9"], "E6": ["D7", "E11", "D8", "H9", "C11"], "F6": ["D7", "E11", "D8", "D10", "C11"], "G6": ["D7", "E11", "D9", "H9", "C11"], "H6": ["D7", "E11", "D9", "D10", "C11"], "A7": ["D7", "F11", "A8", "D8", "G8"], "B7": ["D7", "F11", "A8", "D8", "A9"], "C7": ["D7", "F11", "A8", "D9", "G8"], "D7": ["D7", "F11", "A8", "D9", "A9"], "E7": ["D7", "G11", "D8", "H9", "B10"], "F7": ["D7", "G11", "D8", "D10", "B10"], "G7": ["D7", "G11", "D9", "H9", "B10"], "H7": ["D7", "G11", "D9", "D10", "B10"], "A8": ["E7", "F11", "F10", "E8", "H8"], "B8": ["E7", "F11", "F10", "E8", "B9"], "C8": ["E7", "F11", "F10", "E9", "H8"], "D8": ["E7", "F11", "F10", "E9", "B9"], "E8": ["E7", "G11", "E8", "A10", "H10"], "F8": ["E7", "G11", "E8", "E10", "H10"], "G8": ["E7", "G11", "E9", "A10", "H10"], "H8": ["E7", "G11", "E9", "E10", "H10"], "A9": ["E7", "F11", "B11", "E8", "H8"], "B9": ["E7", "F11", "B11", "E8", "B9"], "C9": ["E7", "F11", "B11", "E9", "H8"], "D9": ["E7", "F11", "B11", "E9", "B9"], "E9": ["E7", "G11", "E8", "A10", "C11"], "F9": ["E7", "G11", "E8", "E10", "C11"], "G9": ["F7", "G11", "E9", "A10", "C11"], "H9": ["F7", "G11", "E9", "E10", "C11"], "A10": ["F7", "H11", "A8", "E8", "H8"], "B10": ["F7", "H11", "A8", "E8", "B9"], "C10": ["F7", "H11", "B8", "E9", "H8"], "D10": ["F7", "H11", "B8", "E9", "B9"], "E10": ["F7", "A12", "E8", "A10", "B10"], "F10": ["F7", "A12", "E8", "E10", "B10"], "G10": ["F7", "A12", "E9", "A10", "C10"], "H10": ["F7", "A12", "E9", "E10", "C10"], "A11": ["F7", "H11", "F10", "E8", "H8"], "B11": ["F7", "H11", "F10", "E8", "B9"], "C11": ["F7", "H11", "G10", "E9", "H8"], "D11": ["F7", "H11", "G10", "E9", "B9"], "E11": ["G7", "A12", "F8", "A10", "H10"], "F11": ["G7", "A12", "F8", "E10", "H10"], "G11": ["G7", "A12", "F9", "A10", "A11"], "H11": ["G7", "A12", "F9", "E10", "A11"]}
tiprack_num=5


def final_assembly(final_assembly_dict, tiprack_num):
    """Implements final assembly reactions using an opentrons OT-2.

    Args:
    final_assembly_dict (dict): Dictionary with keys and values corresponding to destination and associated linker-ligated part wells, respectively.
    tiprack_num (int): Number of tipracks required during run.

    """

    # Constants
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5', '8', '11']
    TIPRACK_TYPE = 'opentrons-tiprack-10ul'
    PIPETTE_MOUNT = 'right'
    MAG_PLATE_TYPE = '4ti-0960_FrameStar'
    MAG_PLATE_POSITION = '4'
    TUBE_RACK_TYPE = 'tube-rack_E1415-1500'
    TUBE_RACK_POSITION = '7'
    DESTINATION_PLATE_TYPE = '4ti-0960_FrameStar'
    DESTINATION_PLATE_POSITION = '1'
    MIX_SETTINGS = (1, 3)

    # Errors
    sample_number = len(final_assembly_dict.keys())
    if sample_number > 96:
        raise ValueError('Final assembly nummber cannot exceed 96.')

    # Tips and pipette
    slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
    tipracks = [labware.load(TIPRACK_TYPE, slot)
                for slot in slots]
    pipette = instruments.P10_Single(mount=PIPETTE_MOUNT, tip_racks=tipracks)

    # Define Labware
    magbead_plate = labware.load(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
    tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
    destination_plate = labware.load(
        DESTINATION_PLATE_TYPE, DESTINATION_PLATE_POSITION)

    # Master mix transfers
    final_assembly_lens = []
    for values in final_assembly_dict.values():
        final_assembly_lens.append(len(values))
    unique_assemblies_lens = list(set(final_assembly_lens))
    master_mix_well_letters = ['A', 'B', 'C', 'D']
    for x in unique_assemblies_lens:
        master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
        destination_inds = [i for i, lens in enumerate(
            final_assembly_lens) if lens == x]
        destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
        destination_wells = list(destination_wells[destination_inds])
        pipette.pick_up_tip()
        pipette.transfer(10 - x, tube_rack.wells(master_mix_well),
                         destination_plate.wells(destination_wells), new_tip='never')
        pipette.drop_tip()

    # Part transfers
    for key, values in list(final_assembly_dict.items()):
        pipette.transfer(1, magbead_plate.wells(values),
                         destination_plate.wells(key), mix_after=MIX_SETTINGS, new_tip='always')


final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

for c in robot.commands():
    print(c)
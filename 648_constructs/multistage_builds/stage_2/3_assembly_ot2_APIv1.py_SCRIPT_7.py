from opentrons import labware, instruments, modules, robot
import numpy as np


final_assembly_dict={"A1": [["G5", "C9", "B11"], [1, 2, 1]], "B1": [["G5", "C9", "C11"], [1, 2, 1]], "C1": [["G5", "C9", "D11"], [1, 2, 1]], "D1": [["G5", "C9", "E11"], [1, 2, 1]], "E1": [["G5", "C9", "F11"], [1, 2, 1]], "F1": [["G5", "C9", "G11"], [1, 2, 1]], "G1": [["G5", "C9", "H11"], [1, 2, 1]], "H1": [["G5", "C9", "A12"], [1, 2, 1]], "A2": [["G5", "C9", "B12"], [1, 2, 1]], "B2": [["H5", "D9", "B11"], [1, 2, 1]], "C2": [["H5", "D9", "C11"], [1, 2, 1]], "D2": [["H5", "D9", "D11"], [1, 2, 1]], "E2": [["H5", "D9", "E11"], [1, 2, 1]], "F2": [["H5", "D9", "F11"], [1, 2, 1]], "G2": [["H5", "D9", "G11"], [1, 2, 1]], "H2": [["H5", "D9", "H11"], [1, 2, 1]], "A3": [["H5", "D9", "A12"], [1, 2, 1]], "B3": [["H5", "D9", "B12"], [1, 2, 1]], "C3": [["H5", "E9", "F12"], [1, 2, 1]], "D3": [["H5", "E9", "G12"], [1, 2, 1]], "E3": [["H5", "E9", "H12"], [1, 2, 1]], "F3": [["H5", "E9", "A1"], [1, 2, 2]], "G3": [["H5", "E9", "B1"], [1, 2, 2]], "H3": [["H5", "E9", "C1"], [1, 2, 2]], "A4": [["A6", "E9", "D1"], [1, 2, 2]], "B4": [["A6", "E9", "E1"], [1, 2, 2]], "C4": [["A6", "E9", "F1"], [1, 2, 2]], "D4": [["A6", "F9", "F12"], [1, 2, 1]], "E4": [["A6", "F9", "G12"], [1, 2, 1]], "F4": [["A6", "F9", "H12"], [1, 2, 1]], "G4": [["A6", "F9", "A1"], [1, 2, 2]], "H4": [["A6", "F9", "B1"], [1, 2, 2]], "A5": [["A6", "F9", "C1"], [1, 2, 2]], "B5": [["A6", "F9", "D1"], [1, 2, 2]], "C5": [["A6", "F9", "E1"], [1, 2, 2]], "D5": [["A6", "F9", "F1"], [1, 2, 2]], "E5": [["A6", "G9", "F12"], [1, 2, 1]], "F5": [["A6", "G9", "G12"], [1, 2, 1]], "G5": [["A6", "G9", "H12"], [1, 2, 1]], "H5": [["B6", "G9", "A1"], [1, 2, 2]], "A6": [["B6", "G9", "B1"], [1, 2, 2]], "B6": [["B6", "G9", "C1"], [1, 2, 2]], "C6": [["B6", "G9", "D1"], [1, 2, 2]], "D6": [["B6", "G9", "E1"], [1, 2, 2]], "E6": [["B6", "G9", "F1"], [1, 2, 2]], "F6": [["B6", "H9", "B2"], [1, 2, 2]], "G6": [["B6", "H9", "C2"], [1, 2, 2]], "H6": [["B6", "H9", "D2"], [1, 2, 2]], "A7": [["B6", "H9", "E2"], [1, 2, 2]], "B7": [["B6", "H9", "F2"], [1, 2, 2]], "C7": [["B6", "H9", "G2"], [1, 2, 2]], "D7": [["B6", "H9", "H2"], [1, 2, 2]], "E7": [["B6", "H9", "A3"], [1, 2, 2]], "F7": [["B6", "H9", "B3"], [1, 2, 2]], "G7": [["C6", "A10", "B2"], [1, 2, 2]], "H7": [["C6", "A10", "C2"], [1, 2, 2]], "A8": [["C6", "A10", "D2"], [1, 2, 2]], "B8": [["C6", "A10", "E2"], [1, 2, 2]], "C8": [["C6", "A10", "F2"], [1, 2, 2]], "D8": [["C6", "A10", "G2"], [1, 2, 2]], "E8": [["C6", "A10", "H2"], [1, 2, 2]], "F8": [["C6", "A10", "A3"], [1, 2, 2]], "G8": [["C6", "A10", "B3"], [1, 2, 2]], "H8": [["C6", "B10", "B2"], [1, 2, 2]], "A9": [["C6", "B10", "C2"], [1, 2, 2]], "B9": [["C6", "B10", "D2"], [1, 2, 2]], "C9": [["C6", "B10", "E2"], [1, 2, 2]], "D9": [["C6", "B10", "F2"], [1, 2, 2]], "E9": [["C6", "B10", "G2"], [1, 2, 2]], "F9": [["D6", "B10", "H2"], [1, 2, 2]], "G9": [["D6", "B10", "A3"], [1, 2, 2]], "H9": [["D6", "B10", "B3"], [1, 2, 2]]}
tiprack_num=3


def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="tiprack-10ul"):
    """Implements final assembly reactions using an opentrons OT-2.

    Args:
    final_assembly_dict (dict): Dictionary with keys and values corresponding to destination and associated linker-ligated part wells, respectively.
    tiprack_num (int): Number of tipracks required during run.

    """

    # Constants
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5', '8', '11']
    PIPETTE_MOUNT = 'right'
    MAG_PLATE_TYPE = '4ti-0960_FrameStar'
    MAG_PLATE_POSITION = '1'
    TUBE_RACK_TYPE = 'tube-rack_E1415-1500'
    TUBE_RACK_POSITION = '7'
    DESTINATION_PLATE_TYPE = 'aluminium-block_4ti-0960_FrameStar'
    TEMPDECK_SLOT = '4'
    TEMP = 20
    TOTAL_VOL = 15
    PART_VOL = 1.5
    MIX_SETTINGS = (1, 3)

    # Errors
    sample_number = len(final_assembly_dict.keys())
    if sample_number > 96:
        raise ValueError('Final assembly nummber cannot exceed 96.')

    # Tips and pipette
    slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
    tipracks = [labware.load(tiprack_type, slot)
                for slot in slots]
    pipette = instruments.P10_Single(mount=PIPETTE_MOUNT, tip_racks=tipracks)

    # Define Labware and set temperature
    magbead_plate = labware.load(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
    tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
    tempdeck = modules.load('tempdeck', TEMPDECK_SLOT)
    destination_plate = labware.load(
        DESTINATION_PLATE_TYPE, TEMPDECK_SLOT, share=True)
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()

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
        destination_wells = np.array(
            [key for key, value in list(final_assembly_dict.items())])
        destination_wells = list(destination_wells[destination_inds])
        pipette.pick_up_tip()
        pipette.transfer(TOTAL_VOL - x * PART_VOL, tube_rack.wells(master_mix_well),
                         destination_plate.wells(destination_wells),
                         new_tip='never')
        pipette.drop_tip()

    # Part transfers
    for key, values in list(final_assembly_dict.items()):
        pipette.transfer(PART_VOL, magbead_plate.wells(values),
                         destination_plate.wells(key), mix_after=MIX_SETTINGS,
                         new_tip='always')

    tempdeck.deactivate()


final_assembly(final_assembly_dict=final_assembly_dict,
               tiprack_num=tiprack_num)

for c in robot.commands():
    print(c)

from opentrons import labware, instruments, modules, robot
import numpy as np


final_assembly_dict={"A1": [["A1", "H7", "E2"], [1, 2, 2]], "B1": [["A1", "H7", "F2"], [1, 2, 2]], "C1": [["A1", "H7", "G2"], [1, 2, 2]], "D1": [["A1", "H7", "H2"], [1, 2, 2]], "E1": [["A1", "H7", "A3"], [1, 2, 2]], "F1": [["A1", "H7", "B3"], [1, 2, 2]], "G1": [["A1", "A8", "F6"], [1, 2, 1]], "H1": [["A1", "A8", "G6"], [1, 2, 1]], "A2": [["A1", "A8", "H6"], [1, 2, 1]], "B2": [["A1", "A8", "A7"], [1, 2, 1]], "C2": [["A1", "A8", "B7"], [1, 2, 1]], "D2": [["A1", "A8", "C7"], [1, 2, 1]], "E2": [["A1", "A8", "D7"], [1, 2, 1]], "F2": [["A1", "A8", "E7"], [1, 2, 1]], "G2": [["A1", "A8", "F7"], [1, 2, 1]], "H2": [["B1", "B8", "F6"], [1, 2, 1]], "A3": [["B1", "B8", "G6"], [1, 2, 1]], "B3": [["B1", "B8", "H6"], [1, 2, 1]], "C3": [["B1", "B8", "A7"], [1, 2, 1]], "D3": [["B1", "B8", "B7"], [1, 2, 1]], "E3": [["B1", "B8", "C7"], [1, 2, 1]], "F3": [["B1", "B8", "D7"], [1, 2, 1]], "G3": [["B1", "B8", "E7"], [1, 2, 1]], "H3": [["B1", "B8", "F7"], [1, 2, 1]], "A4": [["B1", "C8", "F6"], [1, 2, 1]], "B4": [["B1", "C8", "G6"], [1, 2, 1]], "C4": [["B1", "C8", "H6"], [1, 2, 1]], "D4": [["B1", "C8", "A7"], [1, 2, 1]], "E4": [["B1", "C8", "B7"], [1, 2, 1]], "F4": [["B1", "C8", "C7"], [1, 2, 1]], "G4": [["C1", "C8", "D7"], [1, 2, 1]], "H4": [["C1", "C8", "E7"], [1, 2, 1]], "A5": [["C1", "C8", "F7"], [1, 2, 1]], "B5": [["C1", "D8", "B8"], [1, 2, 1]], "C5": [["C1", "D8", "C8"], [1, 2, 1]], "D5": [["C1", "D8", "D8"], [1, 2, 1]], "E5": [["C1", "D8", "E8"], [1, 2, 1]], "F5": [["C1", "D8", "F8"], [1, 2, 1]], "G5": [["C1", "D8", "G8"], [1, 2, 1]], "H5": [["C1", "D8", "H8"], [1, 2, 1]], "A6": [["C1", "D8", "A9"], [1, 2, 1]], "B6": [["C1", "D8", "B9"], [1, 2, 1]], "C6": [["C1", "E8", "B8"], [1, 2, 1]], "D6": [["C1", "E8", "C8"], [1, 2, 1]], "E6": [["C1", "E8", "D8"], [1, 2, 1]], "F6": [["D1", "E8", "E8"], [1, 2, 1]], "G6": [["D1", "E8", "F8"], [1, 2, 1]], "H6": [["D1", "E8", "G8"], [1, 2, 1]], "A7": [["D1", "E8", "H8"], [1, 2, 1]], "B7": [["D1", "E8", "A9"], [1, 2, 1]], "C7": [["D1", "E8", "B9"], [1, 2, 1]], "D7": [["D1", "F8", "B8"], [1, 2, 1]], "E7": [["D1", "F8", "C8"], [1, 2, 1]], "F7": [["D1", "F8", "D8"], [1, 2, 1]], "G7": [["D1", "F8", "E8"], [1, 2, 1]], "H7": [["D1", "F8", "F8"], [1, 2, 1]], "A8": [["D1", "F8", "G8"], [1, 2, 1]], "B8": [["D1", "F8", "H8"], [1, 2, 1]], "C8": [["D1", "F8", "A9"], [1, 2, 1]], "D8": [["D1", "F8", "B9"], [1, 2, 1]], "E8": [["E1", "G8", "F9"], [1, 2, 1]], "F8": [["E1", "G8", "G9"], [1, 2, 1]], "G8": [["E1", "G8", "H9"], [1, 2, 1]], "H8": [["E1", "G8", "A10"], [1, 2, 1]], "A9": [["E1", "G8", "B10"], [1, 2, 1]], "B9": [["E1", "G8", "C10"], [1, 2, 1]], "C9": [["E1", "G8", "D10"], [1, 2, 1]], "D9": [["E1", "G8", "E10"], [1, 2, 1]], "E9": [["E1", "G8", "F10"], [1, 2, 1]], "F9": [["E1", "H8", "F9"], [1, 2, 1]], "G9": [["E1", "H8", "G9"], [1, 2, 1]], "H9": [["E1", "H8", "H9"], [1, 2, 1]], "A10": [["E1", "H8", "A10"], [1, 2, 1]], "B10": [["E1", "H8", "B10"], [1, 2, 1]], "C10": [["E1", "H8", "C10"], [1, 2, 1]], "D10": [["F1", "H8", "D10"], [1, 2, 1]], "E10": [["F1", "H8", "E10"], [1, 2, 1]], "F10": [["F1", "H8", "F10"], [1, 2, 1]], "G10": [["F1", "A9", "F9"], [1, 2, 1]], "H10": [["F1", "A9", "G9"], [1, 2, 1]], "A11": [["F1", "A9", "H9"], [1, 2, 1]], "B11": [["F1", "A9", "A10"], [1, 2, 1]], "C11": [["F1", "A9", "B10"], [1, 2, 1]], "D11": [["F1", "A9", "C10"], [1, 2, 1]], "E11": [["F1", "A9", "D10"], [1, 2, 1]], "F11": [["F1", "A9", "E10"], [1, 2, 1]], "G11": [["F1", "A9", "F10"], [1, 2, 1]], "H11": [["F1", "B9", "B11"], [1, 2, 1]], "A12": [["F1", "B9", "C11"], [1, 2, 1]], "B12": [["F1", "B9", "D11"], [1, 2, 1]], "C12": [["G1", "B9", "E11"], [1, 2, 1]], "D12": [["G1", "B9", "F11"], [1, 2, 1]], "E12": [["G1", "B9", "G11"], [1, 2, 1]], "F12": [["G1", "B9", "H11"], [1, 2, 1]], "G12": [["G1", "B9", "A12"], [1, 2, 1]], "H12": [["G1", "B9", "B12"], [1, 2, 1]]}
tiprack_num=4


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

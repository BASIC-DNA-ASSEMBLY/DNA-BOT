from opentrons import labware, instruments, modules, robot
import numpy as np


final_assembly_dict={"A1": [["B4", "E6", "D10"], [1, 2, 1]], "B1": [["B4", "E6", "E10"], [1, 2, 1]], "C1": [["B4", "E6", "F10"], [1, 2, 1]], "D1": [["B4", "F6", "F9"], [1, 2, 1]], "E1": [["B4", "F6", "G9"], [1, 2, 1]], "F1": [["B4", "F6", "H9"], [1, 2, 1]], "G1": [["C4", "F6", "A10"], [1, 2, 1]], "H1": [["C4", "F6", "B10"], [1, 2, 1]], "A2": [["C4", "F6", "C10"], [1, 2, 1]], "B2": [["C4", "F6", "D10"], [1, 2, 1]], "C2": [["C4", "F6", "E10"], [1, 2, 1]], "D2": [["C4", "F6", "F10"], [1, 2, 1]], "E2": [["C4", "G6", "F9"], [1, 2, 1]], "F2": [["C4", "G6", "G9"], [1, 2, 1]], "G2": [["C4", "G6", "H9"], [1, 2, 1]], "H2": [["C4", "G6", "A10"], [1, 2, 1]], "A3": [["C4", "G6", "B10"], [1, 2, 1]], "B3": [["C4", "G6", "C10"], [1, 2, 1]], "C3": [["C4", "G6", "D10"], [1, 2, 1]], "D3": [["C4", "G6", "E10"], [1, 2, 1]], "E3": [["C4", "G6", "F10"], [1, 2, 1]], "F3": [["D4", "H6", "B11"], [1, 2, 1]], "G3": [["D4", "H6", "C11"], [1, 2, 1]], "H3": [["D4", "H6", "D11"], [1, 2, 1]], "A4": [["D4", "H6", "E11"], [1, 2, 1]], "B4": [["D4", "H6", "F11"], [1, 2, 1]], "C4": [["D4", "H6", "G11"], [1, 2, 1]], "D4": [["D4", "H6", "H11"], [1, 2, 1]], "E4": [["D4", "H6", "A12"], [1, 2, 1]], "F4": [["D4", "H6", "B12"], [1, 2, 1]], "G4": [["D4", "A7", "B11"], [1, 2, 1]], "H4": [["D4", "A7", "C11"], [1, 2, 1]], "A5": [["D4", "A7", "D11"], [1, 2, 1]], "B5": [["D4", "A7", "E11"], [1, 2, 1]], "C5": [["D4", "A7", "F11"], [1, 2, 1]], "D5": [["D4", "A7", "G11"], [1, 2, 1]], "E5": [["E4", "A7", "H11"], [1, 2, 1]], "F5": [["E4", "A7", "A12"], [1, 2, 1]], "G5": [["E4", "A7", "B12"], [1, 2, 1]], "H5": [["E4", "B7", "B11"], [1, 2, 1]], "A6": [["E4", "B7", "C11"], [1, 2, 1]], "B6": [["E4", "B7", "D11"], [1, 2, 1]], "C6": [["E4", "B7", "E11"], [1, 2, 1]], "D6": [["E4", "B7", "F11"], [1, 2, 1]], "E6": [["E4", "B7", "G11"], [1, 2, 1]], "F6": [["E4", "B7", "H11"], [1, 2, 1]], "G6": [["E4", "B7", "A12"], [1, 2, 1]], "H6": [["E4", "B7", "B12"], [1, 2, 1]], "A7": [["E4", "C7", "F12"], [1, 2, 1]], "B7": [["E4", "C7", "G12"], [1, 2, 1]], "C7": [["E4", "C7", "H12"], [1, 2, 1]], "D7": [["F4", "C7", "A1"], [1, 2, 2]], "E7": [["F4", "C7", "B1"], [1, 2, 2]], "F7": [["F4", "C7", "C1"], [1, 2, 2]], "G7": [["F4", "C7", "D1"], [1, 2, 2]], "H7": [["F4", "C7", "E1"], [1, 2, 2]], "A8": [["F4", "C7", "F1"], [1, 2, 2]], "B8": [["F4", "D7", "F12"], [1, 2, 1]], "C8": [["F4", "D7", "G12"], [1, 2, 1]], "D8": [["F4", "D7", "H12"], [1, 2, 1]], "E8": [["F4", "D7", "A1"], [1, 2, 2]], "F8": [["F4", "D7", "B1"], [1, 2, 2]], "G8": [["F4", "D7", "C1"], [1, 2, 2]], "H8": [["F4", "D7", "D1"], [1, 2, 2]], "A9": [["F4", "D7", "E1"], [1, 2, 2]], "B9": [["F4", "D7", "F1"], [1, 2, 2]], "C9": [["G4", "E7", "F12"], [1, 2, 1]], "D9": [["G4", "E7", "G12"], [1, 2, 1]], "E9": [["G4", "E7", "H12"], [1, 2, 1]], "F9": [["G4", "E7", "A1"], [1, 2, 2]], "G9": [["G4", "E7", "B1"], [1, 2, 2]], "H9": [["G4", "E7", "C1"], [1, 2, 2]], "A10": [["G4", "E7", "D1"], [1, 2, 2]], "B10": [["G4", "E7", "E1"], [1, 2, 2]], "C10": [["G4", "E7", "F1"], [1, 2, 2]], "D10": [["G4", "F7", "B2"], [1, 2, 2]], "E10": [["G4", "F7", "C2"], [1, 2, 2]], "F10": [["G4", "F7", "D2"], [1, 2, 2]], "G10": [["G4", "F7", "E2"], [1, 2, 2]], "H10": [["G4", "F7", "F2"], [1, 2, 2]], "A11": [["G4", "F7", "G2"], [1, 2, 2]], "B11": [["H4", "F7", "H2"], [1, 2, 2]], "C11": [["H4", "F7", "A3"], [1, 2, 2]], "D11": [["H4", "F7", "B3"], [1, 2, 2]], "E11": [["H4", "G7", "B2"], [1, 2, 2]], "F11": [["H4", "G7", "C2"], [1, 2, 2]], "G11": [["H4", "G7", "D2"], [1, 2, 2]], "H11": [["H4", "G7", "E2"], [1, 2, 2]], "A12": [["H4", "G7", "F2"], [1, 2, 2]], "B12": [["H4", "G7", "G2"], [1, 2, 2]], "C12": [["H4", "G7", "H2"], [1, 2, 2]], "D12": [["H4", "G7", "A3"], [1, 2, 2]], "E12": [["H4", "G7", "B3"], [1, 2, 2]], "F12": [["H4", "H7", "B2"], [1, 2, 2]], "G12": [["H4", "H7", "C2"], [1, 2, 2]], "H12": [["H4", "H7", "D2"], [1, 2, 2]]}
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

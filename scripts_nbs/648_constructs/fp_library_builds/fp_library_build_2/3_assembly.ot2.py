from opentrons import labware, instruments, modules, robot
import numpy as np


final_assembly_dict={"A1": ["A7", "G7", "A8", "B8", "D8"], "B1": ["A7", "G7", "A8", "B8", "F8"], "C1": ["A7", "G7", "H8", "A9", "C9"], "D1": ["A7", "G7", "H8", "A9", "D8"], "E1": ["A7", "G7", "H8", "A9", "F8"], "F1": ["A7", "G7", "H8", "E9", "C9"], "G1": ["A7", "G7", "H8", "E9", "D8"], "H1": ["A7", "G7", "H8", "E9", "F8"], "A2": ["A7", "G7", "H8", "B8", "C9"], "B2": ["A7", "G7", "H8", "B8", "D8"], "C2": ["A7", "G7", "H8", "B8", "F8"], "D2": ["A7", "G7", "G9", "A9", "C9"], "E2": ["A7", "G7", "G9", "A9", "D8"], "F2": ["A7", "G7", "G9", "A9", "F8"], "G2": ["A7", "G7", "G9", "E9", "C9"], "H2": ["B7", "H7", "G9", "E9", "D8"], "A3": ["B7", "H7", "G9", "E9", "F8"], "B3": ["B7", "H7", "G9", "B8", "C9"], "C3": ["B7", "H7", "G9", "B8", "D8"], "D3": ["B7", "H7", "G9", "B8", "F8"], "E3": ["B7", "H9", "A9", "D10", "F10"], "F3": ["B7", "H9", "A9", "D10", "G10"], "G3": ["B7", "H9", "A9", "D10", "H10"], "H3": ["B7", "H9", "A9", "A11", "F10"], "A4": ["B7", "H9", "A9", "A11", "G10"], "B4": ["B7", "H9", "A9", "A11", "H10"], "C4": ["B7", "H9", "A9", "B11", "F10"], "D4": ["B7", "H9", "A9", "B11", "G10"], "E4": ["B7", "H9", "A9", "B11", "H10"], "F4": ["B7", "H9", "E9", "D10", "F10"], "G4": ["C7", "H9", "E9", "D10", "G10"], "H4": ["C7", "H9", "E9", "D10", "H10"], "A5": ["C7", "H9", "E9", "A11", "F10"], "B5": ["C7", "H9", "E9", "A11", "G10"], "C5": ["C7", "H9", "E9", "A11", "H10"], "D5": ["C7", "A10", "E9", "B11", "F10"], "E5": ["C7", "A10", "E9", "B11", "G10"], "F5": ["C7", "A10", "E9", "B11", "H10"], "G5": ["C7", "A10", "B8", "D10", "F10"], "H5": ["C7", "A10", "B8", "D10", "G10"], "A6": ["C7", "A10", "B8", "D10", "H10"], "B6": ["C7", "A10", "B8", "A11", "F10"], "C6": ["C7", "A10", "B8", "A11", "G10"], "D6": ["C7", "A10", "B8", "A11", "H10"], "E6": ["C7", "A10", "B8", "B11", "F10"], "F6": ["D7", "A10", "C8", "B11", "G10"], "G6": ["D7", "A10", "C8", "B11", "H10"], "H6": ["D7", "A10", "C11", "D11", "C9"], "A7": ["D7", "A10", "C11", "D11", "D8"], "B7": ["D7", "A10", "C11", "D11", "F8"], "C7": ["D7", "B10", "C11", "E11", "C9"], "D7": ["D7", "B10", "C11", "E11", "D8"], "E7": ["D7", "B10", "C11", "E11", "F8"], "F7": ["D7", "B10", "C11", "F11", "C9"], "G7": ["D7", "B10", "C11", "F11", "D8"], "H7": ["D7", "B10", "C11", "F11", "F8"], "A8": ["D7", "B10", "G11", "D11", "C9"], "B8": ["D7", "B10", "G11", "D11", "D8"], "C8": ["D7", "B10", "G11", "D11", "F8"], "D8": ["D7", "B10", "G11", "E11", "C9"], "E8": ["E7", "B10", "G11", "E11", "D8"], "F8": ["E7", "B10", "G11", "E11", "F8"], "G8": ["E7", "B10", "G11", "F11", "C9"], "H8": ["E7", "B10", "G11", "F11", "D8"], "A9": ["E7", "B10", "G11", "F11", "F8"], "B9": ["E7", "C10", "H11", "D11", "C9"], "C9": ["E7", "C10", "H11", "D11", "D8"], "D9": ["E7", "C10", "H11", "D11", "F8"], "E9": ["E7", "C10", "H11", "E11", "C9"], "F9": ["E7", "C10", "H11", "E11", "D8"], "G9": ["E7", "C10", "H11", "E11", "F8"], "H9": ["E7", "C10", "H11", "F11", "C9"], "A10": ["E7", "C10", "H11", "F11", "E8"], "B10": ["E7", "C10", "H11", "F11", "G8"], "C10": ["E7", "A12", "D10", "A8", "B12"], "D10": ["F7", "A12", "D10", "A8", "C12"], "E10": ["F7", "A12", "D10", "A8", "D12"], "F10": ["F7", "A12", "D10", "H8", "B12"], "G10": ["F7", "A12", "D10", "H8", "C12"], "H10": ["F7", "A12", "D10", "H8", "D12"], "A11": ["F7", "A12", "E10", "G9", "B12"], "B11": ["F7", "A12", "E10", "G9", "C12"], "C11": ["F7", "A12", "E10", "G9", "D12"], "D11": ["F7", "A12", "A11", "A8", "B12"], "E11": ["F7", "A12", "A11", "A8", "C12"], "F11": ["F7", "A12", "A11", "A8", "D12"], "G11": ["F7", "A12", "A11", "H8", "B12"], "H11": ["F7", "A12", "A11", "H8", "C12"]}
tiprack_num=5


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

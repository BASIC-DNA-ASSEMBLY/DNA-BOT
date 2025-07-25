from opentrons import protocol_api
import numpy as np
import json

# metadata
metadata = {
'protocolName': 'My Protocol',
'description': 'Simple protocol to get started using OT2',
'apiLevel': '2.8'
}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["B1", "B2", "C2"], [2, 2, 2]],
    "B1": [["B1", "D2", "B7"], [2, 2, 2]],
    "C1": [["C1", "H4", "G2"], [2, 2, 2]],
    "D1": [["C1", "H2", "A3"], [2, 2, 2]],
    "E1": [["C1", "C7", "C5"], [2, 2, 2]],
    "F1": [["C1", "D3", "H5"], [2, 2, 2]],
    "G1": [["C1", "D7", "H5"], [2, 2, 2]],
    "H1": [["C1", "H3", "H5"], [2, 2, 2]],
    "A2": [["C1", "B4", "C4"], [2, 2, 2]],
    "B2": [["C1", "D4", "E4"], [2, 2, 2]],
    "C2": [["C1", "F4", "E7"], [2, 2, 2]],
    "D2": [["C1", "H4", "A5"], [2, 2, 2]],
    "E2": [["C1", "H2", "A3"], [2, 2, 2]],
    "F2": [["C1", "B5", "C5"], [2, 2, 2]],
    "G2": [["C1", "F7", "G7"], [2, 2, 2]],
    "H2": [["D1", "F4", "G4"], [2, 2, 2]],
    "A3": [["D1", "H4", "A5"], [2, 2, 2]],
    "B3": [["D1", "H7", "A3"], [2, 2, 2]],
    "C3": [["D1", "B5", "C5"], [2, 2, 2]],
    "D3": [["D1", "A8", "E5"], [2, 2, 2]],
    "E3": [["D1", "D3", "F5"], [2, 2, 2]],
    "F3": [["D1", "D7", "H5"], [2, 2, 2]],
    "G3": [["D1", "H3", "F5"], [2, 2, 2]],
    "H3": [["D1", "B4", "E5"], [2, 2, 2]]
}
tiprack_num = 1

# protocol run function. the part after the colon lets your editor know


# test dict can be used for simulation
#final_assembly_dict={ "A1": ['A7', 'B7', 'C7', 'F7'], "B1": ['A7', 'B7', 'D7', 'G7'], "C1": ['A7', 'B7', 'E7', 'H7']}
#tiprack_num=1
def run(protocol: protocol_api.ProtocolContext):
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_tiprack_20ul"):
            # Constants, we update all the labware name in version 2
            #Tiprack
            CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5', '8', '11']
            PIPETTE_MOUNT = 'right'
            #Plate of sample after  purification
            MAG_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            MAG_PLATE_POSITION = '1'
            #Tuberack
            TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            TUBE_RACK_POSITION = '7'
            #Destination plate
            DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            #Temperature control plate
            TEMPDECK_SLOT = '4'
            TEMP = 20
            TOTAL_VOL = 15
            PART_VOL = 1.5
            MIX_SETTINGS = (1, 3)
            tiprack_num=tiprack_num+1
            # Errors
            sample_number = len(final_assembly_dict.keys())
            if sample_number > 96:
                raise ValueError('Final assembly nummber cannot exceed 96.')

            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            pipette = protocol.load_instrument('p20_single_gen2', PIPETTE_MOUNT, tip_racks=tipracks)


            # Define Labware and set temperature
            magbead_plate = protocol.load_labware(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
            tempdeck = protocol.load_module('tempdeck', TEMPDECK_SLOT)
            destination_plate = tempdeck.load_labware(
            DESTINATION_PLATE_TYPE, TEMPDECK_SLOT)
            tempdeck.set_temperature(TEMP)

             # Master mix transfers
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))
            master_mix_well_letters = ['A', 'B', 'C', 'D']
            for x in unique_assemblies_lens:
                master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
                destination_inds = [i for i, lens in enumerate(final_assembly_lens) if lens == x]
                destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
                destination_wells = list(destination_wells[destination_inds])
                for destination_well in destination_wells:# make tube_rack_wells and destination_plate.wells in the same type
                    pipette.pick_up_tip()
                    pipette.transfer(TOTAL_VOL - x * PART_VOL, tube_rack.wells(master_mix_well),
                                     destination_plate.wells(destination_well), new_tip='never')#transfer water and buffer in the pipette

                    pipette.drop_tip()

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for value in values:# magbead_plate.wells and destination_plate.wells in the same type
                    pipette.transfer(PART_VOL, magbead_plate.wells(value),
                                     destination_plate.wells(key), mix_after=MIX_SETTINGS,
                                     new_tip='always')#transfer parts in one tube

            tempdeck.deactivate() #stop increasing the temperature

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

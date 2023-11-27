from opentrons import protocol_api
import numpy as np
# metadata
metadata = {
'protocolName': 'My Protocol',
'description': 'Simple protocol to get started using OT2',
'apiLevel': '2.8'
}

# protocol run function. the part after the colon lets your editor know


# test dict can be used for simulation
#final_assembly_dict={ "A1": ['A7', 'B7', 'C7', 'F7'], "B1": ['A7', 'B7', 'D7', 'G7'], "C1": ['A7', 'B7', 'E7', 'H7']}
#tiprack_num=1
final_assembly_dict={"A1": ["A7", "G7", "B8", "C8", "D8"], "B1": ["A7", "G7", "B8", "E8", "F8"], "C1": ["A7", "G7", "B8", "E8", "G8"], "D1": ["A7", "G7", "B8", "E8", "D8"], "E1": ["A7", "G7", "B8", "H8", "F8"], "F1": ["A7", "G7", "B8", "H8", "G8"], "G1": ["A7", "G7", "B8", "H8", "D8"], "H1": ["A7", "G7", "A9", "C9", "D9"], "A2": ["A7", "G7", "A9", "C9", "E9"], "B2": ["A7", "G7", "A9", "C9", "F9"], "C2": ["A7", "G7", "A9", "G9", "D9"], "D2": ["A7", "G7", "A9", "G9", "E9"], "E2": ["A7", "G7", "A9", "G9", "F9"], "F2": ["A7", "G7", "A9", "H9", "D9"], "G2": ["A7", "G7", "A9", "H9", "E9"], "H2": ["B7", "H7", "A9", "H9", "F9"], "A3": ["B7", "H7", "A10", "C9", "D9"], "B3": ["B7", "H7", "A10", "C9", "E9"], "C3": ["B7", "H7", "A10", "C9", "F9"], "D3": ["B7", "H7", "A10", "G9", "D9"], "E3": ["B7", "H7", "A10", "G9", "E9"], "F3": ["B7", "H7", "A10", "G9", "F9"], "G3": ["B7", "H7", "A10", "H9", "D9"], "H3": ["B7", "H7", "A10", "H9", "E9"], "A4": ["B7", "H7", "A10", "H9", "F9"], "B4": ["B7", "H7", "C10", "C9", "D9"], "C4": ["B7", "H7", "C10", "C9", "E9"], "D4": ["B7", "H7", "C10", "C9", "F9"], "E4": ["B7", "H7", "C10", "G9", "D9"], "F4": ["B7", "H7", "C10", "G9", "E9"], "G4": ["C7", "A8", "C10", "G9", "F9"], "H4": ["C7", "A8", "C10", "H9", "D9"], "A5": ["C7", "A8", "C10", "H9", "E9"], "B5": ["C7", "A8", "C10", "H9", "F9"], "C5": ["C7", "D10", "C8", "G10", "H10"], "D5": ["C7", "D10", "C8", "G10", "A11"], "E5": ["C7", "D10", "C8", "G10", "B11"], "F5": ["C7", "D10", "C8", "C11", "H10"], "G5": ["C7", "D10", "C8", "C11", "A11"], "H5": ["C7", "D10", "C8", "C11", "B11"], "A6": ["C7", "D10", "C8", "D11", "H10"], "B6": ["C7", "D10", "C8", "D11", "A11"], "C6": ["C7", "D10", "C8", "D11", "B11"], "D6": ["C7", "D10", "E8", "G10", "H10"], "E6": ["C7", "D10", "E8", "G10", "A11"], "F6": ["D7", "D10", "E8", "G10", "B11"], "G6": ["D7", "D10", "E8", "C11", "H10"], "H6": ["D7", "D10", "E8", "C11", "A11"], "A7": ["D7", "D10", "E8", "C11", "B11"], "B7": ["D7", "E10", "E8", "D11", "H10"], "C7": ["D7", "E10", "E8", "D11", "A11"], "D7": ["D7", "E10", "E8", "D11", "B11"], "E7": ["D7", "E10", "H8", "G10", "H10"], "F7": ["D7", "E10", "H8", "G10", "A11"], "G7": ["D7", "E10", "H8", "G10", "B11"], "H7": ["D7", "E10", "H8", "C11", "H10"], "A8": ["D7", "E10", "H8", "C11", "A11"], "B8": ["D7", "E10", "H8", "C11", "B11"], "C8": ["D7", "E10", "H8", "D11", "H10"], "D8": ["D7", "E10", "H8", "D11", "A11"], "E8": ["E7", "E10", "H8", "D11", "B11"], "F8": ["E7", "E10", "E11", "A9", "F8"], "G8": ["E7", "E10", "E11", "A9", "G8"], "H8": ["E7", "E10", "E11", "A9", "D8"], "A9": ["E7", "F10", "E11", "A10", "F8"], "B9": ["E7", "F10", "E11", "A10", "G8"], "C9": ["E7", "F10", "E11", "A10", "D8"], "D9": ["E7", "F10", "E11", "C10", "F8"], "E9": ["E7", "F10", "E11", "C10", "G8"], "F9": ["E7", "F10", "E11", "C10", "D8"], "G9": ["E7", "F10", "F11", "A9", "F8"], "H9": ["E7", "F10", "F11", "A9", "G8"], "A10": ["E7", "F10", "F11", "A9", "D8"], "B10": ["E7", "F10", "F11", "A10", "F8"], "C10": ["E7", "F10", "F11", "A10", "G8"], "D10": ["F7", "F10", "F11", "A10", "D8"]}
tiprack_num=4

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

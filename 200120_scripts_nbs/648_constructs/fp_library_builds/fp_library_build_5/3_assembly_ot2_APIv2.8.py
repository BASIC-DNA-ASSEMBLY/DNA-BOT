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
final_assembly_dict={"A1": ["A7", "G7", "A8", "C8", "D8"], "B1": ["A7", "G7", "A8", "C8", "E8"], "C1": ["A7", "G7", "A8", "F8", "G8"], "D1": ["A7", "G7", "A8", "F8", "D8"], "E1": ["A7", "G7", "A8", "F8", "E8"], "F1": ["A7", "G7", "A8", "H8", "G8"], "G1": ["A7", "G7", "A8", "H8", "D8"], "H1": ["A7", "G7", "A8", "H8", "E8"], "A2": ["A7", "G7", "A9", "C8", "G8"], "B2": ["A7", "G7", "A9", "C8", "D8"], "C2": ["A7", "G7", "A9", "C8", "E8"], "D2": ["A7", "G7", "A9", "F8", "G8"], "E2": ["A7", "G7", "A9", "F8", "D8"], "F2": ["A7", "G7", "A9", "F8", "E8"], "G2": ["A7", "G7", "A9", "H8", "G8"], "H2": ["B7", "H7", "A9", "H8", "D8"], "A3": ["B7", "H7", "A9", "H8", "E8"], "B3": ["B7", "H7", "C9", "C8", "G8"], "C3": ["B7", "H7", "C9", "C8", "D8"], "D3": ["B7", "H7", "C9", "C8", "E8"], "E3": ["B7", "H7", "C9", "F8", "G8"], "F3": ["B7", "H7", "C9", "F8", "D8"], "G3": ["B7", "H7", "C9", "F8", "E8"], "H3": ["B7", "H7", "C9", "H8", "G8"], "A4": ["B7", "H7", "C9", "H8", "D8"], "B4": ["B7", "H7", "C9", "H8", "E8"], "C4": ["B7", "E9", "A10", "A8", "B10"], "D4": ["B7", "E9", "A10", "A8", "C10"], "E4": ["B7", "E9", "A10", "A8", "D10"], "F4": ["B7", "E9", "A10", "A9", "B10"], "G4": ["C7", "E9", "A10", "A9", "C10"], "H4": ["C7", "E9", "A10", "A9", "D10"], "A5": ["C7", "E9", "A10", "C9", "B10"], "B5": ["C7", "E9", "A10", "C9", "C10"], "C5": ["C7", "E9", "A10", "C9", "D10"], "D5": ["C7", "E9", "E10", "A8", "B10"], "E5": ["C7", "E9", "E10", "A8", "C10"], "F5": ["C7", "E9", "E10", "A8", "D10"], "G5": ["C7", "E9", "E10", "A9", "B10"], "H5": ["C7", "E9", "E10", "A9", "C10"], "A6": ["C7", "E9", "E10", "A9", "D10"], "B6": ["C7", "F9", "E10", "C9", "B10"], "C6": ["C7", "F9", "E10", "C9", "C10"], "D6": ["C7", "F9", "E10", "C9", "D10"], "E6": ["C7", "F9", "F10", "A8", "B10"], "F6": ["D7", "F9", "F10", "B8", "C10"], "G6": ["D7", "F9", "F10", "B8", "D10"], "H6": ["D7", "F9", "F10", "B9", "B10"], "A7": ["D7", "F9", "F10", "B9", "C10"], "B7": ["D7", "F9", "F10", "B9", "D10"], "C7": ["D7", "F9", "F10", "D9", "B10"], "D7": ["D7", "F9", "F10", "D9", "C10"], "E7": ["D7", "F9", "F10", "D9", "D10"], "F7": ["D7", "F9", "G10", "H10", "B11"], "G7": ["D7", "F9", "G10", "H10", "C11"], "H7": ["D7", "F9", "G10", "H10", "D11"], "A8": ["D7", "G9", "G10", "E11", "B11"], "B8": ["D7", "G9", "G10", "E11", "C11"], "C8": ["D7", "G9", "G10", "E11", "D11"], "D8": ["D7", "G9", "G10", "F11", "B11"], "E8": ["E7", "G9", "G10", "F11", "C11"], "F8": ["E7", "G9", "G10", "F11", "D11"], "G8": ["E7", "G9", "G11", "H10", "B11"], "H8": ["E7", "G9", "G11", "H10", "C11"], "A9": ["E7", "G9", "G11", "H10", "D11"], "B9": ["E7", "G9", "G11", "E11", "B11"], "C9": ["E7", "G9", "G11", "E11", "C11"], "D9": ["E7", "G9", "G11", "E11", "D11"], "E9": ["E7", "G9", "G11", "F11", "B11"], "F9": ["E7", "G9", "G11", "F11", "C11"], "G9": ["E7", "G9", "G11", "F11", "D11"], "H9": ["E7", "H9", "H11", "H10", "B11"], "A10": ["E7", "H9", "H11", "H10", "C11"], "B10": ["E7", "H9", "H11", "H10", "D11"], "C10": ["E7", "H9", "H11", "E11", "B11"], "D10": ["F7", "H9", "H11", "E11", "C11"], "E10": ["F7", "H9", "H11", "E11", "D11"], "F10": ["F7", "H9", "H11", "F11", "B11"], "G10": ["F7", "H9", "H11", "F11", "C11"], "H10": ["F7", "H9", "H11", "F11", "D11"], "A11": ["F7", "A12", "H10", "B12", "G8"], "B11": ["F7", "A12", "H10", "B12", "D8"], "C11": ["F7", "A12", "H10", "B12", "E8"], "D11": ["F7", "A12", "H10", "C12", "G8"], "E11": ["F7", "A12", "H10", "C12", "D8"], "F11": ["F7", "A12", "H10", "C12", "E8"], "G11": ["F7", "A12", "A11", "D12", "G8"], "H11": ["F7", "A12", "A11", "D12", "D8"]}
tiprack_num=5

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

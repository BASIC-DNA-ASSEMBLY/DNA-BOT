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
final_assembly_dict={"A73": [["E6", "A10", "H11"], [1, 2, 1]], "B73": [["E6", "A10", "A12"], [1, 2, 1]], "C73": [["E6", "A10", "B12"], [1, 2, 1]], "D73": [["E6", "A10", "C12"], [1, 2, 1]], "E73": [["E6", "A10", "D12"], [1, 2, 1]], "F73": [["E6", "A10", "E12"], [1, 2, 1]], "G73": [["E6", "A10", "F12"], [1, 2, 1]], "H73": [["E6", "A10", "G12"], [1, 2, 1]], "A74": [["E6", "A10", "H12"], [1, 2, 1]], "B74": [["F6", "B10", "H11"], [1, 2, 1]], "C74": [["F6", "B10", "A12"], [1, 2, 1]], "D74": [["F6", "B10", "B12"], [1, 2, 1]], "E74": [["F6", "B10", "C12"], [1, 2, 1]], "F74": [["F6", "B10", "D12"], [1, 2, 1]], "G74": [["F6", "B10", "E12"], [1, 2, 1]], "H74": [["F6", "B10", "F12"], [1, 2, 1]], "A75": [["F6", "B10", "G12"], [1, 2, 1]], "B75": [["F6", "B10", "H12"], [1, 2, 1]], "C75": [["F6", "C10", "D1"], [1, 2, 2]], "D75": [["F6", "C10", "E1"], [1, 2, 2]], "E75": [["F6", "C10", "F1"], [1, 2, 2]], "F75": [["F6", "C10", "G1"], [1, 2, 2]], "G75": [["G6", "C10", "H1"], [1, 2, 2]], "H75": [["G6", "C10", "A2"], [1, 2, 2]], "A76": [["G6", "C10", "B2"], [1, 2, 2]], "B76": [["G6", "C10", "C2"], [1, 2, 2]], "C76": [["G6", "C10", "D2"], [1, 2, 2]], "D76": [["G6", "D10", "D1"], [1, 2, 2]], "E76": [["G6", "D10", "E1"], [1, 2, 2]], "F76": [["G6", "D10", "F1"], [1, 2, 2]], "G76": [["G6", "D10", "G1"], [1, 2, 2]], "H76": [["G6", "D10", "H1"], [1, 2, 2]], "A77": [["G6", "D10", "A2"], [1, 2, 2]], "B77": [["G6", "D10", "B2"], [1, 2, 2]], "C77": [["G6", "D10", "C2"], [1, 2, 2]], "D77": [["H6", "D10", "D2"], [1, 2, 2]], "E77": [["H6", "E10", "D1"], [1, 2, 2]], "F77": [["H6", "E10", "E1"], [1, 2, 2]], "G77": [["H6", "E10", "F1"], [1, 2, 2]], "H77": [["H6", "E10", "G1"], [1, 2, 2]], "A78": [["H6", "E10", "H1"], [1, 2, 2]], "B78": [["H6", "E10", "A2"], [1, 2, 2]], "C78": [["H6", "E10", "B2"], [1, 2, 2]], "D78": [["H6", "E10", "C2"], [1, 2, 2]], "E78": [["H6", "E10", "D2"], [1, 2, 2]], "F78": [["H6", "F10", "H2"], [1, 2, 2]], "G78": [["H6", "F10", "A3"], [1, 2, 2]], "H78": [["H6", "F10", "B3"], [1, 2, 2]], "A79": [["A7", "F10", "C3"], [1, 2, 2]], "B79": [["A7", "F10", "D3"], [1, 2, 2]], "C79": [["A7", "F10", "E3"], [1, 2, 2]], "D79": [["A7", "F10", "F3"], [1, 2, 2]], "E79": [["A7", "F10", "G3"], [1, 2, 2]], "F79": [["A7", "F10", "H3"], [1, 2, 2]], "G79": [["A7", "G10", "H2"], [1, 2, 2]], "H79": [["A7", "G10", "A3"], [1, 2, 2]], "A80": [["A7", "G10", "B3"], [1, 2, 2]], "B80": [["A7", "G10", "C3"], [1, 2, 2]], "C80": [["A7", "G10", "D3"], [1, 2, 2]], "D80": [["A7", "G10", "E3"], [1, 2, 2]], "E80": [["A7", "G10", "F3"], [1, 2, 2]], "F80": [["B7", "G10", "G3"], [1, 2, 2]], "G80": [["B7", "G10", "H3"], [1, 2, 2]], "H80": [["B7", "H10", "H2"], [1, 2, 2]], "A81": [["B7", "H10", "A3"], [1, 2, 2]], "B81": [["B7", "H10", "B3"], [1, 2, 2]], "C81": [["B7", "H10", "C3"], [1, 2, 2]], "D81": [["B7", "H10", "D3"], [1, 2, 2]], "E81": [["B7", "H10", "E3"], [1, 2, 2]], "F81": [["B7", "H10", "F3"], [1, 2, 2]], "G81": [["B7", "H10", "G3"], [1, 2, 2]], "H81": [["B7", "H10", "H3"], [1, 2, 2]]}
tiprack_num=3

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

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
final_assembly_dict={"F61": [["F5", "F8", "H3"], [1, 2, 2]], "G61": [["F5", "G8", "D7"], [1, 2, 1]], "H61": [["F5", "G8", "E7"], [1, 2, 1]], "A62": [["F5", "G8", "F7"], [1, 2, 1]], "B62": [["F5", "G8", "G7"], [1, 2, 1]], "C62": [["F5", "G8", "H7"], [1, 2, 1]], "D62": [["F5", "G8", "A8"], [1, 2, 1]], "E62": [["F5", "G8", "B8"], [1, 2, 1]], "F62": [["F5", "G8", "C8"], [1, 2, 1]], "G62": [["G5", "G8", "D8"], [1, 2, 1]], "H62": [["G5", "H8", "D7"], [1, 2, 1]], "A63": [["G5", "H8", "E7"], [1, 2, 1]], "B63": [["G5", "H8", "F7"], [1, 2, 1]], "C63": [["G5", "H8", "G7"], [1, 2, 1]], "D63": [["G5", "H8", "H7"], [1, 2, 1]], "E63": [["G5", "H8", "A8"], [1, 2, 1]], "F63": [["G5", "H8", "B8"], [1, 2, 1]], "G63": [["G5", "H8", "C8"], [1, 2, 1]], "H63": [["G5", "H8", "D8"], [1, 2, 1]], "A64": [["G5", "A9", "D7"], [1, 2, 1]], "B64": [["G5", "A9", "E7"], [1, 2, 1]], "C64": [["G5", "A9", "F7"], [1, 2, 1]], "D64": [["H5", "A9", "G7"], [1, 2, 1]], "E64": [["H5", "A9", "H7"], [1, 2, 1]], "F64": [["H5", "A9", "A8"], [1, 2, 1]], "G64": [["H5", "A9", "B8"], [1, 2, 1]], "H64": [["H5", "A9", "C8"], [1, 2, 1]], "A65": [["H5", "A9", "D8"], [1, 2, 1]], "B65": [["H5", "B9", "H8"], [1, 2, 1]], "C65": [["H5", "B9", "A9"], [1, 2, 1]], "D65": [["H5", "B9", "B9"], [1, 2, 1]], "E65": [["H5", "B9", "C9"], [1, 2, 1]], "F65": [["H5", "B9", "D9"], [1, 2, 1]], "G65": [["H5", "B9", "E9"], [1, 2, 1]], "H65": [["H5", "B9", "F9"], [1, 2, 1]], "A66": [["A6", "B9", "G9"], [1, 2, 1]], "B66": [["A6", "B9", "H9"], [1, 2, 1]], "C66": [["A6", "C9", "H8"], [1, 2, 1]], "D66": [["A6", "C9", "A9"], [1, 2, 1]], "E66": [["A6", "C9", "B9"], [1, 2, 1]], "F66": [["A6", "C9", "C9"], [1, 2, 1]], "G66": [["A6", "C9", "D9"], [1, 2, 1]], "H66": [["A6", "C9", "E9"], [1, 2, 1]], "A67": [["A6", "C9", "F9"], [1, 2, 1]], "B67": [["A6", "C9", "G9"], [1, 2, 1]], "C67": [["A6", "C9", "H9"], [1, 2, 1]], "D67": [["A6", "D9", "H8"], [1, 2, 1]], "E67": [["A6", "D9", "A9"], [1, 2, 1]], "F67": [["B6", "D9", "B9"], [1, 2, 1]], "G67": [["B6", "D9", "C9"], [1, 2, 1]], "H67": [["B6", "D9", "D9"], [1, 2, 1]], "A68": [["B6", "D9", "E9"], [1, 2, 1]], "B68": [["B6", "D9", "F9"], [1, 2, 1]], "C68": [["B6", "D9", "G9"], [1, 2, 1]], "D68": [["B6", "D9", "H9"], [1, 2, 1]], "E68": [["B6", "E9", "D10"], [1, 2, 1]], "F68": [["B6", "E9", "E10"], [1, 2, 1]], "G68": [["B6", "E9", "F10"], [1, 2, 1]], "H68": [["B6", "E9", "G10"], [1, 2, 1]], "A69": [["B6", "E9", "H10"], [1, 2, 1]], "B69": [["B6", "E9", "A11"], [1, 2, 1]], "C69": [["C6", "E9", "B11"], [1, 2, 1]], "D69": [["C6", "E9", "C11"], [1, 2, 1]], "E69": [["C6", "E9", "D11"], [1, 2, 1]], "F69": [["C6", "F9", "D10"], [1, 2, 1]], "G69": [["C6", "F9", "E10"], [1, 2, 1]], "H69": [["C6", "F9", "F10"], [1, 2, 1]], "A70": [["C6", "F9", "G10"], [1, 2, 1]], "B70": [["C6", "F9", "H10"], [1, 2, 1]], "C70": [["C6", "F9", "A11"], [1, 2, 1]], "D70": [["C6", "F9", "B11"], [1, 2, 1]], "E70": [["C6", "F9", "C11"], [1, 2, 1]], "F70": [["C6", "F9", "D11"], [1, 2, 1]], "G70": [["C6", "G9", "D10"], [1, 2, 1]], "H70": [["D6", "G9", "E10"], [1, 2, 1]], "A71": [["D6", "G9", "F10"], [1, 2, 1]], "B71": [["D6", "G9", "G10"], [1, 2, 1]], "C71": [["D6", "G9", "H10"], [1, 2, 1]], "D71": [["D6", "G9", "A11"], [1, 2, 1]], "E71": [["D6", "G9", "B11"], [1, 2, 1]], "F71": [["D6", "G9", "C11"], [1, 2, 1]], "G71": [["D6", "G9", "D11"], [1, 2, 1]], "H71": [["D6", "H9", "H11"], [1, 2, 1]], "A72": [["D6", "H9", "A12"], [1, 2, 1]], "B72": [["D6", "H9", "B12"], [1, 2, 1]], "C72": [["D6", "H9", "C12"], [1, 2, 1]], "D72": [["D6", "H9", "D12"], [1, 2, 1]], "E72": [["E6", "H9", "E12"], [1, 2, 1]], "F72": [["E6", "H9", "F12"], [1, 2, 1]], "G72": [["E6", "H9", "G12"], [1, 2, 1]], "H72": [["E6", "H9", "H12"], [1, 2, 1]], "A73": [["E6", "A10", "H11"], [1, 2, 1]], "B73": [["E6", "A10", "A12"], [1, 2, 1]], "C73": [["E6", "A10", "B12"], [1, 2, 1]], "D73": [["E6", "A10", "C12"], [1, 2, 1]], "E73": [["E6", "A10", "D12"], [1, 2, 1]], "F73": [["E6", "A10", "E12"], [1, 2, 1]]}
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

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
final_assembly_dict={"D37": [["G3", "A6", "G1"], [1, 2, 2]], "E37": [["G3", "A6", "H1"], [1, 2, 2]], "F37": [["G3", "A6", "A2"], [1, 2, 2]], "G37": [["G3", "A6", "B2"], [1, 2, 2]], "H37": [["G3", "A6", "C2"], [1, 2, 2]], "A38": [["G3", "A6", "D2"], [1, 2, 2]], "B38": [["G3", "B6", "H2"], [1, 2, 2]], "C38": [["G3", "B6", "A3"], [1, 2, 2]], "D38": [["H3", "B6", "B3"], [1, 2, 2]], "E38": [["H3", "B6", "C3"], [1, 2, 2]], "F38": [["H3", "B6", "D3"], [1, 2, 2]], "G38": [["H3", "B6", "E3"], [1, 2, 2]], "H38": [["H3", "B6", "F3"], [1, 2, 2]], "A39": [["H3", "B6", "G3"], [1, 2, 2]], "B39": [["H3", "B6", "H3"], [1, 2, 2]], "C39": [["H3", "C6", "H2"], [1, 2, 2]], "D39": [["H3", "C6", "A3"], [1, 2, 2]], "E39": [["H3", "C6", "B3"], [1, 2, 2]], "F39": [["H3", "C6", "C3"], [1, 2, 2]], "G39": [["H3", "C6", "D3"], [1, 2, 2]], "H39": [["H3", "C6", "E3"], [1, 2, 2]], "A40": [["A4", "C6", "F3"], [1, 2, 2]], "B40": [["A4", "C6", "G3"], [1, 2, 2]], "C40": [["A4", "C6", "H3"], [1, 2, 2]], "D40": [["A4", "D6", "H2"], [1, 2, 2]], "E40": [["A4", "D6", "A3"], [1, 2, 2]], "F40": [["A4", "D6", "B3"], [1, 2, 2]], "G40": [["A4", "D6", "C3"], [1, 2, 2]], "H40": [["A4", "D6", "D3"], [1, 2, 2]], "A41": [["A4", "D6", "E3"], [1, 2, 2]], "B41": [["A4", "D6", "F3"], [1, 2, 2]], "C41": [["A4", "D6", "G3"], [1, 2, 2]], "D41": [["A4", "D6", "H3"], [1, 2, 2]], "E41": [["A4", "E6", "D7"], [1, 2, 1]], "F41": [["B4", "E6", "E7"], [1, 2, 1]], "G41": [["B4", "E6", "F7"], [1, 2, 1]], "H41": [["B4", "E6", "G7"], [1, 2, 1]], "A42": [["B4", "E6", "H7"], [1, 2, 1]], "B42": [["B4", "E6", "A8"], [1, 2, 1]], "C42": [["B4", "E6", "B8"], [1, 2, 1]], "D42": [["B4", "E6", "C8"], [1, 2, 1]], "E42": [["B4", "E6", "D8"], [1, 2, 1]], "F42": [["B4", "F6", "D7"], [1, 2, 1]], "G42": [["B4", "F6", "E7"], [1, 2, 1]], "H42": [["B4", "F6", "F7"], [1, 2, 1]], "A43": [["B4", "F6", "G7"], [1, 2, 1]], "B43": [["B4", "F6", "H7"], [1, 2, 1]], "C43": [["C4", "F6", "A8"], [1, 2, 1]], "D43": [["C4", "F6", "B8"], [1, 2, 1]], "E43": [["C4", "F6", "C8"], [1, 2, 1]], "F43": [["C4", "F6", "D8"], [1, 2, 1]], "G43": [["C4", "G6", "D7"], [1, 2, 1]], "H43": [["C4", "G6", "E7"], [1, 2, 1]], "A44": [["C4", "G6", "F7"], [1, 2, 1]], "B44": [["C4", "G6", "G7"], [1, 2, 1]], "C44": [["C4", "G6", "H7"], [1, 2, 1]], "D44": [["C4", "G6", "A8"], [1, 2, 1]], "E44": [["C4", "G6", "B8"], [1, 2, 1]], "F44": [["C4", "G6", "C8"], [1, 2, 1]], "G44": [["C4", "G6", "D8"], [1, 2, 1]], "H44": [["D4", "H6", "H8"], [1, 2, 1]], "A45": [["D4", "H6", "A9"], [1, 2, 1]], "B45": [["D4", "H6", "B9"], [1, 2, 1]], "C45": [["D4", "H6", "C9"], [1, 2, 1]], "D45": [["D4", "H6", "D9"], [1, 2, 1]], "E45": [["D4", "H6", "E9"], [1, 2, 1]], "F45": [["D4", "H6", "F9"], [1, 2, 1]], "G45": [["D4", "H6", "G9"], [1, 2, 1]], "H45": [["D4", "H6", "H9"], [1, 2, 1]], "A46": [["D4", "A7", "H8"], [1, 2, 1]], "B46": [["D4", "A7", "A9"], [1, 2, 1]], "C46": [["D4", "A7", "B9"], [1, 2, 1]], "D46": [["D4", "A7", "C9"], [1, 2, 1]], "E46": [["E4", "A7", "D9"], [1, 2, 1]], "F46": [["E4", "A7", "E9"], [1, 2, 1]], "G46": [["E4", "A7", "F9"], [1, 2, 1]], "H46": [["E4", "A7", "G9"], [1, 2, 1]], "A47": [["E4", "A7", "H9"], [1, 2, 1]], "B47": [["E4", "B7", "H8"], [1, 2, 1]], "C47": [["E4", "B7", "A9"], [1, 2, 1]], "D47": [["E4", "B7", "B9"], [1, 2, 1]], "E47": [["E4", "B7", "C9"], [1, 2, 1]], "F47": [["E4", "B7", "D9"], [1, 2, 1]], "G47": [["E4", "B7", "E9"], [1, 2, 1]], "H47": [["E4", "B7", "F9"], [1, 2, 1]], "A48": [["E4", "B7", "G9"], [1, 2, 1]], "B48": [["F4", "B7", "H9"], [1, 2, 1]], "C48": [["F4", "C7", "D10"], [1, 2, 1]], "D48": [["F4", "C7", "E10"], [1, 2, 1]], "E48": [["F4", "C7", "F10"], [1, 2, 1]], "F48": [["F4", "C7", "G10"], [1, 2, 1]], "G48": [["F4", "C7", "H10"], [1, 2, 1]], "H48": [["F4", "C7", "A11"], [1, 2, 1]], "A49": [["F4", "C7", "B11"], [1, 2, 1]], "B49": [["F4", "C7", "C11"], [1, 2, 1]], "C49": [["F4", "C7", "D11"], [1, 2, 1]], "D49": [["F4", "D7", "D10"], [1, 2, 1]]}
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

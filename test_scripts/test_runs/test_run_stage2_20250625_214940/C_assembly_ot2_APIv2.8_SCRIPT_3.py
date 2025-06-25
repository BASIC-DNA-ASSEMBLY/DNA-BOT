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
final_assembly_dict={"E25": [["H2", "F4", "G9"], [1, 2, 1]], "F25": [["H2", "F4", "H9"], [1, 2, 1]], "G25": [["H2", "G4", "H8"], [1, 2, 1]], "H25": [["H2", "G4", "A9"], [1, 2, 1]], "A26": [["H2", "G4", "B9"], [1, 2, 1]], "B26": [["H2", "G4", "C9"], [1, 2, 1]], "C26": [["H2", "G4", "D9"], [1, 2, 1]], "D26": [["H2", "G4", "E9"], [1, 2, 1]], "E26": [["H2", "G4", "F9"], [1, 2, 1]], "F26": [["H2", "G4", "G9"], [1, 2, 1]], "G26": [["H2", "G4", "H9"], [1, 2, 1]], "H26": [["H2", "H4", "H8"], [1, 2, 1]], "A27": [["A3", "H4", "A9"], [1, 2, 1]], "B27": [["A3", "H4", "B9"], [1, 2, 1]], "C27": [["A3", "H4", "C9"], [1, 2, 1]], "D27": [["A3", "H4", "D9"], [1, 2, 1]], "E27": [["A3", "H4", "E9"], [1, 2, 1]], "F27": [["A3", "H4", "F9"], [1, 2, 1]], "G27": [["A3", "H4", "G9"], [1, 2, 1]], "H27": [["A3", "H4", "H9"], [1, 2, 1]], "A28": [["A3", "A5", "D10"], [1, 2, 1]], "B28": [["A3", "A5", "E10"], [1, 2, 1]], "C28": [["A3", "A5", "F10"], [1, 2, 1]], "D28": [["A3", "A5", "G10"], [1, 2, 1]], "E28": [["A3", "A5", "H10"], [1, 2, 1]], "F28": [["B3", "A5", "A11"], [1, 2, 1]], "G28": [["B3", "A5", "B11"], [1, 2, 1]], "H28": [["B3", "A5", "C11"], [1, 2, 1]], "A29": [["B3", "A5", "D11"], [1, 2, 1]], "B29": [["B3", "B5", "D10"], [1, 2, 1]], "C29": [["B3", "B5", "E10"], [1, 2, 1]], "D29": [["B3", "B5", "F10"], [1, 2, 1]], "E29": [["B3", "B5", "G10"], [1, 2, 1]], "F29": [["B3", "B5", "H10"], [1, 2, 1]], "G29": [["B3", "B5", "A11"], [1, 2, 1]], "H29": [["B3", "B5", "B11"], [1, 2, 1]], "A30": [["B3", "B5", "C11"], [1, 2, 1]], "B30": [["B3", "B5", "D11"], [1, 2, 1]], "C30": [["C3", "C5", "D10"], [1, 2, 1]], "D30": [["C3", "C5", "E10"], [1, 2, 1]], "E30": [["C3", "C5", "F10"], [1, 2, 1]], "F30": [["C3", "C5", "G10"], [1, 2, 1]], "G30": [["C3", "C5", "H10"], [1, 2, 1]], "H30": [["C3", "C5", "A11"], [1, 2, 1]], "A31": [["C3", "C5", "B11"], [1, 2, 1]], "B31": [["C3", "C5", "C11"], [1, 2, 1]], "C31": [["C3", "C5", "D11"], [1, 2, 1]], "D31": [["C3", "D5", "H11"], [1, 2, 1]], "E31": [["C3", "D5", "A12"], [1, 2, 1]], "F31": [["C3", "D5", "B12"], [1, 2, 1]], "G31": [["C3", "D5", "C12"], [1, 2, 1]], "H31": [["D3", "D5", "D12"], [1, 2, 1]], "A32": [["D3", "D5", "E12"], [1, 2, 1]], "B32": [["D3", "D5", "F12"], [1, 2, 1]], "C32": [["D3", "D5", "G12"], [1, 2, 1]], "D32": [["D3", "D5", "H12"], [1, 2, 1]], "E32": [["D3", "E5", "H11"], [1, 2, 1]], "F32": [["D3", "E5", "A12"], [1, 2, 1]], "G32": [["D3", "E5", "B12"], [1, 2, 1]], "H32": [["D3", "E5", "C12"], [1, 2, 1]], "A33": [["D3", "E5", "D12"], [1, 2, 1]], "B33": [["D3", "E5", "E12"], [1, 2, 1]], "C33": [["D3", "E5", "F12"], [1, 2, 1]], "D33": [["D3", "E5", "G12"], [1, 2, 1]], "E33": [["E3", "E5", "H12"], [1, 2, 1]], "F33": [["E3", "F5", "H11"], [1, 2, 1]], "G33": [["E3", "F5", "A12"], [1, 2, 1]], "H33": [["E3", "F5", "B12"], [1, 2, 1]], "A34": [["E3", "F5", "C12"], [1, 2, 1]], "B34": [["E3", "F5", "D12"], [1, 2, 1]], "C34": [["E3", "F5", "E12"], [1, 2, 1]], "D34": [["E3", "F5", "F12"], [1, 2, 1]], "E34": [["E3", "F5", "G12"], [1, 2, 1]], "F34": [["E3", "F5", "H12"], [1, 2, 1]], "G34": [["E3", "G5", "D1"], [1, 2, 2]], "H34": [["E3", "G5", "E1"], [1, 2, 2]], "A35": [["E3", "G5", "F1"], [1, 2, 2]], "B35": [["F3", "G5", "G1"], [1, 2, 2]], "C35": [["F3", "G5", "H1"], [1, 2, 2]], "D35": [["F3", "G5", "A2"], [1, 2, 2]], "E35": [["F3", "G5", "B2"], [1, 2, 2]], "F35": [["F3", "G5", "C2"], [1, 2, 2]], "G35": [["F3", "G5", "D2"], [1, 2, 2]], "H35": [["F3", "H5", "D1"], [1, 2, 2]], "A36": [["F3", "H5", "E1"], [1, 2, 2]], "B36": [["F3", "H5", "F1"], [1, 2, 2]], "C36": [["F3", "H5", "G1"], [1, 2, 2]], "D36": [["F3", "H5", "H1"], [1, 2, 2]], "E36": [["F3", "H5", "A2"], [1, 2, 2]], "F36": [["F3", "H5", "B2"], [1, 2, 2]], "G36": [["G3", "H5", "C2"], [1, 2, 2]], "H36": [["G3", "H5", "D2"], [1, 2, 2]], "A37": [["G3", "A6", "D1"], [1, 2, 2]], "B37": [["G3", "A6", "E1"], [1, 2, 2]], "C37": [["G3", "A6", "F1"], [1, 2, 2]], "D37": [["G3", "A6", "G1"], [1, 2, 2]], "E37": [["G3", "A6", "H1"], [1, 2, 2]], "F37": [["G3", "A6", "A2"], [1, 2, 2]]}
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

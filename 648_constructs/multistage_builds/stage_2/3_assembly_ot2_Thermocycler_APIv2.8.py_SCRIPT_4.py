from __future__ import unicode_literals
from opentrons import protocol_api
import numpy as np
# metadata
metadata = {
'protocolName': 'DNABOT Assembly Thermocycler',
'description': 'DNABOT Assembly Step3 with Thermocycler',
'apiLevel': '2.8'
}

# It is possible to run 88 assemblies with this new module. The heat block module is removed. 
# Assembly reactions is set up on thermocycler module.


# test dictionary can be used for simulation 3 or 88 assemblies
final_assembly_dict={"A1": ['A7', 'B7', 'C7', 'F7'], "B1": ['A7', 'B7', 'D7', 'G7'], "C1": ['A7', 'B7', 'E7', 'H7']}
tiprack_num=1

#final_assembly_dict={"A1": ["A7", "G7", "H7", "A8", "B8"], "B1": ["A7", "D8", "E8", "F8", "G8"], "C1": ["A7", "D8", "H7", "H8", "B9"], "D1": ["A7", "C9", "E9", "G9", "B8"], "E1": ["A7", "H9", "B10", "E9", "D10"], "F1": ["A7", "C9", "H8", "F10", "D10"], "G1": ["A7", "C9", "H10", "E8", "B9"], "H1": ["A7", "H9", "F8", "H10", "B11"], "A2": ["A7", "G7", "E8", "B10", "G8"], "B2": ["A7", "G7", "D11", "A8", "B9"], "C2": ["A7", "C9", "E9", "G9", "B9"], "D2": ["A7", "G7", "H7", "H8", "B8"], "E2": ["A7", "F11", "H11", "H7", "B12"], "F2": ["A7", "C9", "H8", "H11", "D10"], "G2": ["A7", "G7", "D11", "A8", "B8"], "H2": ["B7", "F11", "B10", "H10", "B11"], "A3": ["B7", "D8", "H7", "H8", "B8"], "B3": ["B7", "C9", "H10", "G9", "B8"], "C3": ["B7", "D12", "H8", "H11", "B11"], "D3": ["B7", "D12", "E9", "E8", "B8"], "E3": ["B7", "D12", "E9", "E8", "B9"], "F3": ["B7", "H9", "B10", "H10", "D10"], "G3": ["B7", "G7", "D11", "H8", "B8"], "H3": ["B7", "D12", "H10", "G9", "B9"], "A4": ["B7", "F11", "F10", "D11", "B12"], "B4": ["B7", "G7", "H7", "A8", "B9"], "C4": ["B7", "G7", "E8", "B10", "B12"], "D4": ["B7", "H9", "H11", "H7", "G8"], "E4": ["B7", "D8", "E8", "F8", "B12"], "F4": ["B7", "D12", "E9", "G9", "B8"], "G4": ["C7", "H9", "B10", "E9", "B11"], "H4": ["C7", "F11", "B10", "H10", "D10"], "A5": ["C7", "H9", "F8", "E9", "B11"], "B5": ["C7", "D12", "H8", "F10", "B11"], "C5": ["C7", "F11", "F8", "H10", "B11"], "D5": ["C7", "F11", "H11", "H7", "G8"], "E5": ["C7", "D8", "D11", "A8", "B9"], "F5": ["C7", "H9", "H11", "H7", "B12"], "G5": ["C7", "C9", "H10", "G9", "B9"], "H5": ["C7", "H9", "F10", "H7", "G8"], "A6": ["C7", "D12", "A8", "H11", "D10"], "B6": ["C7", "C9", "A8", "H11", "B11"], "C6": ["C7", "F11", "H11", "D11", "B12"], "D6": ["C7", "D8", "E8", "B10", "G8"], "E6": ["C7", "C9", "H8", "H11", "B11"], "F6": ["D7", "D8", "G9", "F8", "G8"], "G6": ["D7", "C9", "A8", "F10", "B11"], "H6": ["D7", "F11", "F10", "H7", "B12"], "A7": ["D7", "C9", "A8", "F10", "D10"], "B7": ["D7", "H9", "F8", "E9", "D10"], "C7": ["D7", "G7", "G9", "F8", "B12"], "D7": ["D7", "D12", "A8", "H11", "B11"], "E7": ["D7", "D12", "H10", "G9", "B8"], "F7": ["D7", "H9", "H11", "D11", "B12"], "G7": ["D7", "C9", "H8", "F10", "B11"], "H7": ["D7", "D8", "D11", "H8", "B8"], "A8": ["D7", "C9", "E9", "E8", "B9"], "B8": ["D7", "H9", "F10", "D11", "G8"], "C8": ["D7", "H9", "H11", "D11", "G8"], "D8": ["D7", "D12", "A8", "F10", "D10"], "E8": ["E7", "G7", "G9", "F8", "G8"], "F8": ["E7", "D12", "A8", "F10", "B11"], "G8": ["E7", "H9", "F10", "D11", "B12"], "H8": ["E7", "D8", "E8", "B10", "B12"], "A9": ["E7", "C9", "E9", "E8", "B8"], "B9": ["E7", "F11", "B10", "E9", "D10"], "C9": ["E7", "D12", "H8", "F10", "D10"], "D9": ["E7", "H9", "B10", "H10", "B11"], "E9": ["E7", "D8", "G9", "F8", "B12"], "F9": ["E7", "F11", "B10", "E9", "B11"], "G9": ["E7", "F11", "F8", "E9", "C11"], "H9": ["E7", "G7", "G9", "B10", "B12"], "A10": ["E7", "D8", "G9", "B10", "B12"], "B10": ["E7", "D8", "D11", "A8", "B8"], "C10": ["E7", "F11", "F10", "H7", "G8"], "D10": ["F7", "F11", "F8", "E9", "D10"], "E10": ["F7", "H9", "F10", "H7", "B12"], "F10": ["F7", "D12", "H10", "E8", "B9"], "G10": ["F7", "C9", "H10", "E8", "B8"], "H10": ["F7", "F11", "F8", "H10", "D10"], "A11": ["F7", "D12", "H10", "E8", "B8"], "B11": ["F7", "G7", "H7", "H8", "B9"], "C11": ["F7", "G7", "G9", "B10", "G8"], "D11": ["F7", "D12", "H8", "H11", "D10"], "E11": ["F7", "D9", "A8", "H11", "D10"], "F11": ["F7", "G7", "D11", "H8", "B9"], "G11": ["F7", "F11", "A12", "D11", "G8"], "H11": ["F7", "D8", "D11", "A9", "B9"]}
#tiprack_num=5

# opentrons_simulate.exe dnabot\template_ot2_scripts\assembly_template_Thermocycler_module_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

final_assembly_dict={"A1": [["A1", "C5", "F12"], [1, 2, 1]], "B1": [["A1", "C5", "G12"], [1, 2, 1]], "C1": [["A1", "C5", "H12"], [1, 2, 1]], "D1": [["A1", "C5", "A1"], [1, 2, 2]], "E1": [["A1", "C5", "B1"], [1, 2, 2]], "F1": [["A1", "C5", "C1"], [1, 2, 2]], "G1": [["A1", "C5", "D1"], [1, 2, 2]], "H1": [["A1", "C5", "E1"], [1, 2, 2]], "A2": [["A1", "C5", "F1"], [1, 2, 2]], "B2": [["A1", "D5", "B2"], [1, 2, 2]], "C2": [["A1", "D5", "C2"], [1, 2, 2]], "D2": [["A1", "D5", "D2"], [1, 2, 2]], "E2": [["A1", "D5", "E2"], [1, 2, 2]], "F2": [["A1", "D5", "F2"], [1, 2, 2]], "G2": [["A1", "D5", "G2"], [1, 2, 2]], "H2": [["B1", "D5", "H2"], [1, 2, 2]], "A3": [["B1", "D5", "A3"], [1, 2, 2]], "B3": [["B1", "D5", "B3"], [1, 2, 2]], "C3": [["B1", "E5", "B2"], [1, 2, 2]], "D3": [["B1", "E5", "C2"], [1, 2, 2]], "E3": [["B1", "E5", "D2"], [1, 2, 2]], "F3": [["B1", "E5", "E2"], [1, 2, 2]], "G3": [["B1", "E5", "F2"], [1, 2, 2]], "H3": [["B1", "E5", "G2"], [1, 2, 2]], "A4": [["B1", "E5", "H2"], [1, 2, 2]], "B4": [["B1", "E5", "A3"], [1, 2, 2]], "C4": [["B1", "E5", "B3"], [1, 2, 2]], "D4": [["B1", "F5", "B2"], [1, 2, 2]], "E4": [["B1", "F5", "C2"], [1, 2, 2]], "F4": [["B1", "F5", "D2"], [1, 2, 2]], "G4": [["C1", "F5", "E2"], [1, 2, 2]], "H4": [["C1", "F5", "F2"], [1, 2, 2]], "A5": [["C1", "F5", "G2"], [1, 2, 2]], "B5": [["C1", "F5", "H2"], [1, 2, 2]], "C5": [["C1", "F5", "A3"], [1, 2, 2]], "D5": [["C1", "F5", "B3"], [1, 2, 2]], "E5": [["C1", "G5", "F6"], [1, 2, 1]], "F5": [["C1", "G5", "G6"], [1, 2, 1]], "G5": [["C1", "G5", "H6"], [1, 2, 1]], "H5": [["C1", "G5", "A7"], [1, 2, 1]], "A6": [["C1", "G5", "B7"], [1, 2, 1]], "B6": [["C1", "G5", "C7"], [1, 2, 1]], "C6": [["C1", "G5", "D7"], [1, 2, 1]], "D6": [["C1", "G5", "E7"], [1, 2, 1]], "E6": [["C1", "G5", "F7"], [1, 2, 1]], "F6": [["D1", "H5", "F6"], [1, 2, 1]], "G6": [["D1", "H5", "G6"], [1, 2, 1]], "H6": [["D1", "H5", "H6"], [1, 2, 1]], "A7": [["D1", "H5", "A7"], [1, 2, 1]], "B7": [["D1", "H5", "B7"], [1, 2, 1]], "C7": [["D1", "H5", "C7"], [1, 2, 1]], "D7": [["D1", "H5", "D7"], [1, 2, 1]], "E7": [["D1", "H5", "E7"], [1, 2, 1]], "F7": [["D1", "H5", "F7"], [1, 2, 1]], "G7": [["D1", "A6", "F6"], [1, 2, 1]], "H7": [["D1", "A6", "G6"], [1, 2, 1]], "A8": [["D1", "A6", "H6"], [1, 2, 1]], "B8": [["D1", "A6", "A7"], [1, 2, 1]], "C8": [["D1", "A6", "B7"], [1, 2, 1]], "D8": [["D1", "A6", "C7"], [1, 2, 1]], "E8": [["E1", "A6", "D7"], [1, 2, 1]], "F8": [["E1", "A6", "E7"], [1, 2, 1]], "G8": [["E1", "A6", "F7"], [1, 2, 1]], "H8": [["E1", "B6", "B8"], [1, 2, 1]], "A9": [["E1", "B6", "C8"], [1, 2, 1]], "B9": [["E1", "B6", "D8"], [1, 2, 1]], "C9": [["E1", "B6", "E8"], [1, 2, 1]], "D9": [["E1", "B6", "F8"], [1, 2, 1]], "E9": [["E1", "B6", "G8"], [1, 2, 1]], "F9": [["E1", "B6", "H8"], [1, 2, 1]], "G9": [["E1", "B6", "A9"], [1, 2, 1]], "H9": [["E1", "B6", "B9"], [1, 2, 1]], "A10": [["E1", "C6", "B8"], [1, 2, 1]], "B10": [["E1", "C6", "C8"], [1, 2, 1]], "C10": [["E1", "C6", "D8"], [1, 2, 1]], "D10": [["F1", "C6", "E8"], [1, 2, 1]], "E10": [["F1", "C6", "F8"], [1, 2, 1]], "F10": [["F1", "C6", "G8"], [1, 2, 1]], "G10": [["F1", "C6", "H8"], [1, 2, 1]], "H10": [["F1", "C6", "A9"], [1, 2, 1]], "A11": [["F1", "C6", "B9"], [1, 2, 1]], "B11": [["F1", "D6", "B8"], [1, 2, 1]], "C11": [["F1", "D6", "C8"], [1, 2, 1]], "D11": [["F1", "D6", "D8"], [1, 2, 1]], "E11": [["F1", "D6", "E8"], [1, 2, 1]], "F11": [["F1", "D6", "F8"], [1, 2, 1]], "G11": [["F1", "D6", "G8"], [1, 2, 1]], "H11": [["F1", "D6", "H8"], [1, 2, 1]], "A12": [["F1", "D6", "A9"], [1, 2, 1]], "B12": [["F1", "D6", "B9"], [1, 2, 1]], "C12": [["G1", "E6", "F9"], [1, 2, 1]], "D12": [["G1", "E6", "G9"], [1, 2, 1]], "E12": [["G1", "E6", "H9"], [1, 2, 1]], "F12": [["G1", "E6", "A10"], [1, 2, 1]], "G12": [["G1", "E6", "B10"], [1, 2, 1]], "H12": [["G1", "E6", "C10"], [1, 2, 1]]}
tiprack_num=4


def run(protocol: protocol_api.ProtocolContext):
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_tiprack_20ul"):
            # Constants, we update all the labware name in version 2

            # Plate of sample after  purification
            # MAG_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            # MAG_PLATE_POSITION = '1'
            # magbead_plate = protocol.load_labware(MAG_PLATE_TYPE, MAG_PLATE_POSITION)         #########################

            # Source plate(s)
            SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'

            # Tuberack
            TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            TUBE_RACK_POSITION = '4'
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)

            # Destination plate
            DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            TOTAL_VOL = 15
            PART_VOL = 1.5
            MIX_SETTINGS = (1, 3)
            tiprack_num=tiprack_num+1

            # Error trapping
            sample_number = len(final_assembly_dict.keys())
            if sample_number > 96:
                raise ValueError('Final assembly nummber cannot exceed 96.')

            #Tiprack
            TIPS_REQUIRED = sum([len(value) for key, value in list(final_assembly_dict.items())])

            print(TIPS_REQUIRED)


            CANDIDATE_TIPRACK_SLOTS = ['3', '5', '6', '9']
            TIP_COUNT = [96 for i in CANDIDATE_TIPRACK_SLOTS]
            PIPETTE_MOUNT = 'right'

            
            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            pipette = protocol.load_instrument('p20_single_gen2', PIPETTE_MOUNT, tip_racks=tipracks)
            
            return
            
            #Thermocycler Module
            tc_mod = protocol.load_module('Thermocycler Module')
            destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
            tc_mod.set_block_temperature(20)


             # Master mix transfers
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))
            master_mix_well_letters = ['A', 'B', 'C', 'D']
            destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
            
            for x in unique_assemblies_lens: 
                master_mix_well = master_mix_well_letters[0] + str(x - 1)    # select well in tube rack, for assembly with x number of parts put tube in column x 
                destination_inds = [i for i, lens in enumerate(final_assembly_lens) if lens == x]   # find all assemblies of length x
                destination_wells_for_len = list(destination_wells[destination_inds])

                # Common backbones (up to 3, else complete original way)
                backbones = []
                for components in final_assembly_dict.values():
                    backbone = components[0]
                    backbones.append(backbone) 
                unique_backbones = list(set(backbones))

                if len(unique_backbones) < 3:
                    print('finish this part')
                    ''''Notes
                        The idea here was to produce a master mix for each backbone and dilution combination.
                        I could do away with the user having to premix their own dilution of the assembly buffer
                        and have the robot make a series of master mixes for each unique number of parts (i.e. 
                        dilution) and backbone. This would save a lot of time and tips.

                        This could be done by having the user make a stronger buffer and putting a limit on the 
                        max number of parts for an assembly 
                    '''
                
                pipette.pick_up_tip()
                for destination_well in destination_wells_for_len: # make tube_rack_wells and destination_plate.wells in the same type
                    
                    pipette.transfer(TOTAL_VOL - x * PART_VOL, tube_rack.wells(master_mix_well),
                                     destination_plate.wells(destination_well), new_tip='never')    #transfer water and buffer in the pipette

                pipette.drop_tip()
            return

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for value in values:# magbead_plate.wells and destination_plate.wells in the same type
                    pipette.transfer(PART_VOL, magbead_plate.wells(value),
                                     destination_plate.wells(key), mix_after=MIX_SETTINGS,
                                     new_tip='always')#transfer parts in one tube



            #Thermocycler Module
            tc_mod.close_lid()
            tc_mod.set_lid_temperature(105)
            tc_mod.set_block_temperature(50, hold_time_minutes=45, block_max_volume=15)
            tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=30)
            # Increase the hold time at 4 C if necessary
            tc_mod.set_lid_temperature(37)
            tc_mod.open_lid()

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

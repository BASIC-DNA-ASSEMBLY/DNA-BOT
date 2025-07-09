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
final_assembly_dict={"A1": [["F4", "C7", "B11"], [1, 2, 1]], "B1": [["F4", "C7", "C11"], [1, 2, 1]], "C1": [["F4", "C7", "D11"], [1, 2, 1]], "D1": [["F4", "D7", "D10"], [1, 2, 1]], "E1": [["F4", "D7", "E10"], [1, 2, 1]], "F1": [["F4", "D7", "F10"], [1, 2, 1]], "G1": [["G4", "D7", "G10"], [1, 2, 1]], "H1": [["G4", "D7", "H10"], [1, 2, 1]], "A2": [["G4", "D7", "A11"], [1, 2, 1]], "B2": [["G4", "D7", "B11"], [1, 2, 1]], "C2": [["G4", "D7", "C11"], [1, 2, 1]], "D2": [["G4", "D7", "D11"], [1, 2, 1]], "E2": [["G4", "E7", "D10"], [1, 2, 1]], "F2": [["G4", "E7", "E10"], [1, 2, 1]], "G2": [["G4", "E7", "F10"], [1, 2, 1]], "H2": [["G4", "E7", "G10"], [1, 2, 1]], "A3": [["G4", "E7", "H10"], [1, 2, 1]], "B3": [["G4", "E7", "A11"], [1, 2, 1]], "C3": [["G4", "E7", "B11"], [1, 2, 1]], "D3": [["H4", "E7", "C11"], [1, 2, 1]], "E3": [["H4", "E7", "D11"], [1, 2, 1]], "F3": [["H4", "F7", "H11"], [1, 2, 1]], "G3": [["H4", "F7", "A12"], [1, 2, 1]], "H3": [["H4", "F7", "B12"], [1, 2, 1]], "A4": [["H4", "F7", "C12"], [1, 2, 1]], "B4": [["H4", "F7", "D12"], [1, 2, 1]], "C4": [["H4", "F7", "E12"], [1, 2, 1]], "D4": [["H4", "F7", "F12"], [1, 2, 1]], "E4": [["H4", "F7", "G12"], [1, 2, 1]], "F4": [["H4", "F7", "H12"], [1, 2, 1]], "G4": [["H4", "G7", "H11"], [1, 2, 1]], "H4": [["H4", "G7", "A12"], [1, 2, 1]], "A5": [["A5", "G7", "B12"], [1, 2, 1]], "B5": [["A5", "G7", "C12"], [1, 2, 1]], "C5": [["A5", "G7", "D12"], [1, 2, 1]], "D5": [["A5", "G7", "E12"], [1, 2, 1]], "E5": [["A5", "G7", "F12"], [1, 2, 1]], "F5": [["A5", "G7", "G12"], [1, 2, 1]], "G5": [["A5", "G7", "H12"], [1, 2, 1]], "H5": [["A5", "H7", "H11"], [1, 2, 1]], "A6": [["A5", "H7", "A12"], [1, 2, 1]], "B6": [["A5", "H7", "B12"], [1, 2, 1]], "C6": [["A5", "H7", "C12"], [1, 2, 1]], "D6": [["A5", "H7", "D12"], [1, 2, 1]], "E6": [["A5", "H7", "E12"], [1, 2, 1]], "F6": [["B5", "H7", "F12"], [1, 2, 1]], "G6": [["B5", "H7", "G12"], [1, 2, 1]], "H6": [["B5", "H7", "H12"], [1, 2, 1]], "A7": [["B5", "A8", "D1"], [1, 2, 2]], "B7": [["B5", "A8", "E1"], [1, 2, 2]], "C7": [["B5", "A8", "F1"], [1, 2, 2]], "D7": [["B5", "A8", "G1"], [1, 2, 2]], "E7": [["B5", "A8", "H1"], [1, 2, 2]], "F7": [["B5", "A8", "A2"], [1, 2, 2]], "G7": [["B5", "A8", "B2"], [1, 2, 2]], "H7": [["B5", "A8", "C2"], [1, 2, 2]], "A8": [["B5", "A8", "D2"], [1, 2, 2]], "B8": [["B5", "B8", "D1"], [1, 2, 2]], "C8": [["C5", "B8", "E1"], [1, 2, 2]], "D8": [["C5", "B8", "F1"], [1, 2, 2]], "E8": [["C5", "B8", "G1"], [1, 2, 2]], "F8": [["C5", "B8", "H1"], [1, 2, 2]], "G8": [["C5", "B8", "A2"], [1, 2, 2]], "H8": [["C5", "B8", "B2"], [1, 2, 2]], "A9": [["C5", "B8", "C2"], [1, 2, 2]], "B9": [["C5", "B8", "D2"], [1, 2, 2]], "C9": [["C5", "C8", "D1"], [1, 2, 2]], "D9": [["C5", "C8", "E1"], [1, 2, 2]], "E9": [["C5", "C8", "F1"], [1, 2, 2]], "F9": [["C5", "C8", "G1"], [1, 2, 2]], "G9": [["C5", "C8", "H1"], [1, 2, 2]], "H9": [["D5", "C8", "A2"], [1, 2, 2]], "A10": [["D5", "C8", "B2"], [1, 2, 2]], "B10": [["D5", "C8", "C2"], [1, 2, 2]], "C10": [["D5", "C8", "D2"], [1, 2, 2]], "D10": [["D5", "D8", "H2"], [1, 2, 2]], "E10": [["D5", "D8", "A3"], [1, 2, 2]], "F10": [["D5", "D8", "B3"], [1, 2, 2]], "G10": [["D5", "D8", "C3"], [1, 2, 2]], "H10": [["D5", "D8", "D3"], [1, 2, 2]], "A11": [["D5", "D8", "E3"], [1, 2, 2]], "B11": [["D5", "D8", "F3"], [1, 2, 2]], "C11": [["D5", "D8", "G3"], [1, 2, 2]], "D11": [["D5", "D8", "H3"], [1, 2, 2]], "E11": [["E5", "E8", "H2"], [1, 2, 2]], "F11": [["E5", "E8", "A3"], [1, 2, 2]], "G11": [["E5", "E8", "B3"], [1, 2, 2]], "H11": [["E5", "E8", "C3"], [1, 2, 2]], "A12": [["E5", "E8", "D3"], [1, 2, 2]], "B12": [["E5", "E8", "E3"], [1, 2, 2]], "C12": [["E5", "E8", "F3"], [1, 2, 2]], "D12": [["E5", "E8", "G3"], [1, 2, 2]], "E12": [["E5", "E8", "H3"], [1, 2, 2]], "F12": [["E5", "F8", "H2"], [1, 2, 2]], "G12": [["E5", "F8", "A3"], [1, 2, 2]], "H12": [["E5", "F8", "B3"], [1, 2, 2]]}
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

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
final_assembly_dict={"A1": [["G2", "F4", "C9"], [1, 2, 1]], "B1": [["G2", "F4", "D9"], [1, 2, 1]], "C1": [["G2", "F4", "E9"], [1, 2, 1]], "D1": [["H2", "F4", "F9"], [1, 2, 1]], "E1": [["H2", "F4", "G9"], [1, 2, 1]], "F1": [["H2", "F4", "H9"], [1, 2, 1]], "G1": [["H2", "G4", "H8"], [1, 2, 1]], "H1": [["H2", "G4", "A9"], [1, 2, 1]], "A2": [["H2", "G4", "B9"], [1, 2, 1]], "B2": [["H2", "G4", "C9"], [1, 2, 1]], "C2": [["H2", "G4", "D9"], [1, 2, 1]], "D2": [["H2", "G4", "E9"], [1, 2, 1]], "E2": [["H2", "G4", "F9"], [1, 2, 1]], "F2": [["H2", "G4", "G9"], [1, 2, 1]], "G2": [["H2", "G4", "H9"], [1, 2, 1]], "H2": [["H2", "H4", "H8"], [1, 2, 1]], "A3": [["A3", "H4", "A9"], [1, 2, 1]], "B3": [["A3", "H4", "B9"], [1, 2, 1]], "C3": [["A3", "H4", "C9"], [1, 2, 1]], "D3": [["A3", "H4", "D9"], [1, 2, 1]], "E3": [["A3", "H4", "E9"], [1, 2, 1]], "F3": [["A3", "H4", "F9"], [1, 2, 1]], "G3": [["A3", "H4", "G9"], [1, 2, 1]], "H3": [["A3", "H4", "H9"], [1, 2, 1]], "A4": [["A3", "A5", "D10"], [1, 2, 1]], "B4": [["A3", "A5", "E10"], [1, 2, 1]], "C4": [["A3", "A5", "F10"], [1, 2, 1]], "D4": [["A3", "A5", "G10"], [1, 2, 1]], "E4": [["A3", "A5", "H10"], [1, 2, 1]], "F4": [["B3", "A5", "A11"], [1, 2, 1]], "G4": [["B3", "A5", "B11"], [1, 2, 1]], "H4": [["B3", "A5", "C11"], [1, 2, 1]], "A5": [["B3", "A5", "D11"], [1, 2, 1]], "B5": [["B3", "B5", "D10"], [1, 2, 1]], "C5": [["B3", "B5", "E10"], [1, 2, 1]], "D5": [["B3", "B5", "F10"], [1, 2, 1]], "E5": [["B3", "B5", "G10"], [1, 2, 1]], "F5": [["B3", "B5", "H10"], [1, 2, 1]], "G5": [["B3", "B5", "A11"], [1, 2, 1]], "H5": [["B3", "B5", "B11"], [1, 2, 1]], "A6": [["B3", "B5", "C11"], [1, 2, 1]], "B6": [["B3", "B5", "D11"], [1, 2, 1]], "C6": [["C3", "C5", "D10"], [1, 2, 1]], "D6": [["C3", "C5", "E10"], [1, 2, 1]], "E6": [["C3", "C5", "F10"], [1, 2, 1]], "F6": [["C3", "C5", "G10"], [1, 2, 1]], "G6": [["C3", "C5", "H10"], [1, 2, 1]], "H6": [["C3", "C5", "A11"], [1, 2, 1]], "A7": [["C3", "C5", "B11"], [1, 2, 1]], "B7": [["C3", "C5", "C11"], [1, 2, 1]], "C7": [["C3", "C5", "D11"], [1, 2, 1]], "D7": [["C3", "D5", "H11"], [1, 2, 1]], "E7": [["C3", "D5", "A12"], [1, 2, 1]], "F7": [["C3", "D5", "B12"], [1, 2, 1]], "G7": [["C3", "D5", "C12"], [1, 2, 1]], "H7": [["D3", "D5", "D12"], [1, 2, 1]], "A8": [["D3", "D5", "E12"], [1, 2, 1]], "B8": [["D3", "D5", "F12"], [1, 2, 1]], "C8": [["D3", "D5", "G12"], [1, 2, 1]], "D8": [["D3", "D5", "H12"], [1, 2, 1]], "E8": [["D3", "E5", "H11"], [1, 2, 1]], "F8": [["D3", "E5", "A12"], [1, 2, 1]], "G8": [["D3", "E5", "B12"], [1, 2, 1]], "H8": [["D3", "E5", "C12"], [1, 2, 1]], "A9": [["D3", "E5", "D12"], [1, 2, 1]], "B9": [["D3", "E5", "E12"], [1, 2, 1]], "C9": [["D3", "E5", "F12"], [1, 2, 1]], "D9": [["D3", "E5", "G12"], [1, 2, 1]], "E9": [["E3", "E5", "H12"], [1, 2, 1]], "F9": [["E3", "F5", "H11"], [1, 2, 1]], "G9": [["E3", "F5", "A12"], [1, 2, 1]], "H9": [["E3", "F5", "B12"], [1, 2, 1]], "A10": [["E3", "F5", "C12"], [1, 2, 1]], "B10": [["E3", "F5", "D12"], [1, 2, 1]], "C10": [["E3", "F5", "E12"], [1, 2, 1]], "D10": [["E3", "F5", "F12"], [1, 2, 1]], "E10": [["E3", "F5", "G12"], [1, 2, 1]], "F10": [["E3", "F5", "H12"], [1, 2, 1]], "G10": [["E3", "G5", "D1"], [1, 2, 2]], "H10": [["E3", "G5", "E1"], [1, 2, 2]], "A11": [["E3", "G5", "F1"], [1, 2, 2]], "B11": [["F3", "G5", "G1"], [1, 2, 2]], "C11": [["F3", "G5", "H1"], [1, 2, 2]], "D11": [["F3", "G5", "A2"], [1, 2, 2]], "E11": [["F3", "G5", "B2"], [1, 2, 2]], "F11": [["F3", "G5", "C2"], [1, 2, 2]], "G11": [["F3", "G5", "D2"], [1, 2, 2]], "H11": [["F3", "H5", "D1"], [1, 2, 2]], "A12": [["F3", "H5", "E1"], [1, 2, 2]], "B12": [["F3", "H5", "F1"], [1, 2, 2]], "C12": [["F3", "H5", "G1"], [1, 2, 2]], "D12": [["F3", "H5", "H1"], [1, 2, 2]], "E12": [["F3", "H5", "A2"], [1, 2, 2]], "F12": [["F3", "H5", "B2"], [1, 2, 2]], "G12": [["G3", "H5", "C2"], [1, 2, 2]], "H12": [["G3", "H5", "D2"], [1, 2, 2]]}
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

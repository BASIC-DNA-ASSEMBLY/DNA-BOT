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
final_assembly_dict={"A1": [["H1", "A1", "F12"], [1, 2, 1]], "B1": [["H1", "A1", "G12"], [1, 2, 1]], "C1": [["H1", "A1", "H12"], [1, 2, 1]], "D1": [["H1", "B1", "H11"], [1, 2, 1]], "E1": [["H1", "B1", "A12"], [1, 2, 1]], "F1": [["H1", "B1", "B12"], [1, 2, 1]], "G1": [["H1", "B1", "C12"], [1, 2, 1]], "H1": [["H1", "B1", "D12"], [1, 2, 1]], "A2": [["A2", "B1", "E12"], [1, 2, 1]], "B2": [["A2", "B1", "F12"], [1, 2, 1]], "C2": [["A2", "B1", "G12"], [1, 2, 1]], "D2": [["A2", "B1", "H12"], [1, 2, 1]], "E2": [["A2", "C1", "D1"], [1, 2, 2]], "F2": [["A2", "C1", "E1"], [1, 2, 2]], "G2": [["A2", "C1", "F1"], [1, 2, 2]], "H2": [["A2", "C1", "G1"], [1, 2, 2]], "A3": [["A2", "C1", "H1"], [1, 2, 2]], "B3": [["A2", "C1", "A2"], [1, 2, 2]], "C3": [["A2", "C1", "B2"], [1, 2, 2]], "D3": [["A2", "C1", "C2"], [1, 2, 2]], "E3": [["A2", "C1", "D2"], [1, 2, 2]], "F3": [["B2", "E2", "D1"], [1, 2, 2]], "G3": [["B2", "E2", "E1"], [1, 2, 2]], "H3": [["B2", "E2", "F1"], [1, 2, 2]], "A4": [["B2", "E2", "G1"], [1, 2, 2]], "B4": [["B2", "E2", "H1"], [1, 2, 2]], "C4": [["B2", "E2", "A2"], [1, 2, 2]], "D4": [["B2", "E2", "B2"], [1, 2, 2]], "E4": [["B2", "E2", "C2"], [1, 2, 2]], "F4": [["B2", "E2", "D2"], [1, 2, 2]], "G4": [["B2", "F2", "D1"], [1, 2, 2]], "H4": [["B2", "F2", "E1"], [1, 2, 2]], "A5": [["B2", "F2", "F1"], [1, 2, 2]], "B5": [["B2", "F2", "G1"], [1, 2, 2]], "C5": [["C2", "F2", "H1"], [1, 2, 2]], "D5": [["C2", "F2", "A2"], [1, 2, 2]], "E5": [["C2", "F2", "B2"], [1, 2, 2]], "F5": [["C2", "F2", "C2"], [1, 2, 2]], "G5": [["C2", "F2", "D2"], [1, 2, 2]], "H5": [["C2", "G2", "H2"], [1, 2, 2]], "A6": [["C2", "G2", "A3"], [1, 2, 2]], "B6": [["C2", "G2", "B3"], [1, 2, 2]], "C6": [["C2", "G2", "C3"], [1, 2, 2]], "D6": [["C2", "G2", "D3"], [1, 2, 2]], "E6": [["C2", "G2", "E3"], [1, 2, 2]], "F6": [["C2", "G2", "F3"], [1, 2, 2]], "G6": [["C2", "G2", "G3"], [1, 2, 2]], "H6": [["D2", "G2", "H3"], [1, 2, 2]], "A7": [["D2", "A4", "H2"], [1, 2, 2]], "B7": [["D2", "A4", "A3"], [1, 2, 2]], "C7": [["D2", "A4", "B3"], [1, 2, 2]], "D7": [["D2", "A4", "C3"], [1, 2, 2]], "E7": [["D2", "A4", "D3"], [1, 2, 2]], "F7": [["D2", "A4", "E3"], [1, 2, 2]], "G7": [["D2", "A4", "F3"], [1, 2, 2]], "H7": [["D2", "A4", "G3"], [1, 2, 2]], "A8": [["D2", "A4", "H3"], [1, 2, 2]], "B8": [["D2", "B4", "H2"], [1, 2, 2]], "C8": [["D2", "B4", "A3"], [1, 2, 2]], "D8": [["D2", "B4", "B3"], [1, 2, 2]], "E8": [["E2", "B4", "C3"], [1, 2, 2]], "F8": [["E2", "B4", "D3"], [1, 2, 2]], "G8": [["E2", "B4", "E3"], [1, 2, 2]], "H8": [["E2", "B4", "F3"], [1, 2, 2]], "A9": [["E2", "B4", "G3"], [1, 2, 2]], "B9": [["E2", "B4", "H3"], [1, 2, 2]], "C9": [["E2", "C4", "D7"], [1, 2, 1]], "D9": [["E2", "C4", "E7"], [1, 2, 1]], "E9": [["E2", "C4", "F7"], [1, 2, 1]], "F9": [["E2", "C4", "G7"], [1, 2, 1]], "G9": [["E2", "C4", "H7"], [1, 2, 1]], "H9": [["E2", "C4", "A8"], [1, 2, 1]], "A10": [["E2", "C4", "B8"], [1, 2, 1]], "B10": [["F2", "C4", "C8"], [1, 2, 1]], "C10": [["F2", "C4", "D8"], [1, 2, 1]], "D10": [["F2", "D4", "D7"], [1, 2, 1]], "E10": [["F2", "D4", "E7"], [1, 2, 1]], "F10": [["F2", "D4", "F7"], [1, 2, 1]], "G10": [["F2", "D4", "G7"], [1, 2, 1]], "H10": [["F2", "D4", "H7"], [1, 2, 1]], "A11": [["F2", "D4", "A8"], [1, 2, 1]], "B11": [["F2", "D4", "B8"], [1, 2, 1]], "C11": [["F2", "D4", "C8"], [1, 2, 1]], "D11": [["F2", "D4", "D8"], [1, 2, 1]], "E11": [["F2", "E4", "D7"], [1, 2, 1]], "F11": [["F2", "E4", "E7"], [1, 2, 1]], "G11": [["G2", "E4", "F7"], [1, 2, 1]], "H11": [["G2", "E4", "G7"], [1, 2, 1]], "A12": [["G2", "E4", "H7"], [1, 2, 1]], "B12": [["G2", "E4", "A8"], [1, 2, 1]], "C12": [["G2", "E4", "B8"], [1, 2, 1]], "D12": [["G2", "E4", "C8"], [1, 2, 1]], "E12": [["G2", "E4", "D8"], [1, 2, 1]], "F12": [["G2", "F4", "H8"], [1, 2, 1]], "G12": [["G2", "F4", "A9"], [1, 2, 1]], "H12": [["G2", "F4", "B9"], [1, 2, 1]]}
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

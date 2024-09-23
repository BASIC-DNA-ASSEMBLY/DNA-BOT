from opentrons import protocol_api
import numpy as np
# metadata
metadata = {
    'apiLevel': '2.19',
    'protocolName': 'DNABOT Step 3: Assembly with thermocycler Gen2',
    'description': 'DNABOT Step 3: Assembly with thermocycler Gen2',
    }
# Construct assemblies are set up on thermocycler module gen2 by combining purified clip parts.

# Test dictionary can be used for simulation 3 or 88 assemblies
final_assembly_dict={
 "A1": ['A7', 'B7', 'C7', 'F7','E7'], 
 "B1": ['A7', 'B7', 'D7', 'G7'], 
 "C1": ['A7', 'B7', 'E7', 'H7']
 }
tiprack_num=1

#final_assembly_dict={"A1": ["A7", "G7", "H7", "A8", "B8"], "B1": ["A7", "D8", "E8", "F8", "G8"], "C1": ["A7", "D8", "H7", "H8", "B9"], "D1": ["A7", "C9", "E9", "G9", "B8"], "E1": ["A7", "H9", "B10", "E9", "D10"], "F1": ["A7", "C9", "H8", "F10", "D10"], "G1": ["A7", "C9", "H10", "E8", "B9"], "H1": ["A7", "H9", "F8", "H10", "B11"], "A2": ["A7", "G7", "E8", "B10", "G8"], "B2": ["A7", "G7", "D11", "A8", "B9"], "C2": ["A7", "C9", "E9", "G9", "B9"], "D2": ["A7", "G7", "H7", "H8", "B8"], "E2": ["A7", "F11", "H11", "H7", "B12"], "F2": ["A7", "C9", "H8", "H11", "D10"], "G2": ["A7", "G7", "D11", "A8", "B8"], "H2": ["B7", "F11", "B10", "H10", "B11"], "A3": ["B7", "D8", "H7", "H8", "B8"], "B3": ["B7", "C9", "H10", "G9", "B8"], "C3": ["B7", "D12", "H8", "H11", "B11"], "D3": ["B7", "D12", "E9", "E8", "B8"], "E3": ["B7", "D12", "E9", "E8", "B9"], "F3": ["B7", "H9", "B10", "H10", "D10"], "G3": ["B7", "G7", "D11", "H8", "B8"], "H3": ["B7", "D12", "H10", "G9", "B9"], "A4": ["B7", "F11", "F10", "D11", "B12"], "B4": ["B7", "G7", "H7", "A8", "B9"], "C4": ["B7", "G7", "E8", "B10", "B12"], "D4": ["B7", "H9", "H11", "H7", "G8"], "E4": ["B7", "D8", "E8", "F8", "B12"], "F4": ["B7", "D12", "E9", "G9", "B8"], "G4": ["C7", "H9", "B10", "E9", "B11"], "H4": ["C7", "F11", "B10", "H10", "D10"], "A5": ["C7", "H9", "F8", "E9", "B11"], "B5": ["C7", "D12", "H8", "F10", "B11"], "C5": ["C7", "F11", "F8", "H10", "B11"], "D5": ["C7", "F11", "H11", "H7", "G8"], "E5": ["C7", "D8", "D11", "A8", "B9"], "F5": ["C7", "H9", "H11", "H7", "B12"], "G5": ["C7", "C9", "H10", "G9", "B9"], "H5": ["C7", "H9", "F10", "H7", "G8"], "A6": ["C7", "D12", "A8", "H11", "D10"], "B6": ["C7", "C9", "A8", "H11", "B11"], "C6": ["C7", "F11", "H11", "D11", "B12"], "D6": ["C7", "D8", "E8", "B10", "G8"], "E6": ["C7", "C9", "H8", "H11", "B11"], "F6": ["D7", "D8", "G9", "F8", "G8"], "G6": ["D7", "C9", "A8", "F10", "B11"], "H6": ["D7", "F11", "F10", "H7", "B12"], "A7": ["D7", "C9", "A8", "F10", "D10"], "B7": ["D7", "H9", "F8", "E9", "D10"], "C7": ["D7", "G7", "G9", "F8", "B12"], "D7": ["D7", "D12", "A8", "H11", "B11"], "E7": ["D7", "D12", "H10", "G9", "B8"], "F7": ["D7", "H9", "H11", "D11", "B12"], "G7": ["D7", "C9", "H8", "F10", "B11"], "H7": ["D7", "D8", "D11", "H8", "B8"], "A8": ["D7", "C9", "E9", "E8", "B9"], "B8": ["D7", "H9", "F10", "D11", "G8"], "C8": ["D7", "H9", "H11", "D11", "G8"], "D8": ["D7", "D12", "A8", "F10", "D10"], "E8": ["E7", "G7", "G9", "F8", "G8"], "F8": ["E7", "D12", "A8", "F10", "B11"], "G8": ["E7", "H9", "F10", "D11", "B12"], "H8": ["E7", "D8", "E8", "B10", "B12"], "A9": ["E7", "C9", "E9", "E8", "B8"], "B9": ["E7", "F11", "B10", "E9", "D10"], "C9": ["E7", "D12", "H8", "F10", "D10"], "D9": ["E7", "H9", "B10", "H10", "B11"], "E9": ["E7", "D8", "G9", "F8", "B12"], "F9": ["E7", "F11", "B10", "E9", "B11"], "G9": ["E7", "F11", "F8", "E9", "C11"], "H9": ["E7", "G7", "G9", "B10", "B12"], "A10": ["E7", "D8", "G9", "B10", "B12"], "B10": ["E7", "D8", "D11", "A8", "B8"], "C10": ["E7", "F11", "F10", "H7", "G8"], "D10": ["F7", "F11", "F8", "E9", "D10"], "E10": ["F7", "H9", "F10", "H7", "B12"], "F10": ["F7", "D12", "H10", "E8", "B9"], "G10": ["F7", "C9", "H10", "E8", "B8"], "H10": ["F7", "F11", "F8", "H10", "D10"], "A11": ["F7", "D12", "H10", "E8", "B8"], "B11": ["F7", "G7", "H7", "H8", "B9"], "C11": ["F7", "G7", "G9", "B10", "G8"], "D11": ["F7", "D12", "H8", "H11", "D10"], "E11": ["F7", "D9", "A8", "H11", "D10"], "F11": ["F7", "G7", "D11", "H8", "B9"], "G11": ["F7", "F11", "A12", "D11", "G8"], "H11": ["F7", "D8", "D11", "A9", "B9"]}
#tiprack_num=5

# __LABWARES is expected to be redefined by "generate_ot2_script" method
# Test dict - generic labware for simulation
__LABWARES={
     "p20_single": {"id": "p20_single_gen2"}, 
     #"p300_multi": {"id": "p300_multi_gen2"}, 
     #"mag_deck": {"id": "magdeck"},
     "clip_plate":{"id":"biorad_96_wellplate_200ul_pcr"},
     "final_assembly_plate":{"id":"biorad_96_wellplate_200ul_pcr"},
     "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, 
     #"96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, 
     "24_tuberack_2000ul": {"id": "opentrons_24_tuberack_generic_2ml_screwcap"}, 
     #"96_wellplate_200ul_pcr_step_14": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"96_wellplate_200ul_pcr_step_23": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"agar_plate_step_4": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"12_reservoir_21000ul": {"id": "nest_12_reservoir_15ml"}, 
     #"96_deepwellplate_2ml": {"id": "nest_96_wellplate_2ml_deep"}
     #corning_12_wellplate_6.9ml_flat
     }

def run(protocol: protocol_api.ProtocolContext):

    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type=__LABWARES['96_tiprack_20ul']['id']):
        
            # Constants, we update all the labware name in version 2
            #Tiprack
            CANDIDATE_TIPRACK_SLOTS = ['2', '3', '5', '6', '9']
            PIPETTE_MOUNT = 'right'
            #Plate of sample after  purification
            CLIP_PLATE_TYPE = __LABWARES['clip_plate']['id']
            CLIP_PLATE_POSITION = '1'
            #Tuberack
            TUBE_RACK_TYPE = __LABWARES['24_tuberack_2000ul']['id']
            TUBE_RACK_POSITION = '4'
            #Destination plate
            DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
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
            pipette = protocol.load_instrument(__LABWARES['p20_single']['id'], PIPETTE_MOUNT, tip_racks=tipracks)

            # Define Labware and set temperature
            purified_clip_plate = protocol.load_labware(CLIP_PLATE_TYPE, CLIP_PLATE_POSITION)
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
                  
            #thermocycler module gen2
            tc_mod = protocol.load_module(module_name="thermocyclerModuleV2")
            destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
            tc_mod.open_lid()
            tc_mod.deactivate_lid()
            tc_mod.set_block_temperature(4)

             # Master mix transfers
            final_assembly_lengths = []
            for values in final_assembly_dict.values():
                final_assembly_lengths.append(len(values))
            unique_assemblies_lengths = list(set(final_assembly_lengths))
            master_mix_well_letters = ['A', 'B', 'C', 'D']

            for x in unique_assemblies_lengths:
                master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
                destination_inds = [i for i, lengths in enumerate(final_assembly_lengths) if lengths == x]
                destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
                destination_wells = list(destination_wells[destination_inds])
                
                pipette.flow_rate.aspirate = 6
                pipette.flow_rate.dispense = 6
                pipette.flow_rate.blow_out = 15
                high = 2
                normal = 1
                slow = 0.5
                vslow = 0.2
                pipette.well_bottom_clearance.aspirate = 1 
                pipette.well_bottom_clearance.dispense = 2

                pipette.pick_up_tip()
                for destination_well in destination_wells:# make tube_rack_wells and destination_plate.wells in the same type  
                    pipette.distribute(TOTAL_VOL - x * PART_VOL, tube_rack[master_mix_well], destination_plate[destination_well],blow_out=True, blowout_location="source well", new_tip='never')
                pipette.drop_tip()

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for value in values:# purified_clip_plate.wells and destination_plate.wells in the same type
                    #pipette.transfer(PART_VOL, purified_clip_plate.wells(value), destination_plate.wells(key), mix_after=MIX_SETTINGS, new_tip='always')#transfer parts in one tube
                    pipette.pick_up_tip()
                    pipette.well_bottom_clearance.aspirate = 1  # tip is 2 mm above well bottom
                    pipette.well_bottom_clearance.dispense = 2  # tip is 2 mm above well bottom
                    #Prefix Transfer
                    pipette.aspirate(PART_VOL, purified_clip_plate[value].bottom(1), rate=slow)
                    pipette.dispense(PART_VOL, destination_plate[key].bottom(2), rate=slow)
                    #mix after transfer
                    pipette.aspirate(2, destination_plate[key].bottom(1), rate=normal)
                    pipette.dispense(2, destination_plate[key].bottom(3), rate=high)
                    pipette.aspirate(3, destination_plate[key].bottom(2), rate=normal)
                    pipette.dispense(3, destination_plate[key].bottom(1), rate=normal)
                    pipette.aspirate(4, destination_plate[key].bottom(2), rate=slow)
                    pipette.dispense(4, destination_plate[key].bottom(3), push_out=0.5, rate=vslow)
                    pipette.move_to(destination_plate[key].top(-8))
                    pipette.blow_out()
                    pipette.touch_tip(radius=0.6, v_offset=-8, speed=10)
                    pipette.drop_tip()

            #thermocycler module gen2
            tc_mod.close_lid()
            tc_mod.set_lid_temperature(105)
            tc_mod.set_block_temperature(50, hold_time_minutes=45, block_max_volume=15)
            tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=30)
            # Increase the hold time at 4 C if necessary
            tc_mod.set_lid_temperature(37)
            protocol.delay(seconds=120)
            tc_mod.deactivate_lid()
            tc_mod.open_lid()
            tc_mod.set_block_temperature(4)
            #for line in protocol.commands(): 
                #print(line)

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)
    
    #output command actions in simulate
    for line in protocol.commands(): 
       print(line)
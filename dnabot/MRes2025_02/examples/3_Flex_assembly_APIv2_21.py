from opentrons import protocol_api
import numpy as np

requirements = {"robotType": "Flex"}

final_assembly_dict={"A1": ["A7", "B7", "C7", "D7", "E7"], "B1": ["A7", "B7", "C7", "D7", "F7"], "C1": ["A7", "B7", "C7", "G7", "E7"], "D1": ["A7", "B7", "C7", "G7", "F7"], "E1": ["A7", "H7", "D7", "A8", "B8"], "F1": ["A7", "H7", "D7", "C8", "B8"], "G1": ["A7", "H7", "G7", "A8", "B8"], "H1": ["A7", "H7", "G7", "C8", "B8"], "A2": ["A7", "B7", "D8", "D7", "E7"]}
tiprack_num=1
__HARDWARE={"robot_type": {"id": "Flex"}, "single_pipette": {"id": "Flex_1channel_50"}, "single_pipette_mount": {"id": "right"}, "multi_pipette": {"id": "Flex_8channel_1000"}, "multi_pipette_mount": {"id": "left"}, "thermocycler": {"id": "thermocyclerModuleV2"}, "mag_deck": {"id": "magneticBlockV1"}}
__LABWARES={"tiprack_20ul": {"id": "opentrons_flex_96_tiprack_50ul"}, "tiprack_300ul": {"id": "opentrons_flex_96_tiprack_200ul"}, "flex_96_tiprack_200ul": {"id": "opentrons_flex_96_tiprack_200ul"}, "flex_96_tiprack_1000ul": {"id": "opentrons_flex_96_tiprack_1000ul"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "clip_source_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "clip_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "mix_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "final_assembly_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "transform_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "agar_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "mag_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}, "12_corning_wellplate": {"id": "corning_12_wellplate_6.9ml_flat"}}
__PARAMETERS={"clip_keep_thermo_lid_closed": {"value": "No"}, "premix_linkers": {"value": "Yes"}, "premix_parts": {"value": "Yes"}, "linkers_volume": {"value": 20}, "parts_volume": {"value": 20}, "thermo_temp": {"value": 4}, "purif_magdeck_height": {"value": 10.8}, "purif_wash_time": {"value": 0.5}, "purif_bead_ratio": {"value": 1.8}, "purif_incubation_time": {"value": 5}, "purif_settling_time": {"value": 2}, "purif_drying_time": {"value": 5}, "purif_elution_time": {"value": 2}, "transform_incubation_temp": {"value": 4}, "transform_incubation_time": {"value": 20}}

# metadata
metadata = {'protocolName': 'DNABOT Step 3: Assembly with thermocycler Gen2', "apiLevel": "2.21"}
# Construct assemblies are set up on thermocycler module gen2 by combining purified clip parts.

# Test dictionary can be used for simulation 3 or 88 assemblies
'''final_assembly_dict={
 "A1": ['A7', 'B7', 'C7', 'F7','E7'], 
 "B1": ['A7', 'B7', 'D7', 'G7'], 
 "C1": ['A7', 'B7', 'E7', 'H7']
 }
tiprack_num=1'''

#final_assembly_dict={"A1": ["A7", "G7", "H7", "A8", "B8"], "B1": ["A7", "D8", "E8", "F8", "G8"], "C1": ["A7", "D8", "H7", "H8", "B9"], "D1": ["A7", "C9", "E9", "G9", "B8"], "E1": ["A7", "H9", "B10", "E9", "D10"], "F1": ["A7", "C9", "H8", "F10", "D10"], "G1": ["A7", "C9", "H10", "E8", "B9"], "H1": ["A7", "H9", "F8", "H10", "B11"], "A2": ["A7", "G7", "E8", "B10", "G8"], "B2": ["A7", "G7", "D11", "A8", "B9"], "C2": ["A7", "C9", "E9", "G9", "B9"], "D2": ["A7", "G7", "H7", "H8", "B8"], "E2": ["A7", "F11", "H11", "H7", "B12"], "F2": ["A7", "C9", "H8", "H11", "D10"], "G2": ["A7", "G7", "D11", "A8", "B8"], "H2": ["B7", "F11", "B10", "H10", "B11"], "A3": ["B7", "D8", "H7", "H8", "B8"], "B3": ["B7", "C9", "H10", "G9", "B8"], "C3": ["B7", "D12", "H8", "H11", "B11"], "D3": ["B7", "D12", "E9", "E8", "B8"], "E3": ["B7", "D12", "E9", "E8", "B9"], "F3": ["B7", "H9", "B10", "H10", "D10"], "G3": ["B7", "G7", "D11", "H8", "B8"], "H3": ["B7", "D12", "H10", "G9", "B9"], "A4": ["B7", "F11", "F10", "D11", "B12"], "B4": ["B7", "G7", "H7", "A8", "B9"], "C4": ["B7", "G7", "E8", "B10", "B12"], "D4": ["B7", "H9", "H11", "H7", "G8"], "E4": ["B7", "D8", "E8", "F8", "B12"], "F4": ["B7", "D12", "E9", "G9", "B8"], "G4": ["C7", "H9", "B10", "E9", "B11"], "H4": ["C7", "F11", "B10", "H10", "D10"], "A5": ["C7", "H9", "F8", "E9", "B11"], "B5": ["C7", "D12", "H8", "F10", "B11"], "C5": ["C7", "F11", "F8", "H10", "B11"], "D5": ["C7", "F11", "H11", "H7", "G8"], "E5": ["C7", "D8", "D11", "A8", "B9"], "F5": ["C7", "H9", "H11", "H7", "B12"], "G5": ["C7", "C9", "H10", "G9", "B9"], "H5": ["C7", "H9", "F10", "H7", "G8"], "A6": ["C7", "D12", "A8", "H11", "D10"], "B6": ["C7", "C9", "A8", "H11", "B11"], "C6": ["C7", "F11", "H11", "D11", "B12"], "D6": ["C7", "D8", "E8", "B10", "G8"], "E6": ["C7", "C9", "H8", "H11", "B11"], "F6": ["D7", "D8", "G9", "F8", "G8"], "G6": ["D7", "C9", "A8", "F10", "B11"], "H6": ["D7", "F11", "F10", "H7", "B12"], "A7": ["D7", "C9", "A8", "F10", "D10"], "B7": ["D7", "H9", "F8", "E9", "D10"], "C7": ["D7", "G7", "G9", "F8", "B12"], "D7": ["D7", "D12", "A8", "H11", "B11"], "E7": ["D7", "D12", "H10", "G9", "B8"], "F7": ["D7", "H9", "H11", "D11", "B12"], "G7": ["D7", "C9", "H8", "F10", "B11"], "H7": ["D7", "D8", "D11", "H8", "B8"], "A8": ["D7", "C9", "E9", "E8", "B9"], "B8": ["D7", "H9", "F10", "D11", "G8"], "C8": ["D7", "H9", "H11", "D11", "G8"], "D8": ["D7", "D12", "A8", "F10", "D10"], "E8": ["E7", "G7", "G9", "F8", "G8"], "F8": ["E7", "D12", "A8", "F10", "B11"], "G8": ["E7", "H9", "F10", "D11", "B12"], "H8": ["E7", "D8", "E8", "B10", "B12"], "A9": ["E7", "C9", "E9", "E8", "B8"], "B9": ["E7", "F11", "B10", "E9", "D10"], "C9": ["E7", "D12", "H8", "F10", "D10"], "D9": ["E7", "H9", "B10", "H10", "B11"], "E9": ["E7", "D8", "G9", "F8", "B12"], "F9": ["E7", "F11", "B10", "E9", "B11"], "G9": ["E7", "F11", "F8", "E9", "C11"], "H9": ["E7", "G7", "G9", "B10", "B12"], "A10": ["E7", "D8", "G9", "B10", "B12"], "B10": ["E7", "D8", "D11", "A8", "B8"], "C10": ["E7", "F11", "F10", "H7", "G8"], "D10": ["F7", "F11", "F8", "E9", "D10"], "E10": ["F7", "H9", "F10", "H7", "B12"], "F10": ["F7", "D12", "H10", "E8", "B9"], "G10": ["F7", "C9", "H10", "E8", "B8"], "H10": ["F7", "F11", "F8", "H10", "D10"], "A11": ["F7", "D12", "H10", "E8", "B8"], "B11": ["F7", "G7", "H7", "H8", "B9"], "C11": ["F7", "G7", "G9", "B10", "G8"], "D11": ["F7", "D12", "H8", "H11", "D10"], "E11": ["F7", "D9", "A8", "H11", "D10"], "F11": ["F7", "G7", "D11", "H8", "B9"], "G11": ["F7", "F11", "A12", "D11", "G8"], "H11": ["F7", "D8", "D11", "A9", "B9"]}
#tiprack_num=5

# __LABWARES is expected to be redefined by "generate_ot2_script" method
# Test dict - generic labware for simulation
'''__LABWARES={
     "p50_single": {"id": "flex_1channel_50"}, 
     #"p300_multi": {"id": "p300_multi_gen2"}, 
     #"mag_deck": {"id": "magdeck"},
     "clip_plate":{"id":"biorad_96_wellplate_200ul_pcr"},
     "final_assembly_plate":{"id":"biorad_96_wellplate_200ul_pcr"},
     #"96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, 
     "flex_96_tiprack_50ul": {"id": "opentrons_flex_96_tiprack_50ul"}, 
     #"96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, 
     "24_tuberack_2000ul": {"id": "opentrons_24_tuberack_generic_2ml_screwcap"}, 
     #"96_wellplate_200ul_pcr_step_14": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"96_wellplate_200ul_pcr_step_23": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"agar_plate_step_4": {"id": "biorad_96_wellplate_200ul_pcr"}, 
     #"12_reservoir_21000ul": {"id": "nest_12_reservoir_15ml"}, 
     #"96_deepwellplate_2ml": {"id": "nest_96_wellplate_2ml_deep"}
     #corning_12_wellplate_6.9ml_flat
     }'''

# final_assembly_dict={"A1": ["A7", "B7", "C7", "D7", "E7"], "B1": ["A7", "B7", "C7", "D7", "E7"], "C1": ["A7", "B7", "C7", "F7"], "D1": ["A7", "B7", "C7", "F7"]}
# tiprack_num=1
# __LABWARES={"p50_single": {"id": "flex_1channel_50"},
#             "p50_multi": {"id": "flex_8channel_50"}, 
#             "p1000_single": {"id": "flex_1channel_1000"}, 
#             "p1000_multi": {"id": "flex_8channel_1000"},
#             "mag_block": {"id": "magneticBlockV1"},
#             "mag_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "flex_96_tiprack_50ul": {"id": "opentrons_flex_96_tiprack_50ul"}, 
#             "flex_96_tiprack_200ul": {"id": "opentrons_flex_96_tiprack_200ul"},
#             "flex_96_tiprack_1000ul": {"id": "opentrons_flex_96_tiprack_1000ul"}, 
#             "24_tuberack_1500ul": {"id": "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap"}, 
#             "clip_source_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"},
#             "clip_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "mix_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "final_assembly_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "transfo_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "transfo_plate_wo_thermo": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "agar_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
#             "12_reservoir_21000ul": {"id": "nest_12_reservoir_15ml"}, 
#             "96_deepwellplate_2ml": {"id": "nest_96_wellplate_2ml_deep"}, 
#             "12_corning_wellplate": {"id": "corning_12_wellplate_6.9ml_flat"}}

def run(protocol: protocol_api.ProtocolContext):

    robot_type=__HARDWARE['robot_type']['id']
    if robot_type=='Flex':
        trash = protocol.load_trash_bin("A3")
        tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'], location = "B1")
        tiprack_type = ['opentrons_flex_96_tiprack_50ul']
    elif robot_type=='OT-2':
        tiprack_type = __LABWARES['tiprack_20ul']['id']
    else:
        raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")

       
    
    def final_assembly(robot_type, final_assembly_dict, tiprack_num, tiprack_type):
        
        # Constants, we update all the labware name in version 2
        #Tiprack
        #CANDIDATE_TIPRACK_SLOTS = ['2', '3', '5', '6', '9']
        # CANDIDATE_TIPRACK_SLOTS = ['D2', 'D3', 'C2', 'C3', 'B3']
        # PIPETTE_MOUNT = 'right'
        #Plate of sample after  purification
        
        # OLD VERSION
        # if robot_type=='OT-2':
        #     CLIP_PLATE_TYPE = __LABWARES['clip_plate']['id']
        #     CLIP_PLATE_POSITION = '1'
        # elif robot_type=='Flex':
        #     CLIP_PLATE_TYPE = __LABWARES['clip_plate']['id']
        #     CLIP_PLATE_POSITION = 'D1'
        # #Tuberack
        #     if robot_type=='OT-2':
        #     TUBE_RACK_TYPE = __LABWARES['24_tuberack_1500ul']['id']
        #     TUBE_RACK_POSITION = '4'
        # elif robot_type=='Flex':
        #     TUBE_RACK_TYPE = __LABWARES['24_tuberack_1500ul']['id']
        #     TUBE_RACK_POSITION = 'C1'           
        #     if robot_type=='OT-2':
        #     DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
        #     TUBE_RACK_POSITION = '4'
        # elif robot_type=='Flex':
        #     DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
        #     TUBE_RACK_POSITION = 'C1'  
        # #Destination plate
        # DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
        
        # ADDED VERSION
        if robot_type=='OT-2':
            CLIP_PLATE_TYPE = __LABWARES['clip_plate']['id']
            CLIP_PLATE_POSITION = '1'
        elif robot_type=='Flex':
            CLIP_PLATE_TYPE = __LABWARES['clip_plate']['id']
            CLIP_PLATE_POSITION = 'D1'
        #Tuberack
        if robot_type=='OT-2':
            TUBE_RACK_TYPE = __LABWARES['24_tuberack_1500ul']['id']
            TUBE_RACK_POSITION = '4'
        elif robot_type=='Flex':
            TUBE_RACK_TYPE = __LABWARES['24_tuberack_1500ul']['id']
            TUBE_RACK_POSITION = 'C1'           
        if robot_type=='OT-2':
            DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
            DESTINATION_RACK_POSITION = '2'                                     #originally named TUBE_RACK_POSITION but changed due to clashing on names and positions
        elif robot_type=='Flex':
            DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
            DESTINATION_RACK_POSITION = 'C2'                                    #originally named TUBE_RACK_POSITION but changed due to clashing on names and positions
        
        #Destination plate
        DESTINATION_PLATE_TYPE = __LABWARES['final_assembly_plate']['id']
        # After defining *_TYPE and *_POSITION variables
        clip_plate = protocol.load_labware(CLIP_PLATE_TYPE, CLIP_PLATE_POSITION)
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        destination_plate = protocol.load_labware(DESTINATION_PLATE_TYPE, DESTINATION_RACK_POSITION)  # IS TUBE RACK POSITION THE SAME AS DESTINATION PLATE POSITION??????????
        purified_clip_plate = clip_plate  # or load another if needed


        TOTAL_VOL = 15
        PART_VOL = 1.5
        MIX_SETTINGS = (1, 3)
        tiprack_num=tiprack_num+1
        # Errors
        sample_number = len(final_assembly_dict.keys())
        if sample_number > 96:
            raise ValueError('Final assembly nummber cannot exceed 96.')

        # Constants
        INITIAL_TIP = 'A1'
        # Candidate Tiprack Slots according to robot type
        if robot_type=='OT-2':
            CANDIDATE_TIPRACK_SLOTS = ['3', '5', '6', '9']    #'2' removed as being used as destination plate position
        elif robot_type=='Flex':
            CANDIDATE_TIPRACK_SLOTS = ['D2', 'D3', 'C3', 'B3']    #'C2' see above
        else:
            raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")
        PIPETTE_TYPE = __HARDWARE['single_pipette']['id']
        PIPETTE_MOUNT = __HARDWARE['single_pipette_mount']['id']
        
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        
        pipette = protocol.load_instrument(PIPETTE_TYPE, PIPETTE_MOUNT, tip_racks=tipracks)
        if robot_type=='Flex':
            if(PIPETTE_TYPE)=="flex_1channel_50":
                pipette.flow_rate.aspirate = 20
                pipette.flow_rate.dispense = 20
                pipette.flow_rate.blow_out = 35
            elif robot_type=='OT2':            
                if(PIPETTE_TYPE)=="p20_single_gen2":
                    pipette.flow_rate.aspirate = 8
                    pipette.flow_rate.dispense = 8
                    pipette.flow_rate.blow_out = 15
                else: 
                    raise ValueError("Don't have a single-channel P20 or P50 pipette loaded")
    # relative rates for fine-tuning pipetting steps
        high = 2
        normal = 1
        slow = 0.4
        vslow = 0.2
        # Define thermocycler and set temperature
        if robot_type=='Flex':
            tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'], location = "B1")
        else:
            tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'])
        tc_mod.open_lid()
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(temperature=4)         

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
                pipette.aspirate(10, destination_plate[key].bottom(1), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(3), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(2), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(1), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(3), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(1), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(2), rate=slow)
                pipette.dispense(10, destination_plate[key].bottom(3), push_out=0.5, rate=vslow)
                protocol.delay(seconds=5)     #changed from (5 seconds) to seconds = 5
                pipette.move_to(destination_plate[key].top(-8))
                pipette.blow_out()
                pipette.touch_tip(radius=0.6, v_offset=-8, speed=10)
                pipette.drop_tip()

            #thermocycler module gen2
            tc_mod.close_lid()
            tc_mod.set_lid_temperature(105)
            tc_mod.set_block_temperature(72, hold_time_minutes=1)
            tc_mod.set_block_temperature(70, hold_time_minutes=1)
            tc_mod.set_block_temperature(68, hold_time_minutes=1)
            tc_mod.set_block_temperature(66, hold_time_minutes=1)
            tc_mod.set_block_temperature(64, hold_time_minutes=2)
            tc_mod.set_block_temperature(62, hold_time_minutes=2)
            tc_mod.set_block_temperature(60, hold_time_minutes=2)
            tc_mod.set_block_temperature(58, hold_time_minutes=2)
            tc_mod.set_block_temperature(56, hold_time_minutes=2)
            tc_mod.set_block_temperature(54, hold_time_minutes=2)
            tc_mod.set_block_temperature(52, hold_time_minutes=2)
            tc_mod.set_block_temperature(50, hold_time_minutes=2)
            tc_mod.set_block_temperature(48, hold_time_minutes=2)
            tc_mod.set_block_temperature(46, hold_time_minutes=2)
            tc_mod.set_block_temperature(44, hold_time_minutes=1)
            tc_mod.set_block_temperature(42, hold_time_minutes=1)
            tc_mod.set_block_temperature(40, hold_time_minutes=1)
            tc_mod.set_block_temperature(4)
            # Increase the hold time at 4 C if necessary
            tc_mod.set_lid_temperature(37)
            protocol.delay(minutes=2)
            tc_mod.deactivate_lid()
            tc_mod.open_lid()
            tc_mod.set_block_temperature(4)
            #for line in protocol.commands(): 
                #print(line)

    #final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)
    final_assembly(robot_type=robot_type, final_assembly_dict= final_assembly_dict,tiprack_num= tiprack_num,tiprack_type= tiprack_type)
    
    #output command actions in simulate
    for line in protocol.commands(): 
       print(line)

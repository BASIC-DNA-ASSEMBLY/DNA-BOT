from opentrons import protocol_api
        
# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.8',
     'protocolName': 'CLIP_With_Thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module.'}

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template

clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["5", "5", "5", "5", "5"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}




def run(protocol: protocol_api.ProtocolContext):
# added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity
        
    #Tiprack
    tiprack_type="opentrons_96_tiprack_10ul"
    INITIAL_TIP = 'A1'
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9']
        
    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = 'p20_single_gen2'
             # API 2 supports gen_1 pipettes like the p10_single
    PIPETTE_MOUNT = 'right'
        ### Load Pipette        
        # checks if it's a P10 Single pipette
    if PIPETTE_TYPE != 'p20_single_gen2':
        print('Define labware must be changed to use', PIPETTE_TYPE)
        exit()           
#Thermocycler Module
    tc_mod = protocol.load_module('Thermocycler Module')
# Destination Plates    
    DESTINATION_PLATE_TYPE = '4ti_96_wellplate_200ul'
    # Loads destination plate onto Thermocycler Module
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
        # not supported by API 1
    # DESTINATION_PLATE_POSITION removed, as it is on the thermocylcer module (which is always on the same slot)
    # DESTINATION_PLATE_POSITION = '1'   
    # INITIAL_DESTINATION_WELL constant removed, as destination_plate.wells() automatically starts from A1
          
    # Source Plates
    SOURCE_PLATE_TYPE = '4ti_96_wellplate_200ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
         
    # Tube Rack
    TUBE_RACK_TYPE = 'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
    TUBE_RACK_POSITION = '4'
    MASTER_MIX_WELL = 'A1'
    WATER_WELL = 'A2'
    MASTER_MIX_VOLUME = 20
        
    # Mix settings
    LINKER_MIX_SETTINGS = (1, 3)
    PART_MIX_SETTINGS = (4, 5)
    
    def clip(
            prefixes_wells,
            prefixes_plates,
            suffixes_wells,
            suffixes_plates,
            parts_wells,
            parts_plates,
            parts_vols,
            water_vols):
        
        ### Loading Tiprack        
        # Calculates whether one, two, or three tipracks are needed, which are in slots 3, 6, and 9 respectively
        total_tips = 4 * len(parts_wells)
        letter_dict = {'A': 0, 'B': 1, 'C': 2,
                       'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
        tiprack_1_tips = (
            13 - int(INITIAL_TIP[1:])) * 8 - letter_dict[INITIAL_TIP[0]]
        if total_tips > tiprack_1_tips:
            tiprack_num = 1 + (total_tips - tiprack_1_tips) // 96 + \
            (1 if (total_tips - tiprack_1_tips) % 96 > 0 else 0)
        else:
            tiprack_num = 1
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        
        # loads the correct number of tipracks
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            # changed to protocol.load_labware for API 2.8
        
        # Loads pipette according to constants assigned above
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)
            # changed to protocol.load_labware for API 2.8
            # removed 'pipette.start_at_tip(tipracks[0].well(INITIAL_TIP))'
                # start_at_tip supported by API v1 up to API v2.7, but returns error for 2.8:
                    # 'InstrumentContext' object has no attribute 'start_at_tip'
            # Because of this, ability to specify INITIAL_TIP was removed. Possible future improvement.
        
        # Defines where the destination wells are within the destination plate
        destination_wells = destination_plate.wells()[0:len(parts_wells)]
      
            # old code:
                # destination_wells = destination_plate.wells(INITIAL_DESTINATION_WELL, length=int(len(parts_wells)))
            # For API 2.8 and above, '.wells' will no longer take length arguments
            # Therefore the length arguement replaced by '[0:len(parts_wells)]'
            # Because of this, ability to specify INITIAL_DESTINATION_WELL was removed. Possible future improvement.
        
        ### Load Tube Rack        
        # Loads tube rack according to constants assigned above
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
            # changed to protocol.load_labware for API 2.8
            
        # Defines positions of master mix and water within the tube rack
        master_mix = tube_rack.wells(MASTER_MIX_WELL)
        water = tube_rack.wells(WATER_WELL)
   
         ### Loading Source Plates
        # Makes a source plate key for where prefixes, suffixes, and parts are located, according to the dictionary generated by the DNA-BOT
        source_plates = {}
        source_plates_keys = list(set((prefixes_plates + suffixes_plates + parts_plates)))
        
        # Loads plates according to the source plate key
        for key in source_plates_keys:
            source_plates[key]=protocol.load_labware(SOURCE_PLATE_TYPE, key)
                # changed to protocol.load_labware for API 2.8      
    
        ### Transfers
        
        # transfer master mix into destination wells
            # added blowout into destination wells ('blowout_location' only works for API 2.8 and above)
        pipette.pick_up_tip()
        pipette.transfer(MASTER_MIX_VOLUME, master_mix, destination_wells, blow_out=True, blowout_location='destination well', new_tip='never')
        pipette.drop_tip()
        
        # transfer water into destination wells
            # added blowout into destination wells ('blowout_location' only works for API 2.8 and above)
        pipette.transfer(water_vols, water, destination_wells, blow_out=True, blowout_location='destination well', new_tip='always')
        
        #transfer prefixes, suffixes, and parts into destination wells
            # added blowout into destination wells ('blowout_location' only works for API 2.8 and above)
        for clip_num in range(len(parts_wells)):
            pipette.transfer(1, source_plates[prefixes_plates[clip_num]].wells(prefixes_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(1, source_plates[suffixes_plates[clip_num]].wells(suffixes_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(parts_vols[clip_num], source_plates[parts_plates[clip_num]].wells(parts_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS)
       
    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(**clips_dict)
    ### PCR Reaction in Thermocycler
        
    # close lid and set lid temperature, PCR will not start until lid reaches 37C
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)
        
    # Runs 20 cycles of 37C for 2 minutes and 20C for 1 minute, then holds for 60C for 10 minutes
    profile = [
        {'temperature': 37, 'hold_time_minutes': 2},
        {'temperature': 20, 'hold_time_minutes': 1}]
    tc_mod.execute_profile(steps=profile, repetitions=20, block_max_volume=30)
    tc_mod.set_block_temperature(60, hold_time_minutes=10, block_max_volume=30)
    tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=30)
    #Q Does block_max_volume define total volume in block or individual wells?   
    tc_mod.set_lid_temperature(37)
    tc_mod.open_lid()


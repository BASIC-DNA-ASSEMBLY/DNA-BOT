from opentrons import protocol_api
#from mix_functions import mix_linkers_function, mix_parts_function
import numpy as np
from opentrons.types import Point

# gripper module for Flex is currently not included - this would only impact script 3 purification.

#metadata
metadata = {
     'protocolName': 'DNABOT Step 1: Clip Reaction with thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons Flex, including the thermocycler module gen1 or gen2.', "apiLevel": "2.21",
}

DECK_LAYOUT = {
    'OT-2': {
        'candidate_tiprack_slots': ['3', '6', '9'],
    },
    'Flex': {
        'candidate_tiprack_slots': ['D3', 'C3', 'B3'],
    },
}

def run(protocol: protocol_api.ProtocolContext):
    robot_type=__HARDWARE['robot_type']['id']
    layout = DECK_LAYOUT[robot_type]
    if robot_type=='Flex':
        trash = protocol.load_trash_bin("A3")
        #tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'], location = "B1")
        tiprack_type = __LABWARES['tiprack_20ul']['id'] # CHANGE, CALLED ID AND NOT KEY!!
        #tiprack_1000 = ['Flex_tiprack_1000ul'] 'opentrons_flex_96_tiprack_1000ul'
    elif robot_type=='OT-2':
        tiprack_type = __LABWARES['tiprack_20ul']['id']
    else:
        raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")
    # Constants
    INITIAL_TIP = 'A1'
    # Candidate Tiprack Slots according to robot type
    CANDIDATE_TIPRACK_SLOTS = layout['candidate_tiprack_slots']
    
    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = __HARDWARE['single_pipette']['id']
    #PIPETTE_TYPE = 'flex_1channel_50'
    PIPETTE_MOUNT = __HARDWARE['single_pipette_mount']['id']

    if robot_type=='Flex':
        tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'], location = "B1")
    else:
        tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'])
    tc_mod.open_lid()
    tc_mod.deactivate_lid()
    tc_mod.set_block_temperature(temperature=__PARAMETERS['thermo_temp']['value']) 
    # Destination Plates
    DESTINATION_PLATE_TYPE = __LABWARES['clip_plate']['id']
    # Loads destination plate onto thermocycler module gen2
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)

    # Source Plates
    SOURCE_PLATE_TYPE = __LABWARES['clip_source_plate']['id']
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used

    # Tube Rack
    TUBE_RACK_TYPE = __LABWARES['24_tuberack_1500ul']['id']
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
    #TUBE_RACK_POSITION = '4'
    TUBE_RACK_POSITION = '1'
    MASTER_MIX_WELL = 'A1'
    WATER_WELL = 'A2'
    MASTER_MIX_VOLUME = 20

    #translates the __PARAMETER setting for mix into a boolean for the mix functions
    if __PARAMETERS['premix_linkers']['value']=='Yes':
        Mix_linkers_bool=True
    else:
        Mix_linkers_bool = False
    
    if __PARAMETERS['premix_parts']['value']=='Yes':
        Mix_parts_bool=True
    else:
        Mix_parts_bool = False

    def configure_flex_low_volume_mode(pipette, volume):
        """Configure Flex 50 uL pipettes for the next transfer volume when supported."""
        if robot_type != 'Flex':
            return
        if not hasattr(pipette, 'configure_for_volume'):
            return
        if getattr(pipette, 'max_volume', None) != 50:
            return
        pipette.configure_for_volume(max(1, min(float(volume), 50)))

    def get_low_volume_push_out(volume):
        """Match Flex low-volume defaults explicitly and give OT-2 a defined push-out too."""
        return 7 if float(volume) < 5 else 2

    def get_mix_repetitions(volume):
        """Scale mix repetitions with premix volume."""
        if volume <= 20:
            return 4
        if volume <= 40:
            return 5
        if volume <= 60:
            return 6
        return 8

    def get_premix_dispense_height(volume):
        """Choose a conservative dispense height based on available liquid volume."""
        if volume <= 20:
            return 2
        if volume <= 40:
            return 3
        if volume <= 60:
            return 4
        return 5

    def get_premix_stroke_volume(volume):
        """Derive a premix stroke from the working volume, capped for 20 uL tips."""
        return min(volume * 0.8, 20)

    def perform_premix(pipette, well, volume):
        """Premix with slow aspirates, fast dispenses, then a slower finishing pass."""
        high = 2
        slow = 0.4
        aspirate_height = 1
        dispense_height = get_premix_dispense_height(volume)
        stroke_volume = get_premix_stroke_volume(volume)

        mix_reps = get_mix_repetitions(volume)
        for mix_step in range(mix_reps):
            x_offset = 2 if mix_step % 2 == 0 else -2
            dispense_location = well.bottom(dispense_height).move(Point(x=x_offset, y=0, z=0))
            pipette.aspirate(stroke_volume, well.bottom(aspirate_height), rate=slow)
            pipette.dispense(stroke_volume, dispense_location, rate=high)

        final_dispense_height = max(dispense_height, np.log(volume))
        pipette.aspirate(stroke_volume, well.bottom(aspirate_height), rate=slow)
        pipette.dispense(
            stroke_volume,
            well.bottom(final_dispense_height).move(Point(x=2, y=0, z=0)),
            rate=slow,
            push_out=3,
        )
        pipette.move_to(well.top(-5))
        protocol.delay(seconds=1)
        pipette.blow_out()
        pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
        pipette.drop_tip()

    def perform_clip_mix(pipette, well, mix_volume=20, repetitions=4):
        """Mix a clip reaction well with fixed repeated cycles and a slow finish."""
        high = 2
        normal = 1
        slow = 0.4
        aspirate_volume = mix_volume

        for mix_step in range(repetitions):
            x_offset = 2 if mix_step % 2 == 0 else -2
            dispense_location = well.bottom(3).move(Point(x=x_offset, y=0, z=0))
            pipette.aspirate(aspirate_volume, well.bottom(1), rate=normal)
            pipette.dispense(aspirate_volume, dispense_location, rate=high)

        pipette.aspirate(aspirate_volume, well.bottom(2), rate=slow)
        pipette.dispense(
            aspirate_volume,
            dispense_location,
            rate=slow,
            push_out=3,
        )
        pipette.move_to(well.top(-5))
        protocol.delay(seconds=1)
        pipette.blow_out()
        pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
        pipette.drop_tip()

    def pre_wet(pipette, well, repetitions=1):
        """Pre-wet the tip in the source well before low-volume aspiration."""
        pre_wet_volume = 2
        for _ in range(repetitions):
            pipette.aspirate(pre_wet_volume, well.bottom(1), rate=slow)
            protocol.delay(seconds=0.5)
            pipette.dispense(pre_wet_volume, well.bottom(1).move(Point(x=0, y=0, z=0)), rate=slow)

    def pick_up_tip_with_reload(pipette, tiprack_slots, total_refills_needed=0):
        """Pause for manual refill when the loaded tipracks run out of tips."""
        try:
            pipette.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            refill_note = ""
            if total_refills_needed > 0:
                refill_note = f" This run is expected to need about {total_refills_needed} manual refill(s)."
            protocol.pause(
                f"Refill clip tipracks in slots {tiprack_slots} and resume the protocol.{refill_note}"
            )
            pipette.reset_tipracks()
            pipette.pick_up_tip()

    def mix_linkers_function(Mix_linkers_bool, clips_dict, pipette, source_plates, PIPETTE_TYPE):
            #pipetting speeds - default rates in ul /s
        if robot_type=='Flex':
            if PIPETTE_TYPE=="Flex_1channel_50":
                pipette.flow_rate.aspirate = 20
                pipette.flow_rate.dispense = 20
                pipette.flow_rate.blow_out = 35
        elif robot_type=='OT-2':            
            if PIPETTE_TYPE=="p20_single_gen2":
                pipette.flow_rate.aspirate = 8
                pipette.flow_rate.dispense = 8
                pipette.flow_rate.blow_out = 15
            else: 
                raise ValueError("Don't have a single-channel P20 or P50 pipette loaded")
  

        # Use the user-specified working volume directly; premix stroke size is derived later.
        requested_linkers_volume = __PARAMETERS['linkers_volume']['value']
        if requested_linkers_volume < 20:
            raise ValueError("Linker mixing volume must be at least 20 uL.")
        linkers_volume = requested_linkers_volume
    
        
        #linker_offset=np.log(linker_vol)

        if Mix_linkers_bool:
            #Extracts lists from clips_dict
            prefixes = []
            loop_prefixes_wells = clips_dict["prefixes_wells"]
            loop_prefixes_plates = clips_dict["prefixes_plates"]
            len_prefixes = len(clips_dict["prefixes_wells"])
            #Creates 2d array of wells and plates
            for i in range(len_prefixes):
                prefixes.append([loop_prefixes_plates[i], loop_prefixes_wells[i]])
            #Prunes to unique sets of well/plate so duplicates are removed
            #This means any well/plate combination will only be mixed once
            prefixes_unique = np.unique(np.array(prefixes), axis=0)

            suffixes = []
            loop_suffixes_wells = clips_dict["suffixes_wells"]
            loop_suffixes_plates = clips_dict["suffixes_plates"]
            len_suffixes = len(clips_dict["suffixes_wells"])
            #Creates 2d array of wells and plates
            for i in range(len_suffixes):
                suffixes.append([loop_suffixes_plates[i], loop_suffixes_wells[i]])
            #Prunes to unique sets of well/plate so duplicates are removed
            #This means any well/plate combination will only be mixed once
            suffixes_unique = np.unique(np.array(suffixes), axis=0)

            ##Execute the mix 
            # [clip_num,0] addresses the plate location
            # [clip_num,1] addresses the well location
            for clip_num in range(len(prefixes_unique)):  #high = 2.5, normal = 1, slow = 0.5,  vslow = 0.
                configure_flex_low_volume_mode(pipette, get_premix_stroke_volume(linkers_volume))
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                perform_premix(
                    pipette,
                    source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]],
                    linkers_volume
                )

            for clip_num in range(len(suffixes_unique)):  
                configure_flex_low_volume_mode(pipette, get_premix_stroke_volume(linkers_volume))
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                perform_premix(
                    pipette,
                    source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]],
                    linkers_volume
                )
        else:
            pass

    def mix_parts_function(Mix_parts_bool, clips_dict, pipette_name, source_plates, PIPETTE_TYPE):
        pipette = pipette_name

        if robot_type=='Flex':
            if PIPETTE_TYPE=="Flex_1channel_50":
                pipette.flow_rate.aspirate = 20
                pipette.flow_rate.dispense = 20
                pipette.flow_rate.blow_out = 35
        elif robot_type=='OT-2':            
            if PIPETTE_TYPE=="p20_single_gen2":
                pipette.flow_rate.aspirate = 8
                pipette.flow_rate.dispense = 8
                pipette.flow_rate.blow_out = 15
            else: 
                print("Don't have a single-channel P20 or P50 pipette loaded"),
                protocol.pause()
        # Use the user-specified working volume directly; premix stroke size is derived later.
        requested_parts_volume = __PARAMETERS['parts_volume']['value']
        if requested_parts_volume < 10:
            raise ValueError("Part mixing volume must be at least 10 uL.")
        parts_volume = requested_parts_volume

        if Mix_parts_bool:
            parts = []
            loop_parts_wells = clips_dict["parts_wells"]
            loop_parts_plates = clips_dict["parts_plates"]
            len_parts = len(clips_dict["parts_wells"])

            for i in range(len_parts):
                parts.append([loop_parts_plates[i], loop_parts_wells[i]])

            parts_unique = np.unique(np.array(parts), axis=0)

            for clip_num in range(len(parts_unique)):
                configure_flex_low_volume_mode(pipette, get_premix_stroke_volume(parts_volume))
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
                pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
                perform_premix(
                    pipette,
                    source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]],
                    parts_volume
                )
            else:
                pass

    def clip(
            prefixes_wells,
            prefixes_plates,
            suffixes_wells,
            suffixes_plates,
            parts_wells,
            parts_plates,
            parts_vols,
            water_vols):

        ### Calculating number of unique linkers and parts (again) for tip# calculation
        prefixes = []
        loop_prefixes_wells = clips_dict["prefixes_wells"]
        loop_prefixes_plates = clips_dict["prefixes_plates"]
        len_prefixes = len(clips_dict["prefixes_wells"])
        #Creates 2d array of wells and plates
        for i in range(len_prefixes):
            prefixes.append([loop_prefixes_plates[i], loop_prefixes_wells[i]])
        #Prunes to unique sets of well/plate so duplicates are removed
        #This means any well/plate combination will only be mixed once
        prefixes_unique = np.unique(np.array(prefixes), axis=0)

        suffixes = []
        loop_suffixes_wells = clips_dict["suffixes_wells"]
        loop_suffixes_plates = clips_dict["suffixes_plates"]
        len_suffixes = len(clips_dict["suffixes_wells"])
        for i in range(len_suffixes):
            suffixes.append([loop_suffixes_plates[i], loop_suffixes_wells[i]])
        suffixes_unique = np.unique(np.array(suffixes), axis=0)

        parts = []
        loop_parts_wells = clips_dict["parts_wells"]
        loop_parts_plates = clips_dict["parts_plates"]
        len_parts = len(clips_dict["parts_wells"])
        for i in range(len_parts):
            parts.append([loop_parts_plates[i], loop_parts_wells[i]])
        parts_unique = np.unique(np.array(parts), axis=0)

        clip_count = len(parts_wells)
        purification_script_compatible = clip_count <= 48
        if not purification_script_compatible:
            protocol.comment(
                "Warning: this clip run contains more than 48 clips, so it is "
                "not directly compatible with Step 2 "
                "(2_UNI_purification_template_APIv2_21.py), which is capped at 48."
            )
        
        # Calculates whether one, two, or three tipracks are needed, which are in slots 3, 6, and 9 respectively
        # loads tipracks
        if Mix_linkers_bool: 
            if Mix_parts_bool:             
                total_tips = (4 * clip_count) + len(prefixes_unique) + len(suffixes_unique) + len(parts_unique)
            else: total_tips = (4 * clip_count) + len(prefixes_unique) + len(suffixes_unique)
        else: 
            if Mix_parts_bool:
                total_tips = (4 * clip_count) + len(parts_unique)
            else: total_tips = (4 * clip_count)

        letter_dict = {'A': 0, 'B': 1, 'C': 2,
                       'D': 3, 'E': 4, 'F': 5,
                       'G': 6, 'H': 7
                       }
        tiprack_1_tips = (
            13 - int(INITIAL_TIP[1:])) * 8 - letter_dict[INITIAL_TIP[0]]
        if total_tips > tiprack_1_tips:
            tiprack_num = 1 + (total_tips - tiprack_1_tips) // 96 + \
            (1 if (total_tips - tiprack_1_tips) % 96 > 0 else 0)
        else:
            tiprack_num = 1
        loaded_tiprack_count = min(tiprack_num, len(CANDIDATE_TIPRACK_SLOTS))
        slots = CANDIDATE_TIPRACK_SLOTS[:loaded_tiprack_count]
        loaded_tip_capacity = tiprack_1_tips + max(0, loaded_tiprack_count - 1) * 96
        manual_refills_needed = 0
        if total_tips > loaded_tip_capacity:
            additional_tips = total_tips - loaded_tip_capacity
            manual_refills_needed = 1 + (additional_tips - 1) // (96 * loaded_tiprack_count)

        # loads the correct number of tipracks
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
  
        # Loads pipette according to constants assigned above
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)

        # Defines where the destination wells are within the destination plate
        destination_wells = destination_plate.wells()[0:clip_count]

        ### Load Tube Rack
        # Loads tube rack according to constants assigned above
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)

        # Defines positions of master mix and water within the tube rack
        master_mix = tube_rack[MASTER_MIX_WELL]
        water = tube_rack[WATER_WELL]

         ### Loading Source Plates
        # Makes a source plate key for where prefixes, suffixes, and parts are located, according to the dictionary generated by the DNA-BOT
        source_plates = {}
        source_plates_keys = list(set((prefixes_plates + suffixes_plates + parts_plates)))

        # Loads plates according to the source plate key
        for key in source_plates_keys:
            source_plates[key]=protocol.load_labware(SOURCE_PLATE_TYPE, key)
        
        ###Pre-Mixing of Prefixes and Suffixes or Parts

        mix_linkers_function(Mix_linkers_bool, clips_dict, pipette, source_plates, PIPETTE_TYPE)
        mix_parts_function(Mix_parts_bool, clips_dict, pipette, source_plates, PIPETTE_TYPE)

        ### Use the defaults from the active pipette mode for clip setup transfers.
        
        # transfer master mix into destination wells
        pipette.well_bottom_clearance.aspirate = 1  # tip is x mm above well bottom
        pipette.well_bottom_clearance.dispense = 0  # tip is y mm above well bottom        
        configure_flex_low_volume_mode(pipette, MASTER_MIX_VOLUME)
        pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
        pipette.distribute(MASTER_MIX_VOLUME, master_mix, destination_wells, blow_out=True, blowout_location='source well', new_tip='never')
        pipette.drop_tip()

        # transfer water into destination wells
        pipette.well_bottom_clearance.aspirate = 1  # tip is x mm above well bottom
        pipette.well_bottom_clearance.dispense = 3  # tip is y mm above well bottom
        configure_flex_low_volume_mode(pipette, 1)
        pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
        for destination_well, water_volume in zip(destination_wells, water_vols):
            if water_volume <= 0:
                continue
            pipette.aspirate(water_volume, water.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(
                water_volume,
                destination_well.bottom(3),
                push_out=get_low_volume_push_out(water_volume),
            )
        pipette.drop_tip()

    
        #NEW transfer function for prefix, suffix and parts with custom mix parameters
        for clip_num in range(clip_count):
            pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
            pipette.well_bottom_clearance.dispense = 2  # tip is 2 mm above well bottom
            #Prefix Transfer
            configure_flex_low_volume_mode(pipette, 1)
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            prefix_source = source_plates[prefixes_plates[clip_num]][prefixes_wells[clip_num]]
            pre_wet(pipette, prefix_source)
            pipette.aspirate(1, prefix_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(1, destination_wells[clip_num].bottom(2), push_out=get_low_volume_push_out(1))                                                    #changed from: from_center_cartesian(0.2, 0, -0.9)
            protocol.delay(seconds=0.5)
            #mix after transfer
            perform_clip_mix(pipette, destination_wells[clip_num], mix_volume=20, repetitions=2)
            #Suffix Transfer
            configure_flex_low_volume_mode(pipette, 1)
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            suffix_source = source_plates[suffixes_plates[clip_num]][suffixes_wells[clip_num]]
            pre_wet(pipette, suffix_source)
            pipette.aspirate(1, suffix_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(1, destination_wells[clip_num].bottom(3), push_out=get_low_volume_push_out(1))
            protocol.delay(seconds=0.5)
            #mix after transfer
            perform_clip_mix(pipette, destination_wells[clip_num], mix_volume=20, repetitions=2)
            #Part Transfer
            configure_flex_low_volume_mode(pipette, parts_vols[clip_num])
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            part_source = source_plates[parts_plates[clip_num]][parts_wells[clip_num]]
            pre_wet(pipette, part_source)
            pipette.aspirate(parts_vols[clip_num], part_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(parts_vols[clip_num], destination_wells[clip_num].bottom(3), push_out=get_low_volume_push_out(parts_vols[clip_num]))
            protocol.delay(seconds=0.5)
            #mix after transfer
            perform_clip_mix(pipette, destination_wells[clip_num], mix_volume=20, repetitions=4)

    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(**clips_dict)
    ### PCR Reaction in Thermocycler

    # close lid and set lid temperature, PCR will not start until lid reaches 37C
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)

    # Runs 30 cycles of 37C for 2 minutes and 20C for 1 minute, then holds for 60C for 10 minutes
    #37C for 10 min emphasises final BsaI cleavage at end of reaction
    #65C for 20 min is specified heat deactivation for BsaI HFv2, will also deactivate ligase
    profile = [
        {'temperature': 37, 'hold_time_minutes': 2},
        {'temperature': 20, 'hold_time_minutes': 1}]
    tc_mod.execute_profile(steps=profile, repetitions=30, block_max_volume=30)
    tc_mod.set_block_temperature(37, hold_time_minutes=10, block_max_volume=30)
    tc_mod.set_block_temperature(65, hold_time_minutes=20, block_max_volume=30)
    tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=30)
    
    
    #Q Does block_max_volume define total volume in block or individual wells?
    #Thermo lid at end of reaction
    if __PARAMETERS['clip_keep_thermo_lid_closed']['value']=='Yes':
        Thermo_lid_bool=True
    else:
        Thermo_lid_bool = False

    if Thermo_lid_bool:
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(temperature=4)  # The temperature will be held even after this line
        # Temperature will be maintained even after the end of the script
    else:
        tc_mod.set_lid_temperature(45)
        tc_mod.deactivate_lid()
        tc_mod.open_lid()
         #output command actions in simulate
    for line in protocol.commands(): 
            print(line)

# if __name__ == "__main__":
#     #robot_type = input("Enter robot type (Flex or OT-2): ").strip() or "Flex"
#     robot_type = "Flex"
#     from flex_simulate import FlexibleSimulate
#     # Use the custom FlexSimulate class
#     protocol = FlexibleSimulate.get_protocol_api("2.20", robot_type=robot_type)  # Ensure the correct API level is used

#     # # Debugging: inspect protocol setup
#     print(f"Simulated robot type: {protocol.robot_type}")
#     run(protocol)  # Call the `run` function for the protocol logic

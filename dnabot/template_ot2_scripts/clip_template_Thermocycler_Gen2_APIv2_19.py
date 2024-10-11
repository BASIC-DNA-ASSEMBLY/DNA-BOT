from opentrons import protocol_api
#from mix_functions import mix_linkers_function, mix_parts_function
import numpy as np

# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.19',
     'protocolName': 'DNABOT Step 1: Clip Reaction with thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module gen2.'
}
linker_volume=20

# example dictionary produced by DNA-BOT for a single construct containing 4 parts, un-comment and run to test the template
clips_dict={"prefixes_wells": ["A1", "B1", "C1", "D1"],
            "prefixes_plates": ["2", "2", "2", "2"],
            "suffixes_wells": ["A2", "B2", "C2", "D2"],
            "suffixes_plates": ["2", "2", "2", "2"],
            "parts_wells": ["C1", "C2", "C3", "C4"],
            "parts_plates": ["2", "2", "2", "2"],
            "parts_vols": [1, 1, 1, 1],
            "water_vols": [7.0, 7.0, 7.0, 7.0]}

# __LABWARES is expected to be redefined by "generate_ot2_script" method
# Test dict - values used here for simulation use generic Opentrons definitions to avoid
# specifying custom labware in simulate, which is not straightforward
# custom labware currently commented out
__LABWARES={
    "p20_single": {"id": "p20_single_gen2"}, 
    "p300_multi": {"id": "p300_multi_gen2"}, 
    "mag_deck": {"id": "magdeck"}, 
    "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, 
    "96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, 
    "24_tuberack_1500ul": {"id": "opentrons_24_tuberack_nest_1.5ml_snapcap"}, 
    "96_wellplate_200ul_pcr_step_14": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
    "96_wellplate_200ul_pcr_step_23": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
    "clip_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"},
    "mix_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"},
    #"clip_plate": {"id": "4ti0960rig_96_wellplate_200ul"},
    #"mix_plate": {"id": "4ti0960rig_96_wellplate_200ul"},
    "clip_source_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"},
    "agar_plate_step_4": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
    "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, 
    "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}}
        #BELOW is the code that defines the labwares in the clip script
        # self.user_settings['labwares']['p20_single']['id'] = self.labware_p10_single_entry.get()
        # self.user_settings['labwares']['p300_multi']['id'] = self.labware_p300_multi_entry.get()
        # self.user_settings['labwares']['mag_deck']['id'] = self.labware_mag_deck_entry.get()
        # self.user_settings['labwares']['24_tuberack_1500ul']['id'] = self.labware_24_tuberack_1500ul_entry.get()
        # self.user_settings['labwares']['96_tiprack_20ul']['id'] = self.labware_96_tiprack_20ul_entry.get()
        # self.user_settings['labwares']['96_tiprack_300ul']['id'] = self.labware_96_tiprack_300ul_entry.get()
        # self.user_settings['labwares']['clip_source_plate']['id'] = self.labware_clip_source_plate_entry.get()
        # self.user_settings['labwares']['clip_plate']['id'] = self.labware_clip_plate_entry.get()
        # self.user_settings['labwares']['mix_plate']['id'] = self.labware_mix_plate_entry.get()
        # self.user_settings['labwares']['final_assembly_plate']['id'] = self.labware_final_assembly_plate_entry.get()
        # self.user_settings['labwares']['transfo_plate']['id'] = self.labware_transfo_plate_entry.get()
        # self.user_settings['labwares']['transfo_plate_wo_thermo']['id'] = self.labware_transfo_plate_wo_thermo_entry.get()
        # self.user_settings['labwares']['agar_plate']['id'] = self.agar_plate_entry.get()
        # self.user_settings['labwares']['12_reservoir_21000ul']['id'] = self.labware_12_reservoir_21000ul_entry.get()
        # self.user_settings['labwares']['96_deepwellplate_2ml']['id'] = self.labware_96_deepwellplate_2ml_entry.get()
        # self.user_settings['labwares']['12_corning_wellplate']['id'] = self.labware_12_corning_wellplate_entry.get()

__PARAMETERS={
    "clip_keep_thermo_lid_closed": {"value": "clip_keep_thermo_lid_closed"}
}

# Parameters for the clip reaction step
# self.user_settings["parameters"]["clip_keep_thermo_lid_closed"]["value"] = to_numeric_value(self.param_clip_thermo_lid_closed.get())


def run(protocol: protocol_api.ProtocolContext):

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    tiprack_type=__LABWARES['96_tiprack_20ul']['id']
    INITIAL_TIP = 'A1'
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9']

    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = __LABWARES['p20_single']['id']
    PIPETTE_MOUNT = 'right'
        ### Load Pipette
        # checks if it's a P20 Single pipette
    if PIPETTE_TYPE != 'p20_single_gen2':
        print('Define labware must be changed to use', PIPETTE_TYPE)
        exit()
    #thermocycler module gen2 - turn off lid and cool plate to reduce evaporation
    tc_mod = protocol.load_module(module_name="thermocyclerModuleV2")
    tc_mod.open_lid()
    tc_mod.deactivate_lid()
    tc_mod.set_block_temperature(temperature=4) 
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
    TUBE_RACK_POSITION = '4'
    MASTER_MIX_WELL = 'A1'
    WATER_WELL = 'A2'
    MASTER_MIX_VOLUME = 20

    # Old Mix settings
    LINKER_MIX_SETTINGS = (2, 3)
    PART_MIX_SETTINGS = (2, 5)
    #choose to enable pre-mix for prefixes/suffixes and parts plate
    Mix_linkers_bool = True
    Mix_parts_bool = True

    def mix_linkers_function(Mix_linkers_bool, clips_dict, pipette_name, source_plates):
        pipette = pipette_name
        #pipetting speeds - default rates in ul /s
        pipette.flow_rate.aspirate = 6
        pipette.flow_rate.dispense = 6
        pipette.flow_rate.blow_out = 15
        #pipetting rates below - expressed as multiple of default 
        high = 2.5
        normal = 1
        slow = 0.5
        vslow = 0.2
        #Linker reagent volume - specify minimum volume in linker wells
        #linker_volume=20
        #set maximum volume for mixing calculations as 40 as P20 pipette being used
        #maximum linker mix is set as linker_vol/2
        if linker_volume>40:
            linker_vol=40
        else:
            linker_vol=linker_volume

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
            for clip_num in range(len(prefixes_unique)):
                pipette.pick_up_tip()
                pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(1), rate=high)
                pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(1), rate=normal)
                pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(1.5), rate=slow)
                protocol.delay(seconds=1)
                pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].bottom(2), rate=vslow, push_out=linker_vol/20)
                pipette.move_to(source_plates[prefixes_unique[clip_num, 0]][prefixes_unique[clip_num, 1]].top(-2)) # move to 2mm below the top of current well
                pipette.blow_out()
                pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
                pipette.drop_tip()

            for clip_num in range(len(suffixes_unique)):
                pipette.pick_up_tip()
                pipette.well_bottom_clearance.aspirate = 2  # tip is x mm above well bottom
                pipette.well_bottom_clearance.dispense = 1  # tip is y mm above well bottom
                pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(1), rate=high)
                pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(1), rate=normal)
                pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(1.5), rate=slow)
                protocol.delay(seconds=1)
                pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].bottom(2), rate=vslow, push_out=linker_vol/20)
                pipette.move_to(source_plates[suffixes_unique[clip_num, 0]][suffixes_unique[clip_num, 1]].top(-2)) # move to 2mm below the top of current well
                pipette.blow_out()
                pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
                pipette.drop_tip()
        else:
            pass

    def mix_parts_function(Mix_parts_bool, clips_dict, pipette_name, source_plates):
        pipette = pipette_name
        #pipetting speeds - expressed as multiple of default
        pipette.flow_rate.aspirate = 6
        pipette.flow_rate.dispense = 6
        pipette.flow_rate.blow_out = 15
        high = 2.5
        normal = 1
        slow = 0.5
        vslow = 0.2
        #Linker reagent volume - specify minimum volume in linker wells
        #linker_volume=20
        #set maximum volume for mixing calculations as 40 as P20 pipette being used
        #maximum linker mix is set as linker_vol/2
        if linker_volume>40:
            linker_vol=40
        else:
            linker_vol=linker_volume
        
        if Mix_parts_bool:
            parts = []
            loop_parts_wells = clips_dict["parts_wells"]
            loop_parts_plates = clips_dict["parts_plates"]
            len_parts = len(clips_dict["parts_wells"])

            for i in range(len_parts):
                parts.append([loop_parts_plates[i], loop_parts_wells[i]])

            parts_unique = np.unique(np.array(parts), axis=0)

            for clip_num in range(len(parts_unique)):
                pipette.pick_up_tip()
                pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
                pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
                pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(1), rate=high)
                pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(2), rate=normal)
                pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(1), rate=normal)
                pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(1.5), rate=slow)
                protocol.delay(seconds=1)
                pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].bottom(2), rate=vslow, push_out=linker_vol/20)
                pipette.move_to(source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]].top(-2)) # move to 2mm below the top of current well
                pipette.blow_out()
                pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
                pipette.drop_tip()
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
        
        # Calculates whether one, two, or three tipracks are needed, which are in slots 3, 6, and 9 respectively
        # loads tipracks
        if Mix_linkers_bool: 
            if Mix_parts_bool:             
                total_tips = (4 * len(parts_wells)) + len(prefixes_unique) + len(suffixes_unique) + len(parts_unique)
            else: total_tips = (4 * len(parts_wells)) + len(prefixes_unique) + len(suffixes_unique)
        else: 
            if Mix_parts_bool:
                total_tips = (4 * len(parts_wells)) + len(parts_unique)
            else: total_tips = (4 * len(parts_wells))

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
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]

        # loads the correct number of tipracks
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
  
        # Loads pipette according to constants assigned above
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)

        # Defines where the destination wells are within the destination plate
        destination_wells = destination_plate.wells()[0:len(parts_wells)]

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

        mix_linkers_function(Mix_linkers_bool, clips_dict, pipette, source_plates)
        mix_parts_function(Mix_parts_bool, clips_dict, pipette, source_plates)

        ### Reset pipette clearance for setting up clip reactions - pipetting small volume into larger volume
        pipette.flow_rate.aspirate = 6
        pipette.flow_rate.dispense = 6
        pipette.flow_rate.blow_out = 15
        high = 2.5
        normal = 1
        slow = 0.5
        vslow = 0.2
        pipette.well_bottom_clearance.aspirate = 1  # tip is x mm above well bottom
        pipette.well_bottom_clearance.dispense = 1  # tip is y mm above well bottom
        
        # transfer master mix into destination wells
                
        pipette.pick_up_tip()
        pipette.transfer(MASTER_MIX_VOLUME, master_mix, destination_wells, blow_out=True, blowout_location='destination well', new_tip='never', rate=slow)
        pipette.drop_tip()

        # transfer water into destination wells
        pipette.well_bottom_clearance.aspirate = 1  # tip is x mm above well bottom
        pipette.well_bottom_clearance.dispense = 3  # tip is y mm above well bottom
        
        pipette.pick_up_tip()
        pipette.distribute(water_vols, water, destination_wells, blow_out=True, blowout_location='source well', new_tip='never', rate=slow)
        pipette.drop_tip()

        # OLD transfer prefixes, suffixes, and parts into destination wells     
        #for clip_num in range(len(parts_wells)):
            #pipette.transfer(1, source_plates[prefixes_plates[clip_num]].wells(prefixes_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS, rate=slow)
            # pipette.transfer(1, source_plates[suffixes_plates[clip_num]].wells(suffixes_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS, rate=slow)
            # pipette.transfer(parts_vols[clip_num], source_plates[parts_plates[clip_num]].wells(parts_wells[clip_num]), destination_wells[clip_num], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS, rate=slow)
        
        #NEW transfer function for prefix, suffix and parts with custom mix parameters
        for clip_num in range(len(parts_wells)):
            pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
            pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
            #Prefix Transfer
            pipette.pick_up_tip()
            pipette.aspirate(1, source_plates[prefixes_plates[clip_num]][prefixes_wells[clip_num]].bottom(1), rate=slow)
            pipette.dispense(1, destination_wells[clip_num].bottom(3), rate=slow)
            #mix after transfer
            pipette.aspirate(2, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.dispense(2, destination_wells[clip_num].bottom(3), rate=high)
            pipette.aspirate(3, destination_wells[clip_num].bottom(2), rate=normal)
            pipette.dispense(3, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.aspirate(4, destination_wells[clip_num].bottom(2), rate=slow)
            pipette.dispense(4, destination_wells[clip_num].bottom(3), push_out=1, rate=vslow)
            pipette.move_to(destination_wells[clip_num].top(-5))
            pipette.blow_out()
            pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
            pipette.drop_tip()
            #Suffix Transfer
            pipette.pick_up_tip()
            pipette.aspirate(1, source_plates[suffixes_plates[clip_num]][suffixes_wells[clip_num]].bottom(1), rate=slow)
            pipette.dispense(1, destination_wells[clip_num].bottom(3), rate=slow)
            #mix after transfer
            pipette.aspirate(2, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.dispense(2, destination_wells[clip_num].bottom(3), rate=high)
            pipette.aspirate(3, destination_wells[clip_num].bottom(2), rate=normal)
            pipette.dispense(3, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.aspirate(4, destination_wells[clip_num].bottom(2), rate=slow)
            pipette.dispense(4, destination_wells[clip_num].bottom(3), push_out=1, rate=vslow)
            pipette.move_to(destination_wells[clip_num].top(-5))
            pipette.blow_out()
            pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
            pipette.drop_tip()
            #Part Transfer
            pipette.pick_up_tip()
            pipette.aspirate(parts_vols[clip_num], source_plates[parts_plates[clip_num]][parts_wells[clip_num]].bottom(1), rate=slow)
            pipette.dispense(parts_vols[clip_num], destination_wells[clip_num].bottom(3), rate=slow)
            #mix after transfer
            pipette.aspirate(5, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.dispense(5, destination_wells[clip_num].bottom(3), rate=high)
            pipette.aspirate(10, destination_wells[clip_num].bottom(2), rate=normal)
            pipette.dispense(10, destination_wells[clip_num].bottom(1), rate=normal)
            pipette.aspirate(15, destination_wells[clip_num].bottom(2), rate=slow)
            pipette.dispense(15, destination_wells[clip_num].bottom(3), push_out=1, rate=vslow)
            pipette.move_to(destination_wells[clip_num].top(-5))
            pipette.blow_out()
            pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
            pipette.drop_tip()

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
    if __PARAMETERS["clip_keep_thermo_lid_closed"]["value"] == 1:
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(temperature=4)  # The temperature will be held even after this line
        # Temperature will be maintained even after the end of the script
    else:
        tc_mod.set_lid_temperature(37)
        tc_mod.open_lid()
         #output command actions in simulate
        for line in protocol.commands(): 
            print(line)
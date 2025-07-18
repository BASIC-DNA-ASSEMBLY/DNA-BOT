from opentrons import protocol_api
import json

# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.8',
     'protocolName': 'CLIP_With_Thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module.'}

# Load CLIP data from JSON file
# This will be replaced by the parser with embedded JSON data
with open('clips_data.json') as f:
    clips_dict = json.load(f)

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# opentrons_simulate.exe dnabot\template_ot2_scripts\clip_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["1", "1", "1", "1", "1"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

def run(protocol: protocol_api.ProtocolContext):
# added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    tiprack_type="opentrons_96_tiprack_20ul"
    INITIAL_TIP = 'A1'
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '5']

    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = 'p20_single_gen2'
    PIPETTE_MOUNT = 'right'
        ### Load Pipette
        # checks if it's a P10 Single pipette
    if PIPETTE_TYPE != 'p20_single_gen2':
        print('Define labware must be changed to use', PIPETTE_TYPE)
        exit()

    # Thermocycler Module
    if thermocycler_gen == 'gen1':
        tc_mod = protocol.load_module('Thermocycler Module')
    else:  # gen2
        tc_mod = protocol.load_module('thermocyclerModuleV2')
        
    # Destination Plates
    DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'

    # Loads destination plate onto Thermocycler Module
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
    tc_mod.open_lid()
    tc_mod.set_block_temperature(20)

    # Source Plates
    SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used

    # Tube Rack
    TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
    TUBE_RACK_POSITION = '4'
    MASTER_MIX_WELL = 'A1'
    WATER_WELL = 'A2'
    MASTER_MIX_VOLUME = 20

    # Mix settings
    LINKER_MIX_SETTINGS = (1, 3)
    PART_MIX_SETTINGS = (4, 5)

    def clip(clips_dict):
        ### Loading Tiprack
        total_tips = 4 * len(clips_dict)
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
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)
        # Destination plate on thermocycler
        destination_wells = [destination_plate.wells_by_name()[w] for w in clips_dict.keys()]
        # Tube rack
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        master_mix = tube_rack.wells(MASTER_MIX_WELL)
        water = tube_rack.wells(WATER_WELL)
        # Load source plates
        source_plates = {}
        all_plates = set()
        for well_info in clips_dict.values():
            all_plates.add(well_info['prefix_plate'])
            all_plates.add(well_info['suffix_plate'])
            all_plates.add(well_info['part_plate'])
        for key in all_plates:
            source_plates[key] = protocol.load_labware(SOURCE_PLATE_TYPE, key)
        dest_wells = list(clips_dict.keys())
        # Master mix transfer
        pipette.pick_up_tip()
        pipette.transfer(MASTER_MIX_VOLUME, master_mix, destination_wells, blow_out=True, blowout_location='destination well', new_tip='never')
        pipette.drop_tip()
        # Water transfer (only if needed)
        water_vols = [clips_dict[w]['water_vol'] for w in dest_wells]
        if any([wv > 0 for wv in water_vols]):
            pipette.transfer(water_vols, water, destination_wells, blow_out=True, blowout_location='destination well', new_tip='always')
        # Prefix, suffix, part transfers
        for i, well in enumerate(dest_wells):
            info = clips_dict[well]
            pipette.transfer(1, source_plates[info['prefix_plate']].wells_by_name()[info['prefix_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(1, source_plates[info['suffix_plate']].wells_by_name()[info['suffix_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(info['part_vol'], source_plates[info['part_plate']].wells_by_name()[info['part_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS)
    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(clips_dict)
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
    tc_mod.set_block_temperature(8, block_max_volume=30)
    tc_mod.set_lid_temperature(37)
    # tc_mod.open_lid()                                     # leave lid shut to prevent evaporation

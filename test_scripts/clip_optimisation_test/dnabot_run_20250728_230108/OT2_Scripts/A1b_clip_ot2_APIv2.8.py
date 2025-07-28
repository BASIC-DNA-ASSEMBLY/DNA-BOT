from opentrons import protocol_api
import json

# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use
# Code has been reordered to better group relevant commands and take the constants out of def clip()

#metadata
metadata = {
     'apiLevel': '2.8',
     'protocolName': 'CLIP_No_Thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2. This version does not include the Thermocycler module.'}

# Load CLIP data from JSON file
# This will be replaced by the parser with embedded JSON data
clips_dict = {
    "A7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "sucB_1",
        "part_source_well": "C3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "alaC_1",
        "part_source_well": "B2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "thiS_1",
        "part_source_well": "C5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "serC_1",
        "part_source_well": "E4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "alr_1",
        "part_source_well": "C2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "dadX_1",
        "part_source_well": "G3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tdcB_2",
        "part_source_well": "F6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gadB_1",
        "part_source_well": "G4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ilvA_1",
        "part_source_well": "A8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvH_1",
        "part_source_well": "A2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    }
}

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["5", "5", "5", "5", "5"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

def run(protocol: protocol_api.ProtocolContext):
# added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    tiprack_type="opentrons_96_tiprack_20ul"
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

    # Source Plates
    SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used

    # Destination Plates
    DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
    DESTINATION_PLATE_POSITION = '1'
            # INITIAL_DESTINATION_WELL constant removed, as destination_plate.wells() automatically starts from A1

    # Tube Rack
    TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
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
        destination_plate = protocol.load_labware(DESTINATION_PLATE_TYPE, DESTINATION_PLATE_POSITION)
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        master_mix = tube_rack.wells(MASTER_MIX_WELL)
        water = tube_rack.wells(WATER_WELL)
        # Load source plates
        source_plates = {}
        all_plates = set()
        for well_info in clips_dict.values():
            all_plates.add(well_info['prefix_source_plate'])
            all_plates.add(well_info['suffix_source_plate'])
            all_plates.add(well_info['part_source_plate'])
        for key in all_plates:
            source_plates[key] = protocol.load_labware(SOURCE_PLATE_TYPE, key)
        # Get destination wells in order
        dest_wells = list(clips_dict.keys())
        destination_wells = [destination_plate.wells_by_name()[w] for w in dest_wells]
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
            pipette.transfer(1, source_plates[info['prefix_source_plate']].wells_by_name()[info['prefix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(1, source_plates[info['suffix_source_plate']].wells_by_name()[info['suffix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(info['part_vol'], source_plates[info['part_source_plate']].wells_by_name()[info['part_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS)
    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(clips_dict)

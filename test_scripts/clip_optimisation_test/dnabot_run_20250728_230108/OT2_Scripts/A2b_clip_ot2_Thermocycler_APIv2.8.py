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
clips_dict = {
    "A7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "folA_1",
        "part_source_well": "E1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "cysM_1",
        "part_source_well": "G5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "serC_1",
        "part_source_well": "E4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gadB_1",
        "part_source_well": "G4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "trxA_1",
        "part_source_well": "C8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ubiC_1",
        "part_source_well": "C6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "pabC_1",
        "part_source_well": "A4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "leuB_2",
        "part_source_well": "H1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "carA_1",
        "part_source_well": "A1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "serC_2",
        "part_source_well": "F3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "dadX_1",
        "part_source_well": "G3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "asnA_1",
        "part_source_well": "G7",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvA_1",
        "part_source_well": "A8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tdcB_2",
        "part_source_well": "F6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gltD_1",
        "part_source_well": "A7",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "argI_1",
        "part_source_well": "D2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "asnB_1",
        "part_source_well": "A3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tyrB_1",
        "part_source_well": "E9",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "aspC_1",
        "part_source_well": "E3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "trxC_1",
        "part_source_well": "A6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tdcB_1",
        "part_source_well": "A5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "aspC_1",
        "part_source_well": "E3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "carA_1",
        "part_source_well": "A1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "carB_1",
        "part_source_well": "C1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "carB_1",
        "part_source_well": "C1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "panD_1",
        "part_source_well": "D4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "serC_2",
        "part_source_well": "F3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "aroL_1",
        "part_source_well": "G2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "folM_1",
        "part_source_well": "G10",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "folA_1",
        "part_source_well": "E1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "dmlA_1",
        "part_source_well": "E10",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "leuB_2",
        "part_source_well": "H1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B11": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvN_1",
        "part_source_well": "G11",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ilvH_1",
        "part_source_well": "A2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D11": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ydiB_1",
        "part_source_well": "A11",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "sucB_1",
        "part_source_well": "C3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "folM_1",
        "part_source_well": "G10",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    }
}

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
            all_plates.add(well_info['prefix_source_plate'])
            all_plates.add(well_info['suffix_source_plate'])
            all_plates.add(well_info['part_source_plate'])
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
            pipette.transfer(1, source_plates[info['prefix_source_plate']].wells_by_name()[info['prefix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(1, source_plates[info['suffix_source_plate']].wells_by_name()[info['suffix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(info['part_vol'], source_plates[info['part_source_plate']].wells_by_name()[info['part_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS)
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

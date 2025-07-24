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
    "A1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "A1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "B1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "C1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "D1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "E1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "F1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G1": {
        "prefix_linker": "LMS-P",
        "prefix_source_well": "D12",
        "prefix_source_plate": "1",
        "part": "SV39",
        "part_source_well": "A12",
        "part_source_plate": "1",
        "suffix_linker": "LMP-S",
        "suffix_source_well": "C12",
        "suffix_source_plate": "1",
        "Clip_Well": "G1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H1": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "aroK_1",
        "part_source_well": "E7",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H1",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A2": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "aroL_1",
        "part_source_well": "G2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B2": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "asnA_1",
        "part_source_well": "G7",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C2": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "thiS_2",
        "part_source_well": "D10",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D2": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "argF_1",
        "part_source_well": "E2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E2": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "aspC_2",
        "part_source_well": "H3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F2": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tdcB_1",
        "part_source_well": "A5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G2": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gdhA_2",
        "part_source_well": "B5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H2": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "cysM_1",
        "part_source_well": "G5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H2",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A3": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "cysK_1",
        "part_source_well": "E5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B3": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "panD_1",
        "part_source_well": "D4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C3": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "alr_1",
        "part_source_well": "C2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D3": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tynA_1",
        "part_source_well": "D6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E3": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gltD_1",
        "part_source_well": "A7",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F3": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "pabC_1",
        "part_source_well": "A4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G3": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "thiS_1",
        "part_source_well": "C5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H3": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "alaC_1",
        "part_source_well": "B2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H3",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A4": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gltB_1",
        "part_source_well": "G6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B4": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C4": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "alaC_1",
        "part_source_well": "B2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D4": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "thiS_1",
        "part_source_well": "C5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E4": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "serC_1",
        "part_source_well": "E4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F4": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "alr_1",
        "part_source_well": "C2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G4": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ubiC_1",
        "part_source_well": "C6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H4": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gltB_1",
        "part_source_well": "G6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H4",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A5": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gdhA_1",
        "part_source_well": "C4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B5": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "trxA_1",
        "part_source_well": "C8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C5": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "trxC_1",
        "part_source_well": "A6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D5": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gdhA_2",
        "part_source_well": "B5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E5": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tynA_1",
        "part_source_well": "D6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F5": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gadB_1",
        "part_source_well": "G4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G5": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "cysK_1",
        "part_source_well": "E5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H5": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "aroE_1",
        "part_source_well": "C7",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H5",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A6": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "asnB_1",
        "part_source_well": "A3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B6": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tyrB_1",
        "part_source_well": "E9",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C6": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "argI_1",
        "part_source_well": "D2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D6": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "argI_2",
        "part_source_well": "H9",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E6": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "thiS_2",
        "part_source_well": "D10",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F6": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ilvN_1",
        "part_source_well": "G11",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G6": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "argI_2",
        "part_source_well": "H9",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G6",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H6": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "argF_1",
        "part_source_well": "E2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H6",
        "plate": 2,
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
            all_plates.add(well_info['prefix_plate'])
            all_plates.add(well_info['suffix_plate'])
            all_plates.add(well_info['part_plate'])
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
            pipette.transfer(1, source_plates[info['prefix_plate']].wells_by_name()[info['prefix_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(1, source_plates[info['suffix_plate']].wells_by_name()[info['suffix_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=LINKER_MIX_SETTINGS)
            pipette.transfer(info['part_vol'], source_plates[info['part_plate']].wells_by_name()[info['part_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='always', mix_after=PART_MIX_SETTINGS)
    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(clips_dict)

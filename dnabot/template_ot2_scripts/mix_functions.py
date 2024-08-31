import numpy as np
#from opentrons import protocol_api

# clips_dict={
#     "prefixes_wells": ["A8", "A8", "C7", "C7", "C10"], 
#     "prefixes_plates": ["2", "2", "2", "2", "3"], 
#     "suffixes_wells": ["B7", "C1", "C1", "C1", "B8"], 
#     "suffixes_plates": ["2", "3", "3", "3", "2"], 
#     "parts_wells": ["E2", "F2", "C2", "B2", "D2"], 
#     "parts_plates": ["5", "5", "5", "5", "5"], 
#     "parts_vols": [1, 1, 1, 1, 1], 
#     "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]
#     }

def mix_prefixes_suffixes(clips_dict):
    prefixes = []
    loop_prefixes_wells = clips_dict["prefixes_wells"]
    loop_prefixes_plates = clips_dict["prefixes_plates"]
    len_prefixes = len(clips_dict["prefixes_wells"])

    for i in range(len_prefixes):
        prefixes.append([loop_prefixes_wells[i], loop_prefixes_plates[i]])

    prefixes_unique = np.unique(np.array(prefixes), axis=0)

    suffixes = []
    loop_suffixes_wells = clips_dict["suffixes_wells"]
    loop_suffixes_plates = clips_dict["suffixes_plates"]
    len_suffixes = len(clips_dict["suffixes_wells"])

    for i in range(len_suffixes):
        suffixes.append([loop_suffixes_wells[i], loop_suffixes_plates[i]])

    suffixes_unique = np.unique(np.array(suffixes), axis=0)

    for clip_num in range(len(prefixes_unique)):
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
        pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
        pipette.aspirate(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]))
        pipette.drop_tip()

    for clip_num in range(len(suffixes_unique)):
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
        pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
        pipette.aspirate(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
        pipette.drop_tip()

def mix_parts(clips_dict):
    parts = []
    loop_parts_wells = clips_dict["parts_wells"]
    loop_parts_plates = clips_dict["parts_plates"]
    len_parts = len(clips_dict["parts_wells"])

    for i in range(len_parts):
        parts.append([loop_parts_wells[i], loop_parts_plates[i]])

    parts_unique = np.unique(np.array(parts), axis=0)

    for clip_num in range(len(parts_unique)):
        pipette.pick_up_tip()
        pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
        pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
        pipette.aspirate(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.aspirate(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.dispense(10, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
        pipette.drop_tip()

#-----------------------------------------------------


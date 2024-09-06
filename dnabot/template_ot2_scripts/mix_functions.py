import numpy as np
from opentrons import protocol_api

#metadata
metadata = {
     'apiLevel': '2.19',
     'protocolName': 'DNABOT Step 1: Clip Reaction with thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module gen2.'
}

clips_dict={
    "prefixes_wells": ["A8", "A8", "C7", "C7", "C10"], 
    "prefixes_plates": ["2", "2", "2", "2", "3"], 
    "suffixes_wells": ["B7", "C1", "C1", "C1", "B8"], 
    "suffixes_plates": ["2", "3", "3", "3", "2"], 
    "parts_wells": ["E2", "F2", "C2", "B2", "D2"], 
    "parts_plates": ["5", "5", "5", "5", "5"], 
    "parts_vols": [1, 1, 1, 1, 1], 
    "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]
    }

def mix_prefixes_suffixes_function(Mix_prefix_and_suffix_bool, clips_dict, pipette_name):
    pipette = pipette_name
    #pipetting speeds - expressed as multiple of default
    high = 3
    normal = 1
    slow = 0.4
    #Linker reagent volume - specify minimum volume in linker wells
    linker_vol=20
    if Mix_prefix_and_suffix_bool:
        #Extracts lists from clips_dict
        prefixes = []
        loop_prefixes_wells = clips_dict["prefixes_wells"]
        loop_prefixes_plates = clips_dict["prefixes_plates"]
        len_prefixes = len(clips_dict["prefixes_wells"])
        #Creates 2d array of wells and plates
        for i in range(len_prefixes):
            prefixes.append([loop_prefixes_wells[i], loop_prefixes_plates[i]])
        #Prunes to unique sets of well/plate so duplicates are removed
        #This means any well/plate combination will only be mixed once
        prefixes_unique = np.unique(np.array(prefixes), axis=0)

        suffixes = []
        loop_suffixes_wells = clips_dict["suffixes_wells"]
        loop_suffixes_plates = clips_dict["suffixes_plates"]
        len_suffixes = len(clips_dict["suffixes_wells"])

        for i in range(len_suffixes):
            suffixes.append([loop_suffixes_wells[i], loop_suffixes_plates[i]])

        suffixes_unique = np.unique(np.array(suffixes), axis=0)
        #Execute the mix 
        # [clip_num,0] addresses the plate location
        # [clip_num,1] addresses the well location
        for clip_num in range(len(prefixes_unique)):
            pipette.pick_up_tip()
            pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
            pipette.well_bottom_clearance.dispense = 1  # tip is 1 mm above well bottom
            pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate(normal))
            pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate(high))
            pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate(normal))
            pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate(high))
            pipette.aspirate(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate(slow))
            protocol.delay(seconds=1)
            pipette.dispense(linker_vol/2, source_plates[prefixes_unique[clip_num, 0]].wells(prefixes_unique[clip_num, 1]), rate (slow), push_out=linker_vol/20)
            pipette.blow_out(well.top(-5))
            pipette.touch_tip(speed=10, radius=0.9, v_offset=-5)
            pipette.drop_tip()

        for clip_num in range(len(suffixes_unique)):
            pipette.pick_up_tip()
            pipette.well_bottom_clearance.aspirate = 2  # tip is 2 mm above well bottom
            pipette.well_bottom_clearance.dispense = 1  # tip is 2 mm above well bottom
            pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.aspirate(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[suffixes_unique[clip_num, 0]].wells(suffixes_unique[clip_num, 1]))
            pipette.drop_tip()
    else:
        pass


def mix_parts_function(Mix_parts_plate_bool, clips_dict, pipette_name):
    pipette = pipette_name
    if Mix_parts_plate_bool:
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
            pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.aspirate(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.dispense(linker_vol/2, source_plates[parts_unique[clip_num, 0]].wells(parts_unique[clip_num, 1]))
            pipette.drop_tip()
        else:
            pass

#-----------------------------------------------------


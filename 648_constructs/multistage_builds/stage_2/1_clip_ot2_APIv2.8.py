from opentrons import protocol_api


# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use
# Code has been reordered to better group relevant commands and take the constants out of def clip()

#metadata
metadata = {
     'apiLevel': '2.8',
     'protocolName': 'CLIP_No_Thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2. This version does not include the Thermocycler module.'}

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["5", "5", "5", "5", "5"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

clips_dict={"prefixes_wells": ["G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G3", "H10", "H10", "H10", "H11", "H11", "H11", "H12", "H12", "H12", "G3", "G3", "G3", "H6", "H6", "H6", "H7", "H7", "H7", "H8", "H8", "H8", "G3", "G3", "G3", "H2", "H2", "H2", "H3", "H3", "H3", "H4", "H4", "H4", "G3", "G3", "G3", "H6", "H6", "H6", "H7", "H7", "H7", "H8", "H8", "H8", "G3", "G3", "G3", "H2", "H2", "H2", "H3", "H3", "H3", "H4", "H4", "H4", "G3", "G3", "G3", "H10", "H10", "H10", "H11", "H11", "H11", "H12", "H12", "H12", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3", "G3"], "prefixes_plates": ["2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2"], "suffixes_wells": ["G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "G4", "H9", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H9", "H9", "H5", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H5", "H5", "H1", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H1", "H1", "H5", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H5", "H5", "H1", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H1", "H1", "H9", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "G2", "H9", "H9", "H9", "H9", "H9", "H5", "H5", "H5", "H1", "H1", "H1", "H5", "H5", "H5", "H1", "H1", "H1", "H9", "H9", "H9", "H9", "H9", "H9", "H5", "H5", "H5", "H1", "H1", "H1", "H5", "H5", "H5", "H1", "H1", "H1", "H9", "H9", "H9", "H9", "H9", "H9", "H5", "H5", "H5", "H1", "H1", "H1", "H5", "H5", "H5", "H1", "H1", "H1", "H9", "H9", "H9"], "suffixes_plates": ["2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2"], "parts_wells": ["F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1", "B1", "D1", "D2", "D3", "D1", "D2", "D3", "D1", "D2", "D3", "B2", "B3", "B1", "D4", "D5", "D6", "D4", "D5", "D6", "D4", "D5", "D6", "B2", "B3", "A1", "D7", "D8", "D9", "D7", "D8", "D9", "D7", "D8", "D9", "A2", "A3", "A1", "D10", "D11", "D12", "D10", "D11", "D12", "D10", "D11", "D12", "A2", "A3", "C1", "E1", "E2", "E3", "E1", "E2", "E3", "E1", "E2", "E3", "C2", "C3", "C1", "E4", "E5", "E6", "E4", "E5", "E6", "E4", "E5", "E6", "C2", "C3", "B4", "B5", "B6", "B4", "B5", "B6", "A4", "A5", "A6", "A4", "A5", "A6", "C4", "C5", "C6", "C4", "C5", "C6", "B7", "B8", "B9", "B7", "B8", "B9", "A7", "A8", "A9", "A7", "A8", "A9", "C7", "C8", "C9", "C7", "C8", "C9", "B10", "B11", "B12", "B10", "B11", "B12", "A10", "A11", "A12", "A10", "A11", "A12", "C10", "C11", "C12", "C10", "C11", "C12"], "parts_plates": ["2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2"], "parts_vols": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0]}


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

        # Loads pipette according to constants assigned above
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)

        ### Load Destination Plate
        # Loads destination plate according to constants assigned above
        destination_plate = protocol.load_labware(DESTINATION_PLATE_TYPE, DESTINATION_PLATE_POSITION)

        # Defines where the destination wells are within the destination plate
        destination_wells = destination_plate.wells()[0:len(parts_wells)]

        ### Load Tube Rack
        # Loads tube rack according to constants assigned above
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)

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

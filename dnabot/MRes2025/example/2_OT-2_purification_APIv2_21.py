''''#This code is modified from ot2 version to compatiable with Flex machine manually by Wenhan Hu (CID:02017595),
as the second step in the BASIC workflow.'''

from opentrons import protocol_api

# Rename to 'purification_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

metadata = {
     'protocolName': 'DNABOT Step 2: Purification',
     'description': 'Implements magbead purification reactions for BASIC assembly using an opentrons flex'}


# requirements
sample_number=12
ethanol_well='A11'
__HARDWARE={"robot_type": {"id": "OT-2"}, "single_pipette": {"id": "p20_single_gen2"}, "single_pipette_mount": {"id": "right"}, "multi_pipette": {"id": "p300_multi_gen2"}, "multi_pipette_mount": {"id": "left"}, "thermocycler": {"id": "thermocyclerModuleV2"}, "mag_deck": {"id": "magnetic module gen1"}}
__LABWARES={"OT-2_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, "OT-2_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "clip_source_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "clip_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "mix_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "final_assembly_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "transform_plate": {"id": "4ti0960rig_96_wellplate_200ul"}, "agar_plate": {"id": "thermoomnitrayfor96spots_96_wellplate_50ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}, "12_corning_wellplate": {"id": "corning_12_wellplate_6.9ml_flat"}}
__PARAMETERS={"clip_keep_thermo_lid_closed": {"value": "No", "id": "No"}, "premix_linkers": {"value": "Yes", "id": "Yes"}, "premix_parts": {"value": "Yes", "id": "Yes"}, "linkers_volume": {"value": 20}, "parts_volume": {"value": 20}, "thermo_temp": {"value": 4}, "purif_magdeck_height": {"value": 10.8}, "purif_wash_time": {"value": 0.5}, "purif_bead_ratio": {"value": 1.8}, "purif_incubation_time": {"value": 5}, "purif_settling_time": {"value": 2}, "purif_drying_time": {"value": 5}, "purif_elution_time": {"value": 2}, "transform_incubation_temp": {"value": 4}, "transform_incubation_time": {"value": 20}}

requirements = {"robotType": __HARDWARE['robot_type']['id'], "apiLevel": "2.21"} 


def run(protocol: protocol_api.ProtocolContext):
    robot_type=__HARDWARE['robot_type']['id']
    if robot_type=='Flex':
        trash = protocol.load_trash_bin("A3")
        tc_mod = protocol.load_module(module_name=__HARDWARE['thermocycler']['id'], location = "B1")
        tiprack_type = ['opentrons_flex_96_tiprack_200ul']
        tiprack_1000 = ['opentrons_flex_tiprack_1000ul'] 
    elif robot_type=='OT-2':
        tiprack_type = __LABWARES['OT-2_tiprack_300ul']['id']
    else:
        raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")

    def magbead(sample_number, ethanol_well):
        # sample_number,
        # ethanol_well,
        elution_buffer_well='A1',
        sample_volume=30,
        bead_ratio=__PARAMETERS['purif_bead_ratio']['value'],
        elution_buffer_volume=40,
        incubation_time=__PARAMETERS['purif_incubation_time']['value'],
        settling_time=__PARAMETERS['purif_settling_time']['value'],
            # if using Gen 2 magentic module, need to change time! see: https://docs.opentrons.com/v2/new_modules.html
            # "The GEN2 Magnetic Module uses smaller magnets than the GEN1 version...this means it will take longer for the GEN2 module to attract beads."
            # Recommended Magnetic Module GEN2 bead attraction time:
                # Total liquid volume <= 50 uL: 5 minutes
            # this template was written with the Gen 1 magnetic module, so the settling time is set to 5 minutes
        drying_time=__PARAMETERS['purif_drying_time']['value'],
        elution_time=__PARAMETERS['purif_elution_time']['value'],
        #sample_offset can be used to start the purification at a different column
        sample_offset=0,
        #tiprack_type=__LABWARES['96_tiprack_300ul']['id']):
        #tiprack_200=__LABWARES['flex_96_tiprack_200ul']['id'],            
        #tiprack_1000=__LABWARES['flex_96_tiprack_1000ul']['id']):

        TIPS_PER_SAMPLE = 5
        TIPS_WASH = 2

### Loading Tiprack

        # Calculates whether one/two/three/four/five tipracks are needed, which are in slots 3, 6, 9, 2, and 5 respectively
        total_tips = sample_number * TIPS_PER_SAMPLE
        tiprack_num = total_tips // 96 + (1 if total_tips % 96 > 0 else 0)
        print(str(tiprack_num) + 'of 200ul tipboxes is needed')
            
        wash_tips = sample_number * TIPS_WASH
        tiprack_wash = wash_tips // 96 + (1 if wash_tips % 96 > 0 else 0)
        print(str(tiprack_wash) + 'of 1000ul tipboxes is needed')
        
        # Tiprack
        if robot_type=='OT-2':
            CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5']
            #tiprack_type = __LABWARES['96_tiprack_300ul']['id']
        elif robot_type=='Flex':
            CANDIDATE_TIPRACK_SLOTS = ["D3", "C3", "B3"]
            CANDIDATE_TIPRACK_SLOT_1000 = "C2"
            #tiprack_type = ['opentrons_flex_96_tiprack_200ul']
            #tiprack_1000 = ['opentrons_flex_96_tiprack_1000ul']
        else:
            raise ValueError("Invalid robot type. Must be 'OT-2' or 'Flex'.")
        
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        slots_wash =  CANDIDATE_TIPRACK_SLOT_1000[:tiprack_wash]
        # loads the correct number of tipracks
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        tiprack_wash=[protocol.load_labware(tiprack_1000, slot) for slot in slots_wash]

        #PIPETTE_TYPE = __LABWARES['p300_multi']['id']
        PIPETTE_TYPE = __HARDWARE['multi_pipette']['id']
        PIPETTE_MOUNT = __HARDWARE['multi_pipette_mount']['id']
        ### Loading Pipettes
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT,tip_racks=tipracks)
           #pipetting speeds -
           # Flex 8- channel default speeds:
           # 50 ul tip = 478 ul/s; 200 ul tip = 716 ul/s; 1000 ul tip = 716 ul/s
            # p300_multi_gen2 default speeds:
            # 300 ul tip = 94 ul/s; 1000 ul tip = 94 ul/s
        if PIPETTE_TYPE=="Flex_8channel_1000":
            pipette.flow_rate.aspirate = 200
            pipette.flow_rate.dispense = 200
            pipette.flow_rate.blow_out = 500
            pipette.max_volume = 200
            pipette.max_volume_1000 = 1000
        elif PIPETTE_TYPE=="p300_multi_gen2":
            pipette.flow_rate.aspirate = 60
            pipette.flow_rate.dispense = 60
            pipette.flow_rate.blow_out = 90
            pipette.max_volume = 300
        else: 
            raise ValueError("Don't have a multi-channel pipette loaded"),
 
        # Magnetic Module
        #MAGDECK_POSITION = '1' Magnetic Block Updated
        MAGNETIC_PLATE_TYPE = __LABWARES['mag_plate']['id']
        if robot_type=='OT-2':
            MAGDECK_POSITION = '1'
            MAGDECK = protocol.load_module(__HARDWARE['mag_deck']['id'], location= MAGDECK_POSITION)
            MAGDECK.disengage()
            mag_plate = protocol.load_labware(MAGNETIC_PLATE_TYPE, MAGDECK_POSITION)
        elif robot_type=='Flex':
            MAGDECK_POSITION = 'D1'
            MAGDECK = protocol.load_module(__HARDWARE['mag_deck']['id'], location=MAGDECK_POSITION)
            mag_plate = MAGDECK.load_labware(MAGNETIC_PLATE_TYPE)

        # Mix Plate
        MIX_PLATE_TYPE = __LABWARES['mix_plate']['id']
        if robot_type=='OT-2':
            MIX_PLATE_POSITION = '4'
        elif robot_type=='Flex':
            MIX_PLATE_POSITION = 'C1'
        
        # Reagents
        REAGENT_CONTAINER_TYPE = __LABWARES['12_reservoir_21000ul']['id']
        if robot_type=='OT-2':
            REAGENT_CONTAINER_POSITION = '7'
        elif robot_type=='Flex':
            REAGENT_CONTAINER_POSITION = 'A2'
            # Beads
        BEAD_CONTAINER_TYPE = __LABWARES['96_deepwellplate_2ml']['id']
        if robot_type=='OT-2':
            BEAD_CONTAINER_POSITION = '8'
        elif robot_type=='Flex':
            BEAD_CONTAINER_POSITION = 'B2'
        
        # Settings
        #LIQUID_WASTE_WELL = 'A5'
        BEADS_WELL = 'A1'
        DEAD_TOTAL_VOL = 5
        #SLOW_HEAD_SPEEDS = {'x': 600 // 4, 'y': 400 // 4, 'z': 125 // 10, 'a': 125 // 10}
        #DEFAULT_HEAD_SPEEDS = {'x': 400, 'y': 400, 'z': 125, 'a': 100}
        IMMOBILISE_MIX_REPS = 10
        MAGDECK_HEIGHT = __PARAMETERS['purif_magdeck_height']['value']
        AIR_VOL_COEFF = 0.1
        ETHANOL_VOL = 150
        WASH_TIME = __PARAMETERS['purif_wash_time']['value']
        ETHANOL_DEAD_VOL = 50 
        ELUTION_MIX_REPS = 20
        ELUTANT_SEP_TIME = 1
        ELUTION_DEAD_VOL = 2


        ### Errors
        if sample_number > 48:
            raise ValueError('sample number cannot exceed 48')

        ### Define Labware
        # Mix Plate
        mix_plate = protocol.load_labware(MIX_PLATE_TYPE, MIX_PLATE_POSITION)

        # Reagents
        reagent_container = protocol.load_labware(REAGENT_CONTAINER_TYPE, REAGENT_CONTAINER_POSITION)

        # Beads Container
        bead_container = protocol.load_labware(BEAD_CONTAINER_TYPE, BEAD_CONTAINER_POSITION)

        ### Calculating Columns
        # Total number of columns
        col_num = sample_number // 8 + (1 if sample_number % 8 > 0 else 0)
        print('There will be '+ str(col_num)+' of columns contain sample.')

        # Columns containing samples in location 1 (magentic module)
            # generates a list of lists: [[A1, B1, C1...], [A2, B2, C2...]...]
        samples = [col for col in mag_plate.columns()[sample_offset : col_num + sample_offset]]

        # Columns to mix beads and samples in location 4 (mix plate)
        mixing = [col for col in mix_plate.columns()[sample_offset:col_num + sample_offset]]

        # Columns to dispense output in location 1 (magnetic module)
            # purified parts are dispensed 6 rows to the right of their initial location
            # this is why the number of samples cannot exceed 48

        output = [col for col in mag_plate.columns()[6 + sample_offset:col_num + 6 + sample_offset]]

        ### Defining Wells for Reagents, Liquid Waste, and Beads

        #liquid_waste = reagent_container.wells(LIQUID_WASTE_WELL)
        ethanol = reagent_container.wells(ethanol_well)
        #elution_buffer = reagent_container.wells(elution_buffer_well)
        beads = bead_container[BEADS_WELL]

        ### Define bead and mix volume
        bead_volume = sample_volume * bead_ratio
        if bead_volume / 2 > pipette.max_volume:
            mix_vol = pipette.max_volume
        else:
            mix_vol = bead_volume / 2
        total_vol = bead_volume + sample_volume + DEAD_TOTAL_VOL 


        ### Steps

        # Mix beads and parts
        for target in range(col_num):


            # relative rates for fine-tuning pipetting steps
            high = 1.5
            normal = 1
            slow = 0.4
            vslow = 0.2

            pipette.pick_up_tip()
            # Aspirate beads
            #pipette.pick_up_tip(tiprack_200_1["A1"])
            #for row in samples[target]:
                #pipette.aspirate(bead_volume / len(samples[target]), beads)
            pipette.aspirate(bead_volume, beads.bottom(2), rate=normal)
            #protocol.max_speeds.update(SLOW_HEAD_SPEEDS)

            # Aspirte samples into same tip (saves tips)
            pipette.aspirate(sample_volume + DEAD_TOTAL_VOL, samples[target][0], rate=normal)

            # Transfer and mix on mix_plate
            pipette.dispense(total_vol, mixing[target][0], rate=normal)
                # similar to above, added [0] because samples[target] returned a list of every well in column 1 rather than just one well
            pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, mixing[target][0],rate=high)
                # similar to above, added [0] because samples[target] returned a list of every well in column 1 rather than just one well
            pipette.blow_out()

            # Dispose of tip
            pipette.default_speed = 50
            #protocol.max_speeds.update(DEFAULT_HEAD_SPEEDS)
            pipette.drop_tip(trash)

        # Immobilise sample
        protocol.delay(minutes=incubation_time)

        # Transfer beads+samples back to magblock

        for target in range(len(samples)):
            pipette.pick_up_tip()
            #pipette.pick_up_tip(tiprack_200_1["A7"])
            pipette.aspirate(total_vol, mixing[target][0])  
            pipette.dispense(total_vol, samples[target][0])
            pipette.blow_out()
            pipette.drop_tip()
            #pipette.transfer(total_vol, mixing[target], samples[target], blow_out=True, blowout_location='destination well')
            # added blowout_location=destination well because default location of blowout is waste in API version 2

        # Engagae MagDeck and incubate
        #MAGDECK.engage(height_from_base=MAGDECK_HEIGHT) 
        #modified from 2.14 version MAGDECK.engage(height=MAGDECK_HEIGHT)
        protocol.delay(minutes=settling_time)

        # Remove supernatant from magnetic beads
        for target in range(len(samples)):
            pipette.pick_up_tip()
            #pipette.pick_up_tip(tiprack_200_2["A1"])
            pipette.aspirate(total_vol, samples[target][0])
            pipette.dispense(total_vol, reagent_container['A5'] )
            #protocol.delay(seconds=7)
            pipette.blow_out(reagent_container['A5'] )
            #pipette.transfer(total_vol, target, liquid_waste, blow_out=True)
            pipette.drop_tip()
            
        # Wash beads twice with 70% ethanol
        
        air_vol = pipette.max_volume * AIR_VOL_COEFF
        # for cycle in range(2):
        #     for target in samples:
        #         if robot_type=='Flex'
        #           pipette.pick_up_tip(tiprack_1000['A1'])
        #         else:
        #             pipette.pick_up_tip()
        #         pipette.distribute(ETHANOL_VOL, ethanol.bottom(2), target.bottom(5), air_gap=air_vol, new_tip='never')
        #         pipette.return_tip()
        #         #Reuse the tip since it only comes into contact with ethanol. This approach reduces the number of tips needed for ethanol purification by half, saving a total of 96 tips when processing 48 samples.
                
        #     protocol.delay(minutes=WASH_TIME)
            
        #     for target in range(len(samples)):
        #         #pipette.pick_up_tip()
        #         if robot_type=='Flex'
        #           pipette.pick_up_tip(tiprack_1000['A1'])
        #         else:
        #             pipette.pick_up_tip()
        #         #Tell pipette to restart from A1 column. It will automatically pick up tips from A2 column otherwise.
        #         #pipette.pick_up_tip(tiprack_1000)
        #         pipette.aspirate(ETHANOL_VOL + ETHANOL_DEAD_VOL, samples[target][0]) 
        #         pipette.air_gap(air_vol) 
        #         pipette.dispense(ETHANOL_VOL + ETHANOL_DEAD_VOL + air_vol, reagent_container['A5'].top(5))
        #         pipette.blow_out(reagent_container['A5'] )
        #         #pipette.transfer(ETHANOL_VOL + ETHANOL_DEAD_VOL, target, liquid_waste, air_gap=air_vol)
        #         pipette.drop_tip()
        tip_index = 0  # Track tip positions by row for multi-channel pipette

        def get_flex_tip(tip_index):
            global tip_index
            if tip_index >= len(tiprack_1000.rows()[0]):  # Ensure we don't exceed available tips
                raise IndexError("No more tips available in tiprack_1000")
            tip_pos = tiprack_1000.rows()[0][tip_index]  # Select tips by row for multi-channel pipette
            tip_index += 1  # Move to the next tip for the next cycle
            return tip_pos

        # Dynamic tip handling for OT-2 standard tipracks, ensuring multi-channel row indexing
        def get_ot2_tip(tip_index):
            global tip_index
            if tip_index >= len(pipette.tip_racks[0].rows()[0]):  # Ensure we don't exceed available tips
                raise IndexError("No more tips available in OT-2 tip rack")
            tip_pos = pipette.tip_racks[0].rows()[0][tip_index]  # Select tips by row for multi-channel pipette
            tip_index += 1  # Move to the next tip for the next cycle
            return tip_pos

        for cycle in range(2):
            if robot_type == 'Flex':
                tip = get_flex_tip()  # Get a new tip for each cycle
            else:
                tip = get_ot2_tip()  # Get a new tip for OT-2 with correct row indexing
            
            for target in samples:
                if robot_type == 'Flex':
                    pipette.pick_up_tip(tip)  # Reuse the same tip for this cycle
                else:
                    pipette.pick_up_tip(tip)
                
                pipette.distribute(ETHANOL_VOL, ethanol.bottom(2), target.bottom(5), air_gap=air_vol, new_tip='never')
                pipette.return_tip()  # Return tip for reuse within the cycle
            
            protocol.delay(minutes=WASH_TIME)
            
            for target in range(len(samples)):
                if robot_type == 'Flex':
                    pipette.pick_up_tip(tip)  # Reuse the same tip for this cycle
                else:
                    pipette.pick_up_tip(tip)
                
                pipette.aspirate(ETHANOL_VOL + ETHANOL_DEAD_VOL, samples[target][0]) 
                pipette.air_gap(air_vol) 
                pipette.dispense(ETHANOL_VOL + ETHANOL_DEAD_VOL + air_vol, reagent_container['A5'].top(5))
                pipette.blow_out(reagent_container['A5'])
                pipette.drop_tip()  # Drop the tip at the end of the cycle

        # Dry at room temperature
        protocol.delay(minutes=drying_time)

        # Disengage MagDeck
        #MAGDECK.disengage()

        # Mix beads with elution buffer
        if elution_buffer_volume / 2 > pipette.max_volume:
            mix_vol = pipette.max_volume
        else:   
            mix_vol = elution_buffer_volume / 2
            
        for target in range(len(samples)):
            #pipette.transfer(elution_buffer_volume, elution_buffer, target, mix_after=(ELUTION_MIX_REPS, mix_vol))
            pipette.pick_up_tip()
            pipette.aspirate(elution_buffer_volume, reagent_container[elution_buffer_well] )
            pipette.dispense(elution_buffer_volume, samples[target][0])
            for _ in range(ELUTION_MIX_REPS): 
                pipette.aspirate(mix_vol, samples[target][0])
                pipette.dispense(mix_vol, samples[target][0])
            pipette.blow_out()
            pipette.drop_tip()

        # Incubate at room temperature
        protocol.delay(minutes=elution_time)

        # Engage MagDeck (remains engaged for DNA elution)
        #MAGDECK.engage(height_from_base=MAGDECK_HEIGHT)
        protocol.delay(minutes=ELUTANT_SEP_TIME)

        # Transfer purified parts to a new well        
        for target, dest in zip(samples, output):
            pipette.pick_up_tip()
            #pipette.transfer(elution_buffer_volume - ELUTION_DEAD_VOL, target, dest, blow_out=False)
            pipette.aspirate(elution_buffer_volume - ELUTION_DEAD_VOL, target[0])
            pipette.dispense(elution_buffer_volume - ELUTION_DEAD_VOL, dest[0])
            pipette.drop_tip()
        # Disengage MagDeck
       # MAGDECK.disengage()

    magbead(sample_number=sample_number, ethanol_well=ethanol_well)
    # removed elution buffer well='A1', added that to where the function is defined

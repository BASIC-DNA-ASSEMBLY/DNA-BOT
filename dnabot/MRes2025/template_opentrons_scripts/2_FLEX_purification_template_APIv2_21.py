''''#This code is modified from ot2 version to compatiable with Flex machine manually by Wenhan Hu (CID:02017595),
as the second step in the BASIC workflow.'''

from opentrons import protocol_api

# Rename to 'purification_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

metadata = {
     'protocolName': 'DNABOT Step 2: Purification',
     'description': 'Implements magbead purification reactions for BASIC assembly using an opentrons flex'}


# requirements
requirements = {"robotType": "Flex", "apiLevel": "2.20"}

# example values produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template:
#sample_number=8
#ethanol_well='A3'

# __LABWARES and __PARAMETERS are expected to be redefined by "generate_ot2_script" method
# Test dict
# __LABWARES={"p20_single": {"id": "p20_single_gen2"}, "p300_multi": {"id": "p300_multi_gen2"}, "mag_deck": {"id": "magneticModuleV1"}, "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, "96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "96_wellplate_200ul_pcr_step_14": {"id": "4ti0960rig_96_wellplate_200ul"}, "96_wellplate_200ul_pcr_step_23": {"id": "4ti0960rig_96_wellplate_200ul"}, "agar_plate_step_4": {"id": "4ti0960rig_96_wellplate_200ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}}
# __PARAMETERS={"purif_magdeck_height": {"value": 20.0}, "purif_wash_time": {"value": 0.5}, "purif_bead_ratio": {"value": 1.8}, "purif_incubation_time": {"value": 5.0}, "purif_settling_time": {"value": 2.0}, "purif_drying_time": {"value": 5.0}, "purif_elution_time": {"value": 2.0}, "transfo_incubation_temp": {"value": 4.0}, "transfo_incubation_time": {"value": 20.0}}

sample_number=6 # Vary from 0 ~ 48
ethanol_well='A11'

__LABWARES={"p50_single": {"id": "flex_1channel_50"},
            "p50_multi": {"id": "flex_8channel_50"}, 
            "p1000_single": {"id": "flex_1channel_1000"}, 
            "p1000_multi": {"id": "flex_8channel_1000"},
            "mag_block": {"id": "magneticBlockV1"},
            "mag_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "flex_96_tiprack_50ul": {"id": "opentrons_flex_96_tiprack_50ul"}, 
            "flex_96_tiprack_200ul": {"id": "opentrons_flex_96_tiprack_200ul"},
            "flex_96_tiprack_1000ul": {"id": "opentrons_flex_96_tiprack_1000ul"}, 
            "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, 
            "clip_source_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "clip_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "mix_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "final_assembly_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "transfo_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "transfo_plate_wo_thermo": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "agar_plate": {"id": "nest_96_wellplate_100ul_pcr_full_skirt"}, 
            "12_reservoir_21000ul": {"id": "nest_12_reservoir_15ml"}, 
            "96_deepwellplate_2ml": {"id": "nest_96_wellplate_2ml_deep"}, 
            "12_corning_wellplate": {"id": "corning_12_wellplate_6.9ml_flat"}}

__PARAMETERS={"clip_keep_thermo_lid_closed": {"value": "No", "id": "No"}, 
              "premix_linkers": {"value": "Yes", "id": "No"}, 
              "premix_parts": {"value": "Yes", "id": "Yes"}, 
              "linkers_volume": {"value": 20}, 
              "parts_volume": {"value": 20}, 
              "thermo_temp": {"value": 4}, 
              #"purif_magdeck_height": {"value": 10.8}, 
              "purif_wash_time": {"value": 0.5}, 
              "purif_bead_ratio": {"value": 1.8}, 
              "purif_incubation_time": {"value": 5}, 
              "purif_settling_time": {"value": 2}, 
              "purif_drying_time": {"value": 5}, 
              "purif_elution_time": {"value": 2}, 
              "transfo_incubation_temp": {"value": 4}, 
              "transfo_incubation_time": {"value": 20}}

def run(protocol: protocol_api.ProtocolContext):
    trash = protocol.load_trash_bin(location="A3")# Trash defination for Flex
    
# added run function for API verison 2

    def magbead(
            sample_number,
            ethanol_well,
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
                # this template was written with the Gen 1 magnetic module, as it is compatible with API version 2
            drying_time=__PARAMETERS['purif_drying_time']['value'],
            elution_time=__PARAMETERS['purif_elution_time']['value'],
            sample_offset=0,
            #tiprack_type=__LABWARES['96_tiprack_300ul']['id']):
            tiprack_200=__LABWARES['flex_96_tiprack_200ul']['id'],            
            tiprack_1000=__LABWARES['flex_96_tiprack_1000ul']['id']):
        

        """

        Selected args:
            ethanol_well (str): well in reagent container containing ethanol.
            elution_buffer_well (str): well in reagent container containing elution buffer.
            sample_offset (int): offset the intial sample column by the specified value.

        """


        ### Constants

        # Pipettes
        PIPETTE_ASPIRATE_RATE = 25
        PIPETTE_DISPENSE_RATE = 150
##
        TIPS_PER_SAMPLE = 5
        TIPS_WASH = 2

        #PIPETTE_TYPE = __LABWARES['p300_multi']['id']
        PIPETTE_TYPE = __LABWARES['p1000_multi']['id']

            # new constant for easier swapping between pipette types

        # Tiprack
        #CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5']
        #CANDIDATE_TIPRACK_SLOTS_200 = ['D3', 'C3', 'B3']
        #CANDIDATE_TIPRACK_SLOT_1000 = 'C2'
        tiprack_200_1= protocol.load_labware(tiprack_200, 'D3')
        tiprack_200_2= protocol.load_labware(tiprack_200, 'C3')
        tiprack_200_3= protocol.load_labware(tiprack_200, 'B3')  # 200  ul tip used for Asperation, Transfer, Elusion and so on.
        tiprack_1000 = protocol.load_labware(tiprack_1000, 'C2') # 1000 ul tip used for ethanol wash 
        # Magnetic Module
        #MAGDECK_POSITION = '1' Magnetic Block Updated
        MAGNETIC_BLOCK_TYPE = __LABWARES['mag_block']['id']  # Updated from "magneticBlockV1"
        MAGNETIC_BLOCK_POSITION = 'D1'
        MAGNETIC_PLATE_TYPE = __LABWARES['mag_plate']['id']

        # Mix Plate
        MIX_PLATE_TYPE = __LABWARES['mix_plate']['id']
        MIX_PLATE_POSITION = 'C1' # MIX_PLATE_POSITION = '4' updated for flex now
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
            # also acts as the type of plate loaded onto the magnetic module

        # Reagents
        REAGENT_CONTAINER_TYPE = __LABWARES['12_reservoir_21000ul']['id']
        REAGENT_CONTAINER_POSITION = 'A2' #REAGENT_CONTAINER_POSITION = '7' updated for flex now(Notice!! B1 is not available,so A2 is used)

        # Beads
        BEAD_CONTAINER_TYPE = __LABWARES['96_deepwellplate_2ml']['id']
        BEAD_CONTAINER_POSITION = 'B2' #BEAD_CONTAINER_POSITION = '8' updated for flex now

        # Settings
        #LIQUID_WASTE_WELL = 'A5'
        BEADS_WELL = 'A1'
        DEAD_TOTAL_VOL = 5
        #SLOW_HEAD_SPEEDS = {'x': 600 // 4, 'y': 400 // 4, 'z': 125 // 10, 'a': 125 // 10}
        #DEFAULT_HEAD_SPEEDS = {'x': 400, 'y': 400, 'z': 125, 'a': 100}
        IMMOBILISE_MIX_REPS = 10
        #MAGDECK_HEIGHT = __PARAMETERS['purif_magdeck_height']['value']
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


        ### Loading Tiprack

        # Calculates whether one/two/three/four/five tipracks are needed, which are in slots 3, 6, 9, 2, and 5 respectively
        total_tips = sample_number * TIPS_PER_SAMPLE
        tiprack_num = total_tips // 96 + (1 if total_tips % 96 > 0 else 0)
        print(str(tiprack_num) + 'of 200ul tipboxs is needed')


            # changed to protocol.load_labware for API version 2
            
        wash_tips = sample_number * TIPS_WASH
        tiprack_wash = wash_tips // 96 + (1 if total_tips % 96 > 0 else 0)
        print(str(tiprack_wash) + 'of 1000ul tipboxs is needed')
        
        ### Loading Pipettes

        pipette = protocol.load_instrument(PIPETTE_TYPE, mount="left",tip_racks=[tiprack_200_1,tiprack_200_2,tiprack_200_3,tiprack_1000])
        pipette.aspirate_flow_rate=PIPETTE_ASPIRATE_RATE
        pipette.dispense_flow_rate=PIPETTE_DISPENSE_RATE
            # for reference: default aspirate/dispense flow rate for p300_multi_gen2 is 94 ul/s

        ### Define Labware

        # Magnetic Block used in flex
        #MAGDECK = protocol.load_module(__LABWARES['mag_deck']['id'], location= MAGDECK_POSITION)
        magnetic_block = protocol.load_module(module_name=MAGNETIC_BLOCK_TYPE, location=MAGNETIC_BLOCK_POSITION)
        # Load a 96-well plate on the Magnetic Block
        mag_plate = magnetic_block.load_labware(name=MAGNETIC_PLATE_TYPE)

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

            # Aspirate beads
            #pipette.pick_up_tip(tiprack_200_1["A1"])
            pipette.pick_up_tip()
            pipette.default_speed = 50
            
            #for row in samples[target]:
                #pipette.aspirate(bead_volume / len(samples[target]), beads)
            pipette.aspirate(bead_volume, beads)
            #protocol.max_speeds.update(SLOW_HEAD_SPEEDS)

            # Aspirte samples
            pipette.aspirate(sample_volume + DEAD_TOTAL_VOL, samples[target][0])

            # Transfer and mix on mix_plate
            pipette.dispense(total_vol, mixing[target][0])
                # similar to above, added [0] because samples[target] returned a list of every well in column 1 rather than just one well
            pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, mixing[target][0])
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
            pipette.blow_out(reagent_container['A5'] )
            #pipette.transfer(total_vol, target, liquid_waste, blow_out=True)
            pipette.drop_tip()
            
        # Wash beads twice with 70% ethanol
        
        air_vol = pipette.max_volume * AIR_VOL_COEFF
        for cycle in range(2):
            for target in samples:
                #pipette.pick_up_tip()
                #pipette.pick_up_tip(tiprack_1000)
                pipette.pick_up_tip(tiprack_1000['A1'])
                #Tell p1000_multi to switch mode and pick up p1000 instead of p200
                pipette.distribute(ETHANOL_VOL, ethanol, target, air_gap=air_vol, new_tip='never')
                pipette.return_tip()
                #Reuse the tip since it only comes into contact with ethanol. This approach reduces the number of tips needed for ethanol purification by half, saving a total of 96 tips when processing 48 samples.
                
            protocol.delay(minutes=WASH_TIME)
            
            for target in range(len(samples)):
                #pipette.pick_up_tip()
                pipette.pick_up_tip(tiprack_1000['A1'])
                #Tell pipette to restart from A1 column. It will automatically pick up tips from A2 column otherwise.
                #pipette.pick_up_tip(tiprack_1000)
                pipette.aspirate(ETHANOL_VOL + ETHANOL_DEAD_VOL, samples[target][0]) 
                pipette.air_gap(air_vol) 
                pipette.dispense(ETHANOL_VOL + ETHANOL_DEAD_VOL + air_vol, reagent_container['A5'] )
                pipette.blow_out(reagent_container['A5'] )
                #pipette.transfer(ETHANOL_VOL + ETHANOL_DEAD_VOL, target, liquid_waste, air_gap=air_vol)
                pipette.drop_tip()

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
            pipette.aspirate(elution_buffer_volume, reagent_container['A1'] )
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
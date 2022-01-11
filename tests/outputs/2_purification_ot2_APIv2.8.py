from opentrons import protocol_api

# Rename to 'purification_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

metadata = {
     'apiLevel': '2.8',
     'protocolName': 'purification_template',
     'description': 'Implements magbead purification reactions for BASIC assembly using an opentrons OT-2'}




# example values produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template:
#sample_number=8
#ethanol_well='A3'

# __LABWARES and __PARAMETERS are expected to be redefined by "generate_ot2_script" method
# Test dict
# __LABWARES={"p20_single": {"id": "p20_single_gen2"}, "p300_multi": {"id": "p300_multi_gen2"}, "mag_deck": {"id": "magdeck"}, "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, "96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "96_wellplate_200ul_pcr_step_14": {"id": "4ti0960rig_96_wellplate_200ul"}, "96_wellplate_200ul_pcr_step_23": {"id": "4ti0960rig_96_wellplate_200ul"}, "agar_plate_step_4": {"id": "4ti0960rig_96_wellplate_200ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}}
# __PARAMETERS={"purif_magdeck_height": {"value": 20.0}, "purif_wash_time": {"value": 0.5}, "purif_bead_ratio": {"value": 1.8}, "purif_incubation_time": {"value": 5.0}, "purif_settling_time": {"value": 2.0}, "purif_drying_time": {"value": 5.0}, "purif_elution_time": {"value": 2.0}, "transfo_incubation_temp": {"value": 4.0}, "transfo_incubation_time": {"value": 20.0}}

sample_number=13
ethanol_well='A11'
__LABWARES={"p20_single": {"id": "p20_single_gen2"}, "p300_multi": {"id": "p300_multi_gen2"}, "mag_deck": {"id": "magdeck"}, "96_tiprack_20ul": {"id": "opentrons_96_tiprack_20ul"}, "96_tiprack_300ul": {"id": "opentrons_96_tiprack_300ul"}, "24_tuberack_1500ul": {"id": "e14151500starlab_24_tuberack_1500ul"}, "96_wellplate_200ul_pcr_step_14": {"id": "4ti0960rig_96_wellplate_200ul"}, "96_wellplate_200ul_pcr_step_23": {"id": "4ti0960rig_96_wellplate_200ul"}, "agar_plate_step_4": {"id": "4ti0960rig_96_wellplate_200ul"}, "12_reservoir_21000ul": {"id": "4ti0131_12_reservoir_21000ul"}, "96_deepwellplate_2ml": {"id": "4ti0136_96_wellplate_2200ul"}}
__PARAMETERS={"purif_magdeck_height": {"value": 20}, "purif_wash_time": {"value": 0}, "purif_bead_ratio": {"value": 1}, "purif_incubation_time": {"value": 5}, "purif_settling_time": {"value": 2}, "purif_drying_time": {"value": 5}, "purif_elution_time": {"value": 2}, "transfo_incubation_temp": {"value": 4}, "transfo_incubation_time": {"value": 20}}


def run(protocol: protocol_api.ProtocolContext):
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
            tiprack_type=__LABWARES['96_tiprack_300ul']['id']):

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
        TIPS_PER_SAMPLE = 9
        PIPETTE_TYPE = __LABWARES['p300_multi']['id']
            # new constant for easier swapping between pipette types

        # Tiprack
        CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5']

        # Magnetic Module
        MAGDECK_POSITION = '1'

        # Mix Plate
        MIX_PLATE_TYPE = __LABWARES['96_wellplate_200ul_pcr_step_23']['id']
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
            # also acts as the type of plate loaded onto the magnetic module
        MIX_PLATE_POSITION = '4'

        # Reagents
        REAGENT_CONTAINER_TYPE = __LABWARES['12_reservoir_21000ul']['id']
        REAGENT_CONTAINER_POSITION = '7'

        # Beads
        BEAD_CONTAINER_TYPE = __LABWARES['96_deepwellplate_2ml']['id']
        BEAD_CONTAINER_POSITION = '8'

        # Settings
        LIQUID_WASTE_WELL = 'A5'
        BEADS_WELL = 'A1'
        DEAD_TOTAL_VOL = 5
        SLOW_HEAD_SPEEDS = {'x': 600 // 4, 'y': 400 // 4, 'z': 125 // 10, 'a': 125 // 10}
        DEFAULT_HEAD_SPEEDS = {'x': 400, 'y': 400, 'z': 125, 'a': 100}
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


        ### Loading Tiprack

        # Calculates whether one/two/three/four/five tipracks are needed, which are in slots 3, 6, 9, 2, and 5 respectively
        total_tips = sample_number * TIPS_PER_SAMPLE
        tiprack_num = total_tips // 96 + (1 if total_tips % 96 > 0 else 0)
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            # changed to protocol.load_labware for API version 2


        ### Loading Pipettes

        pipette = protocol.load_instrument(PIPETTE_TYPE, mount="left", tip_racks=tipracks)
        pipette.aspirate_flow_rate=PIPETTE_ASPIRATE_RATE
        pipette.dispense_flow_rate=PIPETTE_DISPENSE_RATE
            # for reference: default aspirate/dispense flow rate for p300_multi_gen2 is 94 ul/s

        ### Define Labware

        # Magnetic Module
        MAGDECK = protocol.load_module(__LABWARES['mag_deck']['id'], MAGDECK_POSITION)
            # 'magdeck' is the gen 1 magnetic module, use 'magnetic module gen2' for the gen 2 magentic module
                # if using gen 2 module, need to change settling time! (see comments under Constants)
        MAGDECK.disengage()
            # disengages the magnets when it is turned on
        mag_plate = MAGDECK.load_labware(MIX_PLATE_TYPE)

        # Mix Plate
        mix_plate = protocol.load_labware(MIX_PLATE_TYPE, MIX_PLATE_POSITION)

        # Reagents
        reagent_container = protocol.load_labware(REAGENT_CONTAINER_TYPE, REAGENT_CONTAINER_POSITION)

        # Beads Container
        bead_container = protocol.load_labware(BEAD_CONTAINER_TYPE, BEAD_CONTAINER_POSITION)


        ### Calculating Columns

        # Total number of columns
        col_num = sample_number // 8 + (1 if sample_number % 8 > 0 else 0)

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

        liquid_waste = reagent_container.wells(LIQUID_WASTE_WELL)
        ethanol = reagent_container.wells(ethanol_well)
        elution_buffer = reagent_container.wells(elution_buffer_well)
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
        for target in range(int(len(samples))):

            # Aspirate beads
            pipette.pick_up_tip()
            pipette.aspirate(bead_volume, beads)
            protocol.max_speeds.update(SLOW_HEAD_SPEEDS)

            # Aspirte samples
            pipette.aspirate(sample_volume + DEAD_TOTAL_VOL, samples[target][0])

            # Transfer and mix on mix_plate
            pipette.dispense(total_vol, mixing[target][0])
                # similar to above, added [0] because samples[target] returned a list of every well in column 1 rather than just one well
            pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, mixing[target][0])
                # similar to above, added [0] because samples[target] returned a list of every well in column 1 rather than just one well
            pipette.blow_out()

            # Dispose of tip
            protocol.max_speeds.update(DEFAULT_HEAD_SPEEDS)
            pipette.drop_tip()

        # Immobilise sample
        protocol.delay(minutes=incubation_time)

        # Transfer beads+samples back to magdeck
        for target in range(int(len(samples))):
            pipette.transfer(total_vol, mixing[target], samples[target], blow_out=True, blowout_location='destination well')
            # added blowout_location=destination well because default location of blowout is waste in API version 2

        # Engagae MagDeck and incubate
        MAGDECK.engage(height=MAGDECK_HEIGHT)
        protocol.delay(minutes=settling_time)

        # Remove supernatant from magnetic beads
        for target in samples:
            pipette.transfer(total_vol, target, liquid_waste, blow_out=True)

        # Wash beads twice with 70% ethanol
        air_vol = pipette.max_volume * AIR_VOL_COEFF
        for cycle in range(2):
            for target in samples:
                pipette.transfer(ETHANOL_VOL, ethanol, target, air_gap=air_vol)
            protocol.delay(minutes=WASH_TIME)
            for target in samples:
                pipette.transfer(ETHANOL_VOL + ETHANOL_DEAD_VOL, target, liquid_waste, air_gap=air_vol)

        # Dry at room temperature
        protocol.delay(minutes=drying_time)

        # Disengage MagDeck
        MAGDECK.disengage()

        # Mix beads with elution buffer
        if elution_buffer_volume / 2 > pipette.max_volume:
            mix_vol = pipette.max_volume
        else:
            mix_vol = elution_buffer_volume / 2
        for target in samples:
            pipette.transfer(elution_buffer_volume, elution_buffer, target, mix_after=(ELUTION_MIX_REPS, mix_vol))

        # Incubate at room temperature
        protocol.delay(minutes=elution_time)

        # Engage MagDeck (remains engaged for DNA elution)
        MAGDECK.engage(height=MAGDECK_HEIGHT)
        protocol.delay(minutes=ELUTANT_SEP_TIME)

        # Transfer purified parts to a new well
        for target, dest in zip(samples, output):
            pipette.transfer(elution_buffer_volume - ELUTION_DEAD_VOL, target,
                             dest, blow_out=False)

        # Disengage MagDeck
        MAGDECK.disengage()

    magbead(sample_number=sample_number, ethanol_well=ethanol_well)
    # removed elution buffer well='A1', added that to where the function is defined

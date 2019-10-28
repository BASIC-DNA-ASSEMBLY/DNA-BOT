from opentrons import labware, instruments, modules, robot


sample_number=38
ethanol_well='A11'


def magbead(
        sample_number,
        ethanol_well,
        elution_buffer_well,
        sample_volume=30,
        bead_ratio=1.8,
        elution_buffer_volume=40,
        incubation_time=5,
        settling_time=2,
        drying_time=5,
        elution_time=2,
        sample_offset=0,
        tiprack_type="opentrons_96_tiprack_300ul"):
    """Implements magbead purification reactions for BASIC assembly using an opentrons OT-2.

    Selected args:
        ethanol_well (str): well in reagent container containing ethanol.
        elution_buffer_well (str): well in reagent container containing elution buffer.
        sample_offset (int): offset the intial sample column by the specified value.

    """

    # Constants
    PIPETTE_ASPIRATE_RATE = 25
    PIPETTE_DISPENSE_RATE = 150
    TIPS_PER_SAMPLE = 9
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5']
    MAGDECK_POSITION = '1'
    MIX_PLATE_TYPE = '4ti-0960_FrameStar'
    MIX_PLATE_POSITION = '4'
    REAGENT_CONTAINER_TYPE = '4ti0131_trough-12'
    REAGENT_CONTAINER_POSITION = '7'
    BEAD_CONTAINER_TYPE = '4ti0136_96_deep-well'
    BEAD_CONTAINER_POSITION = '8'
    LIQUID_WASTE_WELL = 'A12'
    BEADS_WELL = 'A1'
    DEAD_TOTAL_VOL = 5
    SLOW_HEAD_SPEEDS = {'x': 600 // 4, 'y': 400 // 4,
                        'z': 125 // 10, 'a': 125 // 10}
    DEFAULT_HEAD_SPEEDS = {'x': 400, 'y': 400, 'z': 125, 'a': 100}
    IMMOBILISE_MIX_REPS = 10
    MAGDECK_HEIGHT = 20
    AIR_VOL_COEFF = 0.1
    ETHANOL_VOL = 150
    WASH_TIME = 0.5
    ETHANOL_DEAD_VOL = 50
    ELUTION_MIX_REPS = 20
    ELUTANT_SEP_TIME = 1
    ELUTION_DEAD_VOL = 2

    # Errors
    if sample_number > 48:
        raise ValueError('sample number cannot exceed 48')

    # Tips and pipette
    total_tips = sample_number * TIPS_PER_SAMPLE
    tiprack_num = total_tips // 96 + (1 if total_tips % 96 > 0 else 0)
    slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
    tipracks = [labware.load(tiprack_type, slot)
                for slot in slots]
    pipette = instruments.P300_Multi(
        mount="left",
        tip_racks=tipracks,
        aspirate_flow_rate=PIPETTE_ASPIRATE_RATE,
        dispense_flow_rate=PIPETTE_DISPENSE_RATE)

    # Define labware
    MAGDECK = modules.load('magdeck', MAGDECK_POSITION)
    MAGDECK.disengage()
    mag_plate = labware.load(MIX_PLATE_TYPE, MAGDECK_POSITION, share=True)
    mix_plate = labware.load(MIX_PLATE_TYPE, MIX_PLATE_POSITION)
    reagent_container = labware.load(
        REAGENT_CONTAINER_TYPE, REAGENT_CONTAINER_POSITION)
    bead_container = labware.load(BEAD_CONTAINER_TYPE, BEAD_CONTAINER_POSITION)
    col_num = sample_number // 8 + (1 if sample_number % 8 > 0 else 0)
    samples = [col for col in mag_plate.cols(
    )[0 + sample_offset:col_num + sample_offset]]
    output = [col for col in mag_plate.cols(
    )[6 + sample_offset:col_num + 6 + sample_offset]]
    mixing = [col for col in mix_plate.cols(
    )[0 + sample_offset:col_num + sample_offset]]

    # Define reagents and liquid waste
    liquid_waste = reagent_container.wells(LIQUID_WASTE_WELL)
    beads = bead_container.wells(BEADS_WELL)
    ethanol = reagent_container.wells(ethanol_well)
    elution_buffer = reagent_container.wells(elution_buffer_well)

    # Define bead and mix volume
    bead_volume = sample_volume * bead_ratio
    if bead_volume / 2 > pipette.max_volume:
        mix_vol = pipette.max_volume
    else:
        mix_vol = bead_volume / 2
    total_vol = bead_volume + sample_volume + DEAD_TOTAL_VOL

    # Mix beads and PCR samples and incubate
    for target in range(int(len(samples))):
        # Aspirate beads
        pipette.pick_up_tip()
        pipette.aspirate(bead_volume, beads)
        robot.head_speed(**SLOW_HEAD_SPEEDS, combined_speed=max(SLOW_HEAD_SPEEDS.values()))

        # Transfer and mix on  mix_plate
        pipette.aspirate(sample_volume + DEAD_TOTAL_VOL, samples[target])
        pipette.dispense(total_vol, mixing[target])
        pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, mixing[target])
        pipette.blow_out()

        # Dispose of tip
        robot.head_speed(**DEFAULT_HEAD_SPEEDS, combined_speed=max(DEFAULT_HEAD_SPEEDS.values()))
        pipette.drop_tip()

    # Immobilise sample
    pipette.delay(minutes=incubation_time)

    # Transfer sample back to magdeck
    for target in range(int(len(samples))):
        pipette.transfer(total_vol, mixing[target], samples[target],
                         blow_out=True)

    # Engagae MagDeck and incubate
    MAGDECK.engage(height=MAGDECK_HEIGHT)
    pipette.delay(minutes=settling_time)

    # Remove supernatant from magnetic beads
    for target in samples:
        pipette.transfer(total_vol, target, liquid_waste, blow_out=True)

    # Wash beads twice with 70% ethanol
    air_vol = pipette.max_volume * AIR_VOL_COEFF
    for cycle in range(2):
        for target in samples:
            pipette.transfer(ETHANOL_VOL, ethanol, target, air_gap=air_vol)
        pipette.delay(minutes=WASH_TIME)
        for target in samples:
            pipette.transfer(ETHANOL_VOL + ETHANOL_DEAD_VOL, target, liquid_waste,
                             air_gap=air_vol)

    # Dry at RT
    pipette.delay(minutes=drying_time)

    # Disengage MagDeck
    MAGDECK.disengage()

    # Mix beads with elution buffer
    if elution_buffer_volume / 2 > pipette.max_volume:
        mix_vol = pipette.max_volume
    else:
        mix_vol = elution_buffer_volume / 2
    for target in samples:
        pipette.transfer(elution_buffer_volume, elution_buffer,
                         target, mix_after=(ELUTION_MIX_REPS, mix_vol))

    # Incubate at RT for "elution_time" minutes
    pipette.delay(minutes=elution_time)

    # Engagae MagDeck for 1 minute and remain engaged for DNA elution
    MAGDECK.engage(height=MAGDECK_HEIGHT)
    pipette.delay(minutes=ELUTANT_SEP_TIME)

    # Transfer clean PCR product to a new well
    for target, dest in zip(samples, output):
        pipette.transfer(elution_buffer_volume - ELUTION_DEAD_VOL, target,
                         dest, blow_out=False)

    # Disengage MagDeck
    MAGDECK.disengage()


magbead(sample_number=sample_number,
        ethanol_well=ethanol_well, elution_buffer_well='A1')

for c in robot.commands():
    print(c)

from opentrons import labware, instruments, modules, robot
# Robot class added.

tray_1_vols = [[5, 5, 5, 5, 5, 5, 5]]


def sample_num(tray_1_vols):
    """Determines sample number from tray_1_vols input.

    Args:
    tray_1_vols (list): List of lists containing volumes of each transformation reaction to be spotted. Nested lists correspond with columns.

    """
    return len([item for sublist in tray_1_vols for item in sublist])


def final_well(tray_1_vols):
    """Determines well containing the final sample from tray_1_vols input.

    Args:
    tray_1_vols (list): List of lists containing volumes of each transformation reaction to be spotted. Nested lists correspond with columns.

    """
    letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    final_well_column = sample_num(
        tray_1_vols) // 8 + (1 if sample_num(tray_1_vols) % 8 > 0 else 0)
    final_well_row = letter[sample_num(
        tray_1_vols) - (final_well_column - 1) * 8 - 1]
    return final_well_row + str(final_well_column)


def tiprack_slots(
        tray_1_vols,
        tips_per_sample):
    """Calculates tiprack slots required.

    Args:
    tray_1_vols (list): List of lists containing volumes of each transformation reaction to be spotted. Nested lists correspond with columns.
    tips_per_sample (int): Number of tips used per sample.

    """
    total_tips = sample_num(tray_1_vols) * tips_per_sample
    tiprack_slots = total_tips // 96 + 1 if total_tips % 96 > 0 else 0
    return tiprack_slots


def transformation_setup(final_well):
    """Sets up transformation reactions

    Args:
    final_well (str): well containing the final sample. 

    """

    # Constants
    TEMP = 4  # Incubation temperature.
    ASSEMBLY_VOL = 5  # Volume of final assembly added to competent cells.
    INITIAL_WELL = 'A1'  # Initial well in assembly and transformation plates.
    MIX_SETTINGS = (4, 5)  # Mix after setting during final assembly transfers.
    INCUBATION_TIME = 20  # Cells and final assembly incubation time.

    # Set temperature deck to 4 °C and load competent cells
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    robot.pause()
    robot.comment('Load competent cells, uncap and resume run')

    # Transfer final assemblies
    if final_well == 'A1':
        p10_pipette.transfer(ASSEMBLY_VOL, assembly_plate.wells('A1'),
                             transformation_plate.wells('A1'), mix_after=(MIX_SETTINGS))
    else:
        p10_pipette.transfer(ASSEMBLY_VOL,
                             assembly_plate.wells(INITIAL_WELL, to=final_well),
                             transformation_plate.wells(
                                 INITIAL_WELL, to=final_well),
                             new_tip='always', mix_after=(MIX_SETTINGS))

    # Incubate for 20 minutes and remove competent cells for heat shock
    p10_pipette.delay(minutes=INCUBATION_TIME)
    robot.pause()
    robot.comment(
        'Remove transformation reactions for heatshock and resume run')
    tempdeck.deactivate()


def phase_switch(comment='Remove final assembly plate and introduce both agar tray(s) and a deep well plate containing SOC media. Resume run.'):
    """Function pauses run enabling addition/removal of labware.

    Args:
    comment (str): string to be displayed during run following pause.

    """
    robot.pause()
    robot.comment(comment)


def outgrowth(
        final_well,
        soc_well='A1'):
    """Outgrows transformed cells.

    Args:
    final_well (str): Well containing the final sample.
    soc_well (str): Well containing SOC media in relevant plate.

    """

    # Constants
    SOC_VOL = 125  # Volume of SOC media (uL)
    SOC_MIX_SETTINGS = (4, 50)  # mix_after settings when adding SOC
    TEMP = 37  # outgrowth incubation temperature

    # Define wells
    col_num = sample_num(tray_1_vols) // 8 + \
        (1 if sample_num(tray_1_vols) % 8 > 0 else 0)
    transformation_wells = [col for col in transformation_plate.cols()[
        :col_num]]
    soc = soc_plate.wells(soc_well)

    # Add SOC to transformed cells
    p300_pipette.transfer(SOC_VOL, soc, transformation_wells,
                          new_tip='always', mix_after=SOC_MIX_SETTINGS)

    # Incubate for 1 hour at 37 °C
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    p300_pipette.delay(minutes=60)
    tempdeck.deactivate()


def spot_transformations(
        tray_1_vols,
        agar_tray_num=1,
        tray_2_vol=8,
        dead_vol=2,
        spotting_dispense_rate=0.025,
        time_delay=0,
        stabbing_depth=2):
    """Spots transformation reactions.

    Args:
    final_well (str): Well containing the final transformation reaction.
    agar_tray_num (int): Number of agar trays transformations spotted on.
    tray_1_vols (list): List of lists containing volumes of each transformation reaction to be spotted. Nested lists correspond with columns.
    tray_2_vol (float): Fixed volume to be spotted onto agar_trays['2'].
    dead_vol (float): Dead volume aspirated during spotting.
    spotting_dispense_rate (float): Rate p10_pipette dispenses at during spotting.
    time_dealy (float): Seconds following dispensing and prior to touching tip to agar plate.
    stabbing_depth (float): Depth p10_pipette moves into agar during spotting. 

    """

    def spot(
            source,
            target,
            spot_vol):
        """Spots an individual reaction using the p10 pipette.

        Args:
        source (str): Well containing the transformation reaction to be spotted.
        target (str): Well transformation reaction is to be spotted to.
        spot_vol (float): Volume of transformation reaction to be spotted (uL).  

        """

        # Constants
        DEFAULT_HEAD_SPEED = {'x': 400, 'y': 400,
                              'z': 125, 'a': 125}
        SPOT_HEAD_SPEED = {'x': 400, 'y': 400, 'z': 125,
                           'a': 125 // 4}
        DISPENSING_HEIGHT = 5
        SAFE_HEIGHT = 15  # height avoids collision with agar tray.

        # Spot
        p10_pipette.pick_up_tip()
        p10_pipette.aspirate(spot_vol + dead_vol, source)
        p10_pipette.move_to(target.top(SAFE_HEIGHT))
        p10_pipette.move_to(target.top(DISPENSING_HEIGHT))
        p10_pipette.dispense(volume=spot_vol, rate=spotting_dispense_rate)
        p10_pipette.delay(seconds=time_delay)
        robot.head_speed(combined_speed=max(
            SPOT_HEAD_SPEED.values()), **SPOT_HEAD_SPEED)
        p10_pipette.move_to(target.top(-1 * stabbing_depth))
        robot.head_speed(combined_speed=max(
            DEFAULT_HEAD_SPEED.values()), **DEFAULT_HEAD_SPEED)
        p10_pipette.move_to(target.top(SAFE_HEIGHT))

        # Dispose of dead volume and tip
        p10_pipette.dispense(dead_vol, spotting_waste)
        p10_pipette.blow_out()
        p10_pipette.drop_tip()

    # Constants
    TRANSFORMATION_MIX_SETTINGS = [4, 50]

    # Define transformation and agar_tray columns
    col_num = sample_num(tray_1_vols) // 8 + \
        (1 if sample_num(tray_1_vols) % 8 > 0 else 0)
    transformation_cols = [col for col in transformation_plate.cols()[
        :col_num]]
    agar_tray1_cols = [col for col in agar_trays['1'].cols()[
        :col_num]]
    if agar_tray_num == 2:
        agar_tray2_cols = [col for col in agar_trays['2'].cols()[
            :col_num]]

    # Spot transformation reactions column by column
    if agar_tray_num == 1:
        for col in range(col_num):
            p300_pipette.pick_up_tip()
            p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[
                             1], transformation_cols[col])
            p300_pipette.drop_tip()
            for well in range(len(tray_1_vols[col])):
                spot(transformation_cols[col][well],
                     agar_tray1_cols[col][well], tray_1_vols[col][well])
    elif agar_tray_num == 2:
        for col in range(col_num):
            p300_pipette.pick_up_tip()
            p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[
                             1], transformation_cols[col])
            p300_pipette.drop_tip()
            for well in range(len(tray_1_vols[col])):
                spot(transformation_cols[col][well],
                     agar_tray1_cols[col][well], tray_1_vols[col][well])
            for well in range(len(tray_1_vols[col])):
                spot(transformation_cols[col][well],
                     agar_tray1_cols[col][well], tray_2_vol)
    else:
        print('agar_tray number exceeds current maximum')
        exit()


# Run protocol
if __name__ == '__main__':
    # Global variables
    agar_tray_num = 1

    # Constants
    CANDIDATE_P10_SLOTS = ['9', '2', '5']
    CANDIDATE_P300_SLOTS = ['3', '6']
    P300_TIPS_PER_SAMPLE = 2
    P10_MOUNT = 'right'
    P300_MOUNT = 'left'
    ASSEMBLY_PLATE_TYPE = '4ti-0960_FrameStar'
    ASSEMBLY_PLATE_SLOT = '8'
    TEMPDECK_SLOT = '10'
    TRANSFORMATION_PLATE_TYPE = 'aluminium-block_4ti-0753-4_PCR_strips'
    SOC_PLATE_TYPE = '4ti0136_96_deep-well'
    SOC_PLATE_SLOT = '7'
    TUBE_RACK_TYPE = 'tube-rack_E1415-1500'
    TUBE_RACK_SLOT = '11'
    SPOTTING_WASTE_WELL = 'A1'
    AGAR_TRAY_TYPE = 'Nunc_Omnitray'
    AGAR_TRAY_SLOTS = ['1', '4']

    # Tiprack slots
    p10_tips_per_sample = 1 + agar_tray_num
    p10_slots = CANDIDATE_P10_SLOTS[
        :tiprack_slots(tray_1_vols, p10_tips_per_sample)]
    p300_slots = CANDIDATE_P300_SLOTS[
        :tiprack_slots(tray_1_vols, P300_TIPS_PER_SAMPLE)]

    # Define labware
    p10_tipracks = [labware.load('tiprack-10ul', slot) for slot in p10_slots]
    p300_tipracks = [labware.load('opentrons-tiprack-300ul', slot)
                     for slot in p300_slots]
    p10_pipette = instruments.P10_Single(
        mount=P10_MOUNT, tip_racks=p10_tipracks)
    p300_pipette = instruments.P300_Multi(
        mount=P300_MOUNT, tip_racks=p300_tipracks)

    assembly_plate = labware.load(ASSEMBLY_PLATE_TYPE, ASSEMBLY_PLATE_SLOT)
    tempdeck = modules.load('tempdeck', TEMPDECK_SLOT)
    transformation_plate = labware.load(TRANSFORMATION_PLATE_TYPE,
                                        TEMPDECK_SLOT, share=True)
    soc_plate = labware.load(SOC_PLATE_TYPE, SOC_PLATE_SLOT)
    tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_SLOT)
    spotting_waste = tube_rack.wells(SPOTTING_WASTE_WELL)
    agar_trays = {}
    if agar_tray_num == 2:
        for x in range(1, len(AGAR_TRAY_SLOTS) + 1):
            agar_trays[str(x)] = labware.load(
                AGAR_TRAY_TYPE, AGAR_TRAY_SLOTS[x - 1])
    else:
        agar_trays['1'] = labware.load(AGAR_TRAY_TYPE, AGAR_TRAY_SLOTS[0])

    # Register agar_trays for calibration
    if agar_tray_num == 1:
        p10_pipette.transfer(1, agar_trays['1'].wells(
            'A1'), agar_trays['1'].wells('H12'), trash=False)
    else:
        p10_pipette.transfer(1, agar_trays['1'].wells(
            'A1'), agar_trays['2'].wells('H12'), trash=False)
    p10_pipette.start_at_tip(p10_tipracks[0][0])

    # Run functions
    transformation_setup(final_well(tray_1_vols))
    phase_switch()
    outgrowth(final_well=final_well(tray_1_vols))
    spot_transformations(tray_1_vols=tray_1_vols)

for c in robot.commands():
    print(c)

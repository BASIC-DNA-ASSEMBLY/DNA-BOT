from opentrons import labware, instruments, modules, robot


spotting_tuples=((('A1', 'B1', 'C1', 'D1', 'E1'), ('A1', 'B1', 'C1', 'D1', 'E1'), (5, 5, 5, 5, 5)),)


def generate_transformation_wells(spotting_tuples):
    """Evaluates spotting_tuples and returns transformation wells.

    Args:
    spotting_tuples (tuple): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.

    """
    wells = []
    for spotting_tuple in spotting_tuples:
        for source_well in spotting_tuple[0]:
            wells.append(source_well)
    transformation_wells = [well for i, well in enumerate(
        wells) if wells.index(well) == i]
    return transformation_wells


def tiprack_slots(spotting_tuples):
    """Calculates p10 and p300 tiprack slots required.

    Args:
    spotting_tuples (tuple): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.

    """
    # Reactions' number
    transformation_reactions = len(
        generate_transformation_wells(spotting_tuples))
    spotting_reactions = sum(
        len(spotting_tuple[0]) for spotting_tuple in spotting_tuples)

    # p10 tiprack slots
    p10_tips = transformation_reactions + \
        spotting_reactions if agar_tray_num == 1 else int(
            spotting_reactions * 2)
    p10_tiprack_slots = p10_tips // 96 + 1 if p10_tips % 96 > 0 else 0

    # p300 tiprack slots
    p300_tips = transformation_reactions + spotting_reactions
    p300_tiprack_slots = p300_tips // 96 + 1 if p300_tips % 96 > 0 else 0
    return p10_tiprack_slots, p300_tiprack_slots


def transformation_setup(transformation_wells):
    """Sets up transformation reactions

    Args:
    transformation_wells (list). 

    """

    # Constants
    TEMP = 4  # Incubation temperature.
    ASSEMBLY_VOL = 5  # Volume of final assembly added to competent cells.
    MIX_SETTINGS = (4, 5)  # Mix after setting during final assembly transfers.
    INCUBATION_TIME = 20  # Cells and final assembly incubation time.

    # Set temperature deck to 4 °C and load competent cells
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    robot.pause()
    robot.comment('Load competent cells, uncap and resume run')

    # Transfer final assemblies
    p10_pipette.transfer(ASSEMBLY_VOL, assembly_plate.wells(transformation_wells), transformation_plate.wells(
        transformation_wells), new_tip='always', mix_after=(MIX_SETTINGS))

    # Incubate for 20 minutes and remove competent cells for heat shock
    p10_pipette.delay(minutes=INCUBATION_TIME)
    robot.pause()
    robot.comment(
        'Remove transformation reactions for heatshock and resume run')


def phase_switch(comment='Remove final assembly plate. Introduce agar tray(s) and deep well plate containing SOC media. Resume run.'):
    """Function pauses run enabling addition/removal of labware.

    Args:
    comment (str): string to be displayed during run following pause.

    """
    robot.pause()
    robot.comment(comment)


def outgrowth(
        cols,
        soc_well):
    """Outgrows transformed cells.

    Args:
    cols (list of str): list of cols in transformation plate containing samples.
    soc_well (str): Well containing SOC media in relevant plate.

    """

    # Constants
    SOC_VOL = 125
    SOC_MIX_SETTINGS = (4, 50)
    TEMP = 37
    OUTGROWTH_TIME = 60
    SOC_ASPIRATION_RATE = 25
    P300_DEFAULT_ASPIRATION_RATE = 150

    # Define wells
    transformation_cols = transformation_plate.cols(cols)
    soc = soc_plate.wells(soc_well)

    # Add SOC to transformed cells
    p300_pipette.set_flow_rate(aspirate=SOC_ASPIRATION_RATE)
    p300_pipette.transfer(SOC_VOL, soc, transformation_cols,
                          new_tip='always', mix_after=SOC_MIX_SETTINGS)
    p300_pipette.set_flow_rate(aspirate=P300_DEFAULT_ASPIRATION_RATE)

    # Incubate for 1 hour at 37 °C
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    p300_pipette.delay(minutes=OUTGROWTH_TIME)
    tempdeck.deactivate()


def spotting_cols(spotting_tuples):
    """Evaluates spotting_tuples and returns unique cols (str) 
    associated with each spotting_tuple's source wells.

    Args:
    spotting_tuples (tuple): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.

    """
    cols_list = []
    for spotting_tuple in spotting_tuples:
        source_wells_cols = [source_well[1:]
                             for source_well in spotting_tuple[0]]
        unique_cols = [col for i, col in enumerate(
            source_wells_cols) if source_wells_cols.index(col) == i]
        cols_list.append(unique_cols)
    return cols_list


def spot_transformations(
        spotting_tuples,
        agar_tray_num=1,
        tray_2_vol=8,
        dead_vol=2,
        spotting_dispense_rate=0.025,
        stabbing_depth=2):
    """Spots transformation reactions.

    Args:
    spotting_tuples (tuple): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.
    agar_tray_num (int): Number of agar trays transformations spotted on.
    tray_2_vol (float): Fixed volume to be spotted onto agar_trays['2'].
    dead_vol (float): Dead volume aspirated during spotting.
    spotting_dispense_rate (float): Rate p10_pipette dispenses at during spotting.
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

    # Define transformation cols
    transformation_cols = spotting_cols(spotting_tuples)

    # Spot transformation reactions
    for x in range(len(transformation_cols)):
        spotting_tuple = spotting_tuples[x]
        source_wells = spotting_tuple[0]
        target_wells = spotting_tuple[1]
        spot_vols = spotting_tuple[2]
        for col in transformation_cols[x]:
            p300_pipette.pick_up_tip()
            p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[
                1], transformation_plate.cols(col))
            p300_pipette.drop_tip()
        for y in range(len(source_wells)):
            spot(transformation_plate.wells(source_wells[y]),
                 agar_trays['1'].wells(target_wells[y]), spot_vols[y])
        if agar_tray_num == 2:
            for transformation_well, spot_info in spotting_tuples[x].items():
                spot(transformation_plate.wells(source_wells[y]),
                     agar_trays['2'].wells(target_wells[y]), spot_vols[y])


# Run protocol

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
p10_p300_tiprack_slots = tiprack_slots(spotting_tuples)
p10_slots = CANDIDATE_P10_SLOTS[
    :p10_p300_tiprack_slots[0]]
p300_slots = CANDIDATE_P300_SLOTS[
    :p10_p300_tiprack_slots[1]]

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
transformation_setup(generate_transformation_wells(spotting_tuples))
tempdeck.deactivate()
phase_switch()
spotting_tuples_cols = [col for cols in spotting_cols(
    spotting_tuples) for col in cols]
unique_cols = [col for i, col in enumerate(
    spotting_tuples_cols) if spotting_tuples_cols.index(col) == i]
outgrowth(unique_cols, soc_well='A1')
spot_transformations(spotting_tuples)

for c in robot.commands():
    print(c)

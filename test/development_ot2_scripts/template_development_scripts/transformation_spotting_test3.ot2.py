# Robot class added.
from opentrons import labware, instruments, modules, robot


spotting_dicts = (('A1', 'B1', 'C1', 'D1', 'A1', 'B1', 'C1', 'D1'), ('A1', 'B1', 'C1', 'D1', 'A2', 'B2', 'C2', 'D2'), (5, 5, 5, 5, 5, 5, 5, 5)), ((
    'A2', 'B2', 'C2', 'D2', 'A2', 'B2', 'C2', 'D2'), ('A2', 'B2', 'C2', 'D2', 'A3', 'B3', 'C3', 'D3'), (5, 5, 5, 5, 5, 5, 5, 5))


def generate_transformation_wells(spotting_dicts):
    """Evaluates spotting_dicts and returns transformation wells.

    Args:
    spotting_dicts (list): List of dictionaries containing volumes of each transformation reaction to be spotted. Wells in each dictionary are resuspended in unison.

    """
    wells = []
    for dictionary in spotting_dicts:
        for well, _ in dictionary.items():
            wells.append(well)
    transformation_wells = [well for i, well in enumerate(
        wells) if wells.index(well) == i]
    return transformation_wells


def tiprack_slots(spotting_dicts):
    """Calculates p10 and p300 tiprack slots required.

    Args:
    spotting_dicts (list): List of dictionaries containing volumes of each transformation reaction to be spotted. Wells in each dictionary are resuspended in unison.

    """
    # Reactions' number
    transformation_reactions = len(
        generate_transformation_wells(spotting_dicts))
    spotting_reactions = sum([len(element) for element in spotting_dicts])

    # p10 tiprack slots
    p10_tips = transformation_reactions + \
        spotting_reactions if agar_tray_num == 1 else int(
            spotting_reactions * 2)
    p10_tiprack_slots = p10_tips // 96 + 1 if p10_tips % 96 > 0 else 0

    # p300 tiprack slots
    p300_tips = transformation_reactions + spotting_reactions
    p300_tiprack_slots = p300_tips // 96 + 1 if p300_tips % 96 > 0 else 0
    return p10_tiprack_slots, p300_tiprack_slots


def transformation_setup(destination_wells):
    """Sets up transformation reactions

    Args:
    destination_wells (list): Wells containing final assemblies. 

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
    p10_pipette.transfer(ASSEMBLY_VOL, assembly_plate.wells(destination_wells), transformation_plate.wells(
        destination_wells), new_tip='always', mix_after=(MIX_SETTINGS))

    # Incubate for 20 minutes and remove competent cells for heat shock
    p10_pipette.delay(minutes=INCUBATION_TIME)
    robot.pause()
    robot.comment(
        'Remove transformation reactions for heatshock and resume run')
    tempdeck.deactivate()


def phase_switch(comment='Remove final assembly plate. Introduce agar tray(s) and deep well plate containing SOC media. Resume run.'):
    """Function pauses run enabling addition/removal of labware.

    Args:
    comment (str): string to be displayed during run following pause.

    """
    robot.pause()
    robot.comment(comment)


def outgrowth(
        cols,
        soc_well='A1'):
    """Outgrows transformed cells.

    Args:
    cols (list of str): list of cols in transformation plate containing samples.
    soc_well (str): Well containing SOC media in relevant plate.

    """

    # Constants
    SOC_VOL = 125  # Volume of SOC media (uL)
    SOC_MIX_SETTINGS = (4, 50)  # mix_after settings when adding SOC
    TEMP = 37  # outgrowth incubation temperature
    # time transformation reactions are outgrown for (mins)
    OUTGROWTH_TIME = 60

    # Define wells
    transformation_cols = transformation_plate.cols(cols)
    soc = soc_plate.wells(soc_well)

    # Add SOC to transformed cells
    p300_pipette.transfer(SOC_VOL, soc, transformation_cols,
                          new_tip='always', mix_after=SOC_MIX_SETTINGS)

    # Incubate for 1 hour at 37 °C
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    p300_pipette.delay(minutes=OUTGROWTH_TIME)
    tempdeck.deactivate()


def spotting_cols(spotting_dicts):
    """Evaluates spotting_dicts and returns unique cols (str) associated with each dictionary as a list of lists.

    Args:
    spotting_dicts (list): List of dictionaries containing volumes of each transformation reaction to be spotted. Wells in each dictionary are resuspended in unison.

    """
    cols_list = []
    for dictionary in spotting_dicts:
        dict_cols = [well[1:] for well, _ in dictionary.items()]
        unique_dict_cols = [col for i, col in enumerate(
            dict_cols) if dict_cols.index(col) == i]
        cols_list.append(unique_dict_cols)
    return cols_list


def spot_transformations(
        spotting_dicts,
        agar_tray_num=1,
        tray_2_vol=8,
        dead_vol=2,
        spotting_dispense_rate=0.025,
        time_delay=0,
        stabbing_depth=2):
    """Spots transformation reactions.

    Args:
    spotting_dicts (list): List of dictionaries containing volumes of each transformation reaction to be spotted. Wells in each dictionary are resuspended in unison.
    agar_tray_num (int): Number of agar trays transformations spotted on.
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

    # Define transformation cols
    transformation_cols = spotting_cols(spotting_dicts)

    # Spot transformation reactions
    for x in range(len(transformation_cols)):
        for col in transformation_cols[x]:
            p300_pipette.pick_up_tip()
            p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[
                1], transformation_plate.cols(col))
            p300_pipette.drop_tip()
        for transformation_well, spot_info in spotting_dicts[x].items():
            spot(transformation_plate.wells(transformation_well),
                 agar_trays['1'].wells(spot_info[0]), spot_info[1])
        if agar_tray_num == 2:
            for transformation_well, spot_info in spotting_dicts[x].items():
                spot(transformation_plate.wells(transformation_well),
                     agar_trays['2'].wells(spot_info[0]), spot_info[1])


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
p10_p300_tiprack_slots = tiprack_slots(spotting_dicts)
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
transformation_setup(generate_transformation_wells(spotting_dicts))
phase_switch()
spotting_dicts_cols = [col for cols in spotting_cols(
    spotting_dicts) for col in cols]
unique_cols = [col for i, col in enumerate(
    spotting_dicts_cols) if spotting_dicts_cols.index(col) == i]
outgrowth(unique_cols)
spot_transformations(spotting_dicts)

for c in robot.commands():
    print(c)
# python transformation_spotting_test3.ot2.py

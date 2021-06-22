from opentrons import protocol_api, simulate
import numpy as np
protocol = simulate.get_protocol_api('2.8')

# Rename to 'purification_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use
# replace TRANSFORMATION_PLATE_TYPE and AGAR_PLATE_TYPE with custom labware before use!
# email leetingan99@gmail.com if there are issues

metadata = {
     'apiLevel': '2.8',
     'protocolName': 'Transformation',
     'author': 'Baldwin_Lab',
     'description': 'Transformation reactions using an opentrons OT-2 for BASIC assembly.'}

# Example output produced by DNA-BOT for a single construct, uncomment and run to test the template
# spotting_tuples=[(('A1',), ('A1',), (5,))]
# soc_well='A1'

# Example output produced by DNA-BOT for 88 constructs, uncomment and run to test the template
# spotting_tuples=[(('A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1'), ('A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2'), ('A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3'), ('A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4'), ('A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5'), ('A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6'), ('A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7'), ('A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8'), ('A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9'), ('A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10'), ('A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11'), ('A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11'), (5, 5, 5, 5, 5, 5, 5, 5))]
# soc_well='A1'

def run(protocol):
# added run function for API version 2

    def generate_transformation_wells(spotting_tuples):
        """
        Evaluates spotting_tuples and returns transformation wells.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given in the form: ((source wells), (target wells), (spotting volumes)).

        """

        wells = []
        for spotting_tuple in spotting_tuples:
            for source_well in spotting_tuple[0]:
                wells.append(source_well)
        transformation_wells = [well for i, well in enumerate(
            wells) if wells.index(well) == i]
        return transformation_wells


    def tiprack_slots(spotting_tuples, max_spot_vol=5):
        """
        Calculates p10 and p300 tiprack slots required.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given in the form: ((source wells), (target wells), (spotting volumes)).
        max_spot_vol (float): Maximum volume that is spotted per spot reaction.

        """

        # Reactions' number
        transformation_reactions = len(generate_transformation_wells(spotting_tuples))
        spotting_reactions = 0
        for spotting_tuple in spotting_tuples:
            spots = np.array(spotting_tuple[2])/max_spot_vol
            np.ceil(spots)
            spotting_reactions = spotting_reactions + int(np.sum(spots))

        # errrr should be fine lol # REMOVE

        # p10 tiprack slots
        p10_tips = transformation_reactions + spotting_reactions
        p10_tiprack_slots = p10_tips // 96 + 1 if p10_tips % 96 > 0 else p10_tips / 96

        # p300 tiprack slots
        p300_tips = transformation_reactions + spotting_reactions
        p300_tiprack_slots = p300_tips // 96 + \
            1 if p300_tips % 96 > 0 else p300_tips / 96
        return int(p10_tiprack_slots), int(p300_tiprack_slots)


    def transformation_setup(transformation_wells):
        """
        Sets up transformation reactions

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
        # removed: tempdeck.wait_for_temp()
            # API version2 automatically pauses execution until the set temperature is reached
            # thus it no longer uses .wait_for_temp()
        protocol.pause('Load competent cells, uncap and resume run')
        # old code:
            # robot.pause()
            # robot.comment('Load competent cells, uncap and resume run')
        # API version 2 uses 'protocol.' instead of 'robot.' and combines '.pause' and '.comment'

        # Transfer final assemblies
        p10_pipette.transfer(ASSEMBLY_VOL,
                             [assembly_plate.wells_by_name()[well_name] for well_name in transformation_wells],
                             [transformation_plate.wells_by_name()[well_name] for well_name in transformation_wells],
                             new_tip='always',
                             mix_after=(MIX_SETTINGS))
        # old code:
            # p10_pipette.transfer(ASSEMBLY_VOL,
                                 # assembly_plate.wells(transformation_wells),
                                 # transformation_plate.wells(transformation_wells),
                                 # new_tip='always',
                                 # mix_after=(MIX_SETTINGS))
         # .wells() doesn't take lists as arguements, newer wells_by_name() returns a dictionary


        # Incubate for 20 minutes and remove competent cells for heat shock
        protocol.delay(minutes=INCUBATION_TIME)
        # old code:
            # p10_pipette.delay(minutes=INCUBATION_TIME)
        # API version 2 no longer has .delay() for pipettes, it uses protocol.delay() to pause the entire protocol

        protocol.pause('Remove transformation reactions, conduct heatshock and replace.')
        # old code:
            # robot.pause()
            # robot.comment('Remove transformation reactions, conduct heatshock and replace.')
        # API version 2 uses 'protocol.' instead of 'robot.' and combines '.pause' and '.comment'


    def phase_switch():
        """
        Function pauses run enabling addition/removal of labware.

        """
        protocol.pause('Remove final assembly plate. Introduce agar tray and deep well plate containing SOC media. Resume run.')
        # old code:
            # def phase_switch(comment='Remove final assembly plate. Introduce agar tray and deep well plate containing SOC media. Resume run.'):
                # robot.pause()
                # robot.comment(comment)
        # API version 2 uses 'protocol.' instead of 'robot.' and combines '.pause' and '.comment'


    def outgrowth(
            cols,
            soc_well):
        """
        Outgrows transformed cells.

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
        transformation_cols = [transformation_plate.columns_by_name()[column] for column in cols]
        # old code:
            # transformation_cols = transformation_plate.cols(cols)
        # API version 2 labware use .columns attribute, not .cols
        # but .columns() doesn't take lists as arguements, newer columns_by_name() returns a dictionary

        soc = soc_plate.wells(soc_well)

        # Add SOC to transformed cells
        p300_pipette.flow_rate.aspirate = SOC_ASPIRATION_RATE
        # old code:
            # p300_pipette.set_flow_rate(aspirate=SOC_ASPIRATION_RATE)
            # flow rates are set directly in API version 2, brackets not required
        p300_pipette.transfer(SOC_VOL, soc, transformation_cols,
                              new_tip='always', mix_after=SOC_MIX_SETTINGS)
        p300_pipette.flow_rate.aspirate = P300_DEFAULT_ASPIRATION_RATE
        # old code:
            # p300_pipette.set_flow_rate(aspirate=P300_DEFAULT_ASPIRATION_RATE)
            # API version 2 labware use .columns attribute, not .cols

        # Incubate for 1 hour at 37 °C
        tempdeck.set_temperature(TEMP)
        # removed: tempdeck.wait_for_temp()
            # API version2 automatically pauses execution until the set temperature is reached
            # thus it no longer uses .wait_for_temp()

        protocol.delay(minutes=OUTGROWTH_TIME)
        # old code:
            # p300_pipette.delay(minutes=OUTGROWTH_TIME)
        # API version 2 no longer has .delay() for pipettes, it uses protocol.delay() to pause the entire protocol

        tempdeck.deactivate()


    def spotting_cols(spotting_tuples):
        """
        Evaluates spotting_tuples and returns unique cols (str) associated with each spotting_tuple's source wells.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given in the form: ((source wells), (target wells), (spotting volumes)).

        """
        cols_list = []
        for spotting_tuple in spotting_tuples:
            source_wells_cols = [source_well[1:] for source_well in spotting_tuple[0]]
            unique_cols = [col for i, col in enumerate(source_wells_cols) if source_wells_cols.index(col) == i]
            cols_list.append(unique_cols)
        return cols_list


    def spot_transformations(
            spotting_tuples,
            dead_vol=2,
            spotting_dispense_rate=0.025,
            stabbing_depth=2,
            max_spot_vol=5):
        """
        Spots transformation reactions.

        Args:
        spotting_tuples (list): Sets of spotting reactions are given in the form: ((source wells), (target wells), (spotting volumes)).
        dead_vol (float): Dead volume aspirated during spotting.
        spotting_dispense_rate (float): Rate p10_pipette dispenses at during spotting.
        stabbing_depth (float): Depth p10_pipette moves into agar during spotting.
        max_spot_vol (float): Maximum volume that is spotted per spot reaction.

        """

        def spot(
                source,
                target,
                spot_vol):
            """
            Spots an individual reaction using the p10 pipette.

            Args:
            source (str): Well containing the transformation reaction to be spotted.
            target (str): Well transformation reaction is to be spotted to.
            spot_vol (float): Volume of transformation reaction to be spotted (uL).

            """

            # Constants
            DEFAULT_HEAD_SPEED = {'x': 400, 'y': 400,'z': 125, 'a': 125}
            SPOT_HEAD_SPEED = {'x': 400, 'y': 400, 'z': 125,'a': 125 // 4}
            DISPENSING_HEIGHT = 5
            SAFE_HEIGHT = 15  # height avoids collision with agar tray.

            # Spot
            p10_pipette.pick_up_tip()
            p10_pipette.aspirate(spot_vol + dead_vol, source[0])
            # old code:
                # p10_pipette.aspirate(spot_vol + dead_vol, source)
                # returned type error because 'source' was a list containing one item (the well location)
                # source[0] takes the location out of the list

            p10_pipette.move_to(target[0].top(SAFE_HEIGHT))
            p10_pipette.move_to(target[0].top(DISPENSING_HEIGHT))
            # old code:
                # p10_pipette.move_to(target.top(SAFE_HEIGHT))
                # p10_pipette.move_to(target.top(DISPENSING_HEIGHT))
                # returned attribute error because 'target' was a list containing one item (the well location)
                # target[0] takes the location out of the list

            p10_pipette.dispense(volume=spot_vol, rate=spotting_dispense_rate)

            protocol.max_speeds.update(SPOT_HEAD_SPEED)
            # old code:
                # robot.head_speed(combined_speed=max(SPOT_HEAD_SPEED.values()), **SPOT_HEAD_SPEED)
                # robot.head_speed not used in API version 2
                # replaced with protocol.max_speeds
            # new code no longer uses the lower value between combined speed or specified speed
                # just uses each axis' specified speed directly
            p10_pipette.move_to(target[0].top(-1 * stabbing_depth))
            # old code:
                # p10_pipette.move_to(target.top(-1*stabbing_depth))
                # returns attribute error because 'target' was a list containing one item (the well location)
            protocol.max_speeds.update(DEFAULT_HEAD_SPEED)
            # old code:
                # robot.head_speed(combined_speed=max(DEFAULT_HEAD_SPEED.values()), **DEFAULT_HEAD_SPEED)
                # robot.head_speed not used in API version 2
                # replaced with protocol.max_speeds
            # new code no longer uses the lower value between combined speed or specified speed
                # just uses each axis' specified speed directly
            p10_pipette.move_to(target[0].top(SAFE_HEIGHT))
            # old code:
                # p10_pipette.move_to(target[0].top(SAFE_HEIGHT))
                # returns attribute error because 'target' was a list containing one item (the well location)

            # Dispose of dead volume and tip
            p10_pipette.dispense(dead_vol, spotting_waste[0])
            # old code:
                # p10_pipette.dispense(dead_vol, spotting_waste)
                # returns type error because 'target' was a list containing one item (the well location)

            p10_pipette.blow_out()
                # the simple .blow_out command blows out at current position (spotting waste) by defualt
                # unlike blowout=true in complex commands, which by default will blow out in waste

            p10_pipette.drop_tip()

        def spot_tuple(spotting_tuple):
            """
            Spots all reactions defined by the spotting tuple. Requires the function spot.

            Args:
            spotting_tuple (tuple): Spotting reactions given in the form: (source wells), (target wells), (spotting volumes).
            Each unique source well is resuspended once prior to spotting.

            """
            source_wells = spotting_tuple[0]
            target_wells = spotting_tuple[1]
            spot_vols = list(spotting_tuple[2])
            while max(spot_vols) > 0:
                for index, spot_vol in enumerate(spot_vols):
                    if spot_vol == 0:
                        pass
                    else:
                        vol = spot_vol if spot_vol <= max_spot_vol else max_spot_vol
                        spot(source = transformation_plate.wells(source_wells[index]), target = agar_plate.wells(target_wells[index]), spot_vol = vol)
                        spot_vols[index] = spot_vols[index] - vol

        # Constants
        TRANSFORMATION_MIX_SETTINGS = [4, 50]

        # Spot transformation reactions
            # Each unique transformation well is resuspended once prior to spotting.

        for spotting_tuple in spotting_tuples:
            source_wells_cols = [source_well[1:] for source_well in spotting_tuple[0]]
            unique_cols = [col for i, col in enumerate(source_wells_cols) if source_wells_cols.index(col) == i]
            for col in unique_cols:
                p300_pipette.pick_up_tip()
                p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[1],transformation_plate.columns_by_name()[col][0])
                # old code:
                    # p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0], TRANSFORMATION_MIX_SETTINGS[1],transformation_plate.cols(col))
                    # .columns aka .cols doesn't take lists anymore
                        # replaced with .columns_by_name
                    # .mix only takes one location, not several locations
                        # added [0] to specify only the wells in row A
                        # is identical for the protocol, as this is using a multi-channel pipette
                p300_pipette.drop_tip()
            spot_tuple(spotting_tuple)


    ### Run protocol

    # Constants
    CANDIDATE_P10_SLOTS = ['9', '2', '5']
    CANDIDATE_P300_SLOTS = ['3', '6']
    P10_TIPRACK_TYPE = 'opentrons_96_tiprack_10ul'
        # changed from 'tiprack-10ul'
    P300_TIPRACK_TYPE = 'opentrons_96_tiprack_300ul'
    P10_MOUNT = 'right'
    P300_MOUNT = 'left'
    ASSEMBLY_PLATE_TYPE = 'nest_96_wellplate_200ul_flat'
        # changed from '4ti-0960_FrameStar'
        # was previously defined in add.labware.py, API version 2 doesn't support labware.create anymore
    ASSEMBLY_PLATE_SLOT = '8'
    TEMPDECK_SLOT = '10'
    TRANSFORMATION_PLATE_TYPE = 'nest_96_wellplate_200ul_flat'
        # changed from 'Eppendorf_30133366_plate_96'
            # was previously defined in add.labware.py, API version 2 doesn't support labware.create anymore

        # !!! CURRENT PLATE IS A PLACEHOLDER FOR RUNNING SIMMULATION WITHOUT CUSTOM LABWARE!!!
        # name made with Opentron Labware Creator = (not yet made)

        # custom labware not defined: what plate to change to? why is this plate used?
            # according to add_labware.py, it is a 250ul 96 well plate with spacing identical to the 4ti0960-FrameStar 96 wp

    SOC_PLATE_TYPE = 'nest_96_wellplate_2ml_deep'
        # changed from '4ti0136_96_deep-well'
    SOC_PLATE_SLOT = '7'
    TUBE_RACK_TYPE = 'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap'
        # changed from 'tube-rack_E1415-1500'
    TUBE_RACK_SLOT = '11'
    SPOTTING_WASTE_WELL = 'A1'
    AGAR_PLATE_TYPE = 'nest_96_wellplate_200ul_flat'
        # changed from 'Nunc_Omnitray'
            # it is a 1 well plate filled with agar;
            # but for the Opentron to spot in the locations of a 96 wp, it is defined similar to a 96 wp
            # was previously defined in add.labware.py, API version 2 doesn't support labware.create anymore

        # !!! CURRENT PLATE IS A PLACEHOLDER FOR RUNNING SIMMULATION WITHOUT CUSTOM LABWARE!!!
        # !!! name made with Opentron Labware Creator = 'nuncomnitray_96_wellplate_0.01ul' !!!!


        # custom labware made using Opentron's Labware Creator:
            # external dimensions:
                # footprint length = 127.76 mm
                # footrpint width = 85.48 mm
                # footprint height = 15.70 mm
                # taken from Thermofisher's documentation for Nunc Omnitray
                # https://www.thermofisher.com/document-connect/document-connect.html?url=https%3A%2F%2Fassets.thermofisher.com%2FTFS-Assets%2FLSG%2Fmanuals%2FD03023.pdf&title=VGVjaG5pY2FsIERhdGEgU2hlZXQ6IE51bmMgT21uaXRyYXk=
            # well measurements
                # depth = 0.01 mm
                # diameter =  0.01 mm
                # in old add.labware.py, they were defined as 0, but Labware Creator requires a value >0
            # spacing
                # x-offset = 14.38 mm
                # y-offset = 11.24 mm
                # x-spacing = 9.00 mm
                # y-spacing) = 9.00 mm
                # taken from Nest 96 well plates
                # https://labware.opentrons.com/nest_96_wellplate_100ul_pcr_full_skirt/
        # before using protocol, need to upload the 'nuncomnitray_96_wellplate_0.01ul.json' custom labware file into Opentrons app

    AGAR_PLATE_SLOT = '1'

    # Tiprack slots
    p10_p300_tiprack_slots = tiprack_slots(spotting_tuples)
    p10_slots = CANDIDATE_P10_SLOTS[:p10_p300_tiprack_slots[0]]
    p300_slots = CANDIDATE_P300_SLOTS[:p10_p300_tiprack_slots[1]]

    # Define labware
    p10_tipracks = [protocol.load_labware(P10_TIPRACK_TYPE, slot) for slot in p10_slots]
        # changed to protocol.load_labware for API version 2
    p300_tipracks = [protocol.load_labware(P300_TIPRACK_TYPE, slot) for slot in p300_slots]
        # changed to protocol.load_labware for API version 2
    p10_pipette = protocol.load_instrument('p10_single', mount=P10_MOUNT, tip_racks=p10_tipracks)
        # changed to protocol.load_instrument for API version 2
    p300_pipette = protocol.load_instrument('p300_multi', mount=P300_MOUNT, tip_racks=p300_tipracks)
        # changed to protocol.load_instrument for API version 2

    assembly_plate = protocol.load_labware(ASSEMBLY_PLATE_TYPE, ASSEMBLY_PLATE_SLOT)
        # changed to protocol.load_labware for API version 2
    tempdeck = protocol.load_module('tempdeck', TEMPDECK_SLOT)
        # changed to protocol.load_module for API version 2
    transformation_plate = tempdeck.load_labware(TRANSFORMATION_PLATE_TYPE, TEMPDECK_SLOT)
        # changed to protocol.load_labware for API version 2
        # removed share=True, not required in API version 2
        # removed TEMPDECK_SLOT as it is loaded directly onto temperature module
    soc_plate = protocol.load_labware(SOC_PLATE_TYPE, SOC_PLATE_SLOT)
        # changed to protocol.load_labware for API version 2
    tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_SLOT)
        # changed to protocol.load_labware for API version 2
    spotting_waste = tube_rack.wells(SPOTTING_WASTE_WELL)
    agar_plate = protocol.load_labware(AGAR_PLATE_TYPE, AGAR_PLATE_SLOT)
        # changed to protocol.load_labware for API version 2


    # Register agar_plate for calibration
    p10_pipette.transfer(1, agar_plate.wells('A1'), agar_plate.wells('H12'), trash=False)
    # removed:
        # p10_pipette.start_at_tip(p10_tipracks[0][0])
        # pipette automatically starts from 'A1' tiprack location
        # if re-adding, need to use p10.pipette.starting_tip() instead of p10.pipette.start_at_tip()


    # Run functions
    transformation_setup(generate_transformation_wells(spotting_tuples))
    phase_switch()
    spotting_tuples_cols = [col for cols in spotting_cols(spotting_tuples) for col in cols]
    unique_cols = [col for i, col in enumerate(spotting_tuples_cols) if spotting_tuples_cols.index(col) == i]
    outgrowth(cols=unique_cols, soc_well=soc_well)
    spot_transformations(spotting_tuples)
    print(unique_cols)

run(protocol)

for line in protocol.commands():
    print(line)

# replace TRANSFORMATION_PLATE_TYPE and AGAR_PLATE_TYPE with custom labware before use!
# email leetingan99@gmail.com if there are issues

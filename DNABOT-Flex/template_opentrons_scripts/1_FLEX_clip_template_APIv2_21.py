# -*- coding: utf-8 -*-
requirements = {"robotType": "Flex"}
from opentrons import protocol_api
import numpy as np
from opentrons.types import Point

metadata = {
    "protocolName": "DNABOT Step 1: Clip Reaction with thermocycler",
    "description": "Implements linker ligation reactions using an opentrons Flex, including the thermocycler module gen1 or gen2.",
    "apiLevel": "2.21",
}


def run(protocol: protocol_api.ProtocolContext):
    trash = protocol.load_trash_bin("A3")
    tiprack_type = __LABWARES["tiprack_20ul"]["id"]
    INITIAL_TIP = "A1"
    CANDIDATE_TIPRACK_SLOTS = ["D3", "C3", "B3"]
    PIPETTE_TYPE = __HARDWARE["single_pipette"]["id"]
    PIPETTE_MOUNT = __HARDWARE["single_pipette_mount"]["id"]
    tc_mod = protocol.load_module(
        module_name=__HARDWARE["thermocycler"]["id"], location="B1"
    )
    tc_mod.open_lid()
    tc_mod.deactivate_lid()
    tc_mod.set_block_temperature(temperature=__PARAMETERS["thermo_temp"]["value"])
    DESTINATION_PLATE_TYPE = __LABWARES["clip_plate"]["id"]
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
    SOURCE_PLATE_TYPE = __LABWARES["clip_source_plate"]["id"]
    TUBE_RACK_TYPE = __LABWARES["24_tuberack_1500ul"]["id"]
    TUBE_RACK_POSITION = "C1"
    MASTER_MIX_WELL = "A1"
    WATER_WELL = "A2"
    MASTER_MIX_VOLUME = 20
    if __PARAMETERS["premix_linkers"]["value"] == "Yes":
        Mix_linkers_bool = True
    else:
        Mix_linkers_bool = False
    if __PARAMETERS["premix_parts"]["value"] == "Yes":
        Mix_parts_bool = True
    else:
        Mix_parts_bool = False

    def configure_low_volume_mode(pipette, volume):
        """Configure Flex 50 uL pipettes for the next transfer volume when supported."""
        pipette.configure_for_volume(max(1, min(float(volume), 50)))

    def get_low_volume_push_out(volume):
        """Set the dispense push-out for low-volume transfers."""
        return 7 if float(volume) < 5 else 2

    def get_mix_repetitions(volume):
        """Scale mix repetitions with premix volume."""
        if volume <= 20:
            return 4
        if volume <= 40:
            return 5
        if volume <= 60:
            return 6
        return 8

    def get_premix_aspirate_height(volume):
        """Aspirate higher as volume rises to draw liquid toward the tube bottom."""
        if volume <= 20:
            return 2
        if volume <= 40:
            return 3
        if volume <= 60:
            return 4
        return 5

    def get_premix_stroke_volume(volume):
        """Derive a premix stroke from the working volume, capped for 20 uL tips."""
        return min(volume * 0.8, 20)

    def perform_premix(pipette, well, volume):
        """Premix with slow aspirates, fast dispenses, then a slower finishing pass."""
        high = 2
        slow = 0.4
        aspirate_height = get_premix_aspirate_height(volume)
        dispense_height = 2
        stroke_volume = get_premix_stroke_volume(volume)
        mix_reps = get_mix_repetitions(volume)
        for mix_step in range(mix_reps):
            x_offset = 2 if mix_step % 2 == 0 else -2
            dispense_location = well.bottom(dispense_height).move(
                Point(x=x_offset, y=0, z=0)
            )
            pipette.aspirate(stroke_volume, well.bottom(aspirate_height), rate=slow)
            pipette.dispense(stroke_volume, dispense_location, rate=high)
        pipette.aspirate(stroke_volume, well.bottom(aspirate_height), rate=slow)
        pipette.dispense(
            stroke_volume,
            well.bottom(dispense_height).move(Point(x=2, y=0, z=0)),
            rate=slow,
            push_out=3,
        )
        pipette.move_to(well.top(-5))
        protocol.delay(seconds=1)
        pipette.blow_out()
        pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
        pipette.drop_tip()

    def perform_clip_mix(pipette, well, mix_volume=20, repetitions=4):
        configure_low_volume_mode(pipette, mix_volume)
        """Mix a clip reaction well with fixed repeated cycles and a slow finish."""
        high = 2
        normal = 1
        slow = 0.4
        aspirate_volume = mix_volume
        for mix_step in range(repetitions):
            x_offset = 2 if mix_step % 2 == 0 else -2
            dispense_location = well.bottom(3).move(Point(x=x_offset, y=0, z=0))
            pipette.aspirate(aspirate_volume, well.bottom(1), rate=normal)
            pipette.dispense(aspirate_volume, dispense_location, rate=high)
        pipette.aspirate(aspirate_volume, well.bottom(2), rate=slow)
        pipette.dispense(
            aspirate_volume,
            dispense_location,
            rate=high,
            push_out=max(3, aspirate_volume / 5),
        )
        pipette.move_to(well.top(-5))
        protocol.delay(seconds=1)
        pipette.blow_out()
        pipette.touch_tip(radius=0.9, v_offset=-5, speed=10)
        pipette.drop_tip()

    def pre_wet(pipette, well, repetitions=1):
        slow = 0.4
        """Pre-wet the tip in the source well before low-volume aspiration."""
        pre_wet_volume = 2
        for _ in range(repetitions):
            pipette.aspirate(pre_wet_volume, well.bottom(1), rate=slow)
            protocol.delay(seconds=0.5)
            pipette.dispense(
                pre_wet_volume, well.bottom(1).move(Point(x=0, y=0, z=0)), rate=slow
            )

    def pick_up_tip_with_reload(pipette, tiprack_slots, total_refills_needed=0):
        """Pause for manual refill when the loaded tipracks run out of tips."""
        try:
            pipette.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            refill_note = ""
            if total_refills_needed > 0:
                refill_note = f" This run is expected to need about {total_refills_needed} manual refill(s)."
            protocol.pause(
                f"Refill clip tipracks in slots {tiprack_slots} and resume the protocol.{refill_note}"
            )
            pipette.reset_tipracks()
            pipette.pick_up_tip()

    def mix_linkers_function(
        Mix_linkers_bool,
        clips_dict,
        pipette,
        source_plates,
        PIPETTE_TYPE,
        slots,
        manual_refills_needed,
    ):
        if PIPETTE_TYPE == "flex_1channel_50":
            pipette.flow_rate.aspirate = 20
            pipette.flow_rate.dispense = 20
            pipette.flow_rate.blow_out = 35
        requested_linkers_volume = __PARAMETERS["linkers_volume"]["value"]
        if requested_linkers_volume < 20:
            raise ValueError("Linker mixing volume must be at least 20 uL.")
        linkers_volume = requested_linkers_volume
        if Mix_linkers_bool:
            prefixes = []
            loop_prefixes_wells = clips_dict["prefixes_wells"]
            loop_prefixes_plates = clips_dict["prefixes_plates"]
            len_prefixes = len(clips_dict["prefixes_wells"])
            for i in range(len_prefixes):
                prefixes.append([loop_prefixes_plates[i], loop_prefixes_wells[i]])
            prefixes_unique = np.unique(np.array(prefixes), axis=0)
            suffixes = []
            loop_suffixes_wells = clips_dict["suffixes_wells"]
            loop_suffixes_plates = clips_dict["suffixes_plates"]
            len_suffixes = len(clips_dict["suffixes_wells"])
            for i in range(len_suffixes):
                suffixes.append([loop_suffixes_plates[i], loop_suffixes_wells[i]])
            suffixes_unique = np.unique(np.array(suffixes), axis=0)
            for clip_num in range(len(prefixes_unique)):
                configure_low_volume_mode(
                    pipette, get_premix_stroke_volume(linkers_volume)
                )
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                perform_premix(
                    pipette,
                    source_plates[prefixes_unique[clip_num, 0]][
                        prefixes_unique[clip_num, 1]
                    ],
                    linkers_volume,
                )
            for clip_num in range(len(suffixes_unique)):
                configure_low_volume_mode(
                    pipette, get_premix_stroke_volume(linkers_volume)
                )
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                perform_premix(
                    pipette,
                    source_plates[suffixes_unique[clip_num, 0]][
                        suffixes_unique[clip_num, 1]
                    ],
                    linkers_volume,
                )
        else:
            pass

    def mix_parts_function(
        Mix_parts_bool,
        clips_dict,
        pipette_name,
        source_plates,
        PIPETTE_TYPE,
        slots,
        manual_refills_needed,
    ):
        pipette = pipette_name
        if PIPETTE_TYPE == "flex_1channel_50":
            pipette.flow_rate.aspirate = 20
            pipette.flow_rate.dispense = 20
            pipette.flow_rate.blow_out = 35
        requested_parts_volume = __PARAMETERS["parts_volume"]["value"]
        if requested_parts_volume < 10:
            raise ValueError("Part mixing volume must be at least 10 uL.")
        parts_volume = requested_parts_volume
        if Mix_parts_bool:
            parts = []
            loop_parts_wells = clips_dict["parts_wells"]
            loop_parts_plates = clips_dict["parts_plates"]
            len_parts = len(clips_dict["parts_wells"])
            for i in range(len_parts):
                parts.append([loop_parts_plates[i], loop_parts_wells[i]])
            parts_unique = np.unique(np.array(parts), axis=0)
            for clip_num in range(len(parts_unique)):
                configure_low_volume_mode(
                    pipette, get_premix_stroke_volume(parts_volume)
                )
                pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
                pipette.well_bottom_clearance.aspirate = 2
                pipette.well_bottom_clearance.dispense = 1
                perform_premix(
                    pipette,
                    source_plates[parts_unique[clip_num, 0]][parts_unique[clip_num, 1]],
                    parts_volume,
                )
            else:
                pass

    def clip(
        prefixes_wells,
        prefixes_plates,
        suffixes_wells,
        suffixes_plates,
        parts_wells,
        parts_plates,
        parts_vols,
        water_vols,
    ):
        prefixes = []
        loop_prefixes_wells = clips_dict["prefixes_wells"]
        loop_prefixes_plates = clips_dict["prefixes_plates"]
        len_prefixes = len(clips_dict["prefixes_wells"])
        for i in range(len_prefixes):
            prefixes.append([loop_prefixes_plates[i], loop_prefixes_wells[i]])
        prefixes_unique = np.unique(np.array(prefixes), axis=0)
        suffixes = []
        loop_suffixes_wells = clips_dict["suffixes_wells"]
        loop_suffixes_plates = clips_dict["suffixes_plates"]
        len_suffixes = len(clips_dict["suffixes_wells"])
        for i in range(len_suffixes):
            suffixes.append([loop_suffixes_plates[i], loop_suffixes_wells[i]])
        suffixes_unique = np.unique(np.array(suffixes), axis=0)
        parts = []
        loop_parts_wells = clips_dict["parts_wells"]
        loop_parts_plates = clips_dict["parts_plates"]
        len_parts = len(clips_dict["parts_wells"])
        for i in range(len_parts):
            parts.append([loop_parts_plates[i], loop_parts_wells[i]])
        parts_unique = np.unique(np.array(parts), axis=0)
        clip_count = len(parts_wells)
        purification_script_compatible = clip_count <= 48
        if not purification_script_compatible:
            protocol.comment(
                "Warning: this clip run contains more than 48 clips, so it is not directly compatible with Step 2 (2_FLEX_purification_template_APIv2_21.py), which is capped at 48."
            )
        if Mix_linkers_bool:
            if Mix_parts_bool:
                total_tips = (
                    4 * clip_count
                    + len(prefixes_unique)
                    + len(suffixes_unique)
                    + len(parts_unique)
                )
            else:
                total_tips = (
                    4 * clip_count + len(prefixes_unique) + len(suffixes_unique)
                )
        elif Mix_parts_bool:
            total_tips = 4 * clip_count + len(parts_unique)
        else:
            total_tips = 4 * clip_count
        letter_dict = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7}
        tiprack_1_tips = (13 - int(INITIAL_TIP[1:])) * 8 - letter_dict[INITIAL_TIP[0]]
        if total_tips > tiprack_1_tips:
            tiprack_num = (
                1
                + (total_tips - tiprack_1_tips) // 96
                + (1 if (total_tips - tiprack_1_tips) % 96 > 0 else 0)
            )
        else:
            tiprack_num = 1
        loaded_tiprack_count = min(tiprack_num, len(CANDIDATE_TIPRACK_SLOTS))
        slots = CANDIDATE_TIPRACK_SLOTS[:loaded_tiprack_count]
        loaded_tip_capacity = tiprack_1_tips + max(0, loaded_tiprack_count - 1) * 96
        manual_refills_needed = 0
        if total_tips > loaded_tip_capacity:
            additional_tips = total_tips - loaded_tip_capacity
            manual_refills_needed = 1 + (additional_tips - 1) // (
                96 * loaded_tiprack_count
            )
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        pipette = protocol.load_instrument(
            PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks
        )
        destination_wells = destination_plate.wells()[0:clip_count]
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        master_mix = tube_rack[MASTER_MIX_WELL]
        water = tube_rack[WATER_WELL]
        source_plates = {}
        source_plates_keys = list(set(prefixes_plates + suffixes_plates + parts_plates))
        for key in source_plates_keys:
            source_plates[key] = protocol.load_labware(SOURCE_PLATE_TYPE, key)
        mix_linkers_function(
            Mix_linkers_bool,
            clips_dict,
            pipette,
            source_plates,
            PIPETTE_TYPE,
            slots,
            manual_refills_needed,
        )
        mix_parts_function(
            Mix_parts_bool,
            clips_dict,
            pipette,
            source_plates,
            PIPETTE_TYPE,
            slots,
            manual_refills_needed,
        )
        pipette.well_bottom_clearance.aspirate = 1
        pipette.well_bottom_clearance.dispense = 0
        configure_low_volume_mode(pipette, MASTER_MIX_VOLUME)
        pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
        pipette.distribute(
            MASTER_MIX_VOLUME,
            master_mix,
            destination_wells,
            blow_out=True,
            blowout_location="source well",
            new_tip="never",
        )
        pipette.drop_tip()
        pipette.well_bottom_clearance.aspirate = 1
        pipette.well_bottom_clearance.dispense = 3
        configure_low_volume_mode(pipette, 1)
        pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
        for (destination_well, water_volume) in zip(destination_wells, water_vols):
            if water_volume <= 0:
                continue
            configure_low_volume_mode(pipette, water_volume)
            pipette.aspirate(water_volume, water.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(
                water_volume,
                destination_well.bottom(3),
                push_out=get_low_volume_push_out(water_volume),
            )
        pipette.drop_tip()
        for clip_num in range(clip_count):
            pipette.well_bottom_clearance.aspirate = 2
            pipette.well_bottom_clearance.dispense = 2
            configure_low_volume_mode(pipette, 1)
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            prefix_source = source_plates[prefixes_plates[clip_num]][
                prefixes_wells[clip_num]
            ]
            pre_wet(pipette, prefix_source)
            pipette.aspirate(1, prefix_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(
                1,
                destination_wells[clip_num].bottom(2),
                push_out=get_low_volume_push_out(1),
            )
            protocol.delay(seconds=0.5)
            perform_clip_mix(
                pipette, destination_wells[clip_num], mix_volume=20, repetitions=2
            )
            configure_low_volume_mode(pipette, 1)
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            suffix_source = source_plates[suffixes_plates[clip_num]][
                suffixes_wells[clip_num]
            ]
            pre_wet(pipette, suffix_source)
            pipette.aspirate(1, suffix_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(
                1,
                destination_wells[clip_num].bottom(3),
                push_out=get_low_volume_push_out(1),
            )
            protocol.delay(seconds=0.5)
            perform_clip_mix(
                pipette, destination_wells[clip_num], mix_volume=20, repetitions=2
            )
            configure_low_volume_mode(pipette, parts_vols[clip_num])
            pick_up_tip_with_reload(pipette, slots, manual_refills_needed)
            part_source = source_plates[parts_plates[clip_num]][parts_wells[clip_num]]
            pre_wet(pipette, part_source)
            pipette.aspirate(parts_vols[clip_num], part_source.bottom(1))
            protocol.delay(seconds=0.5)
            pipette.dispense(
                parts_vols[clip_num],
                destination_wells[clip_num].bottom(3),
                push_out=get_low_volume_push_out(parts_vols[clip_num]),
            )
            protocol.delay(seconds=0.5)
            perform_clip_mix(
                pipette, destination_wells[clip_num], mix_volume=20, repetitions=4
            )

    clip(**clips_dict)
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)
    profile = [
        {"temperature": 37, "hold_time_minutes": 2},
        {"temperature": 20, "hold_time_minutes": 1},
    ]
    tc_mod.execute_profile(steps=profile, repetitions=30, block_max_volume=30)
    tc_mod.set_block_temperature(37, hold_time_minutes=10, block_max_volume=30)
    tc_mod.set_block_temperature(65, hold_time_minutes=20, block_max_volume=30)
    tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=30)
    if __PARAMETERS["clip_keep_thermo_lid_closed"]["value"] == "Yes":
        Thermo_lid_bool = True
    else:
        Thermo_lid_bool = False
    if Thermo_lid_bool:
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(temperature=4)
    else:
        tc_mod.set_lid_temperature(45)
        tc_mod.deactivate_lid()
        tc_mod.open_lid()
    for line in protocol.commands():
        print(line)

# -*- coding: utf-8 -*-
from opentrons import protocol_api
import numpy as np

metadata = {
    "protocolName": "DNABOT Step 3: Assembly with thermocycler Gen2",
    "apiLevel": "2.21",
}


def run(protocol: protocol_api.ProtocolContext):
    tc_mod = protocol.load_module(__HARDWARE["thermocycler"]["id"])
    tc_mod.open_lid()
    tiprack_type = __LABWARES["tiprack_20ul"]["id"]

    def get_low_volume_push_out(volume):
        """Use a modest 2 uL push-out for OT-2 dispenses."""
        return 2

    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type):
        CLIP_PLATE_TYPE = __LABWARES["clip_plate"]["id"]
        CLIP_PLATE_POSITION = "1"
        TUBE_RACK_TYPE = __LABWARES["24_tuberack_1500ul"]["id"]
        TUBE_RACK_POSITION = "4"
        DESTINATION_PLATE_TYPE = __LABWARES["final_assembly_plate"]["id"]
        DESTINATION_PLATE_TYPE = __LABWARES["final_assembly_plate"]["id"]
        clip_plate = protocol.load_labware(CLIP_PLATE_TYPE, CLIP_PLATE_POSITION)
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
        purified_clip_plate = clip_plate
        TOTAL_VOL = 15
        PART_VOL = 1.5
        MIX_SETTINGS = (1, 3)
        sample_number = len(final_assembly_dict.keys())
        if sample_number > 96:
            raise ValueError("Final assembly nummber cannot exceed 96.")
        INITIAL_TIP = "A1"
        CANDIDATE_TIPRACK_SLOTS = ["3", "5", "6", "9"]
        PIPETTE_TYPE = __HARDWARE["single_pipette"]["id"]
        PIPETTE_MOUNT = __HARDWARE["single_pipette_mount"]["id"]
        if tiprack_num > len(CANDIDATE_TIPRACK_SLOTS):
            raise ValueError("Assembly requires more tipracks than the deck supports.")
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        pipette = protocol.load_instrument(
            PIPETTE_TYPE, PIPETTE_MOUNT, tip_racks=tipracks
        )
        pipette.flow_rate.aspirate = 8
        pipette.flow_rate.dispense = 8
        pipette.flow_rate.blow_out = 15
        high = 2
        normal = 1
        slow = 0.4
        vslow = 0.2
        tc_mod.open_lid()
        tc_mod.deactivate_lid()
        tc_mod.set_block_temperature(temperature=4)
        final_assembly_lengths = []
        for values in final_assembly_dict.values():
            final_assembly_lengths.append(len(values))
        unique_assemblies_lengths = list(set(final_assembly_lengths))
        master_mix_well_letters = ["A", "B", "C", "D"]
        for x in unique_assemblies_lengths:
            master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
            destination_inds = [
                i for (i, lengths) in enumerate(final_assembly_lengths) if lengths == x
            ]
            destination_wells = np.array(
                [key for (key, value) in list(final_assembly_dict.items())]
            )
            destination_wells = list(destination_wells[destination_inds])
            pipette.well_bottom_clearance.aspirate = 1
            pipette.well_bottom_clearance.dispense = 2
            master_mix_volume = TOTAL_VOL - x * PART_VOL
            pipette.pick_up_tip()
            for destination_well in destination_wells:
                pipette.aspirate(
                    master_mix_volume, tube_rack[master_mix_well].bottom(1), rate=slow
                )
                pipette.dispense(
                    master_mix_volume,
                    destination_plate[destination_well].bottom(2),
                    rate=slow,
                    push_out=get_low_volume_push_out(master_mix_volume),
                )
            pipette.drop_tip()
        for (key, values) in list(final_assembly_dict.items()):
            for value in values:
                pipette.pick_up_tip()
                pipette.well_bottom_clearance.aspirate = 1
                pipette.well_bottom_clearance.dispense = 2
                pipette.aspirate(
                    PART_VOL, purified_clip_plate[value].bottom(1), rate=slow
                )
                pipette.dispense(
                    PART_VOL,
                    destination_plate[key].bottom(2),
                    rate=slow,
                    push_out=get_low_volume_push_out(PART_VOL),
                )
                pipette.aspirate(10, destination_plate[key].bottom(1), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(3), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(2), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(1), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(3), rate=normal)
                pipette.dispense(10, destination_plate[key].bottom(1), rate=high)
                pipette.aspirate(10, destination_plate[key].bottom(2), rate=slow)
                pipette.dispense(
                    10, destination_plate[key].bottom(3), push_out=3, rate=vslow
                )
                protocol.delay(seconds=5)
                pipette.move_to(destination_plate[key].top(-8))
                pipette.blow_out()
                pipette.touch_tip(radius=0.6, v_offset=-8, speed=10)
                pipette.drop_tip()
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(105)
        tc_mod.set_block_temperature(72, hold_time_minutes=1)
        tc_mod.set_block_temperature(70, hold_time_minutes=1)
        tc_mod.set_block_temperature(68, hold_time_minutes=1)
        tc_mod.set_block_temperature(66, hold_time_minutes=1)
        tc_mod.set_block_temperature(64, hold_time_minutes=2)
        tc_mod.set_block_temperature(62, hold_time_minutes=2)
        tc_mod.set_block_temperature(60, hold_time_minutes=2)
        tc_mod.set_block_temperature(58, hold_time_minutes=2)
        tc_mod.set_block_temperature(56, hold_time_minutes=2)
        tc_mod.set_block_temperature(54, hold_time_minutes=2)
        tc_mod.set_block_temperature(52, hold_time_minutes=2)
        tc_mod.set_block_temperature(50, hold_time_minutes=2)
        tc_mod.set_block_temperature(48, hold_time_minutes=2)
        tc_mod.set_block_temperature(46, hold_time_minutes=2)
        tc_mod.set_block_temperature(44, hold_time_minutes=1)
        tc_mod.set_block_temperature(42, hold_time_minutes=1)
        tc_mod.set_block_temperature(40, hold_time_minutes=1)
        tc_mod.set_block_temperature(4)
        tc_mod.set_lid_temperature(37)
        protocol.delay(minutes=2)
        tc_mod.deactivate_lid()
        tc_mod.open_lid()
        tc_mod.set_block_temperature(4)

    final_assembly(
        final_assembly_dict=final_assembly_dict,
        tiprack_num=tiprack_num,
        tiprack_type=tiprack_type,
    )
    for line in protocol.commands():
        print(line)

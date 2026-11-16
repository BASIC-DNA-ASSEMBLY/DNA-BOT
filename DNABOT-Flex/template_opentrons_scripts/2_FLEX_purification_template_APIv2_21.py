# -*- coding: utf-8 -*-
requirements = {"robotType": "Flex"}
from opentrons import protocol_api
from opentrons.types import Point

metadata = {
    "protocolName": "DNABOT Step 2: Purification",
    "description": "Implements magbead purification reactions for BASIC assembly.",
    "apiLevel": "2.21",
}


def run(protocol: protocol_api.ProtocolContext):
    protocol.load_trash_bin("A3")
    protocol.load_module(module_name=__HARDWARE["thermocycler"]["id"], location="B1")

    def _dedupe(items):
        seen = set()
        ordered = []
        for item in items:
            if item not in seen:
                seen.add(item)
                ordered.append(item)
        return ordered

    def _tip_positions(racks):
        positions = []
        for rack in racks:
            positions.extend(rack.rows()[0])
        return positions

    def magbead(sample_number, ethanol_well):
        elution_buffer_well = "A1"
        sample_volume = 30
        bead_ratio = __PARAMETERS["purif_bead_ratio"]["value"]
        elution_buffer_volume = 40
        incubation_time = __PARAMETERS["purif_incubation_time"]["value"]
        settling_time = __PARAMETERS["purif_settling_time"]["value"]
        drying_time = __PARAMETERS["purif_drying_time"]["value"]
        elution_time = __PARAMETERS["purif_elution_time"]["value"]
        sample_offset = 0
        SMALL_TIPS_PER_COLUMN = 5
        DEAD_TOTAL_VOL = 5
        IMMOBILISE_MIX_REPS = 10
        AIR_VOL_COEFF = 0.1
        ETHANOL_VOL = 150
        WASH_TIME = __PARAMETERS["purif_wash_time"]["value"]
        ETHANOL_DEAD_VOL = 50
        ELUTION_MIX_REPS = 10
        ELUTANT_SEP_TIME = 1
        ELUTION_DEAD_VOL = 2
        BEADS_WELL = "A1"
        LIQUID_WASTE_WELL = "A5"
        if sample_number > 48:
            raise ValueError("sample number cannot exceed 48")
        col_num = sample_number // 8 + (1 if sample_number % 8 > 0 else 0)
        if col_num == 0:
            return
        LARGE_TIPS_REQUIRED = col_num
        pipette_type = __HARDWARE["multi_pipette"]["id"]
        pipette_mount = __HARDWARE["multi_pipette_mount"]["id"]
        pipette_name = pipette_type.lower()
        small_tiprack_type = __LABWARES["flex_96_tiprack_200ul"]["id"]
        large_tiprack_type = __LABWARES["flex_96_tiprack_1000ul"]["id"]
        small_tiprack_num = SMALL_TIPS_PER_COLUMN * col_num // 12 + (
            1 if SMALL_TIPS_PER_COLUMN * col_num % 12 else 0
        )
        large_tiprack_num = LARGE_TIPS_REQUIRED // 12 + (
            1 if LARGE_TIPS_REQUIRED % 12 else 0
        )
        small_slots = ["D3", "C3", "B3"][:small_tiprack_num]
        large_slots = ["C2"][:large_tiprack_num]
        if len(small_slots) < small_tiprack_num:
            raise ValueError("Not enough 200 uL tiprack slots for Flex purification.")
        if len(large_slots) < large_tiprack_num:
            raise ValueError("Not enough 1000 uL tiprack slots for Flex purification.")
        small_tipracks = [
            protocol.load_labware(small_tiprack_type, slot) for slot in small_slots
        ]
        large_tipracks = [
            protocol.load_labware(large_tiprack_type, slot) for slot in large_slots
        ]
        tipracks = small_tipracks + large_tipracks
        pipette = protocol.load_instrument(
            pipette_type, mount=pipette_mount, tip_racks=tipracks
        )
        PIPETTE_ASPIRATE_RATE = 80
        PIPETTE_DISPENSE_RATE = 100
        PIPETTE_BLOW_OUT_RATE = 300
        pipette.flow_rate.aspirate = PIPETTE_ASPIRATE_RATE
        pipette.flow_rate.dispense = PIPETTE_DISPENSE_RATE
        pipette.flow_rate.blow_out = PIPETTE_BLOW_OUT_RATE
        high = 2
        normal = 1
        slow = 0.4
        vslow = 0.2
        small_tip_positions = _tip_positions(small_tipracks)
        large_tip_positions = _tip_positions(large_tipracks)
        small_tip_index = 0

        def pick_up_small_tip():
            nonlocal small_tip_index
            if small_tip_index >= len(small_tip_positions):
                raise ValueError("Ran out of Flex 200 uL tips during purification.")
            pipette.pick_up_tip(small_tip_positions[small_tip_index])
            small_tip_index += 1

        if len(large_tip_positions) < col_num:
            raise ValueError("Not enough dedicated Flex 1000 uL tips for washes.")

        def pick_up_large_tip_for_column(column_index):
            pipette.pick_up_tip(large_tip_positions[column_index])

        def aspirate_beads(volume, well):
            """Bead aspiration stays low in the reservoir and moves slowly."""
            pipette.aspirate(volume, well.bottom(2), rate=slow)
            protocol.delay(seconds=1)

        def aspirate_sample(volume, well):
            """Sample aspiration uses a shallow height and a short settling pause."""
            pipette.aspirate(volume, well.bottom(1), rate=slow)
            protocol.delay(seconds=1)

        def aspirate_eluate(volume, well):
            """Final eluate aspiration is the gentlest transfer in the protocol."""
            pre_wet_volume = min(max(5, volume / 2), pipette.max_volume / 2)
            pipette.aspirate(pre_wet_volume, well.bottom(1), rate=vslow)
            protocol.delay(seconds=0.5)
            pipette.dispense(pre_wet_volume, well.bottom(2), rate=vslow)
            pipette.aspirate(volume, well.bottom(1), rate=vslow)
            protocol.delay(seconds=1)

        def dispense_gently(volume, well, dispense_height=2, push_out=0):
            """Dispense with a conservative finish to improve small-volume consistency."""
            pipette.dispense(volume, well.bottom(dispense_height), rate=high)
            pipette.dispense(
                0, well.bottom(dispense_height), rate=slow, push_out=push_out
            )

        def elution_mix(well, mix_volume, repetitions):
            """Elution mixing uses normal aspirates and stronger dispenses for resuspension."""
            for mix_step in range(repetitions):
                x_offset = 2 if mix_step % 2 == 0 else -2
                dispense_location = well.bottom(3).move(Point(x=x_offset, y=0, z=0))
                pipette.aspirate(mix_volume, well.bottom(1), rate=normal)
                pipette.dispense(mix_volume, dispense_location, rate=high)
            pipette.aspirate(mix_volume, well.bottom(2), rate=normal)
            pipette.dispense(
                mix_volume,
                well.bottom(3).move(Point(x=2, y=0, z=0)),
                rate=high,
                push_out=max(1, mix_volume / 10),
            )

        mag_module = protocol.load_module(__HARDWARE["mag_deck"]["id"], location="D1")
        mag_plate = mag_module.load_labware(__LABWARES["mag_plate"]["id"])
        mix_plate = protocol.load_labware(__LABWARES["mix_plate"]["id"], "C1")
        reagent_container = protocol.load_labware(
            __LABWARES["12_reservoir_21000ul"]["id"], "A2"
        )
        bead_container = protocol.load_labware(
            __LABWARES["96_deepwellplate_2ml"]["id"], "B2"
        )
        samples = [
            col for col in mag_plate.columns()[sample_offset : col_num + sample_offset]
        ]
        mixing = [
            col for col in mix_plate.columns()[sample_offset : col_num + sample_offset]
        ]
        output = [
            col
            for col in mag_plate.columns()[
                6 + sample_offset : col_num + 6 + sample_offset
            ]
        ]
        ethanol = reagent_container[ethanol_well]
        elution_buffer = reagent_container[elution_buffer_well]
        liquid_waste = reagent_container[LIQUID_WASTE_WELL]
        beads = bead_container[BEADS_WELL]
        bead_volume = sample_volume * bead_ratio
        mix_vol = min(bead_volume / 2, pipette.max_volume)
        total_vol = bead_volume + sample_volume + DEAD_TOTAL_VOL
        for target in range(col_num):
            pick_up_small_tip()
            aspirate_beads(bead_volume, beads)
            aspirate_sample(sample_volume + DEAD_TOTAL_VOL, samples[target][0])
            dispense_gently(total_vol, mixing[target][0], dispense_height=2, push_out=1)
            pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, mixing[target][0])
            pipette.blow_out(mixing[target][0].top())
            pipette.drop_tip()
        protocol.delay(minutes=incubation_time)
        for target in range(col_num):
            pick_up_small_tip()
            aspirate_sample(total_vol, mixing[target][0])
            dispense_gently(
                total_vol, samples[target][0], dispense_height=2, push_out=1
            )
            pipette.blow_out(samples[target][0].top())
            pipette.drop_tip()
        protocol.delay(minutes=settling_time)
        for target in range(col_num):
            pick_up_small_tip()
            aspirate_sample(total_vol, samples[target][0])
            pipette.dispense(total_vol, liquid_waste.top(5), rate=normal)
            pipette.blow_out(liquid_waste.top(5))
            pipette.drop_tip()
        air_vol = min(pipette.max_volume * AIR_VOL_COEFF, 20)
        for wash_cycle in range(2):
            for (target, column) in enumerate(samples):
                well = column[0]
                pick_up_large_tip_for_column(target)
                pipette.distribute(
                    ETHANOL_VOL,
                    ethanol.bottom(2),
                    well.bottom(5),
                    air_gap=air_vol,
                    new_tip="never",
                )
                pipette.return_tip()
            protocol.delay(minutes=WASH_TIME)
            for target in range(col_num):
                pick_up_large_tip_for_column(target)
                pipette.aspirate(
                    ETHANOL_VOL + ETHANOL_DEAD_VOL, samples[target][0], rate=normal
                )
                pipette.air_gap(air_vol)
                pipette.dispense(
                    ETHANOL_VOL + ETHANOL_DEAD_VOL + air_vol, liquid_waste.top(5)
                )
                pipette.blow_out(liquid_waste.top(5))
                if wash_cycle == 0:
                    pipette.return_tip()
                else:
                    pipette.drop_tip()
        protocol.delay(minutes=drying_time)
        protocol.move_labware(mix_plate, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(mag_plate, "C1", use_gripper=False)
        mix_vol = min(elution_buffer_volume / 2, pipette.max_volume)
        for target in range(col_num):
            pick_up_small_tip()
            aspirate_sample(elution_buffer_volume, elution_buffer)
            dispense_gently(
                elution_buffer_volume,
                samples[target][0],
                dispense_height=2,
                push_out=max(1, elution_buffer_volume / 20),
            )
            elution_mix(samples[target][0], mix_vol, ELUTION_MIX_REPS)
            pipette.blow_out(samples[target][0].top())
            pipette.drop_tip()
        protocol.delay(minutes=elution_time)
        protocol.move_labware(mag_plate, mag_module, use_gripper=False)
        protocol.delay(minutes=ELUTANT_SEP_TIME)
        for (target, dest) in zip(samples, output):
            pick_up_small_tip()
            aspirate_eluate(elution_buffer_volume - ELUTION_DEAD_VOL, target[0])
            dispense_gently(
                elution_buffer_volume - ELUTION_DEAD_VOL,
                dest[0],
                dispense_height=2,
                push_out=max(1, (elution_buffer_volume - ELUTION_DEAD_VOL) / 20),
            )
            pipette.drop_tip()

    magbead(sample_number=sample_number, ethanol_well=ethanol_well)

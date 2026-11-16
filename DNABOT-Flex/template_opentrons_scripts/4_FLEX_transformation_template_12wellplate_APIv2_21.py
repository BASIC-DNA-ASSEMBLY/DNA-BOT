"""Transformation and sequential 12-well agar plating for Flex."""
from opentrons import protocol_api

requirements = {"robotType": "Flex"}
metadata = {
    "apiLevel": "2.21",
    "protocolName": "DNABOT Flex Step 4: Transformation and 12-well plating",
}


def run(protocol: protocol_api.ProtocolContext):
    protocol.load_trash_bin("A3")
    tc_mod = protocol.load_module(__HARDWARE["thermocycler"]["id"], "B1")
    tc_mod.open_lid()
    assembly_plate = protocol.load_labware(
        __LABWARES["final_assembly_plate"]["id"], "C2"
    )
    transformation_plate = tc_mod.load_labware(__LABWARES["transform_plate"]["id"])
    soc_plate = protocol.load_labware(__LABWARES["96_deepwellplate_2ml"]["id"], "C1")
    agar_plate = protocol.load_labware(__LABWARES["12_corning_wellplate"]["id"], "D1")
    single_racks = [
        protocol.load_labware(__LABWARES["tiprack_20ul"]["id"], slot)
        for slot in ["D2", "B3"]
    ]
    multi_racks = [
        protocol.load_labware(__LABWARES["tiprack_300ul"]["id"], slot)
        for slot in ["D3", "C3"]
    ]
    single = protocol.load_instrument(
        __HARDWARE["single_pipette"]["id"],
        __HARDWARE["single_pipette_mount"]["id"],
        tip_racks=single_racks,
    )
    multi = protocol.load_instrument(
        __HARDWARE["multi_pipette"]["id"],
        __HARDWARE["multi_pipette_mount"]["id"],
        tip_racks=multi_racks,
    )

    def pickup(pipette):
        try:
            pipette.pick_up_tip()
        except protocol_api.labware.OutOfTipsError:
            protocol.pause(
                "Refill the empty tipracks for " + pipette.name + ", then resume."
            )
            pipette.reset_tipracks()
            pipette.pick_up_tip()

    def volume_mode(volume):
        single.configure_for_volume(max(1, min(float(volume), 50)))

    wells = list(
        dict.fromkeys(
            well for sources, targets, volumes in spotting_tuples for well in sources
        )
    )
    if not wells or len(wells) > 96:
        raise ValueError("Transformation requires between 1 and 96 assemblies.")
    columns = list(dict.fromkeys(well[1:] for well in wells))
    tc_mod.set_block_temperature(
        __PARAMETERS["transform_incubation_temp"]["value"], block_max_volume=50
    )
    protocol.pause(
        "Place competent cells in the transformation plate on the thermocycler. Populate all eight wells of each used column, then resume."
    )
    for well in wells:
        volume_mode(5)
        pickup(single)
        single.transfer(
            5,
            assembly_plate[well].bottom(1),
            transformation_plate[well].bottom(2),
            new_tip="never",
            mix_after=(4, 5),
        )
        single.drop_tip()
    protocol.delay(minutes=__PARAMETERS["transform_incubation_time"]["value"])
    tc_mod.close_lid()
    tc_mod.set_block_temperature(42, hold_time_seconds=30, block_max_volume=50)
    tc_mod.set_block_temperature(4, hold_time_minutes=2, block_max_volume=50)
    tc_mod.open_lid()
    tc_mod.set_block_temperature(20, block_max_volume=50)
    # An eight-channel head takes SOC from all eight wells in the selected column.
    soc_source = soc_plate.columns_by_name()[soc_well[1:]][0]
    protocol.pause(
        "Confirm all eight wells of SOC column "
        + soc_well[1:]
        + " contain sufficient medium, then resume."
    )
    for col in columns:
        pickup(multi)
        multi.transfer(
            100,
            soc_source.bottom(2),
            transformation_plate.columns_by_name()[col][0].bottom(2),
            new_tip="never",
            mix_after=(4, 50),
        )
        multi.drop_tip()
    protocol.pause(
        "Incubate the transformation plate in an external incubator at 37 C for 60 minutes. Return the plate to the open thermocycler, then resume."
    )
    # Each tuple describes a fresh 12-well agar plate, in column-major order.
    for plate_index, (sources, targets, volumes) in enumerate(spotting_tuples, 1):
        protocol.pause(f"Place agar plate {plate_index} in slot D1, then resume.")
        for source, target, total in zip(sources, targets, volumes):
            remaining = float(total)
            if remaining <= 0:
                raise ValueError("Spotting volume must be positive.")
            volume_mode(20)
            pickup(single)
            single.mix(4, 20, transformation_plate[source].bottom(2))
            while remaining > 0:
                amount = min(20, remaining)
                volume_mode(amount)
                single.aspirate(amount, transformation_plate[source].bottom(1))
                single.dispense(amount, agar_plate[target].top(-5), rate=0.2)
                single.move_to(agar_plate[target].top(-13), speed=30)
                single.move_to(agar_plate[target].top(-5))
                single.blow_out(agar_plate[target].top(-5))
                single.move_to(agar_plate[target].top(7))
                remaining -= amount
            single.drop_tip()
    tc_mod.set_block_temperature(4)
    tc_mod.deactivate_lid()

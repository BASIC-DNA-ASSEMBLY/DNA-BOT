from opentrons import protocol_api
import json

metadata = {
    'protocolName': 'Media_Bot: Media Assembly and Subculturing (Flex)',
    'author': 'Liam Hallett, Shishun Liang',
    'description': 'Automated media assembly and subculture from preculture plate based on CSV input for Opentrons Flex'
}
requirements = {"robotType": "Flex", "apiLevel": "2.15"}

# Mix settings
MIX_BEFORE_VOLUME = 150 # µL
MIX_BEFORE_REPS = 3     # Number of mix repetitions before transfer
MIX_AFTER_VOLUME = 50   # µL
MIX_AFTER_REPS = 3      # Number of mix repetitions after transfer
MIX_SPEED = 1.0         # Rate for mixing (1.0 = 100% speed)
DISPENSE_RATE = 0.1     # Rate for dispensing (0.1 = 10% speed)

def custom_transfer(pipette, volume, source, destination, mix_before=None, mix_after=None, mix_speed=1.0, dispense_speed=0.4, new_tip='once'):
    """
    Custom transfer function with fine control over mixing and dispense speeds.
    
    Args:
        pipette: The pipette to use for transfer
        volume: Volume to transfer in µL
        source: Source well
        destination: Destination well
        mix_before: Tuple of (repetitions, volume) for mixing before transfer, or None
        mix_after: Tuple of (repetitions, volume) for mixing after transfer, or None
        mix_speed: Rate for mixing (0-1, default 1.0 = 100% speed)
        dispense_speed: Rate for dispensing (0-1, default 0.4 = 40% speed)
        new_tip: When to get a new tip ('once', 'always', or 'never')
    """
    # Handle tip management
    if new_tip == 'always':
        pipette.pick_up_tip()
    elif new_tip == 'once' and not pipette.has_tip:
        pipette.pick_up_tip()
    
    # Mix before if specified
    if mix_before:
        reps, mix_vol = mix_before
        pipette.mix(reps, mix_vol, source, rate=mix_speed)
    
    # Aspirate at default speed
    pipette.aspirate(volume, source)
    
    # Dispense at specified speed
    pipette.dispense(volume, destination, rate=dispense_speed)
    
    # Mix after if specified
    if mix_after:
        reps, mix_vol = mix_after
        pipette.mix(reps, mix_vol, destination, rate=mix_speed)
    
    # Handle tip disposal
    if new_tip == 'always':
        pipette.drop_tip()

class ReservoirManager:
    def __init__(self, protocol, reservoir, well_capacity=20000):  # 20ml = 20000ul
        self.protocol = protocol
        self.reservoir = reservoir
        self.well_capacity = well_capacity
        self.volumes = {}  # Track volume in each well
        self.well_mapping = {}  # Map media/supplement to list of wells
        self.current_well_indices = {}  # Track which well index we're using for each component
    
    def add_wells(self, name, wells):
        """Add wells for a specific media/supplement"""
        if name not in self.well_mapping:
            self.well_mapping[name] = []
            self.current_well_indices[name] = 0
        self.well_mapping[name].extend(wells)
        for well in wells:
            self.volumes[well] = self.well_capacity
    
    def get_next_well(self, name, required_volume):
        """Get the next available well(s) for a media/supplement that has enough volume.
        
        Returns:
            - A single well if one well has enough volume
            - A dict with 'wells' and 'volumes' if multiple wells are needed
        """
        wells = self.well_mapping[name]
        current_index = self.current_well_indices[name]
        
        # Try to use current well first
        current_well = wells[current_index]
        if self.volumes[current_well] >= required_volume:
            return current_well
        
        # If current well doesn't have enough, try next wells
        remaining_volume = required_volume
        used_wells = []
        used_volumes = []
        
        for i in range(len(wells)):
            well_index = (current_index + i) % len(wells)
            well = wells[well_index]
            
            if self.volumes[well] >= remaining_volume:
                # This well has enough for the remainder
                used_wells.append(well)
                used_volumes.append(remaining_volume)
                self.current_well_indices[name] = well_index  # Update current index
                break
            elif self.volumes[well] > 0:
                # Use what's left in this well
                used_volume = self.volumes[well]
                remaining_volume -= used_volume
                used_wells.append(well)
                used_volumes.append(used_volume)
        
        if remaining_volume > 0:
            raise ValueError(f"No available wells with sufficient volume for {name}. Required: {required_volume}µL")
        
        if len(used_wells) == 1:
            return used_wells[0]
        else:
            return {
                "wells": used_wells,
                "volumes": used_volumes
            }
    
    def use_volume(self, name, volume):
        """Use volume from the appropriate well(s), transitioning if necessary"""
        next_wells = self.get_next_well(name, volume)
        
        if isinstance(next_wells, dict):
            # Multiple wells needed
            for well, vol in zip(next_wells["wells"], next_wells["volumes"]):
                self.volumes[well] -= vol
            return next_wells
        else:
            # Single well
            self.volumes[next_wells] -= volume
            return next_wells

def transfer(experiment_plan, component_name, source_reservoir, mix_plate, single_pipette, media_tips, mix_before=False):
    """Transfer a specific component to multiple wells in the mix plate
    
    Args:
        experiment_plan: List of experiment entries
        component_name: Name of the component to transfer
        source_reservoir: ReservoirManager instance for the source
        mix_plate: Destination plate for mixing
        single_pipette: Single channel pipette
        media_tips: Tip rack for media transfers
        mix_before: Boolean to control whether to mix source before aspirating
    """
    # Create list of wells and volumes for this component
    wells_and_volumes = []
    total_volume_needed = 0
    excess_volume = 1.05
    
    # Make a list of wells and volumes for the component
    for entry in experiment_plan:
        if component_name == entry["media_type"]:
            volume = entry["media_volume"]
        elif component_name == "Water":
            volume = entry["water_volume"]
        elif component_name in entry["supplements"]:
            volume = entry["supplement_volume"]  # Use the calculated supplement volume from experiment plan
        else:
            continue
        
        if volume > 0:
            wells_and_volumes.append({
                "well": mix_plate.wells_by_name()[entry["well"]],
                "volume": volume
            })
            total_volume_needed += volume
    
    # Check if we have enough volume in the reservoir
    try:
        source_reservoir.get_next_well(component_name, total_volume_needed * excess_volume)
    except ValueError as e:
        raise ValueError(f"Insufficient volume for {component_name}. Required: {total_volume_needed}µL") from e
    
    if wells_and_volumes:
        # Pick up the next available tip
        single_pipette.pick_up_tip()
        
        # Get maximum volume from tip specifications
        MAX_VOLUME = media_tips.wells()[0].max_volume  # µl
        
        # Process transfers in batches
        current_batch = []
        current_volume = 0
        
        def process_batch(batch):
            """Helper function to process a batch of transfers.
            
            Handles both single and multiple source wells. When a well is nearly empty,
            the remaining volume is used first, then the required volume is taken from
            the next available well.
            """
            if not batch:
                return
                
            total_volume = sum(t["volume"] for t in batch) * excess_volume
            source_info = source_reservoir.use_volume(component_name, total_volume)
            
            # Handle aspiration from source well(s)
            if isinstance(source_info, dict):
                # Multiple source wells - aspirate from each in sequence
                wells = source_info["wells"]
                volumes = source_info["volumes"]
                source_well = wells[0]  # Use first well for blowout
            else:
                # Single source well - convert to same format for consistent handling
                wells = [source_info]
                volumes = [total_volume]
                source_well = source_info  # Use the single source well for blowout
            
            # Single loop to handle both single and multiple wells
            for well, vol in zip(wells, volumes):
                # Mix before aspirating if mix_before is True
                if mix_before:
                    single_pipette.mix(MIX_BEFORE_REPS, MIX_BEFORE_VOLUME, well)
                # Aspirate the volume (already includes excess)
                single_pipette.aspirate(vol, well)
            
            # Dispense to all wells in the batch
            for t in batch:
                single_pipette.dispense(t["volume"], t["well"].top())
                single_pipette.touch_tip()
            
            # Blow out remaining volume back to source well
            remaining_volume = sum(volumes) * (excess_volume - 1)
            single_pipette.blow_out(source_well)
        
        # Process all transfers in appropriate batch sizes
        for transfer in wells_and_volumes:
            # Check if adding this transfer would exceed max volume
            volume_required = (current_volume + transfer["volume"]) * excess_volume
            if volume_required > MAX_VOLUME:
                # Process current batch if it exists
                process_batch(current_batch)
                
                # Start new batch with current transfer
                current_batch = [transfer]
                current_volume = transfer["volume"]
            else:
                # Add to current batch
                current_batch.append(transfer)
                current_volume += transfer["volume"]
        
        # Process the final batch
        process_batch(current_batch)
        
        # Drop the tip
        single_pipette.drop_tip()

def check_column_match(experiment_plan, strain_map, column_wells):
    """
    Check if all strains in a column match a single column in the preculture plate,
    including order and handling partial columns.
    
    Args:
        experiment_plan (list): List of dictionaries containing experiment details
        strain_map (dict): Mapping of strains to preculture wells
        column_wells (list): List of wells in the current column
        
    Returns:
        tuple: (bool, int) - (whether all strains match a single column in order, current column number)
    """
    # Create a mapping from wells to their experiment entries for quick lookup
    well_to_entry = {entry["well"]: entry for entry in experiment_plan}
    
    # Get the strains in the current column, preserving order
    current_strains = []
    for well in column_wells:
        if well in well_to_entry:
            current_strains.append(well_to_entry[well]["strain"])
        else:
            current_strains.append(None)
    
    # If no strains in this column, no match
    if all(s is None for s in current_strains):
        return False, None
    
    # Create reverse mapping from wells to strains
    well_to_strain = {well: strain for strain, well in strain_map.items()}
    
    # Check each column in the preculture plate
    for col in range(1, 13):  # Assuming 12 columns in preculture plate
        # Get strains in this preculture column
        preculture_strains = []
        for row in range(8):  # Assuming 8 rows
            well = f"{chr(65 + row)}{col}"
            preculture_strains.append(well_to_strain.get(well))
        
        # Directly compare the entire lists, including None values
        if current_strains == preculture_strains:
            return True, col
    
    # If no matching column found
    return False, None

def get_unique_columns(experiment_plan):
    """
    Get a sorted list of unique column numbers from the experiment plan.
    
    Args:
        experiment_plan (list): List of experiment entries
        
    Returns:
        list: Sorted list of unique column numbers
    """
    # Extract column numbers from well positions (e.g., "A1" -> 1)
    columns = set()
    for entry in experiment_plan:
        well = entry["well"]
        column = int(well[1:])  # Get the number part of the well
        columns.add(column)
    
    return sorted(list(columns))

def run(protocol: protocol_api.ProtocolContext):
    ##### Load Data #####
    with open('experiment_plan.json') as f:
        experiment_plan = json.load(f)

    with open('supplement_map.json') as f:
        supplement_map = json.load(f)
        
    with open('strain_map.json') as f:
        strain_map = json.load(f)

    # Turn off lights at start
    protocol.set_rail_lights(False)

    # Number of replicate plates (default 3, can be modified by parser)
    reps = 3
    
    # Verify all strains in experiment plan exist in strain map
    for entry in experiment_plan:
        if entry["strain"] not in strain_map and entry["strain"] not in ["WT", "Blank"]:
            raise ValueError(f"Strain {entry['strain']} not found in strain map")

    ##### Labware Setup #####
    final_plates = [
        protocol.load_labware('corning_96_wellplate_360ul_flat', slot, f'Final Plate {i+1}')
        for i, slot in enumerate(['1', '2', '3'][0:reps])
    ]
    
    # Load tipracks - 1000ul for media transfers, 200ul for subculture transfers
    media_tips = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', '11', 'Media Tips')
    transfer_tips = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', '4', 'Transfer Tips')
    
    # Load Flex pipettes but keep variable names the same
    multi_pipette = protocol.load_instrument('flex_8channel_1000', mount='right', tip_racks=[transfer_tips])
    single_pipette = protocol.load_instrument('flex_1channel_1000', mount='left', tip_racks=[media_tips, transfer_tips])

    mix_plate = protocol.load_labware('4ti0136_96_wellplate_2200ul', '5', 'Mix Plate')
    preculture_plate = protocol.load_labware('4ti0136_96_wellplate_2200ul', '6', 'Preculture Plate')
    
    # Load reservoirs based on supplement map
    reservoirs = {}
    for sup, info in supplement_map.items():
        plate_name = info["plate"]
        if plate_name not in reservoirs:
            if plate_name == "Media reservoir":
                reservoirs[plate_name] = protocol.load_labware('4ti0131_12_reservoir_21000ul', '8', plate_name)
            elif plate_name == "Supplement reservoir":  # Supplement reservoir
                reservoirs[plate_name] = protocol.load_labware('4ti0136_96_wellplate_2200ul', '9', plate_name)
            else:
                ValueError(f"Invalid plate name: {plate_name}")

    # Create reservoir managers
    media_reservoir_manager = ReservoirManager(protocol, reservoirs["Media reservoir"])
    supplement_reservoir_manager = ReservoirManager(protocol, reservoirs["Supplement reservoir"], well_capacity=2000)

    # Set up wells for each media/supplement
    for sup, info in supplement_map.items():
        if "Media" in sup or sup.startswith(("1x_", "2x_")):
            wells = [reservoirs[info["plate"]].wells_by_name()[well] for well in info["wells"]]
            media_reservoir_manager.add_wells(sup, wells)
        elif sup == "Water":
            wells = [reservoirs[info["plate"]].wells_by_name()[well] for well in info["wells"]]
            media_reservoir_manager.add_wells(sup, wells)
        else:
            wells = [reservoirs[info["plate"]].wells_by_name()[well] for well in info["wells"]]
            supplement_reservoir_manager.add_wells(sup, wells)

    ##### Media Manufacturing Stage #####
    # Transfer media components using single channel pipette
    if "Media reservoir" in reservoirs:
        # Transfer media
        for media_type in set(entry["media_type"] for entry in experiment_plan):
            transfer(experiment_plan, media_type, media_reservoir_manager, mix_plate, single_pipette, media_tips, mix_before=True)
        
        # Transfer water
        transfer(experiment_plan, "Water", media_reservoir_manager, mix_plate, single_pipette, media_tips)
    
    # Transfer supplements
    if "Supplement reservoir" in reservoirs:
        for sup in set(sup for entry in experiment_plan for sup in entry["supplements"]):
            transfer(experiment_plan, sup, supplement_reservoir_manager, mix_plate, single_pipette, media_tips)

    ##### Media Distribution Stage #####
    # Get unique columns from experiment plan
    unique_columns = get_unique_columns(experiment_plan)
    
    # Transfer media to final plates using multichannel pipette
    for col in unique_columns:
        # Define well position for this column
        well_pos = f"A{col}"
        
        # Get wells and tips
        tip = transfer_tips.wells_by_name()[well_pos]
        mix_well = mix_plate.wells_by_name()[well_pos]
        
        # Pick up tips from corresponding position
        multi_pipette.pick_up_tip(tip)
        
        # Transfer to each final plate
        for plate in final_plates:
            # Get destination well (top of column)
            dest_well = plate.wells_by_name()[well_pos]
            
            # Transfer 90ul from mix plate to final plate
            custom_transfer(multi_pipette, 90, mix_well, dest_well,
                          mix_before=(MIX_BEFORE_REPS, MIX_BEFORE_VOLUME),
                          mix_speed=MIX_SPEED,
                          dispense_speed=MIX_SPEED * DISPENSE_RATE)
        
        # Return tips attached for subculture stage
        multi_pipette.return_tip()

    ##### Subculture Stage #####
    # First transfer all cultures from preculture to mix plate
    for col in unique_columns:
        # Define well position for this column
        column_well = f"A{col}"
        
        # Get tip for the top of the column
        tip = transfer_tips.wells_by_name()[column_well]
        
        # Get all wells in this column
        column_wells = [f"{chr(65 + i)}{col}" for i in range(8)]
        
        # Check if column matches
        column_matches, current_preculture_col = check_column_match(experiment_plan, strain_map, column_wells)
        
        if column_matches and current_preculture_col is not None:
            # Use multichannel for matching column
            preculture_well = f"A{current_preculture_col}"
            mix_well = column_well
            
            multi_pipette.pick_up_tip(tip)
            
            # Transfer from preculture to mix plate
            custom_transfer(multi_pipette, 10,
                          preculture_plate.wells_by_name()[preculture_well],
                          mix_plate.wells_by_name()[mix_well],
                          mix_before=(MIX_BEFORE_REPS, MIX_BEFORE_VOLUME),
                          mix_after=(MIX_AFTER_REPS, MIX_AFTER_VOLUME),
                          mix_speed=MIX_SPEED,
                          dispense_speed=MIX_SPEED * DISPENSE_RATE)
            
            multi_pipette.return_tip()
        else:
            # Use single channel for non-matching column
            # Transfer each well individually
            for well in column_wells:
                # Find the experiment entry for this well
                entry = next((e for e in experiment_plan if e["well"] == well), None)
                if entry:
                    strain = entry["strain"]
                    preculture_well = strain_map[strain]
                    
                    # Get tip from transfer tips at the same position as the destination well
                    tip = transfer_tips.wells_by_name()[well]
                    single_pipette.pick_up_tip(tip)
                    
                    custom_transfer(single_pipette, 10,
                                  preculture_plate.wells_by_name()[preculture_well],
                                  mix_plate.wells_by_name()[well],
                                  mix_before=(MIX_BEFORE_REPS, MIX_BEFORE_VOLUME),
                                  mix_after=(MIX_AFTER_REPS, MIX_AFTER_VOLUME),
                                  mix_speed=MIX_SPEED,
                                  dispense_speed=MIX_SPEED * DISPENSE_RATE,
                                  new_tip='never')  # We're managing tips manually
                    
                    single_pipette.return_tip()
    
    # Then transfer from mix plate to final plates
    for col in unique_columns:
        # Define well position for this column
        column_well = f"A{col}"
        
        # Get tip for the top of the column
        tip = transfer_tips.wells_by_name()[column_well]
        
        mix_well = mix_plate.wells_by_name()[column_well]
        
        multi_pipette.pick_up_tip(tip)
        
        # Transfer to final plates
        for plate in final_plates:
            # Get destination well (top of column)
            dest_well = plate.wells_by_name()[column_well]
            
            # Transfer 90ul from mix plate to final plate
            custom_transfer(multi_pipette, 10, mix_well, dest_well,
                          mix_before=(MIX_BEFORE_REPS, MIX_BEFORE_VOLUME),
                          mix_after=(MIX_AFTER_REPS, MIX_AFTER_VOLUME),
                          mix_speed=MIX_SPEED,
                          dispense_speed=MIX_SPEED * DISPENSE_RATE)
        
        multi_pipette.drop_tip()

    # Flash lights 3 times at end to indicate completion
    for _ in range(3):
        protocol.set_rail_lights(True)
        protocol.delay(seconds=0.3)
        protocol.set_rail_lights(False)
        protocol.delay(seconds=0.3)


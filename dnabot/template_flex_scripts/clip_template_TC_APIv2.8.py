from opentrons import protocol_api
import json
import time
from typing import Dict, List, Optional, Tuple, Union

# Rename to 'clip_template' and paste into 'template_flex_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.15',
     'robotType': 'Flex',
     'protocolName': 'CLIP_With_Thermocycler_Flex',
     'description': 'Implements linker ligation reactions using an opentrons Flex, including the thermocycler module.'}

# Load CLIP data from JSON file
# This will be replaced by the parser with embedded JSON data
with open('clips_data.json') as f:
    clips_dict = json.load(f)

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# opentrons_simulate.exe dnabot\template_flex_scripts\clip_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["1", "1", "1", "1", "1"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

class TipManager:
    """Manages pipette tips and tracks usage."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 slot: str,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 tip_type: str):
        self.protocol = protocol
        self.pipette = pipette
        self.tip_type = tip_type
        self.rows = 8  # Number of rows in a standard tip rack
        self.cols = 12  # Number of columns in a standard tip rack
        
        # Initialize tip racks array and current rack index
        self.tipracks = []
        self.current_rack = 0
        
        # Load first tip rack using add_tip_rack
        self.add_tip_rack(slot)
        
        self._initialise_tip_arrays()
    
    def add_tip_rack(self, slot: str, tip_type: Optional[str] = None) -> None:
        """Add an additional tip rack to the manager."""
        tip_type = tip_type or self.tip_type
        
        # Check pipette compatibility
        pipette_type = 'p1000' if self.pipette.max_volume >= 1000 else 'p300'
        if tip_type != pipette_type:
            raise ValueError(f"Tip type '{tip_type}' is incompatible with pipette type '{pipette_type}'")
            
        if tip_type == 'p1000':
            self.tipracks.append(self.protocol.load_labware('opentrons_flex_96_tiprack_1000ul', slot))
        elif tip_type == 'p300':
            self.tipracks.append(self.protocol.load_labware('opentrons_flex_96_tiprack_300ul', slot))
        else:
            raise ValueError(f"Unsupported tip type: {tip_type}")
    
    def _initialise_tip_arrays(self) -> None:
        """Initialise tip tracking arrays for normal and inverse tip selection."""
        total_tips = len(self.tipracks[self.current_rack].wells())
        self.normal_tips = list(range(total_tips))  # For normal tip selection
        self.inverse_tips = []  # For inverse tip selection
        
        # Create inverse order array
        for col in range(self.cols):
            for row in range(self.rows-1, -1, -1):  # Start from highest row (H) to lowest (A)
                tip_index = col * self.rows + row
                self.inverse_tips.append(tip_index)
    
    def get_single_tip(self, inverse: bool = False) -> None:
        """Pick up a single tip, optionally using inverse selection to use multichannel pipette for single channel functionality."""
        if not self.normal_tips:  # If either array is empty, we need new tips
            if self.current_rack < len(self.tipracks) - 1:
                # Switch to next tip rack
                self.current_rack += 1
                self._initialise_tip_arrays()
                print(f"Switched to tip rack {self.current_rack + 1}")
            else:
                self._prompt_tip_replacement()
        
        if inverse:
            tip_index = self.inverse_tips[0]  # Get the first tip in inverse order
        else:
            tip_index = min(self.normal_tips)  # Will give us the next tip in normal order
        
        # Remove the tip from both arrays
        self.normal_tips.remove(tip_index)
        self.inverse_tips.remove(tip_index)
        
        self.pipette.pick_up_tip(self.tipracks[self.current_rack].wells()[tip_index])
    
    def _prompt_tip_replacement(self) -> None:
        """Prompt user to replace tip rack and flash lights."""
        # Flash lights 3 times before the prompt
        for _ in range(3):
            self.protocol.set_rail_lights(False)
            time.sleep(0.15)
            self.protocol.set_rail_lights(True)
            time.sleep(0.15)
        
        # Keep lights on and prompt user
        self.protocol.set_rail_lights(True)
        print(f"Please replace the {self.tip_type} tip rack and press 'Resume'")
        self.protocol.pause()
        
        # Reset tip arrays for the new rack
        self._initialise_tip_arrays()

class MasterMixManager:
    """Manages master mix tubes and tracks volume usage."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 labware: protocol_api.labware.Labware,
                 tube_volumes: List[float],
                 transfer_volume: float,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 dead_volume: float = 15.0):
        self.protocol = protocol
        self.labware = labware
        self.tube_volumes = tube_volumes.copy()  # Copy to avoid modifying original
        self.transfer_volume = transfer_volume
        self.pipette = pipette
        self.dead_volume = dead_volume
        self.current_tube = 0
        
    def get_current_tube(self) -> protocol_api.labware.Well:
        """Get the current master mix tube."""
        return self.labware.wells()[self.current_tube]
    
    def can_aspirate(self, volume: float) -> bool:
        """Check if current tube has enough volume for aspiration."""
        return self.tube_volumes[self.current_tube] >= volume + self.dead_volume
    
    def use_volume(self, volume: float) -> None:
        """Use volume from current tube and switch if necessary."""
        if not self.can_aspirate(volume):
            # Switch to next tube
            self.current_tube += 1
            if self.current_tube >= len(self.tube_volumes):
                raise ValueError("No more master mix tubes available")
            print(f"Switched to master mix tube {self.current_tube + 1}")
        
        self.tube_volumes[self.current_tube] -= volume
    
    def distribute_to_wells(self, wells: List[protocol_api.labware.Well], pipette: protocol_api.instrument_context.InstrumentContext) -> None:
        """Distribute master mix to multiple wells efficiently."""
        if not wells:
            return
        
        # Calculate total volume needed
        total_volume = len(wells) * self.transfer_volume
        
        # Check if we have enough volume
        if not self.can_aspirate(total_volume):
            raise ValueError(f"Insufficient master mix volume. Need {total_volume}µL, have {self.tube_volumes[self.current_tube]}µL")
        
        # Pick up tip
        pipette.pick_up_tip()
        
        # Aspirate total volume
        source_well = self.get_current_tube()
        pipette.aspirate(total_volume, source_well)
        
        # Distribute to all wells
        for well in wells:
            pipette.dispense(self.transfer_volume, well)
        
        # Blow out remaining volume
        pipette.blow_out(source_well)
        
        # Drop tip
        pipette.drop_tip()
        
        # Update volume tracking
        self.use_volume(total_volume)

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

def batch_transfer(transfers, pipette, tip_manager, mix_after=None, mix_speed=1.0, dispense_speed=0.4):
    """
    Efficiently transfer multiple volumes using batch processing.
    
    Args:
        transfers: List of dicts with 'source', 'destination', 'volume' keys
        pipette: The pipette to use
        tip_manager: TipManager instance
        mix_after: Optional tuple of (repetitions, volume) for mixing after each transfer
        mix_speed: Rate for mixing (0-1)
        dispense_speed: Rate for dispensing (0-1)
    """
    if not transfers:
        return
    
    # Get maximum volume from pipette
    max_volume = pipette.max_volume
    excess_volume = 1.05  # 5% excess for blowout
    
    def process_batch(batch):
        """Process a batch of transfers."""
        if not batch:
            return
        
        # Calculate total volume needed for batch
        total_volume = sum(t['volume'] for t in batch) * excess_volume
        
        # Pick up tip
        tip_manager.get_single_tip()
        
        # Aspirate total volume from first source
        first_source = batch[0]['source']
        pipette.aspirate(total_volume, first_source)
        
        # Dispense to each destination
        for transfer in batch:
            pipette.dispense(transfer['volume'], transfer['destination'])
            
            # Mix after if specified
            if mix_after:
                reps, mix_vol = mix_after
                pipette.mix(reps, mix_vol, transfer['destination'], rate=mix_speed)
        
        # Blow out remaining volume
        pipette.blow_out(first_source)
        pipette.drop_tip()
    
    # Process transfers in batches
    current_batch = []
    current_volume = 0
    
    for transfer in transfers:
        # Check if adding this transfer would exceed max volume
        volume_required = (current_volume + transfer['volume']) * excess_volume
        if volume_required > max_volume:
            # Process current batch
            process_batch(current_batch)
            
            # Start new batch
            current_batch = [transfer]
            current_volume = transfer['volume']
        else:
            # Add to current batch
            current_batch.append(transfer)
            current_volume += transfer['volume']
    
    # Process final batch
    process_batch(current_batch)

def run(protocol: protocol_api.ProtocolContext):
    # added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    TIPRACK_SLOT = '11'
    TIPRACK_TYPE = 'p1000'
    
    # Master mix tubes
    MASTER_MIX_SLOT = '8'
    MASTER_MIX_TUBES = 8  # Number of master mix tubes
    MASTER_MIX_VOLUME = 2000  # Volume per tube in µL
    MASTER_MIX_TRANSFER_VOLUME = 8  # Volume to transfer per reaction in µL
    
    # Water
    WATER_SLOT = '9'
    WATER_VOLUME = 2000  # Volume per tube in µL
    
    # Parts plates
    PARTS_PLATES_SLOTS = ['5', '6', '7']  # Slots for parts plates
    
    # Destination plate (on thermocycler)
    DESTINATION_PLATE_SLOT = '10'
    
    # Prefix and suffix plates
    PREFIX_SUFFIX_PLATE_SLOT = '2'
    
    # Thermocycler
    THERMOCYCLER_SLOT = '12'
    
    # Mix settings
    LINKER_MIX_SETTINGS = (3, 5)  # (repetitions, volume) for linker mixing
    PART_MIX_SETTINGS = (3, 5)     # (repetitions, volume) for part mixing
    
    def clip(clips_dict):
        ### Loading Tiprack
        tiprack = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', TIPRACK_SLOT, 'Tiprack')
        
        ### Loading Pipettes
        pipette = protocol.load_instrument('flex_1channel_1000', mount='left', tip_racks=[tiprack])
        
        ### Loading Thermocycler
        tc_mod = protocol.load_module('thermocycler module gen2', THERMOCYCLER_SLOT)
        
        ### Loading Labware
        # Master mix tubes
        master_mix_labware = protocol.load_labware('4ti0136_96_wellplate_2200ul', MASTER_MIX_SLOT, 'Master Mix Tubes')
        master_mix_volumes = [MASTER_MIX_VOLUME] * MASTER_MIX_TUBES
        
        # Water tubes
        water_labware = protocol.load_labware('4ti0136_96_wellplate_2200ul', WATER_SLOT, 'Water Tubes')
        
        # Parts plates
        parts_plates = []
        for i, slot in enumerate(PARTS_PLATES_SLOTS):
            parts_plates.append(protocol.load_labware('4ti0136_96_wellplate_2200ul', slot, f'Parts Plate {i+1}'))
        
        # Destination plate (on thermocycler)
        destination_plate = tc_mod.load_labware('4ti0136_96_wellplate_2200ul', 'Destination Plate')
        
        # Prefix and suffix plate
        prefix_suffix_plate = protocol.load_labware('4ti0136_96_wellplate_2200ul', PREFIX_SUFFIX_PLATE_SLOT, 'Prefix and Suffix Plate')
        
        ### Initialize Managers
        tip_manager = TipManager(protocol, TIPRACK_SLOT, pipette, TIPRACK_TYPE)
        master_mix_manager = MasterMixManager(protocol, master_mix_labware, master_mix_volumes, 
                                           MASTER_MIX_TRANSFER_VOLUME, pipette)
        
        ### Extract data from clips_dict
        prefixes_wells = clips_dict["prefixes_wells"]
        prefixes_plates = clips_dict["prefixes_plates"]
        suffixes_wells = clips_dict["suffixes_wells"]
        suffixes_plates = clips_dict["suffixes_plates"]
        parts_wells = clips_dict["parts_wells"]
        parts_plates = clips_dict["parts_plates"]
        parts_vols = clips_dict["parts_vols"]
        water_vols = clips_dict["water_vols"]
        
        ### Master Mix Distribution
        print("Distributing master mix...")
        destination_wells = [destination_plate.wells()[i] for i in range(len(prefixes_wells))]
        master_mix_manager.distribute_to_wells(destination_wells, pipette)
        
        ### Water Distribution
        print("Distributing water...")
        water_transfers = []
        for i, water_vol in enumerate(water_vols):
            if water_vol > 0:
                water_transfers.append({
                    'source': water_labware.wells()[0],  # Use first water tube
                    'destination': destination_plate.wells()[i],
                    'volume': water_vol
                })
        
        batch_transfer(water_transfers, pipette, tip_manager)
        
        ### Prefix Distribution
        print("Distributing prefixes...")
        prefix_transfers = []
        for i, (prefix_well, prefix_plate) in enumerate(zip(prefixes_wells, prefixes_plates)):
            plate_index = int(prefix_plate) - 1
            source_well = prefix_suffix_plate.wells_by_name()[prefix_well]
            prefix_transfers.append({
                'source': source_well,
                'destination': destination_plate.wells()[i],
                'volume': 1.0  # Standard prefix volume
            })
        
        batch_transfer(prefix_transfers, pipette, tip_manager, mix_after=LINKER_MIX_SETTINGS)
        
        ### Parts Distribution
        print("Distributing parts...")
        parts_transfers = []
        for i, (part_well, part_plate, part_vol) in enumerate(zip(parts_wells, parts_plates, parts_vols)):
            plate_index = int(part_plate) - 1
            source_well = parts_plates[plate_index].wells_by_name()[part_well]
            parts_transfers.append({
                'source': source_well,
                'destination': destination_plate.wells()[i],
                'volume': part_vol
            })
        
        batch_transfer(parts_transfers, pipette, tip_manager, mix_after=PART_MIX_SETTINGS)
        
        ### Suffix Distribution
        print("Distributing suffixes...")
        suffix_transfers = []
        for i, (suffix_well, suffix_plate) in enumerate(zip(suffixes_wells, suffixes_plates)):
            plate_index = int(suffix_plate) - 1
            source_well = prefix_suffix_plate.wells_by_name()[suffix_well]
            suffix_transfers.append({
                'source': source_well,
                'destination': destination_plate.wells()[i],
                'volume': 1.0  # Standard suffix volume
            })
        
        batch_transfer(suffix_transfers, pipette, tip_manager, mix_after=LINKER_MIX_SETTINGS)
        
        print("CLIP protocol completed successfully!")
    
    # Execute the CLIP protocol
    clip(clips_dict)
    
    ### PCR Reaction in Thermocycler

    # close lid and set lid temperature, PCR will not start until lid reaches 37C
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)

    # Runs 20 cycles of 37C for 2 minutes and 20C for 1 minute, then holds for 60C for 10 minutes
    profile = [
        {'temperature': 37, 'hold_time_minutes': 2},
        {'temperature': 20, 'hold_time_minutes': 1}]
    tc_mod.execute_profile(steps=profile, repetitions=20, block_max_volume=30)
    tc_mod.set_block_temperature(60, hold_time_minutes=10, block_max_volume=30)
    tc_mod.set_block_temperature(8, block_max_volume=30)
    tc_mod.set_lid_temperature(37)
    # tc_mod.open_lid()                                     # leave lid shut to prevent evaporation 
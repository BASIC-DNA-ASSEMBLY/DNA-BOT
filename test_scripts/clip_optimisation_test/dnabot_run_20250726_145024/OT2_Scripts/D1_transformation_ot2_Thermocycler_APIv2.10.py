from opentrons import protocol_api
import json
import time
from typing import Dict, List, Optional, Tuple, Union

metadata = {
    'apiLevel': '2.10',
    'protocolName': 'Transformation Protocol v2',
    'description': 'Enhanced transformation protocol using thermocycler module',
    'author': 'Liam Hallett'
}

# Load transformation data from JSON file
# This will be replaced by the parser with embedded JSON data
transformation_dict = {
    "transformation_number": 16,
    "source_wells": ["A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3"],
    "source_plates": ["1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"],
    "destination_wells": ["A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3"],
    "dna_volumes": [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    "cell_volumes": [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
    "soc_volumes": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
    "plating_volumes": [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]
}

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# opentrons_simulate.exe dnabot\template_ot2_scripts\transformation_template_TC_APIv2.10.py --custom-labware-path 'labware\Labware definitions'

# Example dictionary structure for transformation data:
# transformation_dict = {
#     'transformation_number': 24,
#     'source_wells': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', ...],
#     'source_plates': ['1', '1', '1', '1', '1', '1', ...],
#     'destination_wells': ['A1', 'A2', 'A3', 'B1', 'B2', 'B3', ...],
#     'dna_volumes': [3, 3, 3, 3, 3, 3, ...],
#     'cell_volumes': [30, 30, 30, 30, 30, 30, ...],
#     'soc_volumes': [100, 100, 100, 100, 100, 100, ...],
#     'plating_volumes': [100, 100, 100, 100, 100, 100, ...]
# }

class TipManager:
    """Manages pipette tips and tracks usage."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 slot: int,
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
    
    def add_tip_rack(self, slot: int, tip_type: Optional[str] = None) -> None:
        """Add an additional tip rack to the manager."""
        tip_type = tip_type or self.tip_type
        
        # Check pipette compatibility
        pipette_type = 'p300' if self.pipette.max_volume >= 300 else 'p20'
        if tip_type != pipette_type:
            raise ValueError(f"Tip type '{tip_type}' is incompatible with pipette type '{pipette_type}'")
            
        if tip_type == 'p300':
            self.tipracks.append(self.protocol.load_labware('opentrons_96_tiprack_300ul', slot))
        elif tip_type == 'p20':
            self.tipracks.append(self.protocol.load_labware('opentrons_96_tiprack_20ul', slot))
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
    
    def get_multi_tip(self, start_col: int = 0) -> None:
        """Pick up multiple tips columnwise for multichannel pipetting."""
        if not self.pipette.channels > 1:
            raise ValueError("Multichannel pipetting requires a multichannel pipette")
            
        if not self.normal_tips:  # If either array is empty, we need new tips
            if self.current_rack < len(self.tipracks) - 1:
                # Switch to next tip rack
                self.current_rack += 1
                self._initialise_tip_arrays()
            else:
                self._prompt_tip_replacement()
        
        # Find a column with enough consecutive tips
        for col in range(start_col, self.cols):
            # Get all tips in this column
            tips_in_col = [t for t in self.normal_tips if t // self.rows == col]
            
            # Only use columns that are full
            if len(tips_in_col) != self.rows:
                continue
                
            # Sort tips in the column
            tips_in_col.sort()
            
            # Check if we have a complete column of consecutive tips
            expected_tips = [col * self.rows + row for row in range(self.rows)]
            if tips_in_col == expected_tips:
                # Remove tips from both arrays
                for tip in tips_in_col:
                    self.normal_tips.remove(tip)
                    self.inverse_tips.remove(tip)
                
                self.pipette.pick_up_tip(self.tipracks[self.current_rack].wells()[tips_in_col[0]])
                return
        
        # If we get here, no suitable column was found in the current rack
        self.normal_tips = []  # Force the next call to switch racks
        self.get_multi_tip(start_col)
    
    def _prompt_tip_replacement(self) -> None:
        """Prompt user to replace tip rack and flash lights."""
        # Flash lights 3 times before the prompt
        for _ in range(3):
            self.protocol.set_rail_lights(False)
            time.sleep(0.15)
            self.protocol.set_rail_lights(True)
            time.sleep(0.15)
        
        self.protocol.pause("Please replace the tip rack")
        
        # Reset current rack and reinitialise tip arrays
        self.current_rack = 0
        self._initialise_tip_arrays()

class WellManager:
    """Manages well positions and column operations."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 transformation_number: int):
        self.protocol = protocol
        self.transformation_number = transformation_number
        self.rows = 8
        self.cols = 12
    
    def get_columns_with_wells(self) -> List[int]:
        """Get list of column indices that contain samples."""
        columns = set()
        for i in range(self.transformation_number):
            col = i % self.cols
            columns.add(col)
        return sorted(list(columns))
    
    def get_source_well(self, index: int, source_plate_a, source_plate_b, source_plates) -> str:
        """Get the source well for a given transformation index."""
        source_plate = source_plate_a if source_plates[index] == '1' else source_plate_b
        return source_plate.wells_by_name()[f'A{index + 1}']
    
    def get_destination_well(self, index: int, tc_plate) -> str:
        """Get the destination well for a given transformation index."""
        return tc_plate.wells_by_name()[f'A{index + 1}']

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
    Perform batch transfers similar to the Media Bot transfer function.
    
    Args:
        transfers: List of dicts with 'source', 'destination', 'volume' keys
        pipette: The pipette to use
        tip_manager: TipManager instance
        mix_after: Tuple of (repetitions, volume) for mixing after transfer, or None
        mix_speed: Rate for mixing (0-1, default 1.0 = 100% speed)
        dispense_speed: Rate for dispensing (0-1, default 0.4 = 40% speed)
    """
    if not transfers:
        return
    
    # Get maximum volume from pipette specifications
    MAX_VOLUME = pipette.max_volume  # µl
    excess_volume = 1.05  # 5% excess for blowout
    
    # Process transfers in batches
    current_batch = []
    current_volume = 0
    
    def process_batch(batch):
        """Helper function to process a batch of transfers."""
        if not batch:
            return
        
        # Pick up tip for this batch
        tip_manager.get_single_tip()
        
        # Calculate total volume needed for this batch
        total_volume = sum(t["volume"] for t in batch) * excess_volume
        
        # Use the first transfer's source for blowout
        source_well = batch[0]["source"]
        
        # Aspirate the total volume
        pipette.aspirate(total_volume, source_well)
        
        # Dispense to all wells in the batch
        for t in batch:
            pipette.dispense(t["volume"], t["destination"].top())
            pipette.touch_tip()
        
        # Blow out remaining volume back to source well
        remaining_volume = total_volume * (excess_volume - 1)
        pipette.blow_out(source_well)
        
        # Mix after if specified
        if mix_after:
            reps, mix_vol = mix_after
            for t in batch:
                pipette.mix(reps, mix_vol, t["destination"], rate=mix_speed)
        
        # Drop the tip
        pipette.drop_tip()
    
    # Process all transfers in appropriate batch sizes
    for transfer in transfers:
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

def run(protocol: protocol_api.ProtocolContext):
    """Main transformation protocol function."""
    
    # Extract data from embedded dictionary
    transformation_number = transformation_dict['transformation_number']
    source_wells = transformation_dict['source_wells']
    source_plates = transformation_dict['source_plates']
    destination_wells = transformation_dict['destination_wells']
    dna_volumes = transformation_dict['dna_volumes']
    cell_volumes = transformation_dict['cell_volumes']
    soc_volumes = transformation_dict['soc_volumes']
    plating_volumes = transformation_dict['plating_volumes']
    
    # Load hardware
    left_pipette = protocol.load_instrument('p300_multi_gen2', 'left')
    right_pipette = protocol.load_instrument('p20_single_gen2', 'right')
    
    # Initialize tip managers
    tip_manager_300 = TipManager(protocol, 6, left_pipette, 'p300')
    tip_manager_20 = TipManager(protocol, 9, right_pipette, 'p20')
    
    # Load plates and vessels
    source_plate_a = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 1)
    source_plate_b = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 2)
    plate_2 = protocol.load_labware('corning_12_wellplate_6.9ml_flat', 5)
    tube_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 4)
    
    # Load thermocycler module
    if thermocycler_gen == 'gen1':
        thermocycler = protocol.load_module('thermocycler', 7)
    else:  # gen2
        thermocycler = protocol.load_module('thermocyclerModuleV2', 7)
    TC_plate = thermocycler.load_labware('biorad_96_wellplate_200ul_pcr')
    
    # Initialize well manager
    well_manager = WellManager(protocol, transformation_number)
    
    # Helper functions
    def flash(n):
        """Flashes light n times to prompt user input."""
        for i in range(n):
            protocol.set_rail_lights(False)
            time.sleep(0.13)
            protocol.set_rail_lights(True)
            time.sleep(0.13)
    
    def tube_counter(n):
        """Returns tube position in tube rack."""
        row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F"}
        row = row_dict[n % 4]
        col = 1 + n // 4
        return row + str(col)
    
    def plate_counter(n):
        """Returns plate position for plating."""
        row_dict = {0: "A", 1: "B", 2: "C"}
        row = row_dict[n % 3]
        col = 1 + n // 3
        return row + str(col)
    
    # Calculate requirements
    cell_tube_vol = 1000  # volume of cells per aliquot
    cell_tube_capacity = cell_tube_vol // 30  # transfers per tube (assuming 30µl per transformation)
    cell_req_tubes = (transformation_number - 1) // cell_tube_capacity + 1
    
    SOC_tube_vol = 1500  # vol per tube of SOC
    SOC_tube_capacity = SOC_tube_vol // 100  # transfers per tube (assuming 100µl per transformation)
    req_SOC_tubes = (transformation_number - 1) // SOC_tube_capacity + 1
    
    plate_capacity = 12  # number of wells available on plate
    req_plates = (transformation_number - 1) // plate_capacity + 1
    
    # Protocol execution
    protocol.set_rail_lights(True)
    thermocycler.set_block_temperature(4)
    
    comment = f"Required:\n {cell_req_tubes} tube(s) of competent cells,\n{req_SOC_tubes} tube(s) of SOC,\n{req_plates} plate(s)"
    protocol.pause(comment)
    flash(5)
    
    # Dispense competent cells
    protocol.comment("Dispensing competent cells")
    cell_transfers = []
    for i in range(transformation_number):
        cell_tube_count = i // cell_tube_capacity
        cell_transfers.append({
            "source": tube_rack[tube_counter(cell_tube_count)],
            "destination": TC_plate[destination_wells[i]],
            "volume": cell_volumes[i]
        })
    
    # Perform batch transfer for cells
    batch_transfer(cell_transfers, left_pipette, tip_manager_300)
    
    # Transfer DNA
    protocol.comment("Transferring DNA")
    dna_transfers = []
    for i in range(transformation_number):
        # Determine source plate
        source_plate = source_plate_a if source_plates[i] == '1' else source_plate_b
        
        dna_transfers.append({
            "source": source_plate[source_wells[i]],
            "destination": TC_plate[destination_wells[i]],
            "volume": dna_volumes[i]
        })
    
    # Perform batch transfer for DNA
    batch_transfer(dna_transfers, right_pipette, tip_manager_20, mix_after=(3, 10))
    
    # Heat shock protocol
    protocol.comment("Starting heat shock protocol")
    thermocycler.close_lid()
    thermocycler.set_block_temperature(4, hold_time_minutes=30, block_max_volume=55)
    thermocycler.set_block_temperature(42, hold_time_seconds=30, block_max_volume=55)
    thermocycler.set_block_temperature(4, hold_time_seconds=90, block_max_volume=55)
    
    # Add SOC
    protocol.comment("Adding SOC")
    thermocycler.open_lid()
    thermocycler.set_block_temperature(37)
    
    soc_transfers = []
    for i in range(transformation_number):
        SOC_tube_count = 4 + (i // SOC_tube_capacity)  # starts at A2 (tube 4)
        soc_transfers.append({
            "source": tube_rack[tube_counter(SOC_tube_count)],
            "destination": TC_plate[destination_wells[i]],
            "volume": soc_volumes[i]
        })
    
    # Perform batch transfer for SOC
    batch_transfer(soc_transfers, left_pipette, tip_manager_300, mix_after=(2, 50))
    
    # Recovery incubation
    protocol.comment("Recovery incubation")
    thermocycler.close_lid()
    thermocycler.set_block_temperature(37, hold_time_minutes=60, block_max_volume=150)
    
    # Plating
    protocol.comment("Plating transformations")
    thermocycler.open_lid()
    thermocycler.set_block_temperature(37)
    
    plating_transfers = []
    for i in range(transformation_number):
        step = i % plate_capacity
        
        if step == 0:
            plate = str(i // plate_capacity + 1)
            comment = f"Place plate {plate} onto robot at deck position 5"
            protocol.pause(comment)
            flash(5)
        
        plating_transfers.append({
            "source": TC_plate[destination_wells[i]],
            "destination": plate_2[plate_counter(step)],
            "volume": plating_volumes[i]
        })
    
    # Perform batch transfer for plating
    batch_transfer(plating_transfers, left_pipette, tip_manager_300)
    
    # Finish protocol
    thermocycler.deactivate()
    protocol.set_rail_lights(False)
    protocol.comment("Transformation protocol complete")

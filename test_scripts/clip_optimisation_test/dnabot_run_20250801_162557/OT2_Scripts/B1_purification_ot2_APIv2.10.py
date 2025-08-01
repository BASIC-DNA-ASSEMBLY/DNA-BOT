from opentrons import protocol_api
import time
from typing import Dict, List, Optional, Tuple, Union

metadata = {
    'apiLevel': '2.10',
    'protocolName': 'DNABOT: B1 Purification v2.10',
    'description': 'Enhanced magnetic bead purification protocol for BASIC assembly using Opentrons OT-2',
    'author': 'Liam Hallett, Matthew Haines'
}

# Protocol Configuration
PROTOCOL_CONFIG = {
    # Sample Configuration
    'clips_number': 88,    # Number of CLIP reactions to process
    'sample_volume': 40,    # Volume of sample in µL
    'bead_ratio': 1.8,      # Ratio of beads to sample volume
    'elution_volume': 40,   # Volume of elution buffer in µL
    
    # Reagent Positions
    'ethanol_well': 'A11',   # Well position for ethanol in reservoir
    'elution_well': 'A10',   # Well position for elution buffer in reservoir
    
    # Timing Parameters
    'incubation_times': {
        'initial': 5,       # Initial incubation time in minutes
        'settling': 5,      # Bead settling time in minutes
        'drying': 5,        # Drying time in minutes
        'elution': 2,       # Elution time in minutes
        'wash': 0.5         # Wash time in minutes
    },
    
    # Hardware Settings
    'magdeck_height': 20,   # Height of magnetic module
    'pipette_rates': {
        'aspirate': 25,     # Aspirate rate in µL/s
        'dispense': 150     # Dispense rate in µL/s
    },
    
    # Protocol Constants
    'dead_volumes': {
        'total': 6,         # Total dead volume
        'ethanol': 50,      # Ethanol dead volume
        'elution': 2        # Elution dead volume
    },
    
    # Mixing Parameters
    'mix_reps': {
        'immobilise': 10,   # Number of mix repetitions for immobilisation
        'elution': 20       # Number of mix repetitions for elution
    }
}

# PLACEHOLDER_PROTOCOL_CONFIG - This will be replaced by the parser

class WellManager:
    """Manages well positions and column operations."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 sample_number: int):
        self.protocol = protocol
        self.sample_number = sample_number
        self.rows = 8
        self.cols = 12
        
        # Calculate number of columns needed
        self.src1_col_num = min((sample_number - 1) // 8 + 1, 6)  # Max 6 columns in first plate
        self.src2_col_num = max(0, (sample_number - 48 - 1) // 8 + 1)  # Remaining columns in second plate
    
    def get_columns_with_wells(self) -> List[int]:
        """Get list of column indices that contain samples."""
        columns = list(range(self.src1_col_num))
        if self.src2_col_num > 0:
            columns.extend(range(self.src1_col_num, self.src1_col_num + self.src2_col_num))
        return columns
    
    def get_source_well(self, col: int, source_plate, source_plate2=None) -> str:
        """Get the source well for a given column index."""
        if col < self.src1_col_num:
            # First plate (columns 0-5)
            return source_plate.wells_by_name()[f'A{col + 1}']
        else:
            # Second plate (columns 6-11)
            if not source_plate2:
                raise ValueError("Source plate 2 is required for columns >= 6")
            return source_plate2.wells_by_name()[f'A{col - self.src1_col_num + 1}']
    
    def get_mag_well(self, col: int, mag_plate) -> str:
        """Get the magnetic plate well for a given column index."""
        return mag_plate.wells_by_name()[f'A{col + 1}']
    
    def get_final_well(self, col: int, final_plate) -> str:
        """Get the final plate well for a given column index."""
        return final_plate.wells_by_name()[f'A{col + 1}']

class TipManager:
    """Manages pipette tips and tracks usage."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 slot: int,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 tip_type: str):
        self.protocol = protocol
        self.pipette = pipette
        self.tip_type = tip_type
        self.rows = 8
        self.cols = 12
        
        # Initialise tip racks array and current rack index
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
    
    def return_tips(self, col: int) -> None:
        """Return tips to their original wells and update tracking arrays."""
        if not self.pipette.has_tip:
            return
            
        # Calculate the tip indices for this column
        tip_indices = [col * self.rows + row for row in range(self.rows)]
        
        # Return tips to their original wells
        self.pipette.return_tip()
        
        # Add the tips back to both arrays
        for tip in tip_indices:
            if tip not in self.normal_tips:
                self.normal_tips.append(tip)
            if tip not in self.inverse_tips:
                self.inverse_tips.append(tip)
        
        # Sort the arrays to maintain order
        self.normal_tips.sort()
        self.inverse_tips.sort()
    
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

def calculate_transfer_time_duration(num_columns: int) -> float:
    """
    Calculate transfer time duration based on number of columns.
    Excludes first column, uses 1/8 minute per additional column.
    
    Args:
        num_columns: Number of columns being processed
        
    Returns:
        Transfer time duration in minutes
    """
    if num_columns <= 1:
        return 0.0
    
    # Exclude first column, calculate time for remaining columns
    additional_columns = num_columns - 1
    return additional_columns * (1/8)  # 1/8 minute per column

def run(protocol: protocol_api.ProtocolContext):
    """Main protocol function."""
    
    # Load hardware
    multi_pipette = protocol.load_instrument('p300_multi_gen2', mount='left')
    multi_pipette.flow_rate.aspirate = PROTOCOL_CONFIG['pipette_rates']['aspirate']
    multi_pipette.flow_rate.dispense = PROTOCOL_CONFIG['pipette_rates']['dispense']
    
    # Load labware
    source_plate = protocol.load_labware('4ti0960rig_96_wellplate_200ul', 1)
    mag_mod = protocol.load_module('magneticModuleV2', 4)
    mag_plate = mag_mod.load_labware('4ti0960rig_96_wellplate_200ul')
    reservoir = protocol.load_labware('4ti0131_12_reservoir_21000ul', 5)
    
    # Load second source plate if needed
    # NB: clips are generated in sets of up to 48 in the clip reaction but the purification protocol can handle up to 96 (2 plates)
    if PROTOCOL_CONFIG['clips_number'] > 48:
        source_plate2 = protocol.load_labware('4ti0960rig_96_wellplate_200ul', 2)
    else:
        source_plate2 = None
    
    # Initialise managers
    well_manager = WellManager(protocol, PROTOCOL_CONFIG['clips_number'])
    
    # Initialise separate tip managers for different stages
    sample_tip_manager = TipManager(protocol, 3, multi_pipette, 'p300')
    wash_tip_manager = TipManager(protocol, 6, multi_pipette, 'p300')
    elution_tip_manager = TipManager(protocol, 9, multi_pipette, 'p300')
    
    # Get columns with wells
    columns_with_wells = well_manager.get_columns_with_wells()
    
    # Calculate volumes
    bead_volume = PROTOCOL_CONFIG['sample_volume'] * PROTOCOL_CONFIG['bead_ratio']
    total_volume = bead_volume + PROTOCOL_CONFIG['sample_volume'] + PROTOCOL_CONFIG['dead_volumes']['total']
    
    # Transfer samples
    protocol.comment("Transferring samples")
    for col in columns_with_wells:
        sample_tip_manager.get_multi_tip(start_col=col)
        source_well = well_manager.get_source_well(col, source_plate, source_plate2)
        mag_well = well_manager.get_mag_well(col, mag_plate)
        
        multi_pipette.transfer(
            PROTOCOL_CONFIG['sample_volume'],
            source_well,
            mag_well,
            mix_after=(PROTOCOL_CONFIG['mix_reps']['immobilise'], total_volume/2),
            new_tip='never'
        )
        sample_tip_manager.return_tips(col)
    
    # Initial incubation
    protocol.comment("Initial incubation")
    transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
    adjusted_delay = max(0, PROTOCOL_CONFIG['incubation_times']['initial'] - transfer_time)
    if adjusted_delay > 0:
        protocol.delay(minutes=adjusted_delay)
    
    # Prompt user to replace source plate 2 with final plate
    if PROTOCOL_CONFIG['clips_number'] > 48:
        protocol.pause("Please remove source plate 2 from slot 2 and add the final plate in its place")
        final_plate = source_plate2
    else:
        protocol.pause("Please add the final plate to slot 2")
        final_plate = protocol.load_labware('4ti0960rig_96_wellplate_200ul', 2)
    
    # Engage magnet and wait
    protocol.comment("Engaging magnet")
    mag_mod.engage(height=PROTOCOL_CONFIG['magdeck_height'])
    transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
    adjusted_delay = max(0, PROTOCOL_CONFIG['incubation_times']['settling'] - transfer_time)
    if adjusted_delay > 0:
        protocol.delay(minutes=adjusted_delay)
    
    # Remove supernatant
    protocol.comment("Removing supernatant")
    for col in columns_with_wells:
        sample_tip_manager.get_multi_tip(start_col=col)
        mag_well = well_manager.get_mag_well(col, mag_plate)
        
        multi_pipette.aspirate(total_volume, mag_well)
        multi_pipette.move_to(reservoir['A1'].top(z=70))
        multi_pipette.blow_out(protocol.fixed_trash["A1"].top())
        sample_tip_manager.return_tips(col)
    
    # Ethanol washes
    for wash in range(2):
        protocol.comment(f"Ethanol wash {wash + 1}")
        
        # Add ethanol
        for col in columns_with_wells:
            wash_tip_manager.get_multi_tip(start_col=col)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            
            multi_pipette.transfer(
                150,  # ETHANOL_VOL
                reservoir[PROTOCOL_CONFIG['ethanol_well']],
                mag_well,
                new_tip='never'
            )
            wash_tip_manager.return_tips(col)
        
        # Calculate dynamic delay based on transfer time
        transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
        adjusted_delay = max(0, PROTOCOL_CONFIG['incubation_times']['wash'] - transfer_time)
        if adjusted_delay > 0:
            protocol.delay(minutes=adjusted_delay)
        
        # Remove ethanol
        for col in columns_with_wells:
            wash_tip_manager.get_multi_tip(start_col=col)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            
            multi_pipette.aspirate(200, mag_well)  # ETHANOL_VOL + ETHANOL_DEAD_VOL
            multi_pipette.move_to(reservoir['A1'].top(z=70))
            multi_pipette.blow_out(protocol.fixed_trash["A1"].top())
            wash_tip_manager.return_tips(col)
    
    # Drying time
    protocol.comment("Drying time")
    transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
    adjusted_delay = max(0, PROTOCOL_CONFIG['incubation_times']['drying'] - transfer_time)
    if adjusted_delay > 0:
        protocol.delay(minutes=adjusted_delay)
    
    # Disengage magnet
    mag_mod.disengage()
    
    # Add elution buffer
    protocol.comment("Adding elution buffer")
    for col in columns_with_wells:
        elution_tip_manager.get_multi_tip(start_col=col)
        mag_well = well_manager.get_mag_well(col, mag_plate)
        
        multi_pipette.transfer(
            PROTOCOL_CONFIG['elution_volume'],
            reservoir[PROTOCOL_CONFIG['elution_well']],
            mag_well,
            mix_after=(PROTOCOL_CONFIG['mix_reps']['elution'], PROTOCOL_CONFIG['elution_volume']/2),
            new_tip='never'
        )
        elution_tip_manager.return_tips(col)
    
    # Elution incubation
    protocol.comment("Elution incubation")
    transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
    adjusted_delay = max(0, PROTOCOL_CONFIG['incubation_times']['elution'] - transfer_time)
    if adjusted_delay > 0:
        protocol.delay(minutes=adjusted_delay)
    
    # Engage magnet
    mag_mod.engage(height=PROTOCOL_CONFIG['magdeck_height'])
    protocol.delay(minutes=1)  # ELUTANT_SEP_TIME
    
    # Transfer to final plate
    protocol.comment("Transferring to final plate")
    for col in columns_with_wells:
        elution_tip_manager.get_multi_tip(start_col=col)
        mag_well = well_manager.get_mag_well(col, mag_plate)
        final_well = well_manager.get_final_well(col, final_plate)
        
        multi_pipette.transfer(
            PROTOCOL_CONFIG['elution_volume'] - PROTOCOL_CONFIG['dead_volumes']['elution'],
            mag_well,
            final_well,
            new_tip='never'
        )
        elution_tip_manager.return_tips(col)
    
    # Disengage magnet
    mag_mod.disengage()
    
    protocol.comment("Protocol complete") 
from opentrons import protocol_api
import time
from typing import Dict, List, Optional, Tuple, Union

metadata = {
    'protocolName': 'OT Magnetic Bead Purification',
    'author': 'Liam Hallett',
    'description': 'Magnetic bead purification protocol for Opentrons OT-2',
    'apiLevel': '2.10'
}

# Protocol Configuration
PROTOCOL_CONFIG = {
    'wells': ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1', 'F3'],  # Example wells
    'sample_volume': 20,    # Example sample volume in µL
    'magbead_ratio': 0.7,   # Ratio of magbeads to sample volume
    'ethanol_volume': 190,  # Volume of ethanol to add (µL)
    'water_volume': 25,     # Volume of water to add (µL)
    'ethanol_well': 'A2',   # Well position for ethanol in reservoir
    'water_well': 'A11',    # Well position for water in reservoir
    'incubation_times': {
        'initial': 300,     # 5 minutes initial incubation
        'mag_bead_separation': 300,  # 5 minutes magnetic bead separation
        'ethanol': 30,      # 30 seconds ethanol incubation
        'drying_time': 60,  # 1 minute drying time
        'water': 300        # 5 minutes water incubation
    }
}

class WellManager:
    """Manages well positions and column operations."""
    
    def __init__(self, wells: List[str]):
        self.wells = wells
        self.rows = 8
        self.cols = 12
    
    def get_columns_with_wells(self) -> List[int]:
        """Get list of columns that contain at least one well from the input list."""
        columns = set()
        for well in self.wells:
            col = int(well[1:]) - 1  # Convert 'A1' to column index 0
            columns.add(col)
        return sorted(list(columns))
    
    def get_wells_in_column(self, col: int) -> List[str]:
        """Get all wells from the input list that are in the specified column."""
        return [well for well in self.wells if int(well[1:]) - 1 == col]

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
    
    def get_multi_tip(self, start_col: int = 0) -> None:
        """Pick up multiple tips for multichannel pipetting."""
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

def run(protocol: protocol_api.ProtocolContext):
    # Load hardware
    multi_pipette = protocol.load_instrument('p300_multi_gen2', mount='left')
    
    # Load labware
    sample_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 1)
    final_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 2)
    reservoir = protocol.load_labware('usascientific_12_reservoir_22ml', 5)  # Reservoir for ethanol and water
    
    # Load magnetic module and plate
    mag_mod = protocol.load_module('magneticModuleV2', 4)
    mag_plate = mag_mod.load_labware('biorad_96_wellplate_200ul_pcr')
    mag_mod.disengage()
    
    # Initialise tip managers
    initial_tip_manager = TipManager(protocol, 9, multi_pipette, 'p300')
    etoh_tip_manager = TipManager(protocol, 6, multi_pipette, 'p300')
    water_tip_manager = TipManager(protocol, 3, multi_pipette, 'p300')
    
    # Calculate magbead volume
    magbead_volume = PROTOCOL_CONFIG['magbead_ratio'] * PROTOCOL_CONFIG['sample_volume']
    
    # Initialise well manager
    well_manager = WellManager(PROTOCOL_CONFIG['wells'])
    
    # Get columns with wells once
    columns_with_wells = well_manager.get_columns_with_wells()
    
    # Transfer samples
    protocol.comment("Transferring samples")
    for col in columns_with_wells:
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        initial_tip_manager.get_multi_tip(start_col=col)
        multi_pipette.aspirate(PROTOCOL_CONFIG['sample_volume'], sample_plate.wells_by_name()[f'A{col+1}'])
        multi_pipette.dispense(PROTOCOL_CONFIG['sample_volume'], mag_plate.wells_by_name()[f'A{col+1}'])
        multi_pipette.mix(3, PROTOCOL_CONFIG['sample_volume'])
        multi_pipette.move_to(reservoir['A1'].top(z=70))
        initial_tip_manager.return_tips(col)
    
    # Initial incubation
    protocol.comment("Initial incubation")
    protocol.delay(seconds=PROTOCOL_CONFIG['incubation_times']['initial'])
    
    # Engage magnet and remove supernatant
    protocol.comment("Engaging magnet and removing supernatant")
    mag_mod.engage()
    protocol.delay(seconds=PROTOCOL_CONFIG['incubation_times']['mag_bead_separation'])
    
    for col in columns_with_wells:
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        initial_tip_manager.get_multi_tip(start_col=col)
        multi_pipette.aspirate(PROTOCOL_CONFIG['sample_volume'] + magbead_volume, mag_plate.wells_by_name()[f'A{col+1}'])
        multi_pipette.move_to(reservoir['A1'].top(z=70))
        multi_pipette.blow_out(protocol.fixed_trash["A1"].top())
        initial_tip_manager.return_tips(col)
    
    # Ethanol washes (2x)
    for wash_num in range(2):
        protocol.comment(f"Ethanol wash {wash_num + 1}")
        for col in columns_with_wells:
            wells_in_col = well_manager.get_wells_in_column(col)
            if not wells_in_col:
                continue
                
            # Add ethanol using ethanol tips
            etoh_tip_manager.get_multi_tip(start_col=col)
            multi_pipette.aspirate(PROTOCOL_CONFIG['ethanol_volume'], reservoir[PROTOCOL_CONFIG['ethanol_well']])
            multi_pipette.dispense(PROTOCOL_CONFIG['ethanol_volume'], mag_plate.wells_by_name()[f'A{col+1}'].top())
            multi_pipette.blow_out(mag_plate.wells_by_name()[f'A{col+1}'].top())
            etoh_tip_manager.return_tips(col)
        
        # Only wait if we have more than 6 columns
        if len(columns_with_wells) < 4: # time taken for 4 columns is > 30 seconds
            protocol.delay(seconds=PROTOCOL_CONFIG['incubation_times']['ethanol'])
        
        # Remove ethanol using initial tips
        for col in columns_with_wells:
            wells_in_col = well_manager.get_wells_in_column(col)
            if not wells_in_col:
                continue
                
            initial_tip_manager.get_multi_tip(start_col=col)
            multi_pipette.aspirate(200, mag_plate.wells_by_name()[f'A{col+1}'])
            multi_pipette.move_to(reservoir['A1'].top(z=70))
            multi_pipette.blow_out(protocol.fixed_trash["A1"].top())
            initial_tip_manager.return_tips(col)
    
    # Drying time
    protocol.comment("Drying time")
    protocol.delay(seconds=PROTOCOL_CONFIG['incubation_times']['drying_time'])
    
    # Disengage magnet and add water
    protocol.comment("Disengaging magnet and adding water")
    mag_mod.disengage()
    
    for col in columns_with_wells:
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        water_tip_manager.get_multi_tip(start_col=col)
        multi_pipette.aspirate(PROTOCOL_CONFIG['water_volume'], reservoir[PROTOCOL_CONFIG['water_well']])
        multi_pipette.dispense(PROTOCOL_CONFIG['water_volume'], mag_plate.wells_by_name()[f'A{col+1}'])
        multi_pipette.mix(3, PROTOCOL_CONFIG['water_volume'])
        water_tip_manager.return_tips(col)
    
    # Final incubation
    protocol.comment("Final water incubation")
    protocol.delay(seconds=PROTOCOL_CONFIG['incubation_times']['water'])
    
    # Re-engage magnet
    mag_mod.engage()
    
    # Wait 5 minutes after magnet engagement
    protocol.comment("Waiting 5 minutes after magnet engagement")
    protocol.delay(seconds=300)
    
    # Transfer purified samples to final plate
    protocol.comment("Transferring purified samples to final plate")
    for col in columns_with_wells:
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        water_tip_manager.get_multi_tip(start_col=col)
        multi_pipette.aspirate(PROTOCOL_CONFIG['water_volume'], mag_plate.wells_by_name()[f'A{col+1}'])
        multi_pipette.dispense(PROTOCOL_CONFIG['water_volume'], final_plate.wells_by_name()[f'A{col+1}'])
        water_tip_manager.return_tips(col)
    
    mag_mod.disengage()
    
    protocol.comment("Protocol complete")

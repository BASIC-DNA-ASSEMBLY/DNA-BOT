from opentrons import protocol_api
import json
import time
from typing import List, Optional

# metadata
metadata = {
'protocolName': 'DNABOT: Transformation and Spotting v2.15 (Flex)',
'description': 'Transformation and spotting protocol for DNA-BOT using Opentrons Flex with thermocycler module'
}
requirements = {"robotType": "Flex", "apiLevel": "2.15"}

# Load transformation data from JSON file
# This will be replaced by the parser with embedded JSON data
with open('transformation_data.json') as f:
    transformation_dict = json.load(f)

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'GEN2'

# opentrons_simulate.exe dnabot\template_flex_scripts\transformation_template_TC_APIv2.10.py --custom-labware-path 'labware\Labware definitions'

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
        """Add an additional tip rack to the manager.
        
        Args:
            slot: The deck slot number for the tip rack
            tip_type: Optional tip type ('flex_1channel_1000' or 'flex_8channel_1000'). If None, uses the manager's default tip type.
            
        Raises:
            ValueError: If the tip type is not supported or incompatible with the pipette
        """
        tip_type = tip_type or self.tip_type
        
        # Check pipette compatibility
        pipette_type = 'flex_8channel_1000' if self.pipette.channels > 1 else 'flex_1channel_1000'
        if tip_type != pipette_type:
            raise ValueError(f"Tip type '{tip_type}' is incompatible with pipette type '{pipette_type}'")
            
        if tip_type == 'flex_8channel_1000':
            self.tipracks.append(self.protocol.load_labware('opentrons_96_tiprack_1000ul', slot))
        elif tip_type == 'flex_1channel_1000':
            self.tipracks.append(self.protocol.load_labware('opentrons_96_tiprack_1000ul', slot))
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
                self.protocol.comment(f"Switched to tip rack {self.current_rack + 1}")
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
        """Pick up multiple tips for multichannel pipetting."""
        self.protocol.comment("Attempting to get multi-channel tips")
        if not self.pipette.channels > 1:
            raise ValueError("Multichannel pipetting requires a multichannel pipette")
            
        if not self.normal_tips:  # If either array is empty, we need new tips
            self.protocol.comment(f"No tips available in current rack: {self.current_rack}")
            if self.current_rack < len(self.tipracks) - 1:
                # Switch to next tip rack
                self.current_rack += 1
                self._initialise_tip_arrays()
                self.protocol.comment(f"Switched to tip rack {self.current_rack + 1}")
            else:
                self._prompt_tip_replacement()
        
        # Find a column with enough consecutive tips
        for col in range(start_col, self.cols):
            # Get all tips in this column
            tips_in_col = [t for t in self.normal_tips if t // self.rows == col]
            self.protocol.comment(f"Tips in column {col + 1}: {tips_in_col}")
            
            # Only use columns that are full
            if len(tips_in_col) != self.rows:
                self.protocol.comment(f"Column {col + 1} not full, only {len(tips_in_col)} tips available")
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
                
                self.protocol.comment(f"Picking up tips from column {col + 1} at {self.tipracks[self.current_rack].wells()[tips_in_col[0]]}")
                self.pipette.pick_up_tip(self.tipracks[self.current_rack].wells()[tips_in_col[0]])
                self.protocol.comment("Successfully picked up multi-channel tips")
                return
        
        # If we get here, no suitable column was found in the current rack
        self.protocol.comment(f"No full columns found in rack {self.current_rack + 1}, switching racks")
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

def run(protocol: protocol_api.ProtocolContext):
    def flash(n):
        """Flash lights n times"""
        for _ in range(n):
            protocol.set_rail_lights(False)
            time.sleep(0.15)
            protocol.set_rail_lights(True)
            time.sleep(0.15)

    def iterate_tip_count(tip_count, pipette):
        """Iterate tip count and get new tip if needed"""
        tip_count += 1
        if tip_count % 96 == 0:  # Every 96 tips, get a new tip rack
            pipette.pick_up_tip()
        return tip_count

    def counter(rows):
        """Counter function for well positions"""
        def inner(n):
            row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}
            if type(n) == int:
                row = row_dict[n // rows]
                col = 1 + n % rows
                return row + f'{col}'
        return inner

    def counter_T(cols):
        """Counter function for column positions"""
        def inner(n):
            col_dict = {0: "1", 1: "2", 2: "3", 3: "4", 4: "5", 5: "6", 6: "7", 7: "8", 8: "9", 9: "10", 10: "11", 11: "12"}
            if type(n) == int:
                col = col_dict[n % cols]
                row = 1 + n // cols
                return f'{row}' + col
        return inner

    def multi2single_tip_counter(n):
        """Convert multi-channel tip count to single-channel tip count"""
        return n * 8

    # Constants
    SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
    DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
    TIPRACK_TYPE = 'opentrons_96_tiprack_1000ul'
    
    # Thermocycler Module
    if thermocycler_gen == 'GEN1':
        tc_mod = protocol.load_module('thermocycler', '7')
    else:  # GEN2
        tc_mod = protocol.load_module('thermocyclerModuleV2', '7')
    
    # Load labware
    source_plate = tc_mod.load_labware(SOURCE_PLATE_TYPE)
    destination_plate = protocol.load_labware(DESTINATION_PLATE_TYPE, '8')
    
    # Load pipettes
    pipette = protocol.load_instrument('flex_1channel_1000', 'right')
    
    # Initialize TipManager for flex_1channel_1000 tips
    # Use available slots for flex_1channel_1000 tip racks (excluding slots 1, 2 for source plates and slots 7, 8, 10, 11 for thermocycler)
    flex_1channel_tip_slots = ['3', '6', '9']  # Available slots for flex_1channel_1000 tip racks
    flex_1channel_tip_manager = TipManager(protocol, int(flex_1channel_tip_slots[0]), pipette, 'flex_1channel_1000')
    
    # Add additional tip racks if needed
    for slot in flex_1channel_tip_slots[1:]:
        flex_1channel_tip_manager.add_tip_rack(int(slot), 'flex_1channel_1000')
    
    # Set up thermocycler
    tc_mod.open_lid()
    tc_mod.set_block_temperature(4)
    
    # Get transformation data
    source_wells = transformation_dict['source_wells']
    destination_wells = transformation_dict['destination_wells']
    volumes = transformation_dict['volumes']
    
    # Transfer samples from thermocycler to destination plate
    for i, (source_well, dest_well, volume) in enumerate(zip(source_wells, destination_wells, volumes)):
        flex_1channel_tip_manager.get_single_tip()
        pipette.transfer(volume, source_plate.wells(source_well), destination_plate.wells(dest_well), new_tip='never')
        pipette.drop_tip()
    
    # Close thermocycler lid
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(37)
    
    # Flash lights to indicate completion
    flash(3) 
from __future__ import unicode_literals
from opentrons import protocol_api
import numpy as np
import json
import time
from typing import List, Optional

# metadata
metadata = {
'protocolName': 'DNABOT: C1 Final Assembly with Thermocycler v2.10',
'description': 'Final assembly protocol for DNA-BOT using Opentrons OT-2 with thermocycler module',
'apiLevel': '2.10'
}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["A1", "A2", "B2"], [1, 1, 1]],
    "B1": [["A1", "C2", "D2"], [1, 1, 1]],
    "C1": [["A1", "E2", "F2"], [1, 1, 1]],
    "D1": [["A1", "G2", "H2"], [1, 1, 1]],
    "E1": [["A1", "A3", "B3"], [1, 1, 1]],
    "F1": [["A1", "C3", "D3"], [1, 1, 1]],
    "G1": [["A1", "E3", "F3"], [1, 1, 1]],
    "H1": [["A1", "G3", "H3"], [1, 1, 1]],
    "A2": [["A1", "A4", "B4"], [1, 1, 1]],
    "B2": [["A1", "C4", "D4"], [1, 1, 1]],
    "C2": [["A1", "E4", "F4"], [1, 1, 1]],
    "D2": [["A1", "G4", "H4"], [1, 1, 1]],
    "E2": [["A1", "A5", "D4"], [1, 1, 1]],
    "F2": [["B1", "B5", "D3"], [1, 1, 1]],
    "G2": [["B1", "E3", "F3"], [1, 1, 1]],
    "H2": [["B1", "G3", "H3"], [1, 1, 1]],
    "A3": [["B1", "A4", "C5"], [1, 1, 1]],
    "B3": [["B1", "D5", "E5"], [1, 1, 1]],
    "C3": [["B1", "F5", "G5"], [1, 1, 1]],
    "D3": [["B1", "H5", "A6"], [1, 1, 1]],
    "E3": [["B1", "B6", "C6"], [1, 1, 1]],
    "F3": [["B1", "D6", "E6"], [1, 1, 1]],
    "G3": [["B1", "F6", "F2"], [1, 1, 1]],
    "H3": [["B1", "G6", "H6"], [1, 1, 1]],
    "A4": [["B1", "A7", "B7"], [1, 1, 1]],
    "B4": [["B1", "C7", "D7"], [1, 1, 1]],
    "C4": [["C1", "E7", "F7"], [1, 1, 1]],
    "D4": [["C1", "G7", "H7"], [1, 1, 1]],
    "E4": [["C1", "A8", "B8"], [1, 1, 1]],
    "F4": [["C1", "C4", "G5"], [1, 1, 1]],
    "G4": [["C1", "E4", "F4"], [1, 1, 1]],
    "H4": [["C1", "E2", "F2"], [1, 1, 1]],
    "A5": [["C1", "D6", "H2"], [1, 1, 1]],
    "B5": [["C1", "C8", "D8"], [1, 1, 1]],
    "C5": [["C1", "A7", "E8"], [1, 1, 1]],
    "D5": [["C1", "F8", "D8"], [1, 1, 1]],
    "E5": [["C1", "G8", "H2"], [1, 1, 1]],
    "F5": [["C1", "G8", "H6"], [1, 1, 1]],
    "G5": [["C1", "A3", "H8"], [1, 1, 1]],
    "H5": [["D1", "C3", "D3"], [1, 1, 1]],
    "A6": [["D1", "E3", "F3"], [1, 1, 1]],
    "B6": [["D1", "G3", "A9"], [1, 1, 1]],
    "C6": [["D1", "B9", "B4"], [1, 1, 1]],
    "D6": [["D1", "C9", "D4"], [1, 1, 1]],
    "E6": [["D1", "E4", "F4"], [1, 1, 1]],
    "F6": [["D1", "D9", "E9"], [1, 1, 1]],
    "G6": [["D1", "G6", "F9"], [1, 1, 1]],
    "H6": [["D1", "G9", "H9"], [1, 1, 1]],
    "A7": [["D1", "C7", "A10"], [1, 1, 1]],
    "B7": [["D1", "G8", "F7"], [1, 1, 1]],
    "C7": [["D1", "G7", "H7"], [1, 1, 1]],
    "D7": [["D1", "A8", "B10"], [1, 1, 1]],
    "E7": [["E1", "C4", "G5"], [1, 1, 1]],
    "F7": [["E1", "E4", "F4"], [1, 1, 1]],
    "G7": [["E1", "E2", "F2"], [1, 1, 1]],
    "H7": [["E1", "C10", "H2"], [1, 1, 1]],
    "A8": [["E1", "G6", "D8"], [1, 1, 1]],
    "B8": [["E1", "D10", "D7"], [1, 1, 1]],
    "C8": [["E1", "E10", "F10"], [1, 1, 1]],
    "D8": [["E1", "G8", "G10"], [1, 1, 1]],
    "E8": [["E1", "G8", "H10"], [1, 1, 1]],
    "F8": [["E1", "A3", "A11"], [1, 1, 1]],
    "G8": [["E1", "C3", "D3"], [1, 1, 1]],
    "H8": [["E1", "F6", "F3"], [1, 1, 1]],
    "A9": [["E1", "G3", "A9"], [1, 1, 1]],
    "B9": [["F1", "B9", "B11"], [1, 1, 1]],
    "C9": [["F1", "C4", "D4"], [1, 1, 1]],
    "D9": [["F1", "E4", "F4"], [1, 1, 1]],
    "E9": [["F1", "D6", "F2"], [1, 1, 1]],
    "F9": [["F1", "G6", "D7"], [1, 1, 1]],
    "G9": [["F1", "A7", "D7"], [1, 1, 1]],
    "H9": [["F1", "C7", "D7"], [1, 1, 1]],
    "A10": [["F1", "G8", "F7"], [1, 1, 1]],
    "B10": [["F1", "G7", "H7"], [1, 1, 1]],
    "C10": [["F1", "A8", "B7"], [1, 1, 1]],
    "D10": [["F1", "C4", "G5"], [1, 1, 1]],
    "E10": [["F1", "E4", "F4"], [1, 1, 1]],
    "F10": [["F1", "E2", "F2"], [1, 1, 1]],
    "G10": [["G1", "C11", "D11"], [1, 1, 1]],
    "H10": [["G1", "A8", "B10"], [1, 1, 1]],
    "A11": [["G1", "C4", "G5"], [1, 1, 1]],
    "B11": [["G1", "E11", "F4"], [1, 1, 1]],
    "C11": [["G1", "E2", "F2"], [1, 1, 1]],
    "D11": [["G1", "H5", "H2"], [1, 1, 1]],
    "E11": [["G1", "G6", "D8"], [1, 1, 1]],
    "F11": [["G1", "A7", "D7"], [1, 1, 1]],
    "G11": [["G1", "C7", "D8"], [1, 1, 1]],
    "H11": [["G1", "G8", "H2"], [1, 1, 1]],
    "A12": [["G1", "G8", "F11"], [1, 1, 1]],
    "B12": [["G1", "A3", "H8"], [1, 1, 1]],
    "C12": [["G1", "C3", "D3"], [1, 1, 1]],
    "D12": [["H1", "E3", "F3"], [1, 1, 1]],
    "E12": [["H1", "G11", "A9"], [1, 1, 1]],
    "F12": [["H1", "B9", "B4"], [1, 1, 1]],
    "G12": [["H1", "C4", "H11"], [1, 1, 1]],
    "H12": [["H1", "E4", "F4"], [1, 1, 1]]
}
tiprack_num = 4

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'GEN2'

# It is possible to run 88 assemblies with this new module. The heat block module is removed. 
# Assembly reactions is set up on thermocycler module.

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
        
        # Initialise tip racks array and current rack index
        self.tipracks = []
        self.current_rack = 0
        
        # Load first tip rack using add_tip_rack
        self.add_tip_rack(slot)
        
        self._initialise_tip_arrays()
    
    def add_tip_rack(self, slot: int, tip_type: Optional[str] = None) -> None:
        """Add an additional tip rack to the manager.
        
        Args:
            slot: The deck slot number for the tip rack
            tip_type: Optional tip type ('p300' or 'p20'). If None, uses the manager's default tip type.
            
        Raises:
            ValueError: If the tip type is not supported or incompatible with the pipette
        """
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

class MasterMixManager:
    """Manages master mix tubes and tracks their volumes.
    
    Args:
        protocol: The protocol context
        labware: The labware containing the master mix tubes
        tube_volumes: List of volumes in each tube (µL)
        transfer_volume: Volume to transfer to each well (µL). This is the amount of master mix dispensed into each destination well during distribution. It is used to calculate how many wells can be filled per aspiration, to track tube depletion, and to ensure the correct amount is dispensed to each well.
        pipette: The pipette to use for transfers
        dead_volume: Minimum volume to leave in each tube (µL)
        tube_position: List of well names (e.g., ["A2"]). Must be user-defined and the same length as tube_volumes.
    """
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 labware: protocol_api.labware.Labware,
                 tube_volumes: list,
                 transfer_volume: float,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 dead_volume: float = 15.0,
                 tube_position: list = None):
        self.protocol = protocol
        self.labware = labware
        self.tube_volumes = tube_volumes.copy()  # Make a copy to avoid modifying the original
        self.transfer_volume = transfer_volume
        self.current_tube = 0
        self.max_volume = pipette.max_volume
        self.dead_volume = dead_volume
        if tube_position is None:
            raise ValueError("tube_position must be provided as a list of well names.")
        if len(tube_position) != len(tube_volumes):
            raise ValueError("tube_position and tube_volumes must have the same length.")
        self.tube_position = tube_position
        
        # Get the deck slot from the labware
        deck_slot = labware.parent
        
        # Print tube locations and volumes
        mm_locations = [f"Tube {i+1}: {tube_position[i]} - {tube_volumes[i]}µL (per well: {transfer_volume}µL)" 
                       for i in range(len(tube_volumes))]
        protocol.comment("Master Mix Tube Locations on slot " + str(deck_slot) + ": " + "; ".join(mm_locations))
    
    def get_current_tube(self) -> protocol_api.labware.Well:
        """Get the current tube's well object."""
        return self.labware.wells()[self.current_tube]
    
    def can_aspirate(self, volume: float) -> bool:
        """Check if we can aspirate the specified volume."""
        if self.current_tube >= len(self.tube_volumes):
            return False
        return self.tube_volumes[self.current_tube] >= volume + self.dead_volume
    
    def use_volume(self, volume: float) -> None:
        """Use a specified volume from the current tube."""
        self.tube_volumes[self.current_tube] -= volume
        
        # If current tube is nearly empty (less than transfer volume + dead volume), switch to next tube
        if self.tube_volumes[self.current_tube] < self.transfer_volume + self.dead_volume:
            self.current_tube += 1
            if self.can_aspirate(self.transfer_volume):
                self.protocol.comment(f"Switching to master mix tube {self.current_tube + 1}")
            else:
                self.protocol.comment("Warning: Not enough master mix to fill all wells!")
    
    def distribute_to_wells(self, wells: List[protocol_api.labware.Well], pipette: protocol_api.instrument_context.InstrumentContext, tip_manager=None) -> None:
        """Distribute master mix to a list of wells.
        
        Args:
            wells: List of wells to distribute to
            pipette: Pipette to use for distribution (can be single or multi-channel)
            tip_manager: Optional TipManager for tip management
        """
        wells_to_fill = wells.copy()
        
        while wells_to_fill:
            # Calculate total volume needed for this aspiration
            total_volume_needed = len(wells_to_fill) * self.transfer_volume
            
            # Calculate maximum volume we can aspirate
            max_aspirate = min(
                self.max_volume,  # Maximum volume of the pipette
                total_volume_needed,  # Total volume needed for all wells
                self.tube_volumes[self.current_tube] - self.dead_volume  # Available volume in current tube
            )
            
            # Safety check: ensure we never exceed pipette max volume
            max_aspirate = min(max_aspirate, self.max_volume)
            
            if not self.can_aspirate(max_aspirate):
                break
            
            # Calculate how many wells we can fill with this aspiration
            wells_per_aspirate = int(max_aspirate // self.transfer_volume)
            current_wells = wells_to_fill[:wells_per_aspirate]
            
            if not current_wells:
                break
                
            # Pick up new tip for each transfer
            if tip_manager:
                tip_manager.get_single_tip(inverse=True)  # Use inverse selection for multichannel
            else:
                pipette.pick_up_tip()
            
            # Aspirate the maximum possible volume
            pipette.aspirate(max_aspirate, self.get_current_tube())
            pipette.touch_tip()
            
            # Dispense into each well
            for well in current_wells:
                pipette.dispense(self.transfer_volume, well)
                pipette.touch_tip()
            
            # Drop tip after transfer
            pipette.drop_tip()
            
            # Update remaining volume and wells to fill
            self.use_volume(max_aspirate)
            wells_to_fill = wells_to_fill[wells_per_aspirate:]

# opentrons_simulate.exe dnabot\template_ot2_scripts\assembly_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

def run(protocol: protocol_api.ProtocolContext):
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_tiprack_20ul"):
        # ============================================================================
        # CONSTANTS
        # ============================================================================
        
        # Pipette settings - pipette instructions in a single location so redefining pipette type is simpler
        PIPETTE_TYPE = 'p20_single_gen2'
        PIPETTE_MOUNT = 'left'
        MULTI_PIPETTE_TYPE = 'p300_multi_gen2'
        MULTI_PIPETTE_MOUNT = 'right'
        
        # Source plates (clip plates) - dynamically loaded based on embeddings
        SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
        
        # Tube rack for master mix
        TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
        TUBE_RACK_POSITION = '4'
        
        # Destination plate
        DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
        DESTINATION_PLATE_SLOT = '7'
        
        # Volume settings
        TOTAL_VOL = 15
        PART_VOL = 1.5
        MIX_SETTINGS = (1, 3)
        
        tiprack_num = tiprack_num + 1
        
        # Error checking
        sample_number = len(final_assembly_dict.keys())
        if sample_number > 96:
            raise ValueError('Assembly number cannot exceed 96.')
        
        # ============================================================================
        # HARDWARE INITIALIZATION
        # ============================================================================
        
        # Load pipettes (without tip_racks argument - TipManager will handle this)
        pipette = protocol.load_instrument(PIPETTE_TYPE, PIPETTE_MOUNT)
        multi_pipette = protocol.load_instrument(MULTI_PIPETTE_TYPE, MULTI_PIPETTE_MOUNT)
        
        # Initialise tip managers
        # Use available slots for p20 tip racks (excluding slots 1, 2 for source plates and slots 7, 8, 10, 11 for thermocycler)
        p20_tip_slots = ['3', '6', '9']  # Available slots for p20 tip racks
        p20_tip_manager = TipManager(protocol, int(p20_tip_slots[0]), pipette, 'p20')  # Initialise with first slot
        
        # Add additional tip racks if needed
        for i in range(1, min(tiprack_num, len(p20_tip_slots))):
            p20_tip_manager.add_tip_rack(int(p20_tip_slots[i]), 'p20')
        
        # Initialise TipManager for p300 tips (slot 5)
        p300_tip_manager = TipManager(protocol, 5, multi_pipette, 'p300')
        
        # Load labware
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        
        # Load source plates dynamically based on embeddings
        source_plates = {}
        source_plate_list = [plate for value in final_assembly_dict.values() for plate in value[1]]  # list of source plates for all clips 
        source_plate_slots = list(set(source_plate_list))  # unique source plates
        for plate in source_plate_slots:
            source_plates[plate] = protocol.load_labware(SOURCE_PLATE_TYPE, plate)

        # Thermocycler Module
        if thermocycler_gen == 'GEN1':
            tc_mod = protocol.load_module('thermocycler', '7')
        else:  # GEN2
            tc_mod = protocol.load_module('thermocyclerModuleV2', '7')

        destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
        tc_mod.open_lid()
        tc_mod.set_block_temperature(20)

        # Error trapping
        sample_number = len(final_assembly_dict.keys())
        if sample_number > 96:
            raise ValueError('Assembly number cannot exceed 96.')

        # Master mix transfers
        final_assembly_lens = [len(values[0]) for values in final_assembly_dict.values()]       # list of assembly lengths (number of clips)
        unique_assemblies_lens = list(set(final_assembly_lens))                                 # unique lengths

        destination_wells = np.array([key for key, value in final_assembly_dict.items()])
        
        # Create separate MasterMixManager for each unique assembly length
        mm_managers_dict = {}
        for assembly_len in unique_assemblies_lens:
            # Calculate master mix volume for this assembly length
            master_mix_volume = TOTAL_VOL - assembly_len * PART_VOL
            # Determine the tube index and well name for this assembly length (A2 = 1, A3 = 2, ...)
            tube_index = assembly_len - 1  # A2 = 1, A3 = 2, etc.
            tube_well_name = f"A{assembly_len}"
            # Create a MasterMixManager for this assembly length, managing only the specific tube
            mm_managers_dict[assembly_len] = MasterMixManager(
                protocol,
                tube_rack,
                [1500],  # Only one tube per manager
                master_mix_volume,
                multi_pipette,  # Use multi-channel pipette for master mix
                dead_volume=15.0,
                tube_position=[tube_well_name]
            )
            # Always use the first tube (index 0) in this manager
            mm_managers_dict[assembly_len].current_tube = 0
            # Store the well index for this manager (for reference/comment)
            protocol.comment(f"Assembly buffer for {assembly_len} parts: {tube_well_name} (well index {tube_index}) - {master_mix_volume}µL per well")
            # Patch get_current_tube to always return the correct well
            def get_current_tube_override(self, idx=tube_index):
                return self.labware.wells()[idx]
            from types import MethodType
            mm_managers_dict[assembly_len].get_current_tube = MethodType(get_current_tube_override, mm_managers_dict[assembly_len])
            
        # Distribute master mix by assembly type using appropriate manager
            for assembly_len in unique_assemblies_lens:
                destination_inds = [i for i, lens in enumerate(final_assembly_lens) if lens == assembly_len]   # find all assemblies of length x
                destination_wells_for_len = list(destination_wells[destination_inds])

                # Create list of destination wells for this assembly type
                wells_to_fill = [destination_plate.wells_by_name()[dest_well] for dest_well in destination_wells_for_len]
                
                # Get the appropriate manager for this assembly length
                mm_manager = mm_managers_dict[assembly_len]
                
                # Distribute master mix to wells for this assembly type using multi-channel pipette
                # Use new tip for each transfer (handled by MasterMixManager)
                mm_manager.distribute_to_wells(wells_to_fill, multi_pipette, tip_manager=p300_tip_manager)

            # Part transfers using TipManager
            for key, values in list(final_assembly_dict.items()):
                for i in range(len(values[0])):                     # find well and plate for every clip in every assembly
                    well  = values[0][i]
                    plate = values[1][i]

                    mix = MIX_SETTINGS
                    # mix = (0,0)
                    # if i == len(values[0])-1:                       # set to mix if on final clip transfer
                    #     mix = MIX_SETTINGS

                    p20_tip_manager.get_single_tip()
                    pipette.transfer(PART_VOL, source_plates[plate].wells(well),
                                     destination_plate.wells(key), mix_after=mix, 
                                     blow_out=True, blowout_location='destination well',
                                     new_tip='never')
                    pipette.drop_tip()

            # Thermocycler Module
            tc_mod.close_lid()
            tc_mod.set_lid_temperature(105)
            tc_mod.set_block_temperature(50, hold_time_minutes=45)
            tc_mod.set_block_temperature(8)
            tc_mod.set_lid_temperature(37)
            # tc_mod.open_lid()                                     # leave lid shut to prevent evaporation
        
        
    def counter(rows):

        def inner(n):
            """ Takes either a value or a well location and converts to other fomat """
            
            row_dict = {0: "A", 1: "B", 2: "C", 3: "D", 4: "E", 5: "F", 6: "G", 7: "H"}

            if type(n) == int:
                row = row_dict[n // rows]
                col = 1 + n % rows
                return row + f'{col}'
                # return row + f'{col:02d}' # for if 2 sf number required (i.e. 'A01' rather than 'A1')

        return inner
    tube_counter = counter(6)


    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

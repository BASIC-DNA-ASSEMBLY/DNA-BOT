from opentrons import protocol_api
import numpy as np
import json
import time
from typing import List, Optional

# metadata
metadata = {
'protocolName': 'DNABOT: Final Assembly v2.15 (Flex)',
'description': 'Final assembly protocol for DNA-BOT using Opentrons OT-2'
}
requirements = {"robotType": "Flex", "apiLevel": "2.15"}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["A1", "C1", "D1"], [2, 2, 2]],
    "B1": [["A1", "E1", "F1"], [2, 2, 2]],
    "C1": [["A1", "G1", "H1"], [2, 2, 2]],
    "D1": [["A1", "A2", "B2"], [2, 2, 2]],
    "E1": [["A1", "C2", "D2"], [2, 2, 2]],
    "F1": [["A1", "E2", "F2"], [2, 2, 2]],
    "G1": [["A1", "G2", "H2"], [2, 2, 2]],
    "H1": [["A1", "A3", "B3"], [2, 2, 2]],
    "A2": [["A1", "C3", "D3"], [2, 2, 2]],
    "B2": [["A1", "E3", "F3"], [2, 2, 2]],
    "C2": [["A1", "G3", "H3"], [2, 2, 2]],
    "D2": [["A1", "A4", "B4"], [2, 2, 2]],
    "E2": [["A1", "C4", "D3"], [2, 2, 2]],
    "F2": [["B1", "D4", "E4"], [2, 2, 2]],
    "G2": [["B1", "F4", "G4"], [2, 2, 2]],
    "H2": [["B1", "H4", "A5"], [2, 2, 2]],
    "A3": [["B1", "B5", "E4"], [2, 2, 2]],
    "B3": [["B1", "C5", "G4"], [2, 2, 2]],
    "C3": [["B1", "D5", "E5"], [2, 2, 2]],
    "D3": [["B1", "F5", "G5"], [2, 2, 2]],
    "E3": [["B1", "H5", "A6"], [2, 2, 2]],
    "F3": [["B1", "B6", "C6"], [2, 2, 2]],
    "G3": [["B1", "D6", "E6"], [2, 2, 2]],
    "H3": [["B1", "F6", "G6"], [2, 2, 2]]
}
tiprack_num = 1

# protocol run function. the part after the colon lets your editor know

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
        pipette_type = 'p1000' if self.pipette.max_volume >= 1000 else 'p50'
        if tip_type != pipette_type:
            raise ValueError(f"Tip type '{tip_type}' is incompatible with pipette type '{pipette_type}'")
            
        if tip_type == 'p1000':
            self.tipracks.append(self.protocol.load_labware('opentrons_flex_96_tiprack_1000ul', slot))
        elif tip_type == 'p50':
            self.tipracks.append(self.protocol.load_labware('opentrons_flex_96_tiprack_50ul', slot))
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

def run(protocol: protocol_api.ProtocolContext):
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_flex_96_tiprack_50ul"):
        # ============================================================================
        # CONSTANTS
        # ============================================================================
        
        # Tiprack settings
        CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '8', '11']  # Exclude slots 1, 2 for source plates
        
        # Pipette settings - pipette instructions in a single location so redefining pipette type is simpler
        PIPETTE_TYPE = 'flex_1channel_50'
        PIPETTE_MOUNT = 'left'
        MULTI_PIPETTE_TYPE = 'flex_8channel_1000'
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
        # Use available slots for p50 tip racks (excluding slots 1, 2 for source plates)
        p50_tip_slots = ['3', '6', '9', '8', '11']  # Available slots for p50 tip racks
        p50_tip_manager = TipManager(protocol, int(p50_tip_slots[0]), pipette, 'p50')  # Initialise with first slot
        
        # Add additional tip racks if needed
        for i in range(1, min(tiprack_num, len(p50_tip_slots))):
            p50_tip_manager.add_tip_rack(int(p50_tip_slots[i]), 'p50')
        
        # Initialise TipManager for p1000 tips (slot 5)
        p1000_tip_manager = TipManager(protocol, 5, multi_pipette, 'p1000')
        
        # Load labware
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        destination_plate = protocol.load_labware(DESTINATION_PLATE_TYPE, DESTINATION_PLATE_SLOT)
        
        # Load source plates dynamically based on embeddings
        source_plates = {}
        source_plate_list = [plate for value in final_assembly_dict.values() for plate in value[1]]  # list of source plates for all clips 
        source_plate_slots = list(set(source_plate_list))  # unique source plates
        for plate in source_plate_slots:
            source_plates[plate] = protocol.load_labware(SOURCE_PLATE_TYPE, plate)

             # Master mix transfers using separate managers for each assembly length
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))

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
                mm_manager.distribute_to_wells(wells_to_fill, multi_pipette, tip_manager=p1000_tip_manager)

            # Part transfers using TipManager
            for key, values in list(final_assembly_dict.items()):
                for i in range(len(values[0])):  # find well and plate for every clip in every assembly
                    well = values[0][i]
                    plate = values[1][i]
                    
                    p50_tip_manager.get_single_tip()
                    pipette.transfer(PART_VOL, source_plates[plate].wells(well),
                                     destination_plate.wells(key), mix_after=MIX_SETTINGS,
                                     new_tip='never')  # transfer parts in one tube
                    pipette.drop_tip()

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

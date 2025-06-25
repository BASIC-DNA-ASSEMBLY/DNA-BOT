from opentrons import protocol_api
import time
from typing import Dict, List, Optional, Tuple, Union

metadata = {
    'protocolName': 'OT Colony PCR',
    'author': 'Liam Hallett',
    'description': 'Colony PCR protocol for Opentrons OT-2',
    'apiLevel': '2.10'
}

# PCR Configuration
PCR_CONFIG = {
    # PCR Wells
    'wells': [
        'A1', 'B1', 'C1', 'D1', 'E1', 'F1',  # Column 1
        'A2', 'B2', 'C2', 'D2', 'E2', 'F2',  # Column 2
        'A3', 'B3', 'C3', 'D3', 'E3', 'F3',  # Column 3
        'A4', 'B4', 'C4', 'D4', 'E4', 'F4',  # Column 4
        'A5', 'B5', 'C5', 'D5', 'E5', 'F5'   # Column 5
    ],
    
    # Blank wells (to be populated by parser)
    'blank_wells': [],
    
    # Transfer volumes (µL)
    'transfer': {
        'master_mix': 10,      # per well
        'primer': 10,          # per well
        'culture': 1,          # per well
        'master_mix_tubes': [1500, 1500, 1500]  # volumes of each master mix tube in order
    },
    
    # PCR Protocol
    'protocol': {
        'initial_denaturation': {
            'temperature': 95,
            'hold_time': 120
        },
        'cycles': 30,
        'denaturation': {
            'temperature': 95,
            'hold_time': 20
        },
        'annealing': {
            'temperature': 55,
            'hold_time': 30
        },
        'extension': {
            'temperature': 72,
            'hold_time': 90
        },
        'final_extension': {
            'temperature': 72,
            'hold_time': 600
        },
        'hold': {
            'temperature': 4
        }
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
        transfer_volume: Volume to transfer to each well (µL)
        pipette: The pipette to use for transfers
        dead_volume: Minimum volume to leave in each tube (µL)
    """
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 labware: protocol_api.labware.Labware,
                 tube_volumes: List[float],
                 transfer_volume: float,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 dead_volume: float = 15.0):
        self.protocol = protocol
        self.labware = labware
        self.tube_volumes = tube_volumes.copy()  # Make a copy to avoid modifying the original
        self.transfer_volume = transfer_volume
        self.current_tube = 0
        self.max_volume = pipette.max_volume
        self.dead_volume = dead_volume
        
        # Get the deck slot from the labware
        deck_slot = labware.parent
        
        # Print tube locations and volumes
        mm_locations = [f"Tube {i+1}: {labware.wells()[i].display_name} - {vol}µL (per well: {transfer_volume}µL)" 
                       for i, vol in enumerate(tube_volumes)]
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
    
    def distribute_to_wells(self, wells: List[protocol_api.labware.Well], pipette: protocol_api.instrument_context.InstrumentContext) -> None:
        """Distribute master mix to a list of wells.
        
        Args:
            wells: List of wells to distribute to
            pipette: Pipette to use for distribution
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
            
            if not self.can_aspirate(max_aspirate):
                break
            
            # Calculate how many wells we can fill with this aspiration
            wells_per_aspirate = max_aspirate // self.transfer_volume
            current_wells = wells_to_fill[:wells_per_aspirate]
            
            if not current_wells:
                break
                
            # Aspirate the maximum possible volume
            pipette.aspirate(max_aspirate, self.get_current_tube())
            pipette.touch_tip()
            
            # Dispense into each well
            for well in current_wells:
                pipette.dispense(self.transfer_volume, well)
                pipette.touch_tip()
            
            # Update remaining volume and wells to fill
            self.use_volume(max_aspirate)
            wells_to_fill = wells_to_fill[wells_per_aspirate:]

def run(protocol: protocol_api.ProtocolContext):

    #### Hardware ####    
    # Load pipettes
    multi_pipette = protocol.load_instrument('p300_multi_gen2', mount='right')
    single_pipette = protocol.load_instrument('p1000_single_gen2', mount='left')
    
    # Load labware
    master_mix_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 1)
    culture_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 5)
    primer_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', 6)
    
    # Load thermocycler and PCR plate
    thermocycler = protocol.load_module('thermocyclerModuleV2', 7)
    PCR_plate = thermocycler.load_labware('biorad_96_wellplate_200ul_pcr')
    
    # Initialise tip managers
    tip_manager_1000 = TipManager(protocol, 4, single_pipette, 'p1000')
    tip_manager_300 = TipManager(protocol, 8, multi_pipette, 'p300')
    tip_manager_300.add_tip_rack(9)  # Add second p300 tip rack
    
    # Initialise well manager
    well_manager = WellManager(PCR_CONFIG['wells'])

    #### Protocol ####
    # Calculate number of master mix tubes needed
    total_master_mix_volume = (len(PCR_CONFIG['wells']) + 3) * PCR_CONFIG['transfer']['master_mix']
    num_tubes = (total_master_mix_volume + PCR_CONFIG['transfer']['master_mix_tubes'][0] - 1) // PCR_CONFIG['transfer']['master_mix_tubes'][0]
    master_mix_wells = [f"A{i+1}" for i in range(num_tubes)]
    protocol.pause(f"Please load {num_tubes} 2x master mix tubes in wells {', '.join(master_mix_wells)} of the master mix rack")
    
    # Set up thermocycler
    protocol.set_rail_lights(True)
    thermocycler.set_block_temperature(8)
    thermocycler.open_lid()

    #### Master Mix Distribution ####
    # Initialise master mix manager
    mm_manager = MasterMixManager(
        protocol,
        master_mix_rack,
        PCR_CONFIG['transfer']['master_mix_tubes'],
        PCR_CONFIG['transfer']['master_mix'],
        single_pipette,
        dead_volume=15.0
    )
    
    # Pick up a single tip for all master mix distributions
    tip_manager_300.get_single_tip(inverse=True)
    
    # Filter out blank wells before distributing master mix
    active_wells = [well for well in PCR_CONFIG['wells'] if well not in PCR_CONFIG['blank_wells']]
    protocol.comment(f"Distributing master mix to {len(active_wells)} wells (excluding {len(PCR_CONFIG['blank_wells'])} blank wells)")
    
    # Distribute master mix to active wells
    mm_manager.distribute_to_wells([PCR_plate[well] for well in active_wells], single_pipette)
    
    # Drop the tip after all distributions are complete
    single_pipette.drop_tip()

    #### Primer Transfer ####    
    for col in well_manager.get_columns_with_wells():
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        # Check if all wells in the column are required
        if len(wells_in_col) == 8:  # All wells in column are used
            protocol.comment(f"Transferring primers to column {col + 1}")
            tip_manager_300.get_multi_tip()
            multi_pipette.aspirate(PCR_CONFIG['transfer']['primer'], primer_plate.wells_by_name()[wells_in_col[0]])
            multi_pipette.dispense(PCR_CONFIG['transfer']['primer'], PCR_plate.wells_by_name()[wells_in_col[0]])
            multi_pipette.touch_tip()
            multi_pipette.drop_tip()
        else:
            # Use single channel for partial columns
            for well in wells_in_col:
                tip_manager_1000.get_single_tip()
                single_pipette.aspirate(PCR_CONFIG['transfer']['primer'], primer_plate[well])
                single_pipette.dispense(PCR_CONFIG['transfer']['primer'], PCR_plate[well])
                single_pipette.touch_tip()
                single_pipette.drop_tip()
    
    #### Culture Transfer ####    
    for col in well_manager.get_columns_with_wells():
        wells_in_col = well_manager.get_wells_in_column(col)
        if not wells_in_col:
            continue
            
        # Check if all wells in the column are required
        if len(wells_in_col) == 8:  # All wells in column are used
            tip_manager_300.get_multi_tip()
            multi_pipette.aspirate(PCR_CONFIG['transfer']['culture'], culture_plate.wells_by_name()[wells_in_col[0]])
            multi_pipette.dispense(PCR_CONFIG['transfer']['culture'], PCR_plate.wells_by_name()[wells_in_col[0]])
            multi_pipette.mix(3, 20)
            multi_pipette.touch_tip()
            multi_pipette.drop_tip()
        else:
            # Use single channel for partial columns
            for well in wells_in_col:
                tip_manager_300.get_single_tip(inverse=True)
                multi_pipette.aspirate(PCR_CONFIG['transfer']['culture'], culture_plate[well])
                multi_pipette.dispense(PCR_CONFIG['transfer']['culture'], PCR_plate[well])
                multi_pipette.mix(3, 20)
                multi_pipette.touch_tip()
                multi_pipette.drop_tip()
    
    #### Run PCR ####    
    protocol.comment("Starting PCR protocol")
    thermocycler.close_lid()
    thermocycler.set_lid_temperature(110)
    
    # Initial denaturation
    protocol.comment("Performing initial denaturation")
    thermocycler.set_block_temperature(
        PCR_CONFIG['protocol']['initial_denaturation']['temperature'],
        hold_time_seconds=PCR_CONFIG['protocol']['initial_denaturation']['hold_time']
    )
    
    # PCR cycles
    for _ in range(PCR_CONFIG['protocol']['cycles']):
        thermocycler.set_block_temperature(
            PCR_CONFIG['protocol']['denaturation']['temperature'],
            hold_time_seconds=PCR_CONFIG['protocol']['denaturation']['hold_time']
        )
        thermocycler.set_block_temperature(
            PCR_CONFIG['protocol']['annealing']['temperature'],
            hold_time_seconds=PCR_CONFIG['protocol']['annealing']['hold_time']
        )
        thermocycler.set_block_temperature(
            PCR_CONFIG['protocol']['extension']['temperature'],
            hold_time_seconds=PCR_CONFIG['protocol']['extension']['hold_time']
        )
    
    # Final extension
    thermocycler.set_block_temperature(
        PCR_CONFIG['protocol']['final_extension']['temperature'],
        hold_time_seconds=PCR_CONFIG['protocol']['final_extension']['hold_time']
    )
    
    # Hold at 4°C
    thermocycler.deactivate_lid()
    thermocycler.set_block_temperature(PCR_CONFIG['protocol']['hold']['temperature'])

    # Protocol complete
    protocol.set_rail_lights(False)

from opentrons import protocol_api
import numpy as np
import json
import time
from typing import Dict, List, Optional, Tuple, Union

# metadata
metadata = {
'protocolName': 'DNABOT Assembly Flex',
'description': 'DNABOT Assembly Step3 without Thermocycler for Opentrons Flex',
'apiLevel': '2.15',
'robotType': 'Flex'
}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["B1", "B2", "C2"], [2, 2, 2]],
    "B1": [["B1", "D2", "B7"], [2, 2, 2]],
    "C1": [["C1", "H4", "G2"], [2, 2, 2]],
    "D1": [["C1", "H2", "A3"], [2, 2, 2]],
    "E1": [["C1", "C7", "C5"], [2, 2, 2]],
    "F1": [["C1", "D3", "H5"], [2, 2, 2]],
    "G1": [["C1", "D7", "H5"], [2, 2, 2]],
    "H1": [["C1", "H3", "H5"], [2, 2, 2]],
    "A2": [["C1", "B4", "C4"], [2, 2, 2]],
    "B2": [["C1", "D4", "E4"], [2, 2, 2]],
    "C2": [["C1", "F4", "E7"], [2, 2, 2]],
    "D2": [["C1", "H4", "A5"], [2, 2, 2]],
    "E2": [["C1", "H2", "A3"], [2, 2, 2]],
    "F2": [["C1", "B5", "C5"], [2, 2, 2]],
    "G2": [["C1", "F7", "G7"], [2, 2, 2]],
    "H2": [["D1", "F4", "G4"], [2, 2, 2]],
    "A3": [["D1", "H4", "A5"], [2, 2, 2]],
    "B3": [["D1", "H7", "A3"], [2, 2, 2]],
    "C3": [["D1", "B5", "C5"], [2, 2, 2]],
    "D3": [["D1", "A8", "E5"], [2, 2, 2]],
    "E3": [["D1", "D3", "F5"], [2, 2, 2]],
    "F3": [["D1", "D7", "H5"], [2, 2, 2]],
    "G3": [["D1", "H3", "F5"], [2, 2, 2]],
    "H3": [["D1", "B4", "E5"], [2, 2, 2]]
}
tiprack_num = 1

# test dict can be used for simulation
#final_assembly_dict={ "A1": ['A7', 'B7', 'C7', 'F7'], "B1": ['A7', 'B7', 'D7', 'G7'], "C1": ['A7', 'B7', 'E7', 'H7']}
#tiprack_num=1

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
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_flex_96_tiprack_300ul"):
            # Constants, we update all the labware name in version 2
            #Tiprack
            CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '2', '5', '8', '11']
            PIPETTE_MOUNT = 'right'
            #Plate of sample after  purification
            MAG_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            MAG_PLATE_POSITION = '1'
            #Tuberack
            TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            TUBE_RACK_POSITION = '7'
            #Destination plate
            DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            #Temperature control plate
            TEMPDECK_SLOT = '4'
            TEMP = 20
            TOTAL_VOL = 15
            PART_VOL = 1.5
            MIX_SETTINGS = (1, 3)
            tiprack_num=tiprack_num+1
            # Errors
            sample_number = len(final_assembly_dict.keys())
            if sample_number > 96:
                raise ValueError('Final assembly nummber cannot exceed 96.')

            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            pipette = protocol.load_instrument('flex_1channel_300', PIPETTE_MOUNT, tip_racks=tipracks)
            
            # Initialize tip manager
            tip_manager = TipManager(protocol, slots[0], pipette, 'p300')
            for slot in slots[1:]:
                tip_manager.add_tip_rack(slot)

            # Define Labware and set temperature
            magbead_plate = protocol.load_labware(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
            tempdeck = protocol.load_module('tempdeck', TEMPDECK_SLOT)
            destination_plate = tempdeck.load_labware(
            DESTINATION_PLATE_TYPE, TEMPDECK_SLOT)
            tempdeck.set_temperature(TEMP)

            # Initialize master mix manager
            mm_manager = MasterMixManager(
                protocol,
                tube_rack,
                [1500],  # Single master mix tube with 1500µL
                TOTAL_VOL,
                pipette,
                dead_volume=15.0
            )

            # Master mix transfers
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))
            master_mix_well_letters = ['A', 'B', 'C', 'D']
            
            # Group assemblies by length for efficient master mix distribution
            for x in unique_assemblies_lens:
                master_mix_well = master_mix_well_letters[(x - 1) // 6] + str(x - 1)
                destination_inds = [i for i, lens in enumerate(final_assembly_lens) if lens == x]
                destination_wells = np.array([key for key, value in list(final_assembly_dict.items())])
                destination_wells = list(destination_wells[destination_inds])
                
                # Calculate master mix volume for this assembly length
                master_mix_volume = TOTAL_VOL - x * PART_VOL
                
                # Create master mix transfers for this assembly length
                master_mix_transfers = []
                for destination_well in destination_wells:
                    master_mix_transfers.append({
                        "source": tube_rack.wells(master_mix_well),
                        "destination": destination_plate.wells(destination_well),
                        "volume": master_mix_volume
                    })
                
                # Perform batch transfer for master mix
                protocol.comment(f"Transferring master mix for {len(destination_wells)} assemblies with {x} parts")
                batch_transfer(master_mix_transfers, pipette, tip_manager)

            # Part transfers
            protocol.comment("Transferring parts")
            part_transfers = []
            for key, values in list(final_assembly_dict.items()):
                for value in values:
                    part_transfers.append({
                        "source": magbead_plate.wells(value),
                        "destination": destination_plate.wells(key),
                        "volume": PART_VOL
                    })
            
            # Perform batch transfer for parts
            batch_transfer(part_transfers, pipette, tip_manager, mix_after=MIX_SETTINGS)

            tempdeck.deactivate() #stop increasing the temperature

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num) 
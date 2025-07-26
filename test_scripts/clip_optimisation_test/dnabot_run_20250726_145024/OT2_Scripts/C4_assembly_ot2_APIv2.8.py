from opentrons import protocol_api
import numpy as np
import json
import time
from typing import Dict, List, Optional, Tuple, Union

# metadata
metadata = {
'protocolName': 'DNABOT Assembly',
'description': 'DNABOT Assembly Step3 without Thermocycler',
'apiLevel': '2.8'
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
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_tiprack_20ul"):
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
            pipette = protocol.load_instrument('p20_single_gen2', PIPETTE_MOUNT, tip_racks=tipracks)
            
            # Initialize tip manager
            tip_manager = TipManager(protocol, int(slots[0]), pipette, 'p20')
            for slot in slots[1:]:
                tip_manager.add_tip_rack(int(slot))

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

from opentrons import protocol_api
import numpy as np
import json
from typing import List

# metadata
metadata = {
'protocolName': 'DNABOT: C3 Final Assembly v2.10',
'description': 'Final assembly protocol for DNA-BOT using Opentrons OT-2',
'apiLevel': '2.10'
}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["A1", "H1", "A2"], [2, 2, 2]],
    "B1": [["A1", "B2", "C2"], [2, 2, 2]],
    "C1": [["A1", "D2", "E2"], [2, 2, 2]],
    "D1": [["A1", "F2", "G2"], [2, 2, 2]],
    "E1": [["A1", "H2", "A3"], [2, 2, 2]],
    "F1": [["A1", "B3", "C3"], [2, 2, 2]],
    "G1": [["A1", "D3", "E3"], [2, 2, 2]],
    "H1": [["A1", "F3", "G3"], [2, 2, 2]],
    "A2": [["A1", "H3", "A4"], [2, 2, 2]],
    "B2": [["A1", "B4", "C4"], [2, 2, 2]],
    "C2": [["A1", "D4", "E4"], [2, 2, 2]],
    "D2": [["A1", "F4", "G4"], [2, 2, 2]],
    "E2": [["A1", "H4", "A5"], [2, 2, 2]],
    "F2": [["B1", "H2", "A3"], [2, 2, 2]],
    "G2": [["B1", "B5", "C5"], [2, 2, 2]],
    "H2": [["B1", "D5", "E5"], [2, 2, 2]],
    "A3": [["B1", "D3", "F5"], [2, 2, 2]],
    "B3": [["B1", "G5", "H5"], [2, 2, 2]],
    "C3": [["B1", "A6", "B6"], [2, 2, 2]],
    "D3": [["B1", "B4", "C6"], [2, 2, 2]],
    "E3": [["B1", "B4", "D6"], [2, 2, 2]],
    "F3": [["B1", "E6", "F6"], [2, 2, 2]],
    "G3": [["B1", "G6", "H6"], [2, 2, 2]],
    "H3": [["B1", "A7", "A2"], [2, 2, 2]]
}
tiprack_num = 1

# protocol run function. the part after the colon lets your editor know

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
    
    def distribute_to_wells(self, wells: List[protocol_api.labware.Well], pipette: protocol_api.instrument_context.InstrumentContext) -> None:
        """Distribute master mix to a list of wells.
        
        Args:
            wells: List of wells to distribute to
            pipette: Pipette to use for distribution (can be single or multi-channel)
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

# test dict can be used for simulation
#final_assembly_dict={ "A1": ['A7', 'B7', 'C7', 'F7'], "B1": ['A7', 'B7', 'D7', 'G7'], "C1": ['A7', 'B7', 'E7', 'H7']}
#tiprack_num=1
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
                raise ValueError('Assembly number cannot exceed 96.')

            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
            # Pipette
            PIPETTE_MOUNT = 'right'      
            pipette = protocol.load_instrument('p20_single_gen2', PIPETTE_MOUNT, tip_racks=tipracks)
            
            # Multi-channel pipette for master mix distribution
            MULTI_PIPETTE_MOUNT = 'left'
            multi_pipette = protocol.load_instrument('p300_multi_gen2', MULTI_PIPETTE_MOUNT, tip_racks=[protocol.load_labware('opentrons_96_tiprack_300ul', '5')])


            # Define Labware and set temperature
            magbead_plate = protocol.load_labware(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
            tempdeck = protocol.load_module('tempdeck', TEMPDECK_SLOT)
            destination_plate = tempdeck.load_labware(
            DESTINATION_PLATE_TYPE, TEMPDECK_SLOT)
            tempdeck.set_temperature(TEMP)

             # Master mix transfers using separate managers for each assembly length
            final_assembly_lens = []
            for values in final_assembly_dict.values():
                final_assembly_lens.append(len(values))
            unique_assemblies_lens = list(set(final_assembly_lens))
            
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
                mm_manager.distribute_to_wells(wells_to_fill, multi_pipette)

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for value in values:# magbead_plate.wells and destination_plate.wells in the same type
                    pipette.transfer(PART_VOL, magbead_plate.wells(value),
                                     destination_plate.wells(key), mix_after=MIX_SETTINGS,
                                     new_tip='always')#transfer parts in one tube

            tempdeck.deactivate() #stop increasing the temperature

    final_assembly(final_assembly_dict=final_assembly_dict, tiprack_num=tiprack_num)

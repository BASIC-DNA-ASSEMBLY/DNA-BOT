from __future__ import unicode_literals
from opentrons import protocol_api
import numpy as np
import json
from typing import List

# metadata
metadata = {
'protocolName': 'DNABOT: C2 Final Assembly with Thermocycler v2.10',
'description': 'Final assembly protocol for DNA-BOT using Opentrons OT-2 with thermocycler module',
'apiLevel': '2.10'
}

# Load assembly data from JSON file
# This will be replaced by the parser with embedded JSON data
final_assembly_dict = {
    "A1": [["B1", "E3", "G4"], [1, 1, 1]],
    "B1": [["B1", "H4", "A5"], [1, 1, 1]],
    "C1": [["B1", "B5", "C5"], [1, 1, 1]],
    "D1": [["B1", "D5", "E5"], [1, 1, 1]],
    "E1": [["B1", "F5", "G5"], [1, 1, 1]],
    "F1": [["B1", "H5", "A6"], [1, 1, 1]],
    "G1": [["B1", "B6", "B2"], [1, 1, 1]],
    "H1": [["B1", "C6", "D6"], [1, 1, 1]],
    "A2": [["B1", "E6", "F6"], [1, 1, 1]],
    "B2": [["B1", "G6", "H6"], [1, 1, 1]],
    "C2": [["C1", "A7", "B7"], [1, 1, 1]],
    "D2": [["C1", "C7", "D7"], [1, 1, 1]],
    "E2": [["C1", "E7", "F7"], [1, 1, 1]],
    "F2": [["C1", "G3", "C5"], [1, 1, 1]],
    "G2": [["C1", "A4", "B4"], [1, 1, 1]],
    "H2": [["C1", "A2", "B2"], [1, 1, 1]],
    "A3": [["C1", "H5", "D2"], [1, 1, 1]],
    "B3": [["C1", "G7", "H7"], [1, 1, 1]],
    "C3": [["C1", "E6", "A8"], [1, 1, 1]],
    "D3": [["C1", "B8", "H7"], [1, 1, 1]],
    "E3": [["C1", "C8", "D2"], [1, 1, 1]],
    "F3": [["C1", "C8", "D6"], [1, 1, 1]],
    "G3": [["C1", "E2", "D8"], [1, 1, 1]],
    "H3": [["D1", "G2", "H2"], [1, 1, 1]]
}
tiprack_num = 1

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# It is possible to run 88 assemblies with this new module. The heat block module is removed. 
# Assembly reactions is set up on thermocycler module.

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

# test dictionary can be used for simulation 3 or 88 assemblies
#final_assembly_dict={"A1": [['A7', 'B7', 'C7', 'F7'], [1, 2, 1, 1]], "B1": [['A7', 'B7', 'D7', 'G7'], [1, 2, 1, 1]], "C1": [['A7', 'E7', 'H7'], [1, 2, 1]]}
#tiprack_num=1

# final_assembly_dict={"A1": [["A1", "C9", "B11"], [1, 2, 1]], "B1": [["A1", "C9", "C11"], [1, 2, 1]], "C1": [["A1", "C9", "D11"], [1, 2, 1]], "D1": [["A1", "C9", "E11"], [1, 2, 1]], "E1": [["A1", "C9", "F11"], [1, 2, 1]], "F1": [["A1", "C9", "G11"], [1, 2, 1]], "G1": [["A1", "C9", "H11"], [1, 2, 1]], "H1": [["A1", "C9", "A12"], [1, 2, 1]], "A2": [["A1", "C9", "B12"], [1, 2, 1]], "B2": [["A1", "D9", "B11"], [1, 2, 1]], "C2": [["A1", "D9", "C11"], [1, 2, 1]], "D2": [["A1", "D9", "D11"], [1, 2, 1]], "E2": [["A1", "D9", "E11"], [1, 2, 1]], "F2": [["A1", "D9", "F11"], [1, 2, 1]], "G2": [["A1", "D9", "G11"], [1, 2, 1]], "H2": [["B1", "D9", "H11"], [1, 2, 1]], "A3": [["B1", "D9", "A12"], [1, 2, 1]], "B3": [["B1", "D9", "B12"], [1, 2, 1]], "C3": [["B1", "E9", "F12"], [1, 2, 1]], "D3": [["B1", "E9", "G12"], [1, 2, 1]], "E3": [["B1", "E9", "H12"], [1, 2, 1]], "F3": [["B1", "E9", "A1"], [1, 2, 2]], "G3": [["B1", "E9", "B1"], [1, 2, 2]], "H3": [["B1", "E9", "C1"], [1, 2, 2]], "A4": [["B1", "E9", "D1"], [1, 2, 2]], "B4": [["B1", "E9", "E1"], [1, 2, 2]], "C4": [["B1", "E9", "F1"], [1, 2, 2]], "D4": [["B1", "F9", "F12"], [1, 2, 1]], "E4": [["B1", "F9", "G12"], [1, 2, 1]], "F4": [["B1", "F9", "H12"], [1, 2, 1]], "G4": [["C1", "F9", "A1"], [1, 2, 2]], "H4": [["C1", "F9", "B1"], [1, 2, 2]], "A5": [["C1", "F9", "C1"], [1, 2, 2]], "B5": [["C1", "F9", "D1"], [1, 2, 2]], "C5": [["C1", "F9", "E1"], [1, 2, 2]], "D5": [["C1", "F9", "F1"], [1, 2, 2]], "E5": [["C1", "G9", "F12"], [1, 2, 1]], "F5": [["C1", "G9", "G12"], [1, 2, 1]], "G5": [["C1", "G9", "H12"], [1, 2, 1]], "H5": [["C1", "G9", "A1"], [1, 2, 2]], "A6": [["C1", "G9", "B1"], [1, 2, 2]], "B6": [["C1", "G9", "C1"], [1, 2, 2]], "C6": [["C1", "G9", "D1"], [1, 2, 2]], "D6": [["C1", "G9", "E1"], [1, 2, 2]], "E6": [["C1", "G9", "F1"], [1, 2, 2]], "F6": [["D1", "H9", "B2"], [1, 2, 2]], "G6": [["D1", "H9", "C2"], [1, 2, 2]], "H6": [["D1", "H9", "D2"], [1, 2, 2]], "A7": [["D1", "H9", "E2"], [1, 2, 2]], "B7": [["D1", "H9", "F2"], [1, 2, 2]], "C7": [["D1", "H9", "G2"], [1, 2, 2]], "D7": [["D1", "H9", "H2"], [1, 2, 2]], "E7": [["D1", "H9", "A3"], [1, 2, 2]], "F7": [["D1", "H9", "B3"], [1, 2, 2]], "G7": [["D1", "A10", "B2"], [1, 2, 2]], "H7": [["D1", "A10", "C2"], [1, 2, 2]], "A8": [["D1", "A10", "D2"], [1, 2, 2]], "B8": [["D1", "A10", "E2"], [1, 2, 2]], "C8": [["D1", "A10", "F2"], [1, 2, 2]], "D8": [["D1", "A10", "G2"], [1, 2, 2]], "E8": [["E1", "A10", "H2"], [1, 2, 2]], "F8": [["E1", "A10", "A3"], [1, 2, 2]], "G8": [["E1", "A10", "B3"], [1, 2, 2]], "H8": [["E1", "B10", "B2"], [1, 2, 2]], "A9": [["E1", "B10", "C2"], [1, 2, 2]], "B9": [["E1", "B10", "D2"], [1, 2, 2]], "C9": [["E1", "B10", "E2"], [1, 2, 2]], "D9": [["E1", "B10", "F2"], [1, 2, 2]], "E9": [["E1", "B10", "G2"], [1, 2, 2]], "F9": [["E1", "B10", "H2"], [1, 2, 2]], "G9": [["E1", "B10", "A3"], [1, 2, 2]], "H9": [["E1", "B10", "B3"], [1, 2, 2]]}
# tiprack_num=3

# opentrons_simulate.exe dnabot\template_ot2_scripts\assembly_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

def run(protocol: protocol_api.ProtocolContext):
    def final_assembly(final_assembly_dict, tiprack_num, tiprack_type="opentrons_96_tiprack_20ul"):
            ### Constants

            # Source plate(s)
            SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            source_plate_list = [plate for value in final_assembly_dict.values() for plate in value[1]]     # list of source plates for all clips 
            source_plate_slots = list(set(source_plate_list))                                               # unique source plates
            source_plates = {plate: protocol.load_labware(SOURCE_PLATE_TYPE, plate) for plate in source_plate_slots}

            # Tuberack
            TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            TUBE_RACK_POSITION = '4'
            tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)

            # Destination plate
            DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            TOTAL_VOL = 15
            PART_VOL = 1.5
            MIX_SETTINGS = (1, 3)
            # tiprack_num += 1                    # + 1 for one index ############################### I think(?)

            # Thermocycler Module
            if thermocycler_gen == 'gen1':
                tc_mod = protocol.load_module('Thermocycler Module')
            else:  # gen2
                tc_mod = protocol.load_module('thermocyclerModuleV2')

            destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
            tc_mod.open_lid()
            tc_mod.set_block_temperature(20)

            # Error trapping
            sample_number = len(final_assembly_dict.keys())
            if sample_number > 96:
                raise ValueError('Assembly number cannot exceed 96.')

            # Tiprack(s)
            CANDIDATE_TIPRACK_SLOTS = ['3', '5', '6', '9']

            if 2 not in source_plate_slots:                  # if only one source plate used, deck slot 2 can be used for a tip rack
                CANDIDATE_TIPRACK_SLOTS.append('2')

            if tiprack_num > len(CANDIDATE_TIPRACK_SLOTS):
                raise ValueError('Not enough tipracks available on deck to satisfy tip requirements. Consider either splitting into multiple builds each with fewer constructs or iterative rounds of building. ')
                  
            slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
            tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]

            # Pipette
            PIPETTE_MOUNT = 'right'      
            pipette = protocol.load_instrument('p20_single_gen2', PIPETTE_MOUNT, tip_racks=tipracks)
            
            # Multi-channel pipette for master mix distribution
            MULTI_PIPETTE_MOUNT = 'left'
            multi_pipette = protocol.load_instrument('p300_multi_gen2', MULTI_PIPETTE_MOUNT, tip_racks=[protocol.load_labware('opentrons_96_tiprack_300ul', '5')])

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
                mm_manager.distribute_to_wells(wells_to_fill, multi_pipette)

            # Part transfers
            for key, values in list(final_assembly_dict.items()):
                for i in range(len(values[0])):                     # find well and plate for every clip in every assembly
                    well  = values[0][i]
                    plate = values[1][i]

                    mix = MIX_SETTINGS
                    # mix = (0,0)
                    # if i == len(values[0])-1:                       # set to mix if on final clip transfer
                    #     mix = MIX_SETTINGS

                    pipette.transfer(PART_VOL, source_plates[plate].wells(well),
                                     destination_plate.wells(key), mix_after=mix, 
                                     blow_out=True, blowout_location='destination well',
                                     new_tip='always')

            # Thermocycler Module
            tc_mod.close_lid()
            tc_mod.set_lid_temperature(105)
            tc_mod.set_block_temperature(50, hold_time_minutes=45, block_max_volume=15)
            tc_mod.set_block_temperature(8, block_max_volume=30)
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

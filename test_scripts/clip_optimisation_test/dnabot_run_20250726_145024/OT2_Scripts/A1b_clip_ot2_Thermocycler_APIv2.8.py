from opentrons import protocol_api
import json
import time
from typing import Dict, List, Optional, Tuple, Union

# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.8',
     'protocolName': 'CLIP_With_Thermocycler',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module.'}

# Load CLIP data from JSON file
# This will be replaced by the parser with embedded JSON data
clips_dict = {
    "A7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "sucB_1",
        "part_source_well": "C3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "alaC_1",
        "part_source_well": "B2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "thiS_1",
        "part_source_well": "C5",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "serC_1",
        "part_source_well": "E4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "alr_1",
        "part_source_well": "C2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "dadX_1",
        "part_source_well": "G3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tdcB_2",
        "part_source_well": "F6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "gadB_1",
        "part_source_well": "G4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H7",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ilvA_1",
        "part_source_well": "A8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvH_1",
        "part_source_well": "A2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D8",
        "plate": 1,
        "part_vol": 1.0,
        "water_vol": 0.0
    }
}

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'gen2'

# opentrons_simulate.exe dnabot\template_ot2_scripts\clip_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["1", "1", "1", "1", "1"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

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
    # added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    tiprack_type="opentrons_96_tiprack_20ul"
    INITIAL_TIP = 'A1'
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '5']

    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = 'p20_single_gen2'
    PIPETTE_MOUNT = 'right'
        ### Load Pipette
        # checks if it's a P10 Single pipette
    if PIPETTE_TYPE != 'p20_single_gen2':
        print('Define labware must be changed to use', PIPETTE_TYPE)
        exit()

    # Thermocycler Module
    if thermocycler_gen == 'gen1':
        tc_mod = protocol.load_module('Thermocycler Module')
    else:  # gen2
        tc_mod = protocol.load_module('thermocyclerModuleV2')
        
    # Destination Plates
    DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'

    # Loads destination plate onto Thermocycler Module
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
    tc_mod.open_lid()
    tc_mod.set_block_temperature(20)

    # Source Plates
    SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used

    # Tube Rack
    TUBE_RACK_TYPE = 'e14151500starlab_24_tuberack_1500ul'
            # modified from custom labware as API 2 doesn't support labware.create anymore, so the old add_labware script can't be used
    TUBE_RACK_POSITION = '4'
    MASTER_MIX_WELL = 'A1'
    WATER_WELL = 'A2'
    MASTER_MIX_VOLUME = 20

    # Mix settings
    LINKER_MIX_SETTINGS = (1, 3)
    PART_MIX_SETTINGS = (4, 5)

    def clip(clips_dict):
        ### Loading Tiprack
        total_tips = 4 * len(clips_dict)
        letter_dict = {'A': 0, 'B': 1, 'C': 2,
                       'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
        tiprack_1_tips = (
            13 - int(INITIAL_TIP[1:])) * 8 - letter_dict[INITIAL_TIP[0]]
        if total_tips > tiprack_1_tips:
            tiprack_num = 1 + (total_tips - tiprack_1_tips) // 96 + \
            (1 if (total_tips - tiprack_1_tips) % 96 > 0 else 0)
        else:
            tiprack_num = 1
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT, tip_racks=tipracks)
        
        # Initialize tip manager
        tip_manager = TipManager(protocol, int(slots[0]), pipette, 'p20')
        for slot in slots[1:]:
            tip_manager.add_tip_rack(int(slot))
        
        # Destination plate on thermocycler
        destination_wells = [destination_plate.wells_by_name()[w] for w in clips_dict.keys()]
        
        # Tube rack
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        master_mix = tube_rack.wells(MASTER_MIX_WELL)
        water = tube_rack.wells(WATER_WELL)
        
        # Load source plates
        source_plates = {}
        all_plates = set()
        for well_info in clips_dict.values():
            all_plates.add(well_info['prefix_plate'])
            all_plates.add(well_info['suffix_plate'])
            all_plates.add(well_info['part_plate'])
        for key in all_plates:
            source_plates[key] = protocol.load_labware(SOURCE_PLATE_TYPE, key)
        
        dest_wells = list(clips_dict.keys())
        
        # Initialize master mix manager
        mm_manager = MasterMixManager(
            protocol,
            tube_rack,
            [1500],  # Single master mix tube with 1500µL
            MASTER_MIX_VOLUME,
            pipette,
            dead_volume=15.0
        )
        
        # Master mix transfer using batch distribution
        protocol.comment("Transferring master mix")
        mm_manager.distribute_to_wells(destination_wells, pipette)
        
        # Water transfer (only if needed)
        water_vols = [clips_dict[w]['water_vol'] for w in dest_wells]
        if any([wv > 0 for wv in water_vols]):
            protocol.comment("Transferring water")
            water_transfers = []
            for i, well in enumerate(dest_wells):
                if water_vols[i] > 0:
                    water_transfers.append({
                        "source": water,
                        "destination": destination_wells[i],
                        "volume": water_vols[i]
                    })
            batch_transfer(water_transfers, pipette, tip_manager)
        
        # Prefix transfers
        protocol.comment("Transferring prefixes")
        prefix_transfers = []
        for i, well in enumerate(dest_wells):
            info = clips_dict[well]
            prefix_transfers.append({
                "source": source_plates[info['prefix_plate']].wells_by_name()[info['prefix_well']],
                "destination": destination_wells[i],
                "volume": 1
            })
        batch_transfer(prefix_transfers, pipette, tip_manager, mix_after=LINKER_MIX_SETTINGS)
        
        # Suffix transfers
        protocol.comment("Transferring suffixes")
        suffix_transfers = []
        for i, well in enumerate(dest_wells):
            info = clips_dict[well]
            suffix_transfers.append({
                "source": source_plates[info['suffix_plate']].wells_by_name()[info['suffix_well']],
                "destination": destination_wells[i],
                "volume": 1
            })
        batch_transfer(suffix_transfers, pipette, tip_manager, mix_after=LINKER_MIX_SETTINGS)
        
        # Part transfers
        protocol.comment("Transferring parts")
        part_transfers = []
        for i, well in enumerate(dest_wells):
            info = clips_dict[well]
            part_transfers.append({
                "source": source_plates[info['part_plate']].wells_by_name()[info['part_well']],
                "destination": destination_wells[i],
                "volume": info['part_vol']
            })
        batch_transfer(part_transfers, pipette, tip_manager, mix_after=PART_MIX_SETTINGS)
    
    # the run function will first define the CLIP function, and then run the CLIP function with the dictionary produced by DNA-BOT
    clip(clips_dict)
    
    ### PCR Reaction in Thermocycler

    # close lid and set lid temperature, PCR will not start until lid reaches 37C
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)

    # Runs 20 cycles of 37C for 2 minutes and 20C for 1 minute, then holds for 60C for 10 minutes
    profile = [
        {'temperature': 37, 'hold_time_minutes': 2},
        {'temperature': 20, 'hold_time_minutes': 1}]
    tc_mod.execute_profile(steps=profile, repetitions=20, block_max_volume=30)
    tc_mod.set_block_temperature(60, hold_time_minutes=10, block_max_volume=30)
    tc_mod.set_block_temperature(8, block_max_volume=30)
    tc_mod.set_lid_temperature(37)
    # tc_mod.open_lid()                                     # leave lid shut to prevent evaporation

from opentrons import protocol_api
import time
from typing import Dict, List, Optional, Tuple, Union

''' To Do
     - explicitly code for tip change
'''

metadata = {
     'apiLevel': '2.8',
     'protocolName': 'purification_template',
     'description': 'Implements magbead purification reactions for BASIC assembly using an opentrons OT-2'}

# example values produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template:
sample_number=24
ethanol_well='A3'

# opentrons_simulate.exe dnabot\template_ot2_scripts\purification_template_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

class WellManager:
    """Manages well positions and column operations."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 sample_number: int):
        self.protocol = protocol
        self.sample_number = sample_number
        self.rows = 8
        self.cols = 12
        
        # Calculate number of columns needed
        self.src1_col_num = min((sample_number - 1) // 8 + 1, 6)  # Max 6 columns in first plate
        self.src2_col_num = max(0, (sample_number - 48 - 1) // 8 + 1)  # Remaining columns in second plate
    
    def get_columns_with_wells(self) -> List[int]:
        """Get list of column indices that contain samples."""
        columns = list(range(self.src1_col_num))
        if self.src2_col_num > 0:
            columns.extend(range(self.src1_col_num, self.src1_col_num + self.src2_col_num))
        return columns
    
    def get_source_well(self, col: int, source_plate, source_plate2=None) -> str:
        """Get the source well for a given column index."""
        if col < self.src1_col_num:
            # First plate (columns 0-5)
            return source_plate.wells_by_name()[f'A{col + 1}']
        else:
            # Second plate (columns 6-11)
            if not source_plate2:
                raise ValueError("Source plate 2 is required for columns >= 6")
            return source_plate2.wells_by_name()[f'A{col - self.src1_col_num + 1}']
    
    def get_mag_well(self, col: int, mag_plate) -> str:
        """Get the magnetic plate well for a given column index."""
        return mag_plate.wells_by_name()[f'A{col + 1}']
    
    def get_final_well(self, col: int, final_plate) -> str:
        """Get the final plate well for a given column index."""
        return final_plate.wells_by_name()[f'A{col + 1}']

class TipManager:
    """Manages pipette tips and tracks usage."""
    
    def __init__(self, protocol: protocol_api.ProtocolContext, 
                 slot: int,
                 pipette: protocol_api.instrument_context.InstrumentContext,
                 tip_type: str):
        self.protocol = protocol
        self.pipette = pipette
        self.tip_type = tip_type
        self.rows = 8
        self.cols = 12
        
        # Initialise tip racks array and current rack index
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
        """Pick up multiple tips columnwise for multichannel pipetting."""
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

def calculate_transfer_time_duration(num_columns: int) -> float:
    """
    Calculate transfer time duration based on number of columns.
    Excludes first column, uses 1/8 minute per additional column.
    
    Args:
        num_columns: Number of columns being processed
        
    Returns:
        Transfer time duration in minutes
    """
    if num_columns <= 1:
        return 0.0
    
    # Exclude first column, calculate time for remaining columns
    additional_columns = num_columns - 1
    return additional_columns * (1/8)  # 1/8 minute per column

def run(protocol: protocol_api.ProtocolContext):
    # added run function for API verison 2

    def magbead(
            sample_number,
            ethanol_well,
            elution_buffer_well='A1',
            sample_volume=30,
            bead_ratio=1.8,
            elution_buffer_volume=40,
            incubation_time=5,
            settling_time=5,                    # 2
                # if using Gen 2 magentic module, need to change time! see: https://docs.opentrons.com/v2/new_modules.html
                # "The GEN2 Magnetic Module uses smaller magnets than the GEN1 version...this means it will take longer for the GEN2 module to attract beads."
                # Recommended Magnetic Module GEN2 bead attraction time:
                    # Total liquid volume <= 50 uL: 5 minutes
                # this template was written with the Gen 1 magnetic module, as it is compatible with API version 2
            drying_time=5,
            elution_time=2,
            sample_offset=0,
            tiprack_type="opentrons_96_tiprack_300ul"):

        """
        Selected args:
            ethanol_well (str): well in reagent container containing ethanol.
            elution_buffer_well (str): well in reagent container containing elution buffer.
            sample_offset (int): offset the intial sample column by the specified value.
        """

        ### Load Labware
        # Pipetting constants
        CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9', '10', '11']      # Calculates whether one/two/three/four/five tipracks are needed, which are in slots 3, 6, 9, 2, and 5 respectively
        TIPS_PER_SAMPLE = 9
        PIPETTE_ASPIRATE_RATE = 25
        PIPETTE_DISPENSE_RATE = 150

        # Tipracks
        total_tips = sample_number * TIPS_PER_SAMPLE
        tiprack_num = (total_tips - 1) // 96 + 1
        slots = CANDIDATE_TIPRACK_SLOTS[:tiprack_num]
        tipracks = [protocol.load_labware(tiprack_type, slot) for slot in slots]

        # Pipettes
        PIPETTE_TYPE = 'p300_multi_gen2'
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount="left", tip_racks=tipracks)
        pipette.aspirate_flow_rate = PIPETTE_ASPIRATE_RATE
        pipette.dispense_flow_rate = PIPETTE_DISPENSE_RATE  # for reference: default aspirate/dispense flow rate for p300_multi_gen2 is 94 ul/s

        # Initialize tip manager
        tip_manager = TipManager(protocol, int(slots[0]), pipette, 'p300')
        for slot in slots[1:]:
            tip_manager.add_tip_rack(int(slot))

        # Source plate(s)
        SOURCE_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
        SRC1_PLATE_POSITION = '1'
        src1_plate = protocol.load_labware(SOURCE_PLATE_TYPE, SRC1_PLATE_POSITION)

        # sample_number = 56
        if sample_number > 48:                              # initialise source plate 2 if more than 48 clips
            SRC2_PLATE_POSITION = '2'
            src2_plate = protocol.load_labware(SOURCE_PLATE_TYPE, SRC2_PLATE_POSITION) 
            src1_sample_number = 48
            src2_sample_number = sample_number - 48
        else:
            src1_sample_number = sample_number
            src2_sample_number = 0

        # Magnetic Module
        MAG_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
        MAGDECK_POSITION = '4'
        MAGDECK = protocol.load_module('magdeck', MAGDECK_POSITION)
        MAGDECK.disengage()                                 # disengages the magnets when it is turned on
        mag_plate = MAGDECK.load_labware(MAG_PLATE_TYPE)

        # Dest Plate
        DEST_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'
        DEST_PLATE_POSITION = '5'
        dest_plate = protocol.load_labware(DEST_PLATE_TYPE, DEST_PLATE_POSITION)

        # Reagents
        REAGENT_CONTAINER_TYPE = '4ti0131_12_reservoir_21000ul'
        REAGENT_CONTAINER_POSITION = '7'
        reagent_container = protocol.load_labware(REAGENT_CONTAINER_TYPE, REAGENT_CONTAINER_POSITION)

        # Beads
        BEAD_CONTAINER_TYPE = '4ti0136_96_wellplate_2200ul'
        BEAD_CONTAINER_POSITION = '8'
        bead_container = protocol.load_labware(BEAD_CONTAINER_TYPE, BEAD_CONTAINER_POSITION)

        # Initialize well manager
        well_manager = WellManager(protocol, sample_number)

        ### Settings
        LIQUID_WASTE_WELL = 'A5'
        BEADS_WELL = 'A1'
        DEAD_TOTAL_VOL = 6
        SLOW_HEAD_SPEEDS = {'x': 600 // 4, 'y': 400 // 4, 'z': 125 // 10, 'a': 125 // 10}
        DEFAULT_HEAD_SPEEDS = {'x': 400, 'y': 400, 'z': 125, 'a': 100}
        IMMOBILISE_MIX_REPS = 10
        MAGDECK_HEIGHT = 20
        AIR_VOL_COEFF = 0.1
        ETHANOL_VOL = 150
        WASH_TIME = 0.5
        ETHANOL_DEAD_VOL = 50
        ELUTION_MIX_REPS = 20
        ELUTANT_SEP_TIME = 1
        ELUTION_DEAD_VOL = 2

        ### Protocol set up
        # Total columns across source plates
        src1_col_num = (src1_sample_number - 1) // 8 + 1
        src2_col_num = (src2_sample_number - 1) // 8 + 1

        # Source plates columns (i.e. position 1 and 2)
            # generates a list of lists: [[A1, B1, C1...], [A2, B2, C2...]...]
        samples = [col for col in src1_plate.columns()[sample_offset : src1_col_num + sample_offset]]

        if sample_number > 48:      # if 2 source plates are required initialise source plate 2 samples and add them to the samples list
            src2_samples = [col for col in src2_plate.columns()[sample_offset : src2_col_num + sample_offset]]
            samples = samples + src2_samples
        
        # Mag plate 
        col_num = src1_col_num + src2_col_num
        mag_cols = [col for col in mag_plate.columns()[sample_offset:col_num + sample_offset]]

        # Output
        output  = [col for col in dest_plate.columns()[sample_offset:col_num + sample_offset]]

        ### Defining Wells for Reagents, Liquid Waste, and Beads
        liquid_waste = reagent_container.wells(LIQUID_WASTE_WELL)
        ethanol = reagent_container.wells(ethanol_well)
        elution_buffer = reagent_container.wells(elution_buffer_well)
        beads = bead_container[BEADS_WELL]

        ### Define bead and mix volume
        bead_volume = sample_volume * bead_ratio
        if bead_volume / 2 > pipette.max_volume:
            mix_vol = pipette.max_volume
        else:
            mix_vol = bead_volume / 2
        total_vol = bead_volume + sample_volume + DEAD_TOTAL_VOL

        # Get columns with wells
        columns_with_wells = well_manager.get_columns_with_wells()

        ### Steps
        # Mix beads and parts
        protocol.comment("Mixing beads and parts")
        for col in columns_with_wells:
            tip_manager.get_multi_tip(start_col=col)
            source_well = well_manager.get_source_well(col, src1_plate, src2_plate if sample_number > 48 else None)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            
            # Aspirate beads
            pipette.aspirate(bead_volume, beads)
            protocol.max_speeds.update(SLOW_HEAD_SPEEDS)

            # Transfer and mix on mag plate
            pipette.mix(IMMOBILISE_MIX_REPS, mix_vol, source_well)
            pipette.transfer(total_vol, source_well, mag_well, new_tip='never', blow_out=True, blowout_location='destination well')

            # Dispose of tip
            protocol.max_speeds.update(DEFAULT_HEAD_SPEEDS)
            tip_manager.return_tips(col)

        # Initial mix and incubation sample
        protocol.comment("Initial incubation")
        transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
        adjusted_delay = max(0, incubation_time - transfer_time)
        if adjusted_delay > 0:
            protocol.delay(minutes=adjusted_delay)

        # Engagae MagDeck and incubate
        protocol.comment("Engaging magnet")
        MAGDECK.engage(height=MAGDECK_HEIGHT)
        transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
        adjusted_delay = max(0, settling_time - transfer_time)
        if adjusted_delay > 0:
            protocol.delay(minutes=adjusted_delay)

        # Remove supernatant from magnetic beads
        protocol.comment("Removing supernatant")
        for col in columns_with_wells:
            tip_manager.get_multi_tip(start_col=col)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            
            pipette.aspirate(total_vol, mag_well)
            pipette.move_to(reagent_container['A1'].top(z=70))
            pipette.blow_out(protocol.fixed_trash["A1"].top())
            tip_manager.return_tips(col)

        # Wash beads twice with 70% ethanol
        air_vol = pipette.max_volume * AIR_VOL_COEFF
        for cycle in range(2):
            protocol.comment(f"Ethanol wash {cycle + 1}")
            
            # Add ethanol
            for col in columns_with_wells:
                tip_manager.get_multi_tip(start_col=col)
                mag_well = well_manager.get_mag_well(col, mag_plate)
                
                pipette.transfer(ETHANOL_VOL, ethanol, mag_well, trash=False, air_gap=air_vol)
                tip_manager.return_tips(col)
            
            # Calculate dynamic delay based on transfer time
            transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
            adjusted_delay = max(0, WASH_TIME - transfer_time)
            if adjusted_delay > 0:
                protocol.delay(minutes=adjusted_delay)
            
            # Remove ethanol
            for col in columns_with_wells:
                tip_manager.get_multi_tip(start_col=col)
                mag_well = well_manager.get_mag_well(col, mag_plate)
                
                pipette.aspirate(ETHANOL_VOL + ETHANOL_DEAD_VOL, mag_well)
                pipette.move_to(reagent_container['A1'].top(z=70))
                pipette.blow_out(protocol.fixed_trash["A1"].top())
                tip_manager.return_tips(col)
        
        # Dry at room temperature
        protocol.comment("Drying time")
        transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
        adjusted_delay = max(0, drying_time - transfer_time)
        if adjusted_delay > 0:
            protocol.delay(minutes=adjusted_delay)

        # Disengage MagDeck
        MAGDECK.disengage()

        # Mix beads with elution buffer (select mix vol to work with either p20 or p300)
        if elution_buffer_volume / 2 > pipette.max_volume:
            mix_vol = pipette.max_volume
        else:
            mix_vol = elution_buffer_volume / 2

        protocol.comment("Adding elution buffer")
        for col in columns_with_wells:
            tip_manager.get_multi_tip(start_col=col)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            
            pipette.transfer(elution_buffer_volume, elution_buffer, mag_well, mix_after=(ELUTION_MIX_REPS, mix_vol))
            tip_manager.return_tips(col)

        # Incubate at room temperature
        protocol.comment("Elution incubation")
        transfer_time = calculate_transfer_time_duration(len(columns_with_wells))
        adjusted_delay = max(0, elution_time - transfer_time)
        if adjusted_delay > 0:
            protocol.delay(minutes=adjusted_delay)

        # Engage MagDeck (remains engaged for DNA elution)
        protocol.comment("Engaging magnet for elution")
        MAGDECK.engage(height=MAGDECK_HEIGHT)
        protocol.delay(minutes=ELUTANT_SEP_TIME)

        # Transfer purified parts to a new well
        protocol.comment("Transferring purified parts")
        for col in columns_with_wells:
            tip_manager.get_multi_tip(start_col=col)
            mag_well = well_manager.get_mag_well(col, mag_plate)
            final_well = well_manager.get_final_well(col, dest_plate)
            
            pipette.transfer(elution_buffer_volume - ELUTION_DEAD_VOL, mag_well, final_well, blow_out=False)
            tip_manager.return_tips(col)

        # Disengage MagDeck
        MAGDECK.disengage()
        
        protocol.comment("Protocol complete")

    # for i in range(96*2):
    #     print(i)
        
    #     magbead(sample_number=i, ethanol_well=ethanol_well)
    
    magbead(sample_number=sample_number, ethanol_well=ethanol_well)
    # removed elution buffer well='A1', added that to where the function is defined
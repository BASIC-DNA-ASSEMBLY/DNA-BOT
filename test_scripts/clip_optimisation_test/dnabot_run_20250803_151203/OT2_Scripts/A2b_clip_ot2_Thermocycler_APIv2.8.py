from opentrons import protocol_api
import json
import time
from typing import List, Optional

# Rename to 'clip_template' and paste into 'template_ot2_scripts' folder in DNA-BOT to use

#metadata
metadata = {
     'apiLevel': '2.10',
     'protocolName': 'DNABOT: A2b CLIP Assembly with Thermocycler v2.10',
     'description': 'Implements linker ligation reactions using an opentrons OT-2, including the thermocycler module.'}

# Load CLIP data from JSON file
# This will be replaced by the parser with embedded JSON data
clips_dict = {
    "A7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "folA_1",
        "part_source_well": "E1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "cysM_1",
        "part_source_well": "G5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "serC_1",
        "part_source_well": "E4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gadB_1",
        "part_source_well": "G4",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "trxA_1",
        "part_source_well": "C8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ubiC_1",
        "part_source_well": "C6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G7": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "pabC_1",
        "part_source_well": "A4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H7": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "leuB_2",
        "part_source_well": "H1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H7",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "carA_1",
        "part_source_well": "A1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "serC_2",
        "part_source_well": "F3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "menA_1",
        "part_source_well": "G8",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "dadX_1",
        "part_source_well": "G3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "asnA_1",
        "part_source_well": "G7",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvA_1",
        "part_source_well": "A8",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G8": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tdcB_2",
        "part_source_well": "F6",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H8": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "gltD_1",
        "part_source_well": "A7",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H8",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "argI_1",
        "part_source_well": "D2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "asnB_1",
        "part_source_well": "A3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "tyrB_1",
        "part_source_well": "E9",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "C9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "aspC_1",
        "part_source_well": "E3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "D9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "trxC_1",
        "part_source_well": "A6",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "tdcB_1",
        "part_source_well": "A5",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G9": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "aspC_1",
        "part_source_well": "E3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "G9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H9": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "carA_1",
        "part_source_well": "A1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "H9",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "carB_1",
        "part_source_well": "C1",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "A10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "carB_1",
        "part_source_well": "C1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "B10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "panD_1",
        "part_source_well": "D4",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "serC_2",
        "part_source_well": "F3",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "aroL_1",
        "part_source_well": "G2",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "E10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "folM_1",
        "part_source_well": "G10",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "F10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "G10": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "folA_1",
        "part_source_well": "E1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "G10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "H10": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "dmlA_1",
        "part_source_well": "E10",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "H10",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "A11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "leuB_2",
        "part_source_well": "H1",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "A11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "B11": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ilvN_1",
        "part_source_well": "G11",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "B11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "C11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "ilvH_1",
        "part_source_well": "A2",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "C11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "D11": {
        "prefix_linker": "LMP-P",
        "prefix_source_well": "B12",
        "prefix_source_plate": "1",
        "part": "ydiB_1",
        "part_source_well": "A11",
        "part_source_plate": "1",
        "suffix_linker": "L1-S",
        "suffix_source_well": "G12",
        "suffix_source_plate": "1",
        "Clip_Well": "D11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "E11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "sucB_1",
        "part_source_well": "C3",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "E11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    },
    "F11": {
        "prefix_linker": "L1-P",
        "prefix_source_well": "F12",
        "prefix_source_plate": "1",
        "part": "folM_1",
        "part_source_well": "G10",
        "part_source_plate": "1",
        "suffix_linker": "LMS-S",
        "suffix_source_well": "E12",
        "suffix_source_plate": "1",
        "Clip_Well": "F11",
        "plate": 2,
        "part_vol": 1.0,
        "water_vol": 0.0
    }
}
all_default_conc = True

# Thermocycler generation setting
# This will be replaced by the parser with embedded thermocycler generation
thermocycler_gen = 'GEN2'

# opentrons_simulate.exe dnabot\template_ot2_scripts\clip_template_TC_APIv2.8.py --custom-labware-path 'labware\Labware definitions'

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
    
    def distribute_to_wells(self, wells: List[protocol_api.labware.Well], pipette: protocol_api.instrument_context.InstrumentContext, tip_manager=None) -> None:
        """Distribute master mix to a list of wells.
        
        Args:
            wells: List of wells to distribute to
            pipette: Pipette to use for distribution
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

# example dictionary produced by DNA-BOT for a single construct containing 5 parts, un-comment and run to test the template
#clips_dict={"prefixes_wells": ["A8", "A7", "C5", "C7", "C10"], "prefixes_plates": ["2", "2", "2", "2", "2"], "suffixes_wells": ["B7", "C1", "C2", "C3", "B8"], "suffixes_plates": ["2", "2", "2", "2", "2"], "parts_wells": ["E2", "F2", "C2", "B2", "D2"], "parts_plates": ["1", "1", "1", "1", "1"], "parts_vols": [1, 1, 1, 1, 1], "water_vols": [7.0, 7.0, 7.0, 7.0, 7.0]}

def run(protocol: protocol_api.ProtocolContext):
# added run function for API 2.8

    ### Constants - these have been moved out of the def clip() for clarity

    #Tiprack
    tiprack_type="opentrons_96_tiprack_20ul"
    INITIAL_TIP = 'A1'
    CANDIDATE_TIPRACK_SLOTS = ['3', '6', '9']  # Exclude slots 1, 2 for source plates and slots 7, 8, 10, 11 for thermocycler

    # Pipettes - pipette instructions in a single location so redefining pipette type is simpler
    PIPETTE_TYPE = 'p20_single_gen2'
    PIPETTE_MOUNT = 'right'
        ### Load Pipette
        # checks if it's a P10 Single pipette
    if PIPETTE_TYPE != 'p20_single_gen2':
        print('Define labware must be changed to use', PIPETTE_TYPE)
        exit()

    # Thermocycler Module
    if thermocycler_gen == 'GEN1':
        tc_mod = protocol.load_module('thermocycler', '7')
    else:  # GEN2
        tc_mod = protocol.load_module('thermocyclerModuleV2', '7')
        
    # Destination Plates
    DESTINATION_PLATE_TYPE = '4ti0960rig_96_wellplate_200ul'

    # Loads destination plate onto Thermocycler Module
    destination_plate = tc_mod.load_labware(DESTINATION_PLATE_TYPE)
    tc_mod.open_lid()
    tc_mod.set_block_temperature(20)

    # Source Plates - dynamically loaded based on embeddings
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
        
        # Load pipettes (without tip_racks argument - TipManager will handle this)
        pipette = protocol.load_instrument(PIPETTE_TYPE, mount=PIPETTE_MOUNT)
        
        # Multi-channel pipette for master mix distribution
        multi_pipette = protocol.load_instrument('p300_multi_gen2', 'left')
        
        # Initialize TipManager for p20 tips
        # Use available slots for p20 tip racks (excluding slots 1, 2 for source plates and slots 7, 8, 10, 11 for thermocycler)
        p20_tip_slots = ['3', '6', '9']  # Available slots for p20 tip racks
        p20_tip_manager = TipManager(protocol, int(p20_tip_slots[0]), pipette, 'p20')  # Initialize with first slot
        
        # Add additional tip racks if needed
        for i in range(1, min(tiprack_num, len(p20_tip_slots))):
            p20_tip_manager.add_tip_rack(int(p20_tip_slots[i]), 'p20')
        
        # Initialize TipManager for p300 tips (slot 5)
        p300_tip_manager = TipManager(protocol, 5, multi_pipette, 'p300')
        
        # Destination plate on thermocycler
        destination_wells = [destination_plate.wells_by_name()[w] for w in clips_dict.keys()]
        # Tube rack
        tube_rack = protocol.load_labware(TUBE_RACK_TYPE, TUBE_RACK_POSITION)
        water = tube_rack.wells(WATER_WELL)
        
        # Load source plates dynamically based on embeddings
        source_plates = {}
        all_plates = set()
        for well_info in clips_dict.values():
            all_plates.add(well_info['prefix_source_plate'])
            all_plates.add(well_info['suffix_source_plate'])
            all_plates.add(well_info['part_source_plate'])
        for key in all_plates:
            source_plates[key] = protocol.load_labware(SOURCE_PLATE_TYPE, key)
        dest_wells = list(clips_dict.keys())
        
        # Set master mix volume based on all_default_conc
        if all_default_conc:
            MASTER_MIX_VOLUME = 27  # Optimised master mix for default concentration
        
        # Master mix distribution using MasterMixManager
        master_mix_tube_volumes = [1500]  # Single tube is sufficient for CLIP reactions
        
        # Initialize MasterMixManager
        mm_manager = MasterMixManager(
            protocol,
            tube_rack,
            master_mix_tube_volumes,
            MASTER_MIX_VOLUME,
            multi_pipette,  # Use multi-channel pipette for master mix
            dead_volume=15.0
        )
        
        # Distribute master mix to all destination wells using multi-channel pipette
        mm_manager.distribute_to_wells(destination_wells, multi_pipette, tip_manager=p300_tip_manager)
        
        # Water transfer (only if needed and not all_default_conc)
        if not all_default_conc:
            water_vols = [clips_dict[w]['water_vol'] for w in dest_wells]
            if any([wv > 0 for wv in water_vols]):
                p20_tip_manager.get_single_tip()
                pipette.transfer(water_vols, water, destination_wells, blow_out=True, blowout_location='destination well', new_tip='never')
                pipette.drop_tip()
        
        # Prefix, suffix, part transfers using TipManager
        for i, well in enumerate(dest_wells):
            info = clips_dict[well]
            
            # Prefix transfer
            p20_tip_manager.get_single_tip()
            pipette.transfer(1, source_plates[info['prefix_source_plate']].wells_by_name()[info['prefix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='never', mix_after=LINKER_MIX_SETTINGS)
            pipette.drop_tip()
            
            # Suffix transfer
            p20_tip_manager.get_single_tip()
            pipette.transfer(1, source_plates[info['suffix_source_plate']].wells_by_name()[info['suffix_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='never', mix_after=LINKER_MIX_SETTINGS)
            pipette.drop_tip()
            
            # Part transfer
            p20_tip_manager.get_single_tip()
            pipette.transfer(info['part_vol'], source_plates[info['part_source_plate']].wells_by_name()[info['part_source_well']], destination_wells[i], blow_out=True, blowout_location='destination well', new_tip='never', mix_after=PART_MIX_SETTINGS)
            pipette.drop_tip()
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

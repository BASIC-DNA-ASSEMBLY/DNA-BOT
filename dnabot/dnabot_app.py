# -*- coding: utf-8 -*-
"""
DNA-BOT: DNA assembly using BASIC on OpenTrons

@author: mh2210

TO DO
    - add in new transformation protocol
    - ammend meta information
    - add robot comments 
    - add different file locations for different protocol versions
    - add clip input functionality
    - ammend instruction manuel
    - new GUI??
 
"""
from typing import List, Dict, Tuple, Union, Optional, Any
import os
import sys
from dataclasses import dataclass
from pathlib import Path

#add dnabot module to syspath
abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, abs_path)

import csv
import argparse
import pandas as pd
import numpy as np
import json
import tkinter as tk
import dnabot_gui as gui
import regex as re
import datetime


@dataclass
class ProtocolConfig:
    """Configuration for DNA assembly protocol parameters."""
    
    # CLIP reaction parameters
    CLIP_DEAD_VOL: int = 60
    CLIP_VOL: int = 30
    CLIP_T4_BUFF_VOL: float = 3.0
    CLIP_BSAI_VOL: float = 1.0
    CLIP_T4_LIG_VOL: float = 0.5
    CLIP_MAST_WATER: float = 15.5
    CLIP_PART_PER_CLIP: int = 200
    CLIP_MIN_VOL: float = 1.0
    CLIP_DEFAULT_PART_VOL: float = 1.0
    
    # Assembly parameters
    ASSEMBLY_MAX_CLIPS_PER_PLATE: int = 48
    ASSEMBLY_MAX_CLIPS_TOTAL: int = 96 * 2
    ASSEMBLY_MAX_ASSEMBLIES_PER_PLATE: int = 96
    ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP: int = 13
    ASSEMBLY_MAX_SOURCE_PLATES: int = 6
    ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS: int = 4
    ASSEMBLY_TIPS_PER_BOX: int = 96


@dataclass
class FileConfig:
    """Configuration for file paths and naming conventions."""
    
    # Directory structure
    TEMPLATE_DIR_NAME: str = 'template_ot2_scripts'
    V2_8_PATH: str = "\\APIv2.8\\"
    V2_8_TC_PATH: str = "\\Thermocycler_APIv2.8\\"
    
    # Template files
    TEMPLATE_FILES: Dict[str, Dict[str, str]] = None
    
    # Output files
    OUTPUT_FILES: Dict[str, Dict[str, str]] = None
    
    def __post_init__(self):
        if self.TEMPLATE_FILES is None:
            self.TEMPLATE_FILES = {
    'CLIP': {
        'V2_8': 'clip_template_APIv2.8.py',
        'V2_8_TC': 'clip_template_Thermocycler_module_APIv2.8.py'
    },
    'MAGBEAD': {
        'V2_8': 'purification_template_APIv2.8.py',
        'V2_10': 'purification_template_APIv2.10.py'
    },
    'F_ASSEMBLY': {
        'V2_8': 'assembly_template_APIv2.8.py',
        'V2_8_TC': 'assembly_template_Thermocycler_module_APIv2.8.py'
    },
    'TRANS_SPOT': {
        'V2_8': 'transformation_template_APIv2.8.py',
        'V2_8_TC': 'transformation_template_Thermocycler_module_APIv2.8.py'
    }
}

        if self.OUTPUT_FILES is None:
            self.OUTPUT_FILES = {
    'CLIP': {
        'V2_8': 'A_clip_ot2_APIv2.8',
        'V2_8_TC': 'A_clip_ot2_Thermocycler_APIv2.8'
    },
    'MAGBEAD': {
        'V2_8': 'B_purification_ot2_APIv2.8',
        'V2_10': 'B_purification_ot2_APIv2.10'
    },
    'F_ASSEMBLY': {
        'V2_8': 'C_assembly_ot2_APIv2.8',
        'V2_8_TC': 'C_assembly_ot2_Thermocycler_APIv2.8'
    },
    'TRANS_SPOT': {
        'V2_8': 'D_transformation_ot2_APIv2.8.py',
        'V2_8_TC': 'D_transformation_ot2_Thermocycler_APIv2.8.py'
    },
    'INFO': {
        'CLIPS': 'clip_run_info.csv',
        'FINAL_ASSEMBLIES': 'final_assembly_run_info.csv',
        'WELL_OUTPUT': 'wells.txt',
        'NEW_CONSTRUCTS': 'new_construct_list.csv',
        'ASSEMBLY_TO_CLIP_MAPPING': 'assembly_to_clip_mapping.csv'
    }
}


@dataclass
class DeckConfig:
    """Configuration for OT-2 deck layout."""
    
    SOURCE_POSITIONS: List[str] = None
    SPOTTING_VOLS: Dict[int, int] = None
    
    def __post_init__(self):
        if self.SOURCE_POSITIONS is None:
            self.SOURCE_POSITIONS = ['1', '2']  # NB for thermocycler protocols, the thermocycler takes up slots 7, 8, 10, 11
        if self.SPOTTING_VOLS is None:
            self.SPOTTING_VOLS = {2: 5, 3: 5, 4: 5, 5: 5, 6: 5, 7: 5}


@dataclass
class ClipReaction:
    """Represents a single CLIP reaction with its components."""
    prefix: str
    part: str
    suffix: str
    number: int = 1
    mag_well: Union[str, Tuple[str, ...]] = None
    plate: Union[int, Tuple[int, ...]] = None
    
    def __post_init__(self):
        if self.mag_well is None:
            self.mag_well = '0'
        if self.plate is None:
            self.plate = 0


@dataclass
class Construct:
    """Represents a DNA construct with its component parts."""
    name: str
    clips: List[ClipReaction]
    assembly_well: str = None
    assembly_plate: int = None


@dataclass
class SourceLocation:
    """Represents the location of a part/linker in a source plate."""
    well: str
    concentration: Optional[float] = None
    deck_position: str = None
    additional_info: List[str] = None
    
    def __post_init__(self):
        if self.additional_info is None:
            self.additional_info = []


@dataclass
class AssemblyPlan:
    """Represents the plan for final assembly of constructs."""
    destination_well: str
    source_clips: List[Tuple[str, int]]  # (well, plate) pairs
    construct_name: str = None


# Global configuration instances
PROTOCOL_CONFIG = ProtocolConfig()
FILE_CONFIG = FileConfig()
DECK_CONFIG = DeckConfig()


def __cli() -> argparse.Namespace:
    """Command line interface.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    desc = "DNA assembly using BASIC on OpenTrons."
    parser = argparse.ArgumentParser(description=desc)

    # Specific options for collecting settings from command line
    subparsers = parser.add_subparsers(help='Optional, to define settings from the terminal instead of the graphical '
                                            'interface. Type "python dnabot_app.py nogui -h" for more info.')
    parser_nogui = subparsers.add_parser('nogui')
    parser_nogui.add_argument('--construct_path', help='Construct CSV file.', required=True)
    parser_nogui.add_argument('--source_paths', help='Source CSV files.', nargs='+', required=True)
    parser_nogui.add_argument('--etoh_well', help='Well coordinate for Ethanol. Default: A3', default='A3', type=str)
    parser_nogui.add_argument('--soc_column', help='Column coordinate for SOC. Default: 1', default=1, type=int)
    parser_nogui.add_argument('--output_dir',
                              help='Output directory. Default: same directory than the one containing the '
                                   '"construct_path" file',
                              default=None, type=str)
    parser_nogui.add_argument('--template_dir',
                              help='Template directory. Default: "template_ot2_scripts" located next to the present '
                                   'script.',
                              default=None, type=str)
    parser_nogui.add_argument('--keep_layout',
                              help='Keep original CSV layout including empty rows. Default: True',
                              type=str, default='True', choices=['True', 'False'])
    
    parser.set_defaults(nogui=False)
    parser_nogui.set_defaults(nogui=True)
    return parser.parse_args()


def __info_from_gui() -> Dict[str, Union[str, List[str], int, bool]]:
    """Pop GUI to collect user inputs using a single persistent window.

    Returns:
        Dict[str, Union[str, List[str], int, bool]]: Dictionary containing user inputs
    """
    user_inputs = {
        'construct_path': None,
        'sources_paths': None,
        'etoh_well': None,
        'soc_column': None,
        'keep_layout': True  # Default to True
    }

    # Create a single persistent window
    print("Requesting user input, if not visible check minimised windows.")
    root = tk.Tk()
    root.title("DNA-BOT Configuration")
    
    # Configure window to stay on top and be modal
    root.lift()
    root.attributes('-topmost', True)
    root.after_idle(root.attributes, '-topmost', False)
    
    # Collect all inputs in the same window
    try:
        # First, get the main configuration (etoh_well, soc_column)
        dnabotinst = gui.DnabotApp(root)
        root.mainloop()
        
        if dnabotinst.quit_status:
            root.destroy()
            sys.exit("User specified 'QUIT' during app.")
        
        # Store the configuration values
        user_inputs['etoh_well'] = dnabotinst.etoh_well
        user_inputs['soc_column'] = dnabotinst.soc_column
        user_inputs['keep_layout'] = dnabotinst.keep_layout
        
        # Now get the construct file path
        user_inputs['construct_path'] = gui.UserDefinedPaths(root, 'Construct csv file').output
        
        # Finally get the source file paths
        user_inputs['sources_paths'] = gui.UserDefinedPaths(root, 'Sources csv files', multiple_files=True).output
        
    finally:
        # Clean up the window
        root.destroy()

    return user_inputs


def _resolve_template_dir(template_dir: str, construct_path: str) -> str:
    """Resolve the template directory path based on user input or defaults."""
    if template_dir is not None:
        print(f"Using provided template directory: {template_dir}")
        return template_dir
    print("Automatically deducing template directory from script location")
    script_path = os.path.abspath(__file__)
    return os.path.abspath(os.path.join(script_path, '..', FILE_CONFIG.TEMPLATE_DIR_NAME))


def _ensure_output_dir(output_dir: str) -> None:
    """Ensure the output directory exists and change to it."""
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    os.chdir(output_dir)


def _create_timestamped_output_dir(base_output_dir: str) -> str:
    """
    Create a timestamped subdirectory within the base output directory.
    
    Args:
        base_output_dir: Base directory where the timestamped folder will be created
        
    Returns:
        Path to the newly created timestamped directory
    """
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    timestamped_dir = os.path.join(base_output_dir, f"dnabot_run_{timestamp}")
    
    if not os.path.exists(timestamped_dir):
        print(f"Creating timestamped output directory: {timestamped_dir}")
        os.makedirs(timestamped_dir)
    
    return timestamped_dir


def _generate_clip_script_name(base_name: str, clip_plate: int) -> str:
    """
    Generate CLIP script name with new naming convention (1a, 1b, 2a, 2b, etc.).
    
    Args:
        base_name: Base name from FILE_CONFIG (e.g., 'A_clip_ot2_APIv2.8')
        clip_plate: 0-based clip plate index
        
    Returns:
        Script name with new convention (e.g., 'A1a_clip_ot2_APIv2.8.py')
    """
    # Calculate plate number and half (a or b)
    plate_num = (clip_plate // 2) + 1
    half_letter = 'a' if clip_plate % 2 == 0 else 'b'
    
    # Insert the number and letter after the stage letter
    # e.g., 'A_clip_ot2_APIv2.8' -> 'A1a_clip_ot2_APIv2.8'
    stage_letter = base_name[0]  # Extract the stage letter (A, B, C, D)
    rest_of_name = base_name[1:]  # Get the rest of the name
    
    return f"{stage_letter}{plate_num}{half_letter}{rest_of_name}.py"


def _generate_script_name_with_number(base_name: str, script_number: int) -> str:
    """
    Generate script name with number moved after the stage letter.
    
    Args:
        base_name: Base name from FILE_CONFIG (e.g., 'C_assembly_ot2_APIv2.8')
        script_number: Script number (1-based)
        
    Returns:
        Script name with new convention (e.g., 'C1_assembly_ot2_APIv2.8.py')
    """
    # Extract the stage letter and rest of the name
    stage_letter = base_name[0]  # Extract the stage letter (A, B, C, D)
    rest_of_name = base_name[1:]  # Get the rest of the name
    
    # Check if the name already ends with .py
    if rest_of_name.endswith('.py'):
        return f"{stage_letter}{script_number}{rest_of_name}"
    else:
        return f"{stage_letter}{script_number}{rest_of_name}.py"


def main() -> None:
    """
    Main function to run the DNA-BOT application.

    This function orchestrates the entire DNA assembly workflow:
    1. Collects user input (via GUI or CLI)
    2. Validates input files and directories
    3. Processes input files to generate intermediate data structures
    4. Generates OT-2 scripts for each protocol stage
    5. Outputs metadata and summary files

    The workflow follows this sequence:
    - CLIP reactions: Create DNA fragments with compatible ends
    - Purification: Clean up CLIP reactions using magnetic beads
    - Final Assembly: Combine purified fragments into final constructs
    - Transformation: Prepare constructs for bacterial transformation (optional)

    Raises:
        FileNotFoundError: If required input files or directories are missing
        ValueError: If input data is malformed or exceeds protocol limits
        Exception: For other processing errors
    """
    try:
        print("Starting DNA-BOT...")
        print("=" * 50)

        # -----------------------------
        # 1. Collect and validate user input
        # -----------------------------
        user_config = _collect_user_input()
        _validate_input_files(user_config)
        
        # -----------------------------
        # 2. Set up directories and paths
        # -----------------------------
        paths = _setup_directories(user_config)
        
        # -----------------------------
        # 3. Process input files and generate data structures
        # -----------------------------
        data_structures = _process_input_files(user_config, paths)
        
        # -----------------------------
        # 4. Generate OT-2 scripts for each protocol stage
        # -----------------------------
        _generate_ot2_scripts(data_structures, paths, user_config)
        
        # -----------------------------
        # 5. Output metadata and summary files
        # -----------------------------
        _output_metadata_files(data_structures, paths, user_config)

        print("=" * 50)
        print('DNA-BOT completed successfully!')

    except Exception as e:
        print(f"\nError: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)


def _collect_user_input() -> Dict[str, Union[str, List[str], int]]:
    """
    Collect user input either via CLI or GUI.
    
    Returns:
        Dictionary containing all user configuration parameters
    """
    args = __cli()
    if args.nogui:
        print("Running in CLI mode...")
        return {
            'etoh_well': args.etoh_well,
            'soc_column': args.soc_column,
            'construct_path': args.construct_path,
            'sources_paths': args.source_paths,
            'output_dir': args.output_dir,
            'template_dir': args.template_dir,
            'keep_layout': args.keep_layout == 'True'
        }
    else:
        print("Running in GUI mode...")
        return __info_from_gui()


def _validate_input_files(user_config: Dict[str, Union[str, List[str], int]]) -> None:
    """
    Validate that all required input files exist.
    
    Args:
        user_config: User configuration dictionary
        
    Raises:
        FileNotFoundError: If any required file is missing
    """
    construct_path = user_config['construct_path']
    sources_paths = user_config['sources_paths']
    
    if not os.path.exists(construct_path):
        raise FileNotFoundError(f"Construct file not found: {construct_path}")
    
    for source_path in sources_paths:
        if not os.path.exists(source_path):
            raise FileNotFoundError(f"Source file not found: {source_path}")


def _setup_directories(user_config: Dict[str, Union[str, List[str], int]]) -> Dict[str, str]:
    """
    Set up output and template directories.
    
    Args:
        user_config: User configuration dictionary
        
    Returns:
        Dictionary containing resolved directory paths
    """
    construct_path = user_config['construct_path']
    output_dir = user_config.get('output_dir')
    template_dir = user_config.get('template_dir')
    
    # Resolve base output directory
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(construct_path))
        print(f"Using construct file directory as base output: {output_dir}")
    
    # Create timestamped output directory
    timestamped_output_dir = _create_timestamped_output_dir(output_dir)
    
    # Resolve template directory
    template_dir_path = _resolve_template_dir(template_dir, construct_path)
    print(f"Template directory path: {template_dir_path}")
    
    if not os.path.exists(template_dir_path):
        raise FileNotFoundError(f"Template directory not found at: {template_dir_path}")
    
    # Change to timestamped output directory
    _ensure_output_dir(timestamped_output_dir)

    return {
        'output_dir': timestamped_output_dir,
        'base_output_dir': output_dir,
        'template_dir': template_dir_path,
        'construct_base': os.path.splitext(os.path.basename(construct_path))[0]
    }


def _process_input_files(user_config: Dict[str, Union[str, List[str], int]], 
                        paths: Dict[str, str]) -> Dict[str, Any]:
    """
    Process input files and generate intermediate data structures.
    
    Args:
        user_config: User configuration dictionary
        paths: Directory paths dictionary
        
    Returns:
        Dictionary containing all processed data structures
        
    Raises:
        ValueError: If input data is malformed or exceeds protocol limits
        FileNotFoundError: If required files cannot be found
    """
    print('Processing input files...')
    print("=" * 60)

    try:
        # Process constructs
        print("\n1. Loading and validating constructs...")
        constructs_dict = generate_constructs_list(user_config['construct_path'], user_config.get('keep_layout', True))
        print(f"✓ Generated {len(constructs_dict)} constructs")
        
        # Generate CLIP reactions first
        print("\n2. Generating CLIP reactions...")
        clips_df = generate_clips_df(constructs_dict)
        print(f"✓ Generated {len(clips_df)} unique clips")
        
        # Validate construct data (now with clips_df for accurate counting)
        validate_construct_data(constructs_dict, clips_df)
        print("✓ Construct data validation passed")
        
        # Validate CLIP data
        validate_clips_data(clips_df)
        print("✓ CLIP data validation passed")
        
        # Process source files
        print("\n3. Loading and validating source data...")
        sources_dict = generate_sources_dict(user_config['sources_paths'])
        print(f"✓ Processed {len(sources_dict)} sources")

        # Validate source data
        validate_sources_data(sources_dict)
        print("✓ Source data validation passed")
        
        # Validate component availability
        print("\n4. Validating component availability...")
        validate_components_availability(constructs_dict, sources_dict)
        print("✓ All required components are available")
        
        # Validate tip usage
        print("\n5. Validating tip usage...")
        validate_tip_usage(constructs_dict, clips_df)
        print("✓ Tip usage within limits")
        
        # Log processing summary
        log_processing_summary(constructs_dict, clips_df, sources_dict)

        print('\n6. Calculating OT-2 variables...')
        
        # Generate CLIP dictionaries for OT-2 scripts using optimised assignment
        clips_dict_list, assembly_to_clip_mapping, optimised_clips_df = generate_optimised_clips_dict_list(constructs_dict, sources_dict)
        print(f"✓ Generated {len(clips_dict_list)} clip plate(s)")

        # Calculate magbead sample distribution based on actual CLIP reactions needed for each optimised plate
        # We need to calculate the CLIP count for each subset of constructs that goes into each optimised plate
        magbead_sample_list = []
        
        # Get the assembly plates to understand the split
        assembly_plates = set()
        for (order, plate, well) in constructs_dict.keys():
            assembly_plates.add(plate)
        assembly_plates = sorted(assembly_plates)
        
        if len(assembly_plates) <= 1:
            # Single assembly plate - use total count
            total_clips = count_clips_for_constructs(constructs_dict)
            # Split into 96-sample plates
            full_plates = total_clips // 96
            remaining_samples = total_clips % 96
            magbead_sample_list = [96] * full_plates
            if remaining_samples > 0:
                magbead_sample_list.append(remaining_samples)
        else:
            # Multiple assembly plates - calculate for each half
            # Split constructs by assembly plates
            mid_point = (len(assembly_plates) + 1) // 2
            first_half_plates = assembly_plates[:mid_point]
            second_half_plates = assembly_plates[mid_point:]
            
            # Create construct subsets
            first_half = {}
            second_half = {}
            for (order, plate, well), construct_df in constructs_dict.items():
                if plate in first_half_plates:
                    first_half[(order, plate, well)] = construct_df
                else:
                    second_half[(order, plate, well)] = construct_df
            
            # Calculate CLIP counts for each half
            if first_half:
                first_half_clips = count_clips_for_constructs(first_half)
                # Split first half into 96-sample plates
                first_full_plates = first_half_clips // 96
                first_remaining = first_half_clips % 96
                magbead_sample_list.extend([96] * first_full_plates)
                if first_remaining > 0:
                    magbead_sample_list.append(first_remaining)
            
            if second_half:
                second_half_clips = count_clips_for_constructs(second_half)
                # Split second half into 96-sample plates
                second_full_plates = second_half_clips // 96
                second_remaining = second_half_clips % 96
                magbead_sample_list.extend([96] * second_full_plates)
                if second_remaining > 0:
                    magbead_sample_list.append(second_remaining)
        
        magbead_sample_number_total = sum(magbead_sample_list)
        print(f"✓ Total magbead samples: {magbead_sample_number_total}")
        
        # Generate final assembly plans
        final_assembly_dict_list = generate_final_assembly_dict_list(constructs_dict, clips_df)
        print(f"✓ Generated {len(final_assembly_dict_list)} assembly plate(s)")

        print("=" * 60)
        print("✓ All processing completed successfully!")

        return {
            'constructs_dict': constructs_dict,
            'clips_df': clips_df,
            'optimised_clips_df': optimised_clips_df,
            'sources_dict': sources_dict,
            'clips_dict_list': clips_dict_list,
            'magbead_sample_list': magbead_sample_list,
            'final_assembly_dict_list': final_assembly_dict_list,
            'assembly_to_clip_mapping': assembly_to_clip_mapping
        }
        
    except Exception as e:
        print("\n" + "=" * 60)
        print("❌ ERROR DURING PROCESSING")
        print("=" * 60)
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        raise


def _generate_ot2_scripts(data_structures: Dict[str, Any], 
                         paths: Dict[str, str],
                         user_config: Dict[str, Union[str, List[str], int]]) -> None:
    """
    Generate OT-2 scripts for each protocol stage.
    
    Args:
        data_structures: Processed data structures
        paths: Directory paths
        user_config: User configuration
    """
    print('Writing OT-2 scripts...')
    
    # Generate CLIP scripts
    for clip_plate, sub_clip_dict in enumerate(data_structures['clips_dict_list']):
        print(f"  Clip plate {clip_plate + 1}...")
        _generate_clip_scripts(sub_clip_dict, clip_plate, paths)

    # Generate magbead purification scripts
    for i, magbead_sample_number in enumerate(data_structures['magbead_sample_list']):
        print(f"  Magbead sample {i + 1}...")
        _generate_magbead_scripts(magbead_sample_number, i, paths, user_config)
    
    # Generate final assembly scripts
    for plate_number, final_assembly_dict in data_structures['final_assembly_dict_list'].items():
        print(f"  Assembly plate {plate_number}...")
        _generate_assembly_scripts(final_assembly_dict, plate_number, paths)


def _generate_clip_scripts(sub_clip_dict: Dict[str, List], 
                          clip_plate: int, 
                          paths: Dict[str, str]) -> None:
    """Generate CLIP reaction scripts for a single plate using Media Bot-style parameterisation."""
    template_dir = paths['template_dir']
    
    # Generate standard CLIP script with new naming convention
    clip_script_name = _generate_clip_script_name(FILE_CONFIG.OUTPUT_FILES['CLIP']['V2_8'], clip_plate)
    _generate_clip_script_media_bot_style(
        clip_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['CLIP']['V2_8']),
        sub_clip_dict
    )
    
    # Generate thermocycler CLIP script with new naming convention
    clip_tc_script_name = _generate_clip_script_name(FILE_CONFIG.OUTPUT_FILES['CLIP']['V2_8_TC'], clip_plate)
    _generate_clip_script_media_bot_style(
        clip_tc_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['CLIP']['V2_8_TC']),
        sub_clip_dict
    )


def _generate_magbead_scripts(magbead_sample_number: int, 
                             script_index: int, 
                             paths: Dict[str, str],
                             user_config: Dict[str, Union[str, List[str], int]]) -> None:
    """Generate magnetic bead purification scripts."""
    import re
    template_dir = paths['template_dir']
    
    # Generate script with new naming convention using v2.10 template
    magbead_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['MAGBEAD']['V2_10'], script_index + 1)
    
    # Convert NumPy types to native Python types to avoid JSON serialisation issues
    clips_number = int(magbead_sample_number) if hasattr(magbead_sample_number, 'item') else magbead_sample_number
    print(f"Clips number: {clips_number}")
    
    template_path = os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['MAGBEAD']['V2_10'])
    with open(template_path, 'r') as f:
        template_content = f.read()
    
    # Replace only the specific lines that need to be dynamic
    # Replace clips_number line
    template_content = re.sub(
        r"'clips_number': \d+,",
        f"'clips_number': {clips_number},",
        template_content
    )
    
    # Replace ethanol_well line
    ethanol_well = user_config.get('etoh_well', 'A3')
    template_content = re.sub(
        r"'ethanol_well': '[^']*',",
        f"'ethanol_well': '{ethanol_well}',",
        template_content
    )
    
    # Replace elution_well line (if it needs to be different from default)
    elution_well = user_config.get('elution_well', 'A10')
    template_content = re.sub(
        r"'elution_well': '[^']*',",
        f"'elution_well': '{elution_well}',",
        template_content
    )
    
    # Write the modified template to the output file
    with open(magbead_script_name, 'w') as f:
        f.write(template_content)
    
    print(f"    ✓ {os.path.basename(magbead_script_name)}")


def _generate_assembly_scripts(final_assembly_dict: Dict[str, List], 
                              plate_number: int, 
                              paths: Dict[str, str]) -> None:
    """Generate final assembly scripts using Media Bot-style parameterisation."""
    template_dir = paths['template_dir']
    final_assembly_tipracks = calculate_final_assembly_tipracks(final_assembly_dict)
    
    # Generate standard assembly script with new naming convention
    assembly_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['F_ASSEMBLY']['V2_8'], plate_number)
    _generate_assembly_script_media_bot_style(
        assembly_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['F_ASSEMBLY']['V2_8']),
        final_assembly_dict,
        final_assembly_tipracks
    )
    
    # Generate thermocycler assembly script with new naming convention
    assembly_tc_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['F_ASSEMBLY']['V2_8_TC'], plate_number)
    _generate_assembly_script_media_bot_style(
        assembly_tc_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['F_ASSEMBLY']['V2_8_TC']),
        final_assembly_dict,
        final_assembly_tipracks
    )


def _output_metadata_files(data_structures: Dict[str, Any], 
                          paths: Dict[str, str],
                          user_config: Dict[str, Union[str, List[str], int]]) -> None:
    """
    Output metadata and summary files.
    
    Args:
        data_structures: Processed data structures
        paths: Directory paths
        user_config: User configuration
    """
    print('Writing metadata files...')
    
    # Create metainformation directory
    metainfo_dir = os.path.join(paths['output_dir'], 'metainformation')
    if not os.path.exists(metainfo_dir):
        os.makedirs(metainfo_dir)
    
    # Change to metainformation directory for file creation
    original_dir = os.getcwd()
    os.chdir(metainfo_dir)
    
    try:
        # Generate master mix information
        master_mix_df = generate_master_mix_df(data_structures['clips_df']['number'].sum())
        
        # Generate source plate information
        sources_paths_df = generate_sources_paths_df(
            user_config['sources_paths'], 
            DECK_CONFIG.SOURCE_POSITIONS
        )
        
        # Write CLIP run information
        # Use optimised clip data if available, otherwise use original clips_df
        if 'optimised_clips_df' in data_structures:
            clip_reactions_df = data_structures['optimised_clips_df']
        else:
            clip_reactions_df = data_structures['clips_df']
            
        dfs_to_csv(
            paths['construct_base'] + '_' + FILE_CONFIG.OUTPUT_FILES['INFO']['CLIPS'],
            index=False,
            MASTER_MIX=master_mix_df,
            SOURCE_PLATES=sources_paths_df,
            CLIP_REACTIONS=clip_reactions_df
        )
        
        # Write final assembly information
        _write_final_assembly_info(data_structures, paths)
        
        # Write well output information
        _write_well_output_info(user_config, paths)
        
        # Write assembly to clip mapping information
        _write_assembly_to_clip_mapping(data_structures.get('assembly_to_clip_mapping', {}), paths)
        
        # Write new constructs information
        new_constructs_df = generate_new_constructs_df(
            user_config['construct_path'], 
            data_structures['final_assembly_dict_list'],
            data_structures['constructs_dict']
        )
        new_constructs_df.to_csv(
            paths['construct_base'] + '_' + FILE_CONFIG.OUTPUT_FILES['INFO']['NEW_CONSTRUCTS'],
            index=False
        )
        
    finally:
        # Restore original directory
        os.chdir(original_dir)


def _write_final_assembly_info(data_structures: Dict[str, Any], paths: Dict[str, str]) -> None:
    """Write final assembly information to CSV file."""
    with open(paths['construct_base'] + '_' + FILE_CONFIG.OUTPUT_FILES['INFO']['FINAL_ASSEMBLIES'],
                  'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for plate_number, final_assembly_dict in data_structures['final_assembly_dict_list'].items():
            for well, construct_clips in final_assembly_dict.items():
                csvwriter.writerow([plate_number, well, construct_clips])


def _write_well_output_info(user_config: Dict[str, Union[str, List[str], int]], 
                           paths: Dict[str, str]) -> None:
    """Write well output information to text file."""
    with open(paths['construct_base'] + '_' + FILE_CONFIG.OUTPUT_FILES['INFO']['WELL_OUTPUT'], 'w') as f:
        f.write(f'Magbead ethanol well: {user_config["etoh_well"]}\n')
        f.write(f'SOC column: {user_config["soc_column"]}')


def _write_assembly_to_clip_mapping(assembly_to_clip_mapping: Dict[int, List[int]], 
                                   paths: Dict[str, str]) -> None:
    """Write assembly to clip mapping information to CSV file."""
    if not assembly_to_clip_mapping:
        return
    
    # Create DataFrame for mapping
    mapping_data = []
    for assembly_plate, clip_plates in assembly_to_clip_mapping.items():
        # Extract magbead plate numbers from clip plate names
        # e.g., "1a, 1b" -> "1", "2a, 2b" -> "2"
        magbead_plates = set()
        for clip_plate in clip_plates:
            # Extract the number before the letter (e.g., "1a" -> "1")
            if isinstance(clip_plate, str) and len(clip_plate) > 0:
                # Find the first non-digit character
                plate_num = ""
                for char in clip_plate:
                    if char.isdigit():
                        plate_num += char
                    else:
                        break
                if plate_num:
                    magbead_plates.add(int(plate_num))
            else:
                # Fallback for numeric clip plates
                magbead_plates.add(clip_plate)
        
        mapping_data.append({
            'Assembly_Plate': assembly_plate,
            'Required_Clip_Plates': ', '.join(map(str, clip_plates)),
            'Required_Magbead_Plates': ', '.join(map(str, sorted(magbead_plates))),
            'Number_of_Clip_Plates': len(magbead_plates)
        })
    
    mapping_df = pd.DataFrame(mapping_data)
    mapping_df.to_csv(
        paths['construct_base'] + '_assembly_to_clip_mapping.csv',
        index=False
    )
    
    print(f"✓ Assembly to clip mapping saved to: {paths['construct_base']}_assembly_to_clip_mapping.csv")


def generate_constructs_list(path: str, keep_layout: bool = True) -> Dict[Tuple[int, int, str], pd.DataFrame]:
    """Generates a dictionary mapping construct positions to their dataframes.
    
    Each dataframe lists components of the CLIP reactions required for that construct,
    including prefix linkers, parts, and suffix linkers.

    Args:
        path (str): Path to the constructs CSV file
        keep_layout (bool): If True, preserve empty rows and associate constructs with their well positions.
                           If False, filter out empty rows and process constructs in order.

    Returns:
        Dict[Tuple[int, int, str], pd.DataFrame]: Dictionary where keys are (order, plate, well) tuples and values are
            dataframes containing CLIP reaction components for each construct. Each dataframe has columns:
            - prefixes: Prefix linker identifiers
            - parts: Part identifiers
            - suffixes: Suffix linker identifiers

    Raises:
        FileNotFoundError: If the constructs file cannot be found
        ValueError: If the constructs file is malformed
    """
    print(f"\n... Loading constructs from: {path}")
    
    def process_construct(construct: List[str], construct_index: int) -> pd.DataFrame:
        """Processes an individual construct into a dataframe of CLIP reactions.
        
        Args:
            construct (List[str]): List of parts and linkers for a single construct
                in the format [linker1, part1, linker2, part2, ...]
            construct_index (int): Index of the construct for error reporting

        Returns:
            pd.DataFrame: DataFrame containing CLIP reaction components with columns:
                - prefixes: Prefix linker identifiers
                - parts: Part identifiers
                - suffixes: Suffix linker identifiers
        """
        def get_suffix_linker(linker: str) -> str:
            """Determines the suffix linker identifier from a linker.
            
            If the linker starts with 'U', it's a UTR linker and gets '-S' suffix.
            Otherwise, just appends '-S' to the linker identifier.

            Args:
                linker (str): Linker identifier

            Returns:
                str: Suffix linker identifier
            """
            if linker.startswith('U'):
                return linker.split('-')[0] + '-S'
            return linker + "-S"

        # Validate construct structure
        if len(construct) < 2:
            raise ValueError(f"Construct {construct_index + 1} has insufficient components. "
                           f"Expected at least 2 (linker, part), got {len(construct)}")
        
        if len(construct) % 2 != 0:
            raise ValueError(f"Construct {construct_index + 1} has invalid structure. "
                           f"Expected an even number of components (linker, part, linker, part, ...), got {len(construct)}")
        
        # Check that the number of linkers equals the number of parts
        linkers = construct[::2]
        parts = construct[1::2]
        if len(linkers) != len(parts):
            raise ValueError(f"Construct {construct_index + 1} has mismatched number of linkers and parts. "
                           f"Linkers: {len(linkers)}, Parts: {len(parts)}")

        # Initialize dictionary to store CLIP reaction components
        clips_info = {
            'prefixes': [],
            'parts': [],
            'suffixes': []
        }

        # Process each part and its surrounding linkers
        for i, sequence in enumerate(construct):
            if i % 2 != 0:  # Only process parts (odd indices)
                # Validate part is not empty
                if not sequence.strip():
                    raise ValueError(f"Construct {construct_index + 1}, position {i}: Empty part found")
                
                clips_info['parts'].append(sequence.strip())
                
                # Validate prefix linker
                prefix_linker = construct[i - 1].strip()
                if not prefix_linker:
                    raise ValueError(f"Construct {construct_index + 1}, position {i-1}: Empty prefix linker found")
                clips_info['prefixes'].append(prefix_linker + '-P')
                
                # Determine suffix linker
                if i == len(construct) - 1:
                    # Last part uses first linker as suffix
                    suffix_linker = get_suffix_linker(construct[0].strip())
                else:
                    # Use next linker as suffix
                    next_linker = construct[i + 1].strip()
                    if not next_linker:
                        raise ValueError(f"Construct {construct_index + 1}, position {i+1}: Empty suffix linker found")
                    suffix_linker = get_suffix_linker(next_linker)
                clips_info['suffixes'].append(suffix_linker)

        return pd.DataFrame.from_dict(clips_info)

    # Process each construct in the CSV file
    constructs_dict = {}
    valid_construct_index = 0  # Counter for valid constructs only
    
    try:
        with open(path, 'r') as csvfile:
            csv_reader = csv.reader(csvfile)
            
            # First, determine the format by reading the header
            header_row = None
            csv_reader_list = list(csv.reader(open(path, 'r')))
            if csv_reader_list:
                header_row = csv_reader_list[0]
            
            # Determine if this is the new format (Plate, Well, ...) or old format (Well, ...)
            has_plate_column = header_row and len(header_row) >= 2 and header_row[0].strip().lower() == 'plate'
            
            # Validate CSV structure for new format
            if has_plate_column:
                if len(header_row) < 2:
                    raise ValueError(f"CSV file appears to use new format but has insufficient columns. "
                                   f"Expected at least 2 columns (Plate, Well), got {len(header_row)}")
                
                first_col = header_row[0].strip().lower()
                second_col = header_row[1].strip().lower()
                
                if first_col != 'plate':
                    raise ValueError(f"CSV file appears to use new format but first column is '{header_row[0]}' instead of 'Plate'")
                
                if second_col != 'well':
                    raise ValueError(f"CSV file appears to use new format but second column is '{header_row[1]}' instead of 'Well'")
                
                print(f"... Validated CSV structure: Plate, Well format detected")
            else:
                # Validate old format structure
                if len(header_row) < 1:
                    raise ValueError(f"CSV file has insufficient columns. "
                                   f"Expected at least 1 column (Well), got {len(header_row)}")
                
                first_col = header_row[0].strip().lower()
                if first_col != 'well':
                    raise ValueError(f"CSV file appears to use old format but first column is '{header_row[0]}' instead of 'Well'")
                
                print(f"... Validated CSV structure: Well format detected")
            
            for index, row in enumerate(csv_reader_list):
                if index == 0:  # Skip header row
                    continue
                
                # Use the format determined from header
                if has_plate_column:
                    # New format: Plate, Well, Linker1, Part1, ...
                    plate_str = row[0].strip() if row[0] else "1"
                    well_position = row[1] if len(row) > 1 else ""
                    construct_components = row[2:]  # Start from column 3 (after Plate, Well)
                    
                    # Validate plate number
                    try:
                        plate_number = int(plate_str)
                        if plate_number < 1:
                            raise ValueError(f"Plate number must be a positive integer, got {plate_number} at row {index + 1}")
                    except ValueError as e:
                        if "invalid literal" in str(e):
                            raise ValueError(f"Plate number must be a valid integer, got '{plate_str}' at row {index + 1}")
                        else:
                            raise ValueError(f"Plate number error at row {index + 1}: {str(e)}")
                else:
                    # Old format: Well, Linker1, Part1, ...
                    plate_number = 1  # Default to plate 1
                    well_position = row[0] if row else ""
                    construct_components = row[1:]  # Start from column 2
                
                # Validate well position if not empty
                if well_position.strip():
                    try:
                        validate_well_format(well_position)
                    except ValueError as e:
                        raise ValueError(f"Well position error at row {index + 1}: {str(e)}")
                
                # Filter out empty strings
                construct_components = list(filter(None, construct_components))
                
                # Skip rows with no construct components (empty rows)
                if not construct_components:
                    continue
                
                if keep_layout:
                    # In keep_layout mode, use plate and well from CSV
                    try:
                        construct_df = process_construct(construct_components, index)
                        # Use plate number from CSV and well position from CSV
                        order = valid_construct_index
                        constructs_dict[(order, plate_number, well_position)] = construct_df
                        valid_construct_index += 1
                    except ValueError as e:
                        raise ValueError(f"Error processing construct at Plate {plate_number}, Well {well_position} (row {index + 1}): {str(e)}")
                else:
                    # When keep_layout=False, assign new plate and well positions based on order
                    try:
                        construct_df = process_construct(construct_components, index)
                        # Calculate new plate and well based on order
                        new_plate = (valid_construct_index // 96) + 1  # 96 wells per plate
                        new_well = tip_counter(valid_construct_index % 96)  # Convert index to well coordinate
                        order = valid_construct_index
                        constructs_dict[(order, new_plate, new_well)] = construct_df
                        valid_construct_index += 1
                    except ValueError as e:
                        raise ValueError(f"Error processing construct at row {index + 1}: {str(e)}")
        
        if keep_layout:
            print(f"... Successfully loaded {len(constructs_dict)} constructs with layout preserved")
        else:
            print(f"... Successfully loaded {len(constructs_dict)} constructs")
        
        return constructs_dict
        
    except FileNotFoundError:
        raise FileNotFoundError(f"Constructs file not found: {path}")
    except Exception as e:
        print(f"[ERROR] Failed to load constructs: {str(e)}")
        raise


def generate_sources_dict(paths: List[str]) -> Dict[str, Tuple[str, ...]]:
    """Creates a dictionary mapping parts/linkers to their source locations.
    
    Args:
        paths (List[str]): List of paths to source CSV files containing part/linker
            information. Each file should have columns:
            - Part/linker identifier (first column)
            - Well location (second column)
            - Concentration (optional, third column)
            - Other information (optional, additional columns)

    Returns:
        Dict[str, Tuple[str, ...]]: Dictionary where:
            - Keys are part/linker identifiers
            - Values are tuples containing:
                - Well location
                - Concentration (if provided)
                - Deck position
                - Other information from CSV

    Raises:
        FileNotFoundError: If any source file cannot be found
        ValueError: If source files are malformed
    """
    sources_dict = {}
    
    for deck_index, path in enumerate(paths):
        print(f"\n... Loading source data from: {path}")
        
        try:
            with open(path, 'r') as csvfile:
                csv_reader = csv.reader(csvfile)
                
                for index, source in enumerate(csv_reader):
                    if index == 0:  # Skip header row if present
                        continue
                    
                    # Strip whitespace from part name
                    part_name = str(source[0]).strip()
                    
                    # Skip rows with empty part names (blank wells)
                    if not part_name:
                        continue
                    
                    # Validate row has minimum required data
                    if len(source) < 2:
                        raise ValueError(f"Row {index + 1}: Insufficient data. Need at least part name and well location.")
                    
                    # Extract values and add deck position
                    csv_values = source[1:]
                    csv_values.append(DECK_CONFIG.SOURCE_POSITIONS[deck_index])
                    
                    # Validate well format (only for non-empty wells)
                    well = source[1] if len(source) > 1 else ""
                    if well.strip():  # Only validate if well is not empty
                        try:
                            validate_well_format(well)
                        except ValueError as e:
                            raise ValueError(f"Row {index + 1}, part '{part_name}': {str(e)}")
                    
                    # Check for duplicate parts
                    if part_name in sources_dict:
                        raise ValueError(f"Duplicate part '{part_name}' found in source files. "
                                       f"First occurrence in file {paths.index(path) + 1}, "
                                       f"second occurrence in file {deck_index + 1}")
                    
                    sources_dict[part_name] = tuple(csv_values)
            
            print(f"... Successfully loaded {len([k for k in sources_dict.keys() if k in sources_dict])} parts from file {deck_index + 1}")
            
        except FileNotFoundError:
            raise FileNotFoundError(f"Source file not found: {path}")
        except Exception as e:
            print(f"[ERROR] Failed to load source file {path}: {str(e)}")
            raise
    
    print(f"\n... Total parts loaded: {len(sources_dict)}")
    return sources_dict


def generate_clips_df(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame]) -> pd.DataFrame:
    """Generates a dataframe containing information about all unique CLIP reactions.
    
    For each unique CLIP reaction, calculates:
    - Number of times it's needed
    - Magnetic bead well locations
    - CLIP plate assignments
    
    Args:
        constructs_dict (Dict[Tuple[int, int, str], pd.DataFrame]): Dictionary mapping construct positions 
            to dataframes containing CLIP reaction components for each construct

    Returns:
        pd.DataFrame: DataFrame containing all unique CLIP reactions with columns:
            - prefixes: Prefix linker identifiers
            - parts: Part identifiers
            - suffixes: Suffix linker identifiers
            - number: Number of times this CLIP reaction is needed
            - mag_well: Tuple of magnetic bead well locations
            - plate: Tuple of CLIP plate numbers

    Raises:
        ValueError: If total number of CLIP reactions exceeds maximum allowed
    """
    def count_unique_clips(clips_df: pd.DataFrame, merged_construct_dfs: pd.DataFrame) -> pd.DataFrame:
        """Counts how many times each unique CLIP reaction is needed.
        
        Args:
            clips_df (pd.DataFrame): DataFrame of unique CLIP reactions
            merged_construct_dfs (pd.DataFrame): DataFrame of all CLIP reactions
            
        Returns:
            pd.DataFrame: DataFrame with added 'number' column indicating count
        """
        clip_count = np.zeros(len(clips_df.index))
        for i, unique_clip in clips_df.iterrows():
            for _, clip in merged_construct_dfs.iterrows():
                if unique_clip.equals(clip):
                    clip_count[i] += 1
                    
        # Calculate number of reactions needed based on final assemblies per CLIP
        clip_count = clip_count // PROTOCOL_CONFIG.ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP + 1
        
        # Add the number column to the DataFrame
        clips_df['number'] = [int(i) for i in clip_count.tolist()]

        return clips_df

    # Merge all constructs and find unique CLIP reactions
    # Extract all construct dataframes from the dictionary
    valid_constructs = list(constructs_dict.values())
    merged_construct_dfs = pd.concat(valid_constructs, ignore_index=True)
    unique_clips_df = merged_construct_dfs.drop_duplicates().reset_index(drop=True)
    clips_df = unique_clips_df.copy()
    
    # Count occurrences of each unique CLIP reaction
    clips_df = count_unique_clips(unique_clips_df, merged_construct_dfs)
    
    # Add columns for well and plate locations
    clips_df['mag_well'] = pd.Series(['0'] * len(clips_df.index), index=clips_df.index)
    clips_df['plate'] = pd.Series(['0'] * len(clips_df.index), index=clips_df.index)
    
    # Assign well and plate locations for each CLIP reaction
    clip_count = 0
    for unique_clip_count, clip_number in clips_df['number'].items():
        mag_wells = []
        plates = []         
        for well in range(clip_count, clip_count + clip_number):
            mag_wells.append(tip_counter(well % 96))
            plates.append(1 + well//96)
        clips_df.at[unique_clip_count, 'mag_well'] = tuple(mag_wells)
        clips_df.at[unique_clip_count, 'plate'] = tuple(plates)
        clip_count += clip_number

    return clips_df


def generate_clips_dict(clips_df: pd.DataFrame, sources_dict: Dict[str, Tuple[str, ...]]) -> Dict[str, List]:
    """Generates dictionary of CLIP reaction information for OT-2 script.
    
    Args:
        clips_df (pd.DataFrame): DataFrame containing CLIP reaction information
        sources_dict (Dict[str, Tuple[str, ...]]): Dictionary mapping parts/linkers
            to their source locations

    Returns:
        Dict[str, List]: Dictionary containing:
            - prefixes_wells: List of prefix linker well locations
            - prefixes_plates: List of prefix linker plate numbers
            - suffixes_wells: List of suffix linker well locations
            - suffixes_plates: List of suffix linker plate numbers
            - parts_wells: List of part well locations
            - parts_plates: List of part plate numbers
            - parts_vols: List of part volumes
            - water_vols: List of water volumes

    Raises:
        ValueError: If required parts/linkers are missing from sources_dict
    """
    # Calculate maximum part volume based on total reaction volume
    max_part_vol = PROTOCOL_CONFIG.CLIP_VOL - (
        PROTOCOL_CONFIG.CLIP_T4_BUFF_VOL + 
        PROTOCOL_CONFIG.CLIP_BSAI_VOL + 
        PROTOCOL_CONFIG.CLIP_T4_LIG_VOL + 
        PROTOCOL_CONFIG.CLIP_MAST_WATER + 2
    )

    # Initialize dictionary for CLIP reaction information
    clips_dict = {
        'prefixes_wells': [],
        'prefixes_plates': [],
        'suffixes_wells': [],
        'suffixes_plates': [],
        'parts_wells': [],
        'parts_plates': [],
        'parts_vols': [],
        'water_vols': []
    }

    # Check for missing parts in sources_dict
    missing_parts = []
    for _, clip_info in clips_df.iterrows():
        prefix_linker = clip_info['prefixes'].strip()
        suffix_linker = clip_info['suffixes'].strip()
        part = clip_info['parts'].strip()
        
        if prefix_linker not in sources_dict:
            missing_parts.append(f"Prefix linker: {prefix_linker}")
        if suffix_linker not in sources_dict:
            missing_parts.append(f"Suffix linker: {suffix_linker}")
        if part not in sources_dict:
            missing_parts.append(f"Part: {part}")
    
    if missing_parts:
        error_msg = "The following parts/linkers are in the constructs but not in the parts plate:\n"
        error_msg += "\n".join(missing_parts)
        raise ValueError(error_msg)

    try:
        # Generate CLIP reaction information
        for _, clip_info in clips_df.iterrows():
            # Process prefix linker
            prefix_linker = clip_info['prefixes'].strip()
            clips_dict['prefixes_wells'].append([sources_dict[prefix_linker][0]] * clip_info['number'])
            clips_dict['prefixes_plates'].append(
                [normalize_source_data(sources_dict[prefix_linker])[2]] * clip_info['number'])
            
            # Process suffix linker
            suffix_linker = clip_info['suffixes'].strip()
            clips_dict['suffixes_wells'].append([sources_dict[suffix_linker][0]] * clip_info['number'])
            clips_dict['suffixes_plates'].append(
                [normalize_source_data(sources_dict[suffix_linker])[2]] * clip_info['number'])
            
            # Process part
            part = clip_info['parts'].strip()
            clips_dict['parts_wells'].append([sources_dict[part][0]] * clip_info['number'])
            clips_dict['parts_plates'].append(
                [normalize_source_data(sources_dict[part])[2]] * clip_info['number'])
            
            # Calculate part and water volumes
            if not sources_dict[part][1]:  # No concentration specified
                clips_dict['parts_vols'].append(
                    [float(PROTOCOL_CONFIG.CLIP_DEFAULT_PART_VOL)] * clip_info['number'])
                clips_dict['water_vols'].append(
                    [float(max_part_vol - PROTOCOL_CONFIG.CLIP_DEFAULT_PART_VOL)] * clip_info['number'])
            else:  # Use specified concentration
                part_vol = round(
                    PROTOCOL_CONFIG.CLIP_PART_PER_CLIP / float(sources_dict[part][1]), 1)
                part_vol = max(PROTOCOL_CONFIG.CLIP_MIN_VOL,
                             min(part_vol, max_part_vol))
                water_vol = max_part_vol - part_vol
                clips_dict['parts_vols'].append([float(part_vol)] * clip_info['number'])
                clips_dict['water_vols'].append([float(water_vol)] * clip_info['number'])
                    
        # Flatten nested lists
        for key, value in clips_dict.items():
            clips_dict[key] = [item for sublist in value for item in sublist]

        return clips_dict
        
    except Exception as e:
        print(f"\nError generating clips dictionary:")
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        raise


def generate_clips_dict_list(clips_df, sources_dict):
    '''Subsets the clips df into chunks of 48, runs the generate_clips_dict function 
    for each and returns a list of the resulting sub clips dicts'''

    def clip_df_long_format(clip_df):
        ''' Takes clips df and returns a long format df with one row per mag well 
        rather than one row per unique clip'''

        df = clip_df.copy()
        long_clip_df = pd.DataFrame(columns = df.columns)

        index = 0

        for row, clip in clip_df.iterrows():
            for clip_rep in range(clip.number):                         # copy row and replace mag well and plate lists with well and plate for current rep
                data = clip.copy()
                data.number = 1
                data.mag_well = clip.mag_well[clip_rep]
                data.plate = clip.plate[clip_rep]
                long_clip_df.loc[index] = data
                index += 1

        return long_clip_df
    
    long_clip_df = clip_df_long_format(clips_df)

    CLIP_COUNT = clips_df['number'].sum()
    CLIP_PLATE_COUNT = int(CLIP_COUNT // PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE + 1)       # plus one to include final partially full plate

    # Error 
    if clips_df['number'].sum() > PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL:
        raise ValueError('Number of CLIP reactions exceeds {}. Reduce number of constructs in construct.csv.'.format(PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL))
    
    if clips_df['number'].sum() < 96:                                   # if less clips required than one full plate, the second desk slot can be used for tips
        PROTOCOL_CONFIG.ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS = 5

    clips_dict_list = []

    for plate in range(CLIP_PLATE_COUNT):
        subset_lower = (plate * PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE)                    # set upper and lower bounds for subset of clips for a given plate
        subset_upper = subset_lower + PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE

        if subset_upper > CLIP_COUNT:                                   # set total number number of clips as upper bound if plate incomplete
            subset_upper = CLIP_COUNT
    
        sub_clip_df = long_clip_df.iloc[subset_lower:subset_upper, :]
        sub_clip_dict = generate_clips_dict(sub_clip_df, sources_dict)
        clips_dict_list.append(sub_clip_dict)                           # generate and append sub_clip_dict to list - allows for multiple clip reactions

    return clips_dict_list


def generate_final_assembly_dict(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], clips_df):
    """Using constructs_dict and clips_df, returns keys and values for a 
    dictionary of final assemblies; with keys defining destination plate 
    well positions, and values indicating which clip reaction wells are used.

    """

    final_assembly_dict_keys = []
    final_assembly_dict_values = []

    clips_count = np.zeros(len(clips_df.index))

    # Process constructs in the order they appear in the CSV (using the order field in the tuple)
    # This ensures CLIP assignment order matches the original system
    sorted_positions = sorted(constructs_dict.keys(), key=lambda x: x[0])  # Sort by order (first element)
    
    for position in sorted_positions:
        construct_df = constructs_dict[position]
        construct_well_list = []
        construct_plate_list = []

        for _, clip in construct_df.iterrows():                                     # for each clip in construct
            clip_info = clips_df[(clips_df['prefixes'] == clip['prefixes']) &       # find clips in clips_df with the correct parts required for this clip
                                 (clips_df['parts'] == clip['parts']) &
                                 (clips_df['suffixes'] == clip['suffixes'])]
            
            clip_num = int(clip_info.index[0])                                      # row index of this clip in clips_df
            clip_wells = clip_info.at[clip_num, 'mag_well']                         # list of all mag_wells for this clip
            clip_plates = clip_info.at[clip_num, 'plate']                           # list of all plates for this clip

            chosen_well_index = int(clips_count[clip_num] // PROTOCOL_CONFIG.ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP)     # next viable mag well for this clip (i.e. the nth well in the set of total number of wells for this unique clip)
            if chosen_well_index >= len(clip_wells):
                raise IndexError(f"Chosen well index {chosen_well_index} out of range for clip wells {clip_wells} (clip_num={clip_num})")
            clip_well = clip_wells[chosen_well_index]
            clip_plate = clip_plates[chosen_well_index]
            construct_well_list.append(clip_well)
            construct_plate_list.append(clip_plate)

            clips_count[clip_num] = clips_count[clip_num] + 1

        # Use the plate and well from the position tuple as the destination key
        order, plate, well = position
        destination_key = (plate, well)  # Use tuple of (plate, well) as key
        final_assembly_dict_keys.append(destination_key)
        final_assembly_dict_values.append([construct_well_list, construct_plate_list])

    return final_assembly_dict_keys, final_assembly_dict_values                     # return list of dict keys and values


def generate_final_assembly_dict_list(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                     clips_df: pd.DataFrame) -> Dict[int, Dict[str, List]]:
    """
    Generate a list of assembly dictionaries, each representing a subset of constructs
    that can be assembled on a single plate while respecting tip and well constraints.
    
    This function ensures that:
    1. No more than MAX_ASSEMBLIES_PER_PLATE constructs are assembled per plate
    2. Tip usage doesn't exceed available tiprack capacity
    3. Each construct gets the correct destination well location
    4. CLIP wells are not reused across different assembly plates
    
    Args:
        constructs_dict: Dictionary mapping construct positions (order, plate, well) to DataFrames containing CLIP reactions
        clips_df: DataFrame containing all unique CLIP reactions with their locations
        
    Returns:
        Dictionary where outer keys are plate numbers and inner dictionaries map destination wells to 
        [clip_wells_list, clip_plates_list] for the constructs in that plate
        
    Raises:
        ValueError: If the number of constructs exceeds protocol limits
    """
    # Count total constructs
    total_constructs = len(constructs_dict)
    
    # Generate the complete assembly plan for all constructs
    assembly_keys, assembly_values = generate_final_assembly_dict(constructs_dict, clips_df)
    
    # Calculate the number of unique construct lengths (affects master mix tip usage)
    unique_construct_lengths = {len(assembly_value[0]) for assembly_value in assembly_values}
    master_mix_tip_count = len(unique_construct_lengths)
    
    # Group constructs by their original plate numbers
    assembly_by_plate = {}
    
    for construct_index, (original_key, construct_value) in enumerate(zip(assembly_keys, assembly_values)):
        # Unpack the original key (plate, well)
        plate, well = original_key
        
        # Initialize plate dictionary if it doesn't exist
        if plate not in assembly_by_plate:
            assembly_by_plate[plate] = {}
        
        # Add construct to its original plate using well as key
        assembly_by_plate[plate][well] = construct_value
    
    return assembly_by_plate


def calculate_final_assembly_tipracks(final_assembly_dict: Dict[str, List]) -> int:
    """
    Calculate the number of tipracks required for final assembly operations.
    
    This function determines how many tipracks are needed based on:
    1. Master mix distribution tips (one per unique construct length)
    2. Individual CLIP transfer tips (one per CLIP reaction per construct)
    
    Args:
        final_assembly_dict: Dictionary mapping destination wells to 
            [clip_wells_list, clip_plates_list] for constructs
            
    Returns:
        Number of tipracks required for the assembly
        
    Raises:
        ValueError: If the calculated tipracks exceed the maximum allowed
    """
    # Count CLIP reactions per construct
    clips_per_construct = []
    for construct_data in final_assembly_dict.values():
        clips_per_construct.append(len(construct_data[0]))
    
    # Calculate unique construct lengths (affects master mix tip usage)
    unique_construct_lengths = set(clips_per_construct)
    master_mix_tip_count = len(unique_construct_lengths)
    
    # Calculate total tips needed
    total_clip_tips = sum(clips_per_construct)
    total_tips = master_mix_tip_count + total_clip_tips
    
    # Calculate number of tipracks needed
    tips_per_box = PROTOCOL_CONFIG.ASSEMBLY_TIPS_PER_BOX
    tipracks_needed = (total_tips - 1) // tips_per_box + 1
    
    max_allowed_tipracks = PROTOCOL_CONFIG.ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS
    
    print(f"Tipracks calculated: {tipracks_needed}, Maximum allowed: {max_allowed_tipracks}")
    
    if tipracks_needed > max_allowed_tipracks:
        raise ValueError(
            f'Final assembly tiprack number ({tipracks_needed}) exceeds maximum allowed ({max_allowed_tipracks}). '
            'Reduce number of constructs in constructs.csv.'
        )
    
    return tipracks_needed


def generate_ot2_script(ot2_script_path, template_path, **kwargs):
    """Generates an ot2 script named 'ot2_script_path', where kwargs are
    written as global variables at the top of the script. For each kwarg, the
    keyword defines the variable name while the value defines the name of the
    variable. The remainder of template file is subsequently written below.

    Args:
        ot2_script_path (str): Path where the OT-2 script will be written
        template_path (str): Path to the template file
        **kwargs: Variables to be written at the top of the script

    Raises:
        FileNotFoundError: If template file is not found
        IOError: If there are issues reading/writing files
    """
    def convert_numpy_types(obj):
        """Convert NumPy types to native Python types for JSON serialisation."""
        import numpy as np
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_types(item) for item in obj]
        else:
            return obj
    
    try:
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file not found: {template_path}")
            
        with open(ot2_script_path, 'w') as wf:
            with open(template_path, 'r') as rf:
                # Find the start of the function definition
                function_start = None
                for index, line in enumerate(rf):
                    if line[:3] == 'def':
                        function_start = index
                        break
                    else:
                        wf.write(line)
                
                if function_start is None:
                    raise ValueError(f"No function definition found in template: {template_path}")
                
                # Write the variables
                for key, value in kwargs.items():
                    wf.write('{}='.format(key))
                    if type(value) == dict:
                        # Convert NumPy types before JSON serialisation
                        converted_value = convert_numpy_types(value)
                        wf.write(json.dumps(converted_value))
                    elif type(value) == str:
                        wf.write("'{}'".format(value))
                    else:
                        # Convert NumPy types for other types too
                        converted_value = convert_numpy_types(value)
                        wf.write(str(converted_value))
                    wf.write('\n')
                wf.write('\n')
                
            # Write the rest of the template, but skip any existing variable definitions
            with open(template_path, 'r') as rf:
                lines = rf.readlines()
                skip_until_function = False
                
                for index, line in enumerate(lines):
                    # Skip lines until we find the function definition
                    if index < function_start - 1:
                        continue
                    
                    # Check if this line defines a variable that we're overriding
                    should_skip_line = False
                    for key in kwargs.keys():
                        # Check for both 'key=' and 'key =' patterns
                        stripped_line = line.strip()
                        if (stripped_line.startswith(f'{key}=') or 
                            stripped_line.startswith(f'{key} =')):
                            should_skip_line = True
                            break
                    
                    if should_skip_line:
                        continue
                    
                    wf.write(line)
                        
        print(f"    ✓ {os.path.basename(ot2_script_path)}")
        
    except Exception as e:
        print(f"\nError generating OT-2 script {os.path.basename(ot2_script_path)}:")
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        raise


def _generate_clip_script_media_bot_style(ot2_script_path: str, template_path: str, clips_dict: Dict[str, List]) -> None:
    """Generate CLIP script using Media Bot-style parameterisation.
    
    This function replaces the JSON file loading code in the template with embedded JSON data,
    similar to how Media Bot parameterises its templates.
    
    Args:
        ot2_script_path (str): Path where the OT-2 script will be written
        template_path (str): Path to the template file
        clips_dict (Dict[str, List]): CLIP reaction data dictionary
        
    Raises:
        FileNotFoundError: If template file is not found
        IOError: If there are issues reading/writing files
    """
    def convert_numpy_types(obj):
        """Convert NumPy types to native Python types for JSON serialisation."""
        import numpy as np
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_types(item) for item in obj]
        else:
            return obj
    
    try:
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file not found: {template_path}")
            
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # Convert clips_dict to JSON string
        converted_clips_dict = convert_numpy_types(clips_dict)
        clips_json = json.dumps(converted_clips_dict, indent=4)
        
        # Replace the JSON file loading code with embedded JSON data
        modified_protocol = template_content.replace(
            "with open('clips_data.json') as f:\n    clips_dict = json.load(f)",
            f"clips_dict = {clips_json}"
        )
        
        # Write the modified protocol
        with open(ot2_script_path, 'w') as f:
            f.write(modified_protocol)
            
        print(f"    ✓ {os.path.basename(ot2_script_path)}")
        
    except Exception as e:
        print(f"\nError generating CLIP script {os.path.basename(ot2_script_path)}:")
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        raise


def _generate_assembly_script_media_bot_style(ot2_script_path: str, template_path: str, 
                                             final_assembly_dict: Dict[str, List], tiprack_num: int) -> None:
    """Generate assembly script using Media Bot-style parameterisation.
    
    This function replaces the JSON file loading code in the template with embedded JSON data,
    similar to how Media Bot parameterises its templates.
    
    Args:
        ot2_script_path (str): Path where the OT-2 script will be written
        template_path (str): Path to the template file
        final_assembly_dict (Dict[str, List]): Final assembly data dictionary
        tiprack_num (int): Number of tipracks required
        
    Raises:
        FileNotFoundError: If template file is not found
        IOError: If there are issues reading/writing files
    """
    def convert_numpy_types(obj):
        """Convert NumPy types to native Python types for JSON serialisation."""
        import numpy as np
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy_types(item) for item in obj]
        else:
            return obj
    
    try:
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file not found: {template_path}")
            
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # Create assembly data dictionary
        assembly_data = {
            'final_assembly_dict': final_assembly_dict,
            'tiprack_num': tiprack_num
        }
        
        # Convert to JSON string with compact formatting
        converted_assembly_data = convert_numpy_types(assembly_data)
        
        # Create a more compact and readable format for the assembly data
        final_assembly_dict = converted_assembly_data['final_assembly_dict']
        tiprack_num = converted_assembly_data['tiprack_num']
        
        # Format the final_assembly_dict in a more compact way
        assembly_dict_lines = []
        assembly_dict_lines.append("final_assembly_dict = {")
        
        for i, (well, (clips, plates)) in enumerate(final_assembly_dict.items()):
            # Format clips as individual list items
            clips_list = [f'"{clip}"' for clip in clips]
            clips_str = f'[{", ".join(clips_list)}]'
            plates_str = f'[{", ".join(map(str, plates))}]'
            
            if i == len(final_assembly_dict) - 1:
                # Last item - no comma
                assembly_dict_lines.append(f'    "{well}": [{clips_str}, {plates_str}]')
            else:
                assembly_dict_lines.append(f'    "{well}": [{clips_str}, {plates_str}],')
        
        assembly_dict_lines.append("}")
        
        # Create the compact assembly data string
        compact_assembly_data = "\n".join(assembly_dict_lines)
        
        # Replace the JSON file loading code with embedded compact data
        modified_protocol = template_content.replace(
            "with open('assembly_data.json') as f:\n    assembly_data = json.load(f)\n    final_assembly_dict = assembly_data['final_assembly_dict']\n    tiprack_num = assembly_data['tiprack_num']",
            f"{compact_assembly_data}\ntiprack_num = {tiprack_num}"
        )
        
        # Write the modified protocol
        with open(ot2_script_path, 'w') as f:
            f.write(modified_protocol)
            
        print(f"    ✓ {os.path.basename(ot2_script_path)}")
        
    except Exception as e:
        print(f"\nError generating assembly script {os.path.basename(ot2_script_path)}:")
        print(f"Error: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        raise


def generate_new_constructs_df(construct_path, final_assembly_dict_list, constructs_dict):
    """Takes the user supplied constructs csv and reallocates the constructs 
    according to the new order applied in the 
    'generate_final_assembly_dict_list' function, returning a new constructs
    csv with the new locations of each construct given.
    
    Args:
        construct_path (str): Path to the original constructs CSV file
        final_assembly_dict_list (Dict[int, Dict[str, List]]): Dictionary of assembly dictionaries by plate
        constructs_dict (Dict[Tuple[int, int, str], pd.DataFrame]): Processed constructs dictionary
        
    Returns:
        pd.DataFrame: DataFrame with the original construct data plus new
            Assembly and Well columns showing the final locations
            
    Note:
        The input CSV should have headers, with the first column being the
        construct ID/name and subsequent columns containing the construct
        components (linkers, parts, etc.).
    """
    # Create a mapping from original positions to new assembly positions
    # Sort constructs by their original order to maintain consistency
    sorted_constructs = sorted(constructs_dict.items(), key=lambda x: x[0][0])  # Sort by order (first element of tuple)
    
    # Create lists for assembly and well assignments
    assembly_list = []
    well_list = []
    
    # Map each construct to its new assembly position
    for (order, original_plate, original_well), construct_df in sorted_constructs:
        # Find which assembly plate this construct ended up in
        for assembly_plate, assembly_dict in final_assembly_dict_list.items():
            if original_well in assembly_dict:
                assembly_list.append(assembly_plate)
                well_list.append(original_well)
                break
    
    # Read the original CSV to preserve the layout
    constructs_df = pd.read_csv(construct_path)
    
    # Determine format from the header
    header_row = constructs_df.columns.tolist()
    has_plate_column = len(header_row) >= 2 and header_row[0].strip().lower() == 'plate'
    
    # Create assembly and well lists that match the original layout
    layout_assembly_list = []
    layout_well_list = []
    construct_index = 0
    
    for index, row in constructs_df.iterrows():
        # Use the format determined from header
        if has_plate_column:
            # New format: Plate, Well, Linker1, Part1, ...
            construct_components = row.iloc[2:].tolist()  # Start from column 3 (after Plate, Well)
        else:
            # Old format: Well, Linker1, Part1, ...
            construct_components = row.iloc[1:].tolist()  # Start from column 2
        
        # Filter out empty strings and NaN values
        construct_components = [comp for comp in construct_components if comp and str(comp).strip() and str(comp).strip() != 'nan']
        
        # Check if this row has valid construct data
        if not construct_components:
            # Empty row - use empty values
            layout_assembly_list.append("")
            layout_well_list.append("")
            continue
        
        # Valid construct - use assembly data
        if construct_index < len(assembly_list):
            layout_assembly_list.append(assembly_list[construct_index])
            layout_well_list.append(well_list[construct_index])
            construct_index += 1
        else:
            # This shouldn't happen, but just in case
            layout_assembly_list.append("")
            layout_well_list.append("")
    
    # Handle column insertion based on format
    if has_plate_column:
        # New format: Plate, Well, Linker1, Part1, ...
        # Remove the first two columns (Plate, Well) and add Assembly and Well columns
        constructs_df.drop(constructs_df.columns[0:2], axis=1, inplace=True)
        constructs_df.insert(0, "Assembly", layout_assembly_list)
        constructs_df.insert(1, "Well", layout_well_list)
    else:
        # Old format: Well, Linker1, Part1, ...
        # Remove the first column (Well) and add Assembly and Well columns
        constructs_df.drop(constructs_df.columns[0], axis=1, inplace=True)
        constructs_df.insert(0, "Assembly", layout_assembly_list)
        constructs_df.insert(1, "Well", layout_well_list)

    return constructs_df


def generate_master_mix_df(clip_number):
    """Generates a dataframe detailing the components required in the clip
    reaction master mix.

    """
    COMPONENTS = {'Component': ['Promega T4 DNA Ligase buffer, 10X',
                                'Water', 'NEB BsaI-HFv2',
                                'Promega T4 DNA Ligase']}
    VOL_COLUMN = 'Volume (uL)'
    master_mix_df = pd.DataFrame.from_dict(COMPONENTS)
    
    # Calculate volumes for each component
    clip_vol = PROTOCOL_CONFIG.CLIP_VOL
    dead_vol = PROTOCOL_CONFIG.CLIP_DEAD_VOL
    t4_buff_vol = PROTOCOL_CONFIG.CLIP_T4_BUFF_VOL
    mast_water = PROTOCOL_CONFIG.CLIP_MAST_WATER
    bsai_vol = PROTOCOL_CONFIG.CLIP_BSAI_VOL
    t4_lig_vol = PROTOCOL_CONFIG.CLIP_T4_LIG_VOL
    
    # Ensure float calculation to avoid integer results
    multiplier = float(clip_number + dead_vol/clip_vol)
    master_mix_df[VOL_COLUMN] = multiplier * \
        np.array([t4_buff_vol, mast_water, bsai_vol, t4_lig_vol])
    return master_mix_df


def generate_sources_paths_df(paths, deck_positions):
    """Generates a dataframe detailing source plate information.

    Args:
        paths (list): list of strings specifying paths to source plates.
        deck_positions (list): list of strings specifying candidate deck positions.

    """
    source_plates_dict = {'Deck position': [], 'Source plate': [], 'Path': []}
    for index, path in enumerate(paths):
        source_plates_dict['Deck position'].append(DECK_CONFIG.SOURCE_POSITIONS[index])
        source_plates_dict['Source plate'].append(os.path.basename(path))
        source_plates_dict['Path'].append(path)
    return pd.DataFrame(source_plates_dict)


def dfs_to_csv(path, index=True, **kw_dfs):
    """Generates a csv file defined by path, where kw_dfs are
    written one after another with each key acting as a title. If index=True,
    df indexes are written to the csv file.

    """
    with open(path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for key, value in kw_dfs.items():
            csvwriter.writerow([str(key)])
            value.to_csv(csvfile, index=index)
            csvwriter.writerow('')


def counter(rows: int):
    """
    Returns a function that converts between well index and well coordinate for a plate with the given number of rows.
    - Integer index (0-based) -> Well coordinate (e.g., 0 -> 'A1', 8 -> 'B1')
    - Well coordinate -> Integer index (e.g., 'A1' -> 0, 'B1' -> 8)
    """
    import re
    def convert(well_position):
        if isinstance(well_position, int):
            row = chr(ord('A') + (well_position % rows))
            col = 1 + well_position // rows
            return f"{row}{col}"
        elif isinstance(well_position, str):
            row, col = re.findall(r'\d+|\D+', well_position)
            col = int(col) - 1
            return col * rows + (ord(row) - ord('A'))
        else:
            raise TypeError(f"Expected int or str, got {type(well_position)}")
    return convert


# Default for 96-well plates (8 rows)
tip_counter = counter(8)


def calculate_construct_position(index: int) -> Tuple[int, int, str]:
    """
    Calculate order, plate number and well coordinate for a construct based on its index.
    
    Args:
        index (int): 0-based index of the construct
        
    Returns:
        Tuple[int, int, str]: (order, plate_number, well_coordinate)
        Example: (0, 1, 'A1') for index 0, (95, 1, 'H12') for index 95, (96, 2, 'A1') for index 96
    """
    plate = index // 96 + 1
    well = tip_counter(index % 96)
    return (index, plate, well)


def normalize_source_data(data_tuple: Union[Tuple, List]) -> Tuple[str, str, str]:
    """
    Normalize source data to ensure consistent 3-column format.
    
    This function handles both 2-column and 3-column CSV formats by ensuring
    the concentration field is always present (empty string if not provided).
    
    Args:
        data_tuple: Tuple or list containing source data
        
    Returns:
        Normalized 3-element tuple: (well, concentration, deck_position)
        
    Examples:
        >>> normalize_source_data(('A8', '2'))
        ('A8', '', '2')
        >>> normalize_source_data(('A8', '', '2'))
        ('A8', '', '2')
    """
    if isinstance(data_tuple, list):
        data_tuple = data_tuple[0]
    
    if len(data_tuple) == 2:
        # Insert empty concentration field
        return (data_tuple[0], "", data_tuple[1])
    elif len(data_tuple) >= 3:
        # Return first three elements
        return (data_tuple[0], data_tuple[1], data_tuple[2])
    else:
        raise ValueError(f"Expected 2 or more elements, got {len(data_tuple)}")


def validate_construct_data(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], clips_df: pd.DataFrame) -> None:
    """
    Validate that construct data is properly formatted and within protocol limits.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing construct information
        clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If constructs are malformed or exceed limits
    """
    if not constructs_dict:
        raise ValueError("No constructs found in input file")
    
    # Use clips_df for accurate clip counting
    total_clips = clips_df['number'].sum()
    
    max_clips = PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL
    
    if total_clips > max_clips:
        raise ValueError(
            f"Total number of CLIP reactions ({total_clips}) exceeds maximum allowed ({max_clips}). "
            "Reduce the number of constructs or simplify construct designs."
        )
    
    # Validate individual constructs
    for i, (position, construct) in enumerate(constructs_dict.items()):
        if construct.empty:
            raise ValueError(f"Construct at {position} is empty")
        
        required_columns = ['prefixes', 'parts', 'suffixes']
        missing_columns = [col for col in required_columns if col not in construct.columns]
        if missing_columns:
            raise ValueError(f"Construct at {position} missing required columns: {missing_columns}")
        
        # Check for empty or invalid values
        for col in required_columns:
            empty_values = construct[construct[col].isna() | (construct[col] == '')]
            if not empty_values.empty:
                raise ValueError(f"Construct at {position} has empty values in column '{col}' at rows: {empty_values.index.tolist()}")


def validate_clips_data(clips_df: pd.DataFrame) -> None:
    """
    Validate that CLIP reaction data is properly formatted.
    
    Args:
        clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If CLIP data is malformed
    """
    if clips_df.empty:
        raise ValueError("No CLIP reactions found")
    
    required_columns = ['prefixes', 'parts', 'suffixes', 'number', 'mag_well', 'plate']
    missing_columns = [col for col in required_columns if col not in clips_df.columns]
    if missing_columns:
        raise ValueError(f"CLIP DataFrame missing required columns: {missing_columns}")
    
    # Check for negative or zero reaction numbers
    invalid_numbers = clips_df[clips_df['number'] <= 0]
    if not invalid_numbers.empty:
        raise ValueError(f"Found {len(invalid_numbers)} CLIP reactions with invalid numbers (≤0): {invalid_numbers.index.tolist()}")
    
    # Check for empty values in required columns
    for col in ['prefixes', 'parts', 'suffixes']:
        empty_values = clips_df[clips_df[col].isna() | (clips_df[col] == '')]
        if not empty_values.empty:
            raise ValueError(f"CLIP DataFrame has empty values in column '{col}' at rows: {empty_values.index.tolist()}")


def validate_sources_data(sources_dict: Dict[str, Tuple[str, ...]]) -> None:
    """
    Validate that source data is properly formatted.
    
    Args:
        sources_dict: Dictionary mapping parts/linkers to source locations
        
    Raises:
        ValueError: If source data is malformed
    """
    if not sources_dict:
        raise ValueError("No source data found")
    
    invalid_entries = []
    wells = []
    
    for part_name, source_data in sources_dict.items():
        if not isinstance(source_data, (tuple, list)):
            invalid_entries.append(f"{part_name}: expected tuple/list, got {type(source_data)}")
            continue
        
        if len(source_data) < 2:
            invalid_entries.append(f"{part_name}: insufficient data (need at least well and deck position)")
            continue
        
        # Validate well format
        well = source_data[0]
        try:
            validate_well_format(well)
            wells.append(well)
        except ValueError as e:
            invalid_entries.append(f"{part_name}: {str(e)}")
    
    if invalid_entries:
        raise ValueError("Invalid source data entries:\n" + "\n".join(invalid_entries))
    
    # Check for duplicate wells
    try:
        validate_unique_wells(wells, "source plate wells")
    except ValueError as e:
        raise ValueError(f"Source data validation failed: {str(e)}")


def validate_components_availability(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                   sources_dict: Dict[str, Tuple[str, ...]]) -> None:
    """
    Validate that all required parts and linkers are available in the source data.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing construct information
        sources_dict: Dictionary mapping parts/linkers to source locations
        
    Raises:
        ValueError: If any required components are missing
    """
    # Collect all required components
    required_components = set()
    for construct in constructs_dict.values():
        for _, row in construct.iterrows():
            required_components.add(row['prefixes'].strip())
            required_components.add(row['parts'].strip())
            required_components.add(row['suffixes'].strip())
    
    # Check which components are missing
    missing_components = required_components - set(sources_dict.keys())
    
    if missing_components:
        raise ValueError(
            f"The following parts/linkers are required but not found in the source data:\n"
            f"{', '.join(sorted(missing_components))}\n\n"
            f"Please add these components to your source CSV files."
        )


def validate_tip_usage(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                      clips_df: pd.DataFrame) -> None:
    """
    Validate that the experiment doesn't exceed available tip capacity.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing construct information
        clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If tip usage exceeds limits
    """
    # Calculate total tips needed
    total_clip_tips = clips_df['number'].sum()
    
    # Calculate unique construct lengths for master mix tips
    unique_construct_lengths = set()
    for construct in constructs_dict.values():
        unique_construct_lengths.add(len(construct))
    
    master_mix_tips = len(unique_construct_lengths)
    total_tips = total_clip_tips + master_mix_tips
    
    max_tips = PROTOCOL_CONFIG.ASSEMBLY_TIPS_PER_BOX * PROTOCOL_CONFIG.ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS
    
    if total_tips > max_tips:
        raise ValueError(
            f"Total tip usage ({total_tips}) exceeds maximum available ({max_tips}).\n"
            f"Breakdown:\n"
            f"  - CLIP transfer tips: {total_clip_tips}\n"
            f"  - Master mix tips: {master_mix_tips}\n"
            f"  - Total: {total_tips}\n"
            f"  - Maximum: {max_tips}\n\n"
            f"Consider reducing the number of constructs or simplifying construct designs."
        )


def log_processing_summary(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                          clips_df: pd.DataFrame,
                          sources_dict: Dict[str, Tuple[str, ...]]) -> None:
    """
    Log a summary of the processing results for user verification.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames
        clips_df: CLIP reactions DataFrame
        sources_dict: Source locations dictionary
    """
    print("\n" + "="*60)
    print("PROCESSING SUMMARY")
    print("="*60)
    print(f"Total constructs: {len(constructs_dict)}")
    print(f"Total CLIP reactions: {clips_df['number'].sum()}")
    print(f"Unique CLIP reactions: {len(clips_df)}")
    print(f"Source parts/linkers: {len(sources_dict)}")
    
    # Calculate some useful statistics
    valid_constructs = list(constructs_dict.values())
    avg_clips_per_construct = sum(len(construct) for construct in valid_constructs) / len(valid_constructs)
    print(f"Average CLIP reactions per construct: {avg_clips_per_construct:.1f}")
    
    max_clips_per_construct = max(len(construct) for construct in valid_constructs)
    min_clips_per_construct = min(len(construct) for construct in valid_constructs)
    print(f"CLIP reactions per construct range: {min_clips_per_construct} - {max_clips_per_construct}")
    
    print("="*60)


def validate_well_format(well: str) -> bool:
    """
    Validate that a well identifier is in the correct format (e.g., 'A1' or 'H12').
    
    Args:
        well (str): Well identifier to validate
        
    Returns:
        bool: True if well is valid
        
    Raises:
        ValueError: If well format is invalid
    """
    if not isinstance(well, str):
        raise ValueError(f"Well identifier must be a string, got {type(well)}")
    
    # Handle empty wells (blank wells)
    if not well.strip():
        return True
    
    if len(well) < 2 or len(well) > 3:
        raise ValueError(f"Well identifier must be 2-3 characters long, got '{well}'")
    
    if not well[0].isalpha() or not well[0].isupper():
        raise ValueError(f"Well row must be an uppercase letter, got '{well[0]}' in '{well}'")
    
    if not well[1:].isdigit():
        raise ValueError(f"Well column must be a number, got '{well[1:]}' in '{well}'")
    
    row = ord(well[0]) - ord('A')
    col = int(well[1:])
    
    if row < 0 or row > 7:
        raise ValueError(f"Well row must be A-H, got '{well[0]}' in '{well}'")
    
    if col < 1 or col > 12:
        raise ValueError(f"Well column must be 1-12, got '{col}' in '{well}'")
    
    return True


def validate_unique_wells(wells: List[str], context: str = "wells") -> None:
    """
    Validate that all wells are unique.
    
    Args:
        wells (List[str]): List of well identifiers
        context (str): Context for error message (e.g., "source plate wells")
        
    Raises:
        ValueError: If duplicate wells are found
    """
    seen = set()
    duplicates = set()
    
    for well in wells:
        if well in seen:
            duplicates.add(well)
        seen.add(well)
    
    if duplicates:
        raise ValueError(f"Duplicate {context} found: {', '.join(sorted(duplicates))}")


def validate_csv_columns(reader: csv.DictReader, required_columns: List[str], file_name: str) -> None:
    """
    Validate that a CSV file contains all required columns.
    
    Args:
        reader: CSV DictReader object
        required_columns: List of required column names
        file_name: Name of the file being validated
        
    Raises:
        ValueError: If required columns are missing
    """
    missing_columns = [col for col in required_columns if col not in reader.fieldnames]
    if missing_columns:
        raise ValueError(f"CSV file '{file_name}' is missing required columns: {', '.join(missing_columns)}")


def count_clips_for_constructs(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame]) -> int:
    """Count total number of clips needed for a set of constructs."""
    if not constructs_dict:
        return 0
    
    # Merge all constructs and find unique CLIP reactions
    valid_constructs = list(constructs_dict.values())
    merged_construct_dfs = pd.concat(valid_constructs, ignore_index=True)
    unique_clips_df = merged_construct_dfs.drop_duplicates().reset_index(drop=True)
    
    # Count occurrences of each unique CLIP reaction
    clip_count = np.zeros(len(unique_clips_df.index))
    for i, unique_clip in unique_clips_df.iterrows():
        for _, clip in merged_construct_dfs.iterrows():
            if unique_clip.equals(clip):
                clip_count[i] += 1
    
    # Calculate number of reactions needed based on final assemblies per CLIP
    clip_count = clip_count // PROTOCOL_CONFIG.ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP + 1
    
    return int(clip_count.sum())


def generate_optimised_clips_dict_list(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                      sources_dict: Dict[str, Tuple[str, ...]]) -> Tuple[List[Dict[str, List]], Dict[int, List[str]], pd.DataFrame]:
    """
    Generate optimised clip assignments by splitting constructs by assembly plates.
    
    This function implements an intelligent clip assignment strategy that:
    1. Validates total clip count doesn't exceed maximum (192 clips = 2 full plates)
    2. Splits constructs into two halves based on assembly plates
    3. Validates each half requires ≤96 clips
    4. Generates separate clip lists for each half
    5. Creates a mapping showing which clip scripts are needed for each assembly plate
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing CLIP reactions
        sources_dict: Dictionary mapping parts/linkers to their source locations
        
    Returns:
        Tuple containing:
        - List of clip dictionaries for OT-2 scripts
        - Dictionary mapping assembly plate numbers to lists of required clip script names (e.g., ['1a', '1b'])
    """
    
    def split_constructs_by_assembly_plates(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame]) -> Tuple[Dict, Dict]:
        """
        Split constructs into two halves based on assembly plates.
        
        Returns:
            Tuple of (first_half, second_half) construct dictionaries
        """
        # Get all assembly plates
        assembly_plates = set()
        for (order, plate, well) in constructs_dict.keys():
            assembly_plates.add(plate)
        
        assembly_plates = sorted(assembly_plates)
        
        if len(assembly_plates) <= 1:
            # Only one assembly plate, return original dict and empty dict
            return constructs_dict, {}
        
        # Split assembly plates into two halves
        mid_point = (len(assembly_plates) + 1) // 2  # Ceiling division
        first_half_plates = assembly_plates[:mid_point]
        second_half_plates = assembly_plates[mid_point:]
        
        # Split constructs based on assembly plates
        first_half = {}
        second_half = {}
        
        for (order, plate, well), construct_df in constructs_dict.items():
            if plate in first_half_plates:
                first_half[(order, plate, well)] = construct_df
            else:
                second_half[(order, plate, well)] = construct_df
        
        return first_half, second_half
    
    def validate_construct_from_clips(construct_df: pd.DataFrame, clips_df: pd.DataFrame, 
                                    sources_dict: Dict[str, Tuple[str, ...]], construct_name: str = "test") -> bool:
        """
        Validate that a construct can be reconstructed from its clips.
        
        Args:
            construct_df: DataFrame containing the original construct
            clips_df: DataFrame containing all clips with their wells
            sources_dict: Dictionary mapping parts to source locations
            construct_name: Name of the construct for error reporting
            
        Returns:
            bool: True if construct can be reconstructed correctly
        """
        try:
            reconstructed_clips = []
            
            for _, clip_row in construct_df.iterrows():
                prefix = clip_row['prefixes'].strip()
                part = clip_row['parts'].strip()
                suffix = clip_row['suffixes'].strip()
                
                # Find this clip in the clips_df
                matching_clips = clips_df[
                    (clips_df['prefixes'] == prefix) &
                    (clips_df['parts'] == part) &
                    (clips_df['suffixes'] == suffix)
                ]
                
                if matching_clips.empty:
                    print(f"  ❌ {construct_name}: Clip {prefix}-{part}-{suffix} not found in clips_df")
                    return False
                
                # Check that source locations match
                if prefix not in sources_dict:
                    print(f"  ❌ {construct_name}: Prefix {prefix} not found in sources")
                    return False
                if part not in sources_dict:
                    print(f"  ❌ {construct_name}: Part {part} not found in sources")
                    return False
                if suffix not in sources_dict:
                    print(f"  ❌ {construct_name}: Suffix {suffix} not found in sources")
                    return False
                
                reconstructed_clips.append({
                    'prefix': prefix,
                    'part': part,
                    'suffix': suffix,
                    'prefix_well': sources_dict[prefix][0],
                    'part_well': sources_dict[part][0],
                    'suffix_well': sources_dict[suffix][0]
                })
            
            print(f"  ✓ {construct_name}: Successfully reconstructed {len(reconstructed_clips)} clips")
            return True
            
        except Exception as e:
            print(f"  ❌ {construct_name}: Error during reconstruction: {str(e)}")
            return False
    
    def validate_clip_assignments(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                 clips_df: pd.DataFrame, 
                                 sources_dict: Dict[str, Tuple[str, ...]]) -> None:
        """
        Validate that clip assignments are correct by reconstructing a sample of constructs.
        
        Args:
            constructs_dict: Dictionary of all constructs
            clips_df: DataFrame containing all clips
            sources_dict: Dictionary mapping parts to source locations
        """
        print("\nValidating clip assignments...")
        
        # Get a sample of constructs to validate (first, middle, last)
        construct_positions = list(constructs_dict.keys())
        sample_positions = []
        
        if len(construct_positions) >= 3:
            sample_positions = [
                construct_positions[0],  # First
                construct_positions[len(construct_positions)//2],  # Middle
                construct_positions[-1]  # Last
            ]
        elif len(construct_positions) >= 1:
            sample_positions = construct_positions[:min(3, len(construct_positions))]
        
        validation_passed = True
        for i, position in enumerate(sample_positions):
            order, plate, well = position
            construct_df = constructs_dict[position]
            construct_name = f"Construct_{order}_{plate}_{well}"
            
            if not validate_construct_from_clips(construct_df, clips_df, sources_dict, construct_name):
                validation_passed = False
        
        if validation_passed:
            print("✓ All sampled constructs validated successfully")
        else:
            print("❌ Some constructs failed validation")
            raise ValueError("Clip assignment validation failed")
    
    # Get all assembly plates
    assembly_plates = set()
    for (order, plate, well) in constructs_dict.keys():
        assembly_plates.add(plate)
    assembly_plates = sorted(assembly_plates)
    
    print(f"\nOptimising clip assignment for {len(assembly_plates)} assembly plate(s)")
    print(f"Assembly plates: {assembly_plates}")
    
    # Initial validation: Check total clips don't exceed maximum
    total_clips = count_clips_for_constructs(constructs_dict)
    max_clips = PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL  # 192 clips = 2 full plates
    
    print(f"Total clips required: {total_clips}")
    print(f"Maximum clips allowed: {max_clips}")
    
    if total_clips > max_clips:
        raise ValueError(
            f"Total clips ({total_clips}) exceed maximum allowed ({max_clips}). "
            "Reduce number of constructs or simplify construct designs."
        )
    
    # If only one assembly plate, use original approach
    if len(assembly_plates) <= 1:
        print("Single assembly plate - using standard clip assignment")
        clips_df = generate_clips_df(constructs_dict)
        clips_dict_list = generate_clips_dict_list(clips_df, sources_dict)
        
        # Validate assignments
        validate_clip_assignments(constructs_dict, clips_df, sources_dict)
        
        # Create mapping with new naming scheme
        script_names = []
        for i in range(len(clips_dict_list)):
            plate_num = (i // 2) + 1
            half_letter = 'a' if i % 2 == 0 else 'b'
            script_names.append(f"{plate_num}{half_letter}")
        
        assembly_to_clip_mapping = {1: script_names}
        
        return clips_dict_list, assembly_to_clip_mapping, clips_df
    
    # Split constructs by assembly plates
    first_half, second_half = split_constructs_by_assembly_plates(constructs_dict)
    
    # Count clips for each half
    first_clips = count_clips_for_constructs(first_half)
    second_clips = count_clips_for_constructs(second_half)
    
    print(f"\nInitial split:")
    print(f"  First half (plates {[p for (o,p,w) in first_half.keys()]}): {first_clips} clips")
    print(f"  Second half (plates {[p for (o,p,w) in second_half.keys()]}): {second_clips} clips")
    
    # Validate each half is within limits
    if first_clips > 96 or second_clips > 96:
        print(f"One or both halves exceed 96 clips - reverting to standard approach")
        print(f"  First half: {first_clips} clips")
        print(f"  Second half: {second_clips} clips")
        
        # Revert to original approach - generate clips for all constructs together
        clips_df = generate_clips_df(constructs_dict)
        clips_dict_list = generate_clips_dict_list(clips_df, sources_dict)
        
        # Validate assignments
        validate_clip_assignments(constructs_dict, clips_df, sources_dict)
        
        # Create mapping with new naming scheme for all scripts
        script_names = []
        for i in range(len(clips_dict_list)):
            plate_num = (i // 2) + 1
            half_letter = 'a' if i % 2 == 0 else 'b'
            script_names.append(f"{plate_num}{half_letter}")
        
        # All assembly plates use all clip scripts
        assembly_to_clip_mapping = {}
        for plate in assembly_plates:
            assembly_to_clip_mapping[plate] = script_names
        
        print(f"Standard approach: {len(clips_dict_list)} clip script(s) for all assembly plates")
        
        return clips_dict_list, assembly_to_clip_mapping, clips_df
    
    # Generate clips for each half
    first_clips_df = generate_clips_df(first_half) if first_half else pd.DataFrame()
    second_clips_df = generate_clips_df(second_half) if second_half else pd.DataFrame()
    
    first_clips_dict_list = generate_clips_dict_list(first_clips_df, sources_dict) if not first_clips_df.empty else []
    second_clips_dict_list = generate_clips_dict_list(second_clips_df, sources_dict) if not second_clips_df.empty else []
    
    # Combine clips dict lists
    clips_dict_list = first_clips_dict_list + second_clips_dict_list
    
    # Create mapping with new naming scheme
    assembly_to_clip_mapping = {}
    
    # Map first half assembly plates to first half clip scripts
    first_half_assembly_plates = set()
    for (order, plate, well) in first_half.keys():
        first_half_assembly_plates.add(plate)
    
    first_script_names = []
    for i in range(len(first_clips_dict_list)):
        plate_num = (i // 2) + 1
        half_letter = 'a' if i % 2 == 0 else 'b'
        first_script_names.append(f"{plate_num}{half_letter}")
    
    for plate in first_half_assembly_plates:
        assembly_to_clip_mapping[plate] = first_script_names
    
    # Map second half assembly plates to second half clip scripts
    second_half_assembly_plates = set()
    for (order, plate, well) in second_half.keys():
        second_half_assembly_plates.add(plate)
    
    second_script_names = []
    for i in range(len(second_clips_dict_list)):
        plate_num = ((len(first_clips_dict_list) + i) // 2) + 1
        half_letter = 'a' if (len(first_clips_dict_list) + i) % 2 == 0 else 'b'
        second_script_names.append(f"{plate_num}{half_letter}")
    
    for plate in second_half_assembly_plates:
        assembly_to_clip_mapping[plate] = second_script_names
    
    # Validate assignments
    combined_clips_df = pd.concat([first_clips_df, second_clips_df], ignore_index=True) if not first_clips_df.empty and not second_clips_df.empty else (first_clips_df if not first_clips_df.empty else second_clips_df)
    validate_clip_assignments(constructs_dict, combined_clips_df, sources_dict)
    
    # Create optimised clips DataFrame with correct plate assignments
    if not first_clips_df.empty and not second_clips_df.empty:
        # Create optimised clips DataFrame
        optimised_clips_df = pd.concat([first_clips_df, second_clips_df], ignore_index=True)
        
        # Update plate assignments for optimised structure
        for idx in range(len(optimised_clips_df)):
            if idx < len(first_clips_df):
                # First half clips - assign to plate 1
                optimised_clips_df.at[idx, 'plate'] = (1,)
            else:
                # Second half clips - assign to plate 2
                optimised_clips_df.at[idx, 'plate'] = (2,)
    else:
        # Single half or no optimisation
        optimised_clips_df = combined_clips_df
    
    print(f"\nOptimisation complete:")
    print(f"  First half: {len(first_clips_dict_list)} clip script(s) for assembly plate(s) {sorted(first_half_assembly_plates)}")
    print(f"  Second half: {len(second_clips_dict_list)} clip script(s) for assembly plate(s) {sorted(second_half_assembly_plates)}")
    print(f"  Total clip scripts: {len(clips_dict_list)}")
    
    return clips_dict_list, assembly_to_clip_mapping, optimised_clips_df


def generate_clips_dict_list(clips_df, sources_dict):
    '''Subsets the clips df into chunks of 48, runs the generate_clips_dict function 
    for each and returns a list of the resulting sub clips dicts'''

    def clip_df_long_format(clip_df):
        ''' Takes clips df and returns a long format df with one row per mag well 
        rather than one row per unique clip'''

        df = clip_df.copy()
        long_clip_df = pd.DataFrame(columns = df.columns)

        index = 0

        for row, clip in clip_df.iterrows():
            for clip_rep in range(clip.number):                         # copy row and replace mag well and plate lists with well and plate for current rep
                data = clip.copy()
                data.number = 1
                data.mag_well = clip.mag_well[clip_rep]
                data.plate = clip.plate[clip_rep]
                long_clip_df.loc[index] = data
                index += 1

        return long_clip_df
    
    long_clip_df = clip_df_long_format(clips_df)

    CLIP_COUNT = clips_df['number'].sum()
    CLIP_PLATE_COUNT = int(CLIP_COUNT // PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE + 1)       # plus one to include final partially full plate

    # Error 
    if clips_df['number'].sum() > PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL:
        raise ValueError('Number of CLIP reactions exceeds {}. Reduce number of constructs in construct.csv.'.format(PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL))
    
    if clips_df['number'].sum() < 96:                                   # if less clips required than one full plate, the second desk slot can be used for tips
        PROTOCOL_CONFIG.ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS = 5

    clips_dict_list = []

    for plate in range(CLIP_PLATE_COUNT):
        subset_lower = (plate * PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE)                    # set upper and lower bounds for subset of clips for a given plate
        subset_upper = subset_lower + PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_PER_PLATE

        if subset_upper > CLIP_COUNT:                                   # set total number number of clips as upper bound if plate incomplete
            subset_upper = CLIP_COUNT
    
        sub_clip_df = long_clip_df.iloc[subset_lower:subset_upper, :]
        sub_clip_dict = generate_clips_dict(sub_clip_df, sources_dict)
        clips_dict_list.append(sub_clip_dict)                           # generate and append sub_clip_dict to list - allows for multiple clip reactions

    return clips_dict_list


if __name__ == '__main__':
    main()
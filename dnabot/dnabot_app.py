# -*- coding: utf-8 -*-
"""
DNA-BOT: DNA assembly using BASIC on OpenTrons

@author: mh2210, ljh119

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
import shutil

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
    CLIP_DEFAULT_PART_CONC: float = 200.0  # Default concentration in ng/µl if not specified
    
    # Assembly parameters
    ASSEMBLY_MAX_CLIPS_PER_PLATE: int = 48
    ASSEMBLY_MAX_CLIPS_TOTAL: int = 96 * 2
    ASSEMBLY_MAX_ASSEMBLIES_PER_PLATE: int = 96
    ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP: int = 13
    ASSEMBLY_MAX_SOURCE_PLATES: int = 6
    ASSEMBLY_MAX_FINAL_ASSEMBLY_TIPRACKS: int = 4
    ASSEMBLY_TIPS_PER_BOX: int = 96
    
    # Thermocycler configuration
    THERMOCYCLER_GENERATION: str = 'gen2'  # Default to Gen2 thermocycler


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
        'V2_8_TC': 'clip_template_TC_APIv2.8.py'
    },
    'MAGBEAD': {
        'V2_8': 'purification_template_APIv2.8.py',
        'V2_10': 'purification_template_APIv2.10.py'
    },
    'F_ASSEMBLY': {
        'V2_8': 'assembly_template_APIv2.8.py',
        'V2_8_TC': 'assembly_template_TC_APIv2.8.py'
    },
    'TRANS_SPOT': {
        'V2_8': 'transformation_template_APIv2.8.py',
        'V2_8_TC': 'transformation_template_TC_APIv2.8.py',
        'V2_10_TC': 'transformation_template_TC_APIv2.10.py'
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
        'V2_8_TC': 'D_transformation_ot2_Thermocycler_APIv2.8.py',
        'V2_10_TC': 'D_transformation_ot2_Thermocycler_APIv2.10.py'
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
            self.SOURCE_POSITIONS = ['1', '2']  # Note: thermocycler protocols use slots 7, 8, 10, 11
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
    parser_nogui.add_argument('--thermocycler_gen',
                              help='Thermocycler generation (gen1 or gen2). Default: gen2',
                              type=str, default='gen2', choices=['gen1', 'gen2'])
    
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
        'thermocycler_gen': 'gen2',  # Default to gen2
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
        user_inputs['thermocycler_gen'] = dnabotinst.thermocycler_gen
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
            'thermocycler_gen': args.thermocycler_gen,
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
    # Copy input CSVs to output directory
    try:
        shutil.copy(construct_path, timestamped_output_dir)
        for src in user_config['sources_paths']:
            shutil.copy(src, timestamped_output_dir)
        print("✓ Input CSVs copied to output directory.")
    except Exception as e:
        print(f"Warning: Failed to copy input CSVs to output directory: {e}")
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
        print("\n" + "=" * 60 + "  \n1. Loading and validating constructs...")
        constructs_dict = generate_constructs_list(user_config['construct_path'], user_config.get('keep_layout', True))
        print(f"✓ Generated {len(constructs_dict)} constructs")
        
        # Generate CLIP reactions first
        print("\n" + "=" * 60 + "  \n2. Generating CLIP reactions...")
        unique_clips_df = generate_unique_clips_df(constructs_dict)
        print(f"✓ Unique clips required: {len(unique_clips_df)}")
        
        # Validate construct data (now with unique_clips_df for accurate counting)
        validate_construct_data(constructs_dict, unique_clips_df)
        print("✓ Construct data validation passed")
        
        # Validate Clip data
        validate_clips_data(unique_clips_df)
        print("✓ Clip data validation passed")
        
        # Process source files
        print("\n" + "=" * 60 + "  \n3. Loading and validating source data...")
        sources_dict = generate_sources_dict(user_config['sources_paths'])
        print(f"✓ Processed {len(sources_dict)} components from {len(user_config['sources_paths'])} source files")

        # Check if all parts/linkers are at default concentration
        all_default_conc = True
        for component, data in sources_dict.items():
            conc = data[1]
            if conc and float(conc) != PROTOCOL_CONFIG.CLIP_DEFAULT_PART_CONC:
                all_default_conc = False
                break
        print(f"✓ All parts/linkers at default concentration: {all_default_conc}")

        # Validate source data
        validate_sources_data(sources_dict)
        print("✓ Source data validation passed")
        
        # Validate component availability
        print("\n" + "=" * 60 + "  \n4. Validating component availability...")
        validate_components_availability(constructs_dict, sources_dict)
        print("✓ All required components are available")
        
        # Validate tip usage
        print("\n" + "=" * 60 + "  \n5. Validating tip usage...")
        validate_tip_usage(constructs_dict, unique_clips_df)
        print("✓ Tip usage within limits")
        
        # Log processing summary
        log_processing_summary(constructs_dict, unique_clips_df, sources_dict)

        print('\n' + "=" * 60 + "  \n6. Calculating OT-2 variables...")
        
        # Generate CLIP dictionaries for OT-2 scripts using optimised assignment
        clips_dict_list, assembly_to_clip_mapping, optimised_clips_df = generate_optimised_clips_dict_list(constructs_dict, sources_dict, all_default_conc)
        print(f"✓ Generated {len(clips_dict_list)} clip plate(s)")
        
        # Calculate magbead sample distribution based on actual CLIP plates generated
        magbead_sample_list = [len(plate_dict) for plate_dict in clips_dict_list]
        magbead_sample_number_total = sum(magbead_sample_list)
        print(f"✓ Total magbead samples: {magbead_sample_number_total}")
        
        # Generate final assembly plans
        print(f"\nAssembly to clip mapping: {assembly_to_clip_mapping}")
        final_assembly_dict_list = generate_final_assembly_dict_list(constructs_dict, optimised_clips_df, assembly_to_clip_mapping)
        print(f"✓ Generated {len(final_assembly_dict_list)} assembly plate(s)")

        # Validate assembly reconstruction
        print("\n" + "=" * 60 + "  \n7. Validating assembly reconstruction for 5 random constructs...")
        validate_final_assembly_reconstruction(
            final_assembly_dict_list, 
            clips_dict_list,  # Always use the unified structure
            constructs_dict, 
            sources_dict
        )
        print("✓ Assembly reconstruction validation passed")

        print("=" * 60)
        print("✓ All processing completed successfully!")

        return {
            'constructs_dict': constructs_dict,
            'optimised_clips_df': optimised_clips_df,  # Use the optimised clips DataFrame
            'sources_dict': sources_dict,
            'magbead_sample_list': magbead_sample_list,
            'final_assembly_dict_list': final_assembly_dict_list,
            'assembly_to_clip_mapping': assembly_to_clip_mapping,
            'all_default_conc': all_default_conc,
            'unique_clips_df': unique_clips_df,
            'clips_dict_list': clips_dict_list,
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
    # Use the unified structure for script generation
    # Generate CLIP scripts
    print("\nGenerating CLIP scripts...")
    for clip_plate, sub_clip_dict in enumerate(data_structures['clips_dict_list']):
        # print(f"  Clip plate {clip_plate + 1}...")
        _generate_clip_scripts(sub_clip_dict, clip_plate, paths, data_structures.get('all_default_conc', False), user_config.get('thermocycler_gen', 'gen2'))

    # Generate magbead purification scripts
    print("\nGenerating magbead purification scripts...")
    for i, magbead_sample_number in enumerate(data_structures['magbead_sample_list']):
        # print(f"  Magbead plate {i + 1}...")
        _generate_magbead_scripts(magbead_sample_number, i, paths, user_config)
    
    # Generate final assembly scripts
    print("\nGenerating final assembly scripts...")
    for plate_number, final_assembly_dict in data_structures['final_assembly_dict_list'].items():
        # print(f"  Assembly plate {plate_number}...")
        _generate_assembly_scripts(final_assembly_dict, plate_number, paths, user_config.get('thermocycler_gen', 'gen2'))
    
    # Generate transformation scripts
    print("\nGenerating transformation scripts...")
    for plate_number, final_assembly_dict in data_structures['final_assembly_dict_list'].items():
        # print(f"  Transformation plate {plate_number}...")
        transformation_dict = generate_transformation_dict(final_assembly_dict, plate_number)
        _generate_transformation_scripts(transformation_dict, plate_number, paths, user_config.get('thermocycler_gen', 'gen2'))


def _generate_clip_scripts(sub_clip_dict: dict, clip_plate: int, paths: dict, all_default_conc: bool = False, thermocycler_gen: str = 'gen2') -> None:
    """Generate CLIP reaction scripts for a single plate using embedded parameterisation.
    If there are more than 48 clips, split into two scripts (a and b) for the same plate number.
    """
    template_dir = paths['template_dir']
    # Get all destination wells for this plate
    dest_wells = list(sub_clip_dict.keys())
    n_clips = len(dest_wells)
    # Split into halves if needed
    halves = []
    if n_clips > 48:
        halves.append((dest_wells[:48], 'a'))
        halves.append((dest_wells[48:], 'b'))
    else:
        halves.append((dest_wells, 'a'))
    for wells, half in halves:
        # Build a sub-dictionary for this half
        half_clip_dict = {w: sub_clip_dict[w] for w in wells}
        # Script name: e.g. A1a_clip_ot2_APIv2.8.py
        base_name = FILE_CONFIG.OUTPUT_FILES['CLIP']['V2_8']
        stage_letter = base_name[0]
        script_name = f"{stage_letter}{clip_plate+1}{half}{base_name[1:]}.py"
        _generate_clip_script_embedded(
            script_name,
            os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['CLIP']['V2_8']),
            half_clip_dict,
            all_default_conc
        )
        # Thermocycler version
        tc_base_name = FILE_CONFIG.OUTPUT_FILES['CLIP']['V2_8_TC']
        tc_script_name = f"{tc_base_name[0]}{clip_plate+1}{half}{tc_base_name[1:]}.py"
        _generate_clip_script_embedded(
            tc_script_name,
            os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['CLIP']['V2_8_TC']),
            half_clip_dict,
            all_default_conc,
            thermocycler_gen
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
    # print(f"Clips number: {clips_number}")
    
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
                              paths: Dict[str, str],
                              thermocycler_gen: str = 'gen2') -> None:
    """Generate final assembly scripts using embedded parameterisation."""
    template_dir = paths['template_dir']
    final_assembly_tipracks = calculate_final_assembly_tipracks(final_assembly_dict)
    
    # Generate standard assembly script with new naming convention
    assembly_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['F_ASSEMBLY']['V2_8'], plate_number)
    _generate_assembly_script_embedded(
        assembly_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['F_ASSEMBLY']['V2_8']),
        final_assembly_dict,
        final_assembly_tipracks
    )
    
    # Generate thermocycler assembly script with new naming convention
    assembly_tc_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['F_ASSEMBLY']['V2_8_TC'], plate_number)
    _generate_assembly_script_embedded(
        assembly_tc_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['F_ASSEMBLY']['V2_8_TC']),
        final_assembly_dict,
        final_assembly_tipracks,
        thermocycler_gen
    )


def generate_transformation_dict(final_assembly_dict: Dict[str, List], 
                               assembly_plate: int) -> Dict[str, Any]:
    """Generate transformation data dictionary from final assembly data.
    
    Args:
        final_assembly_dict: Final assembly dictionary with well -> (clips, plates) mapping
        assembly_plate: Assembly plate number
        
    Returns:
        Dictionary containing transformation data for the template
    """
    # Extract wells and create transformation data
    destination_wells = list(final_assembly_dict.keys())
    transformation_number = len(destination_wells)
    
    # Generate source wells and plates (these would come from the assembly plate)
    # For now, we'll use the destination wells as source wells and plate 1 as source plate
    source_wells = destination_wells.copy()
    source_plates = ['1'] * transformation_number  # All from plate 1 for now
    
    # Standard volumes for transformation
    dna_volumes = [3] * transformation_number  # 3µl DNA per transformation
    cell_volumes = [30] * transformation_number  # 30µl cells per transformation
    soc_volumes = [100] * transformation_number  # 100µl SOC per transformation
    plating_volumes = [100] * transformation_number  # 100µl plating volume
    
    return {
        'transformation_number': transformation_number,
        'source_wells': source_wells,
        'source_plates': source_plates,
        'destination_wells': destination_wells,
        'dna_volumes': dna_volumes,
        'cell_volumes': cell_volumes,
        'soc_volumes': soc_volumes,
        'plating_volumes': plating_volumes
    }


def _generate_transformation_scripts(transformation_dict: Dict[str, Any], 
                                   script_index: int, 
                                   paths: Dict[str, str],
                                   thermocycler_gen: str = 'gen2') -> None:
    """Generate transformation scripts using embedded parameterisation."""
    template_dir = paths['template_dir']
    
    # Generate transformation script with new naming convention
    transformation_script_name = _generate_script_name_with_number(FILE_CONFIG.OUTPUT_FILES['TRANS_SPOT']['V2_10_TC'], script_index)
    _generate_transformation_script_embedded(
        transformation_script_name,
        os.path.join(template_dir, FILE_CONFIG.TEMPLATE_FILES['TRANS_SPOT']['V2_10_TC']),
        transformation_dict,
        thermocycler_gen
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
        master_mix_df = generate_master_mix_df(
            sum([len(plate_dict) for plate_dict in data_structures['clips_dict_list']]),
            data_structures.get('all_default_conc', False)
        )
        
        # Generate source plate information
        sources_paths_df = generate_sources_paths_df(
            user_config['sources_paths'], 
            data_structures['sources_dict']
        )
        
        # Write CLIP run information
        # Use clips_dict_list for canonical clip data
        clip_reactions_df = pd.DataFrame([clip for plate_dict in data_structures['clips_dict_list'] for clip in plate_dict.values()])
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
        f.write(f'SOC column: {user_config["soc_column"]}\n')
        f.write(f'Thermocycler generation: {user_config.get("thermocycler_gen", "gen2")}')


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
    Only accepts the new format: Plate, Well, Linker 1, Part 1, ...
    Raises an error if the format is not correct.
    """
    print(f"\n... Loading constructs from: {path}")
    def process_construct(construct: List[str], construct_index: int) -> pd.DataFrame:
        def get_suffix_linker(linker: str) -> str:
            if linker.startswith('U'):
                return linker.split('-')[0] + '-S'
            return linker + "-S"
        if len(construct) < 2:
            raise ValueError(f"Construct {construct_index + 1} has insufficient components. "
                           f"Expected at least 2 (linker, part), got {len(construct)}")
        if len(construct) % 2 != 0:
            raise ValueError(f"Construct {construct_index + 1} has invalid structure. "
                           f"Expected an even number of components (linker, part, linker, part, ...), got {len(construct)}")
        linkers = construct[::2]
        parts = construct[1::2]
        if len(linkers) != len(parts):
            raise ValueError(f"Construct {construct_index + 1} has mismatched number of linkers and parts. "
                           f"Linkers: {len(linkers)}, Parts: {len(parts)}")
        clips_info = {'prefixes': [], 'parts': [], 'suffixes': []}
        for i, sequence in enumerate(construct):
            if i % 2 != 0:
                if not sequence.strip():
                    raise ValueError(f"Construct {construct_index + 1}, position {i}: Empty part found")
                clips_info['parts'].append(sequence.strip())
                prefix_linker = construct[i - 1].strip()
                if not prefix_linker:
                    raise ValueError(f"Construct {construct_index + 1}, position {i-1}: Empty prefix linker found")
                clips_info['prefixes'].append(prefix_linker + '-P')
                if i == len(construct) - 1:
                    suffix_linker = get_suffix_linker(construct[0].strip())
                else:
                    next_linker = construct[i + 1].strip()
                    if not next_linker:
                        raise ValueError(f"Construct {construct_index + 1}, position {i+1}: Empty suffix linker found")
                    suffix_linker = get_suffix_linker(next_linker)
                clips_info['suffixes'].append(suffix_linker)
        return pd.DataFrame.from_dict(clips_info)
    constructs_dict = {}
    valid_construct_index = 0
    try:
        with open(path, 'r') as csvfile:
            csv_reader_list = list(csv.reader(open(path, 'r')))
            if not csv_reader_list:
                raise ValueError("Constructs file is empty.")
            header_row = csv_reader_list[0]
            # Only accept format: Plate, Well, ...
            if not (len(header_row) >= 2 and header_row[0].strip().lower() == 'plate' and header_row[1].strip().lower() == 'well'):
                raise ValueError(f"Constructs file must use the format with header: Plate, Well, ...\nFound header: {header_row}")
            print(f"... Validated CSV structure: Plate, Well format detected")
            for index, row in enumerate(csv_reader_list):
                if index == 0:
                    continue
                plate_str = row[0].strip() if row[0] else None
                well_position = row[1] if len(row) > 1 else None
                construct_components = row[2:]
                
                # Validate plate number is provided and valid
                if not plate_str:
                    raise ValueError(f"Plate number is missing at row {index + 1}. Plate column cannot be empty.")
                try:
                    plate_number = int(plate_str)
                    if plate_number < 1:
                        raise ValueError(f"Plate number must be a positive integer, got {plate_number} at row {index + 1}")
                except ValueError as e:
                    raise ValueError(f"Plate number must be a valid integer, got '{plate_str}' at row {index + 1}")
                
                # Validate well position is provided and valid
                if not well_position:
                    raise ValueError(f"Well position is missing at row {index + 1}. Well column cannot be empty.")
                if not well_position.strip():
                    raise ValueError(f"Well position is empty at row {index + 1}. Well column cannot be blank.")
                try:
                    validate_well_format(well_position)
                except ValueError as e:
                    raise ValueError(f"Well position error at row {index + 1}: {str(e)}")
                construct_components = list(filter(None, construct_components))
                if not construct_components:
                    continue
                if keep_layout:
                    try:
                        construct_df = process_construct(construct_components, index)
                        order = valid_construct_index
                        constructs_dict[(order, plate_number, well_position)] = construct_df
                        valid_construct_index += 1
                    except ValueError as e:
                        raise ValueError(f"Error processing construct at Plate {plate_number}, Well {well_position} (row {index + 1}): {str(e)}")
                else:
                    try:
                        construct_df = process_construct(construct_components, index)
                        new_plate = (valid_construct_index // 96) + 1
                        new_well = tip_counter(valid_construct_index % 96)
                        order = valid_construct_index
                        constructs_dict[(order, new_plate, new_well)] = construct_df
                        valid_construct_index += 1
                    except ValueError as e:
                        raise ValueError(f"Error processing construct at row {index + 1}: {str(e)}")
        print(f"... Successfully loaded {len(constructs_dict)} constructs with layout preserved")
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
            - Deck position (first column) - must be one of the allowed deck positions
            - Well location (second column)
            - Component name (third column)
            - Concentration (optional, fourth column)
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
        ValueError: If source files are malformed or use invalid deck positions
    """
    sources_dict = {}
    
    for path in paths:
        print(f"\n... Loading source data from: {path}")
        
        try:
            with open(path, 'r') as csvfile:
                csv_reader_list = list(csv.reader(csvfile))
                
                if not csv_reader_list:
                    raise ValueError(f"Source file is empty: {path}")
                
                header_row = csv_reader_list[0]
                
                # Validate CSV structure - must have Deck position, Well, Component name columns
                if len(header_row) < 3:
                    raise ValueError(f"Source file must have at least 3 columns: Deck position, Well, Component name. "
                                   f"Found {len(header_row)} columns in {path}")
                
                # Check header format (case-insensitive)
                header_lower = [col.strip().lower() for col in header_row]
                expected_headers = ['deck position', 'well', 'component name']
                
                for i, expected in enumerate(expected_headers):
                    if i < len(header_lower) and expected not in header_lower[i]:
                        raise ValueError(f"Source file header column {i+1} should be '{expected}' but found '{header_row[i]}' in {path}")
                
                print(f"... Validated CSV structure: Deck position, Well, Component name format detected")
                
                for index, row in enumerate(csv_reader_list):
                    if index == 0:  # Skip header row
                        continue
                    
                    # Validate minimum required data
                    if len(row) < 3:
                        raise ValueError(f"Row {index + 1}: Insufficient data. Need at least deck position, well, and component name.")
                    
                    # Extract and validate deck position
                    deck_position = row[0].strip() if row[0] else None
                    if not deck_position:
                        raise ValueError(f"Deck position is missing at row {index + 1}. Deck position column cannot be empty.")
                    
                    # Validate deck position is one of the allowed positions
                    if deck_position not in DECK_CONFIG.SOURCE_POSITIONS:
                        raise ValueError(f"Invalid deck position '{deck_position}' at row {index + 1}. "
                                       f"Allowed positions are: {', '.join(DECK_CONFIG.SOURCE_POSITIONS)}")
                    
                    # Extract and validate well position
                    well_position = row[1].strip() if row[1] else None
                    if not well_position:
                        raise ValueError(f"Well position is missing at row {index + 1}. Well column cannot be empty.")
                    if not well_position.strip():
                        raise ValueError(f"Well position is empty at row {index + 1}. Well column cannot be blank.")
                    try:
                        validate_well_format(well_position)
                    except ValueError as e:
                        raise ValueError(f"Well position error at row {index + 1}: {str(e)}")
                    
                    # Extract component name
                    component_name = row[2].strip() if row[2] else ""
                    
                    # Skip rows with empty component names (blank wells)
                    if not component_name:
                        continue
                    
                    # Extract concentration and additional data
                    concentration = row[3].strip() if len(row) > 3 and row[3] else ""
                    additional_data = row[4:] if len(row) > 4 else []
                    
                    # Build the data tuple: (well, concentration, deck_position, additional_data...)
                    csv_values = [well_position, concentration, deck_position] + additional_data
                    
                    # Check for duplicate components
                    if component_name in sources_dict:
                        raise ValueError(f"Duplicate component '{component_name}' found in source files. "
                                       f"First occurrence in file {paths.index(path) + 1}, "
                                       f"second occurrence in file {paths.index(path) + 1}")
                    
                    sources_dict[component_name] = tuple(csv_values)
            
            print(f"... Successfully loaded {len([k for k in sources_dict.keys() if k in sources_dict])} components from file {paths.index(path) + 1}")
            
        except FileNotFoundError:
            raise FileNotFoundError(f"Source file not found: {path}")
        except Exception as e:
            print(f"[ERROR] Failed to load source file {path}: {str(e)}")
            raise
    
    print(f"\n... Total components loaded: {len(sources_dict)}")
    return sources_dict


def generate_unique_clips_df(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], start_plate_idx: int = 1) -> pd.DataFrame:
    """Generates a dataframe containing information about all unique CLIP reactions."""
    def count_unique_clips(clips_df: pd.DataFrame, merged_construct_dfs: pd.DataFrame) -> pd.DataFrame:
        clip_count = np.zeros(len(clips_df.index))
        for i, unique_clip in clips_df.iterrows():
            for _, clip in merged_construct_dfs.iterrows():
                if unique_clip.equals(clip):
                    clip_count[i] += 1
        clip_count = clip_count // PROTOCOL_CONFIG.ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP + 1
        clips_df['number'] = [int(i) for i in clip_count.tolist()]
        return clips_df
    valid_constructs = list(constructs_dict.values())
    merged_construct_dfs = pd.concat(valid_constructs, ignore_index=True)
    unique_clips_df = merged_construct_dfs.drop_duplicates().reset_index(drop=True)
    unique_clips_df = count_unique_clips(unique_clips_df, merged_construct_dfs)
    # Add columns for Clip_Well and plate locations
    unique_clips_df['Clip_Well'] = pd.Series(['0'] * len(unique_clips_df.index), index=unique_clips_df.index)
    unique_clips_df['plate'] = pd.Series(['0'] * len(unique_clips_df.index), index=unique_clips_df.index)
    clip_count = 0
    for unique_clip_count, clip_number in unique_clips_df['number'].items():
        clip_wells = []
        plates = []         
        for well in range(clip_count, clip_count + clip_number):
            clip_wells.append(tip_counter(well % 96))
            plates.append(start_plate_idx + well//96)
        unique_clips_df.at[unique_clip_count, 'Clip_Well'] = tuple(clip_wells)
        unique_clips_df.at[unique_clip_count, 'plate'] = tuple(plates)
        clip_count += clip_number
    return unique_clips_df


def generate_clips_dict(unique_clips_df: pd.DataFrame, sources_dict: Dict[str, Tuple[str, ...]], all_default_conc: bool = False) -> Dict[str, dict]:
    """Generates a dictionary for CLIP reactions keyed by destination well, with all info for that well as a dict value."""
    # Calculate maximum part volume based on total reaction volume
    max_part_vol = PROTOCOL_CONFIG.CLIP_VOL - (
        PROTOCOL_CONFIG.CLIP_T4_BUFF_VOL + 
        PROTOCOL_CONFIG.CLIP_BSAI_VOL + 
        PROTOCOL_CONFIG.CLIP_T4_LIG_VOL + 
        PROTOCOL_CONFIG.CLIP_MAST_WATER + 2
    )
    clips_dict = {}
    # For each row in unique_clips_df, assign a destination well (A1, A2, ...)

    for idx, clip_info in unique_clips_df.iterrows():
        wells = clip_info['Clip_Well'] if isinstance(clip_info['Clip_Well'], (tuple, list)) else [clip_info['Clip_Well']]
        plates = clip_info['plate'] if isinstance(clip_info['plate'], (tuple, list)) else [clip_info['plate']]
        for rep, dest_well in enumerate(wells):
            plate_val = plates[rep] if rep < len(plates) else plates[0]
            prefix_linker = clip_info['prefixes'].strip()
            suffix_linker = clip_info['suffixes'].strip()
            part = clip_info['parts'].strip()
            prefix_well = sources_dict[prefix_linker][0]
            prefix_plate = normalize_source_data(sources_dict[prefix_linker])[2]
            suffix_well = sources_dict[suffix_linker][0]
            suffix_plate = normalize_source_data(sources_dict[suffix_linker])[2]
            part_well = sources_dict[part][0]
            part_plate = normalize_source_data(sources_dict[part])[2]
            if not sources_dict[part][1]:
                part_conc = PROTOCOL_CONFIG.CLIP_DEFAULT_PART_CONC
            else:
                part_conc = float(sources_dict[part][1])
            part_vol = round(PROTOCOL_CONFIG.CLIP_PART_PER_CLIP / float(part_conc), 1)
            part_vol = max(PROTOCOL_CONFIG.CLIP_MIN_VOL, min(part_vol, max_part_vol))
            if all_default_conc:
                water_vol = 0.0
            else:
                water_vol = max_part_vol - part_vol
            clips_dict[dest_well] = {
                'prefix_linker': prefix_linker,
                'prefix_source_well': prefix_well,
                'prefix_source_plate': prefix_plate,
                'part': part,
                'part_source_well': part_well,
                'part_source_plate': part_plate,
                'suffix_linker': suffix_linker,
                'suffix_source_well': suffix_well,
                'suffix_source_plate': suffix_plate,
                'Clip_Well': dest_well,
                'plate': plate_val,
                'part_vol': float(part_vol),
                'water_vol': float(water_vol)
            }
    return clips_dict


def generate_final_assembly_dict(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], unique_clips_df, assembly_to_clip_mapping=None):
    """Using constructs_dict and unique_clips_df, returns keys and values for a 
    dictionary of final assemblies; with keys defining destination plate 
    well positions, and values indicating which clip reaction wells are used.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing CLIP reactions
        unique_clips_df: DataFrame containing all unique CLIP reactions with their locations
        assembly_to_clip_mapping: Dictionary mapping assembly plate numbers to list of clip plate numbers
                                 If None, uses all available clips (legacy behavior)
    """

    final_assembly_dict_keys = []
    final_assembly_dict_values = []

    clips_count = np.zeros(len(unique_clips_df.index))

    # Process constructs in the order they appear in the CSV (using the order field in the tuple)
    # This ensures CLIP assignment order matches the original system
    sorted_positions = sorted(constructs_dict.keys(), key=lambda x: x[0])  # Sort by order (first element)
    
    for position in sorted_positions:
        construct_df = constructs_dict[position]
        construct_well_list = []
        construct_plate_list = []
        
        # Get the assembly plate number from the position
        order, assembly_plate, well = position
        
        # Filter unique_clips_df to only include clips from the appropriate clip plates for this assembly plate
        if assembly_to_clip_mapping and assembly_plate in assembly_to_clip_mapping:
            allowed_clip_plates = assembly_to_clip_mapping[assembly_plate]
            # Create a mask for clips that are on the allowed clip plates
            clip_plate_mask = unique_clips_df['plate'].apply(lambda x: any(plate in allowed_clip_plates for plate in x))
            filtered_clips_df = unique_clips_df[clip_plate_mask]
        else:
            # Legacy behavior: use all clips
            filtered_clips_df = unique_clips_df

        for _, clip in construct_df.iterrows():                                     # for each clip in construct
            clip_info = filtered_clips_df[(filtered_clips_df['prefixes'] == clip['prefixes']) &       # find clips in filtered_clips_df with the correct parts required for this clip
                                         (filtered_clips_df['parts'] == clip['parts']) &
                                         (filtered_clips_df['suffixes'] == clip['suffixes'])]
            
            if clip_info.empty:
                # If no clip found in filtered set, fall back to full unique_clips_df (for backward compatibility)
                clip_info = unique_clips_df[(unique_clips_df['prefixes'] == clip['prefixes']) &
                                   (unique_clips_df['parts'] == clip['parts']) &
                                   (unique_clips_df['suffixes'] == clip['suffixes'])]
            
            clip_num = int(clip_info.index[0])                                      # row index of this clip in unique_clips_df
            clip_wells = clip_info.at[clip_num, 'Clip_Well']                         # list of all mag_wells for this clip
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
        destination_key = (assembly_plate, well)  # Use tuple of (plate, well) as key
        final_assembly_dict_keys.append(destination_key)
        final_assembly_dict_values.append([construct_well_list, construct_plate_list])

    return final_assembly_dict_keys, final_assembly_dict_values                     # return list of dict keys and values


def generate_final_assembly_dict_list(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                     unique_clips_df: pd.DataFrame, 
                                     assembly_to_clip_mapping: Dict[int, List[int]] = None) -> Dict[int, Dict[str, List]]:
    """
    Generate a list of assembly dictionaries, each representing a subset of constructs
    that can be assembled on a single plate while respecting tip and well constraints.
    
    This function ensures that:
    1. No more than MAX_ASSEMBLIES_PER_PLATE constructs are assembled per plate
    2. Tip usage doesn't exceed available tiprack capacity
    3. Each construct gets the correct destination well location
    4. CLIP wells are not reused across different assembly plates
    5. Each assembly plate only uses clips from its assigned clip plates (if mapping provided)
    
    Args:
        constructs_dict: Dictionary mapping construct positions (order, plate, well) to DataFrames containing CLIP reactions
        unique_clips_df: DataFrame containing all unique CLIP reactions with their locations
        assembly_to_clip_mapping: Dictionary mapping assembly plate numbers to list of clip plate numbers
        
    Returns:
        Dictionary where outer keys are plate numbers and inner dictionaries map destination wells to 
        [clip_wells_list, clip_plates_list] for the constructs in that plate
        
    Raises:
        ValueError: If the number of constructs exceeds protocol limits
    """
    # Count total constructs
    total_constructs = len(constructs_dict)
    
    # Generate the complete assembly plan for all constructs
    assembly_keys, assembly_values = generate_final_assembly_dict(constructs_dict, unique_clips_df, assembly_to_clip_mapping)
    
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
    
    # print(f"Tipracks calculated: {tipracks_needed}, Maximum allowed: {max_allowed_tipracks}")
    
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


def _generate_clip_script_embedded(ot2_script_path: str, template_path: str, clips_dict: dict, all_default_conc: bool = False, thermocycler_gen: str = 'gen2') -> None:
    """Generate CLIP script using embedded parameterisation.
    Embeds the new per-well dictionary structure directly.
    """
    def convert_numpy_types(obj):
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
        
        # Replace thermocycler generation setting (only for thermocycler templates)
        if 'TC' in template_path:
            modified_protocol = modified_protocol.replace(
                "thermocycler_gen = 'gen2'",
                f"thermocycler_gen = '{thermocycler_gen}'"
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


def _generate_assembly_script_embedded(ot2_script_path: str, template_path: str, 
                                             final_assembly_dict: Dict[str, List], tiprack_num: int, thermocycler_gen: str = 'gen2') -> None:
    """Generate assembly script using embedded parameterisation.
    
    This function replaces the JSON file loading code in the template with embedded JSON data,
    similar to how embedded parameterisation works in templates.
    
    Args:
        ot2_script_path (str): Path where the OT-2 script will be written
        template_path (str): Path to the template file
        final_assembly_dict (Dict[str, List]): Final assembly data dictionary
        tiprack_num (int): Number of tipracks needed
        
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
        
        # Replace thermocycler generation setting (only for thermocycler templates)
        if 'TC' in template_path:
            modified_protocol = modified_protocol.replace(
                "thermocycler_gen = 'gen2'",
                f"thermocycler_gen = '{thermocycler_gen}'"
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


def _generate_transformation_script_embedded(ot2_script_path: str, template_path: str, 
                                           transformation_dict: Dict[str, Any], thermocycler_gen: str = 'gen2') -> None:
    """Generate transformation script with embedded parameterisation.
    
    This function replaces the JSON file loading code in the template with embedded JSON data,
    similar to how embedded parameterisation works in templates.
    
    Args:
        ot2_script_path (str): Path where the OT-2 script will be written
        template_path (str): Path to the template file
        transformation_dict (Dict[str, Any]): Transformation data dictionary
        
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
        
        # Convert to JSON string with compact formatting
        converted_transformation_dict = convert_numpy_types(transformation_dict)
        
        # Create a more compact and readable format for the transformation data
        transformation_dict_lines = []
        transformation_dict_lines.append("transformation_dict = {")
        
        for key, value in converted_transformation_dict.items():
            if isinstance(value, list):
                # Format lists in a compact way
                if all(isinstance(item, str) for item in value):
                    # String list - format as ['item1', 'item2', ...]
                    items = [f'"{item}"' for item in value]
                    value_str = f'[{", ".join(items)}]'
                else:
                    # Numeric list - format as [1, 2, 3, ...]
                    value_str = f'[{", ".join(map(str, value))}]'
            else:
                # Simple value
                value_str = f'"{value}"' if isinstance(value, str) else str(value)
            
            transformation_dict_lines.append(f'    "{key}": {value_str},')
        
        # Remove trailing comma from last line
        if transformation_dict_lines:
            transformation_dict_lines[-1] = transformation_dict_lines[-1].rstrip(',')
        
        transformation_dict_lines.append("}")
        
        # Create the compact transformation data string
        compact_transformation_data = "\n".join(transformation_dict_lines)
        
        # Replace the JSON file loading code with embedded compact data
        modified_protocol = template_content.replace(
            "with open('transformation_data.json') as f:\n    transformation_dict = json.load(f)",
            compact_transformation_data
        )
        
        # Replace thermocycler generation setting (only for thermocycler templates)
        if 'TC' in template_path:
            modified_protocol = modified_protocol.replace(
                "thermocycler_gen = 'gen2'",
                f"thermocycler_gen = '{thermocycler_gen}'"
            )
        
        # Write the modified protocol
        with open(ot2_script_path, 'w') as f:
            f.write(modified_protocol)
            
        print(f"    ✓ {os.path.basename(ot2_script_path)}")
        
    except Exception as e:
        print(f"\nError generating transformation script {os.path.basename(ot2_script_path)}:")
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


def generate_master_mix_df(clip_number, all_default_conc):
    """Generates a dataframe detailing the components required in the clip reaction master mix.
    If all_default_conc is True, use the optimised mastermix (no separate water transfer).
    Always outputs both per-reaction (1x) and total (Nx) columns.
    """
    if all_default_conc:
        COMPONENTS = {'Component': [
            'Promega T4 DNA Ligase buffer, 10X',
            'Water',
            'NEB BsaI-HFv2',
            'Promega T4 DNA Ligase'
        ]}
        PER_REACTION = [3.0, 22.5, 1.0, 0.5]
        VOL_COLUMN_1X = 'Volume per reaction (uL)'
        VOL_COLUMN_TOTAL = 'Total volume (uL)'
        master_mix_df = pd.DataFrame.from_dict(COMPONENTS)
        master_mix_df[ VOL_COLUMN_1X ] = PER_REACTION
        # Add dead volume as before
        multiplier = float(clip_number + PROTOCOL_CONFIG.CLIP_DEAD_VOL / PROTOCOL_CONFIG.CLIP_VOL)
        master_mix_df[ VOL_COLUMN_TOTAL ] = [round(x * multiplier, 2) for x in PER_REACTION]
        return master_mix_df
    else:
        COMPONENTS = {'Component': ['Promega T4 DNA Ligase buffer, 10X',
                                    'Water', 'NEB BsaI-HFv2',
                                    'Promega T4 DNA Ligase']}
        VOL_COLUMN_1X = 'Volume per reaction (uL)'
        VOL_COLUMN_TOTAL = 'Total volume (uL)'
        master_mix_df = pd.DataFrame.from_dict(COMPONENTS)
        # Calculate volumes for each component
        clip_vol = PROTOCOL_CONFIG.CLIP_VOL
        dead_vol = PROTOCOL_CONFIG.CLIP_DEAD_VOL
        t4_buff_vol = PROTOCOL_CONFIG.CLIP_T4_BUFF_VOL
        mast_water = PROTOCOL_CONFIG.CLIP_MAST_WATER
        bsai_vol = PROTOCOL_CONFIG.CLIP_BSAI_VOL
        t4_lig_vol = PROTOCOL_CONFIG.CLIP_T4_LIG_VOL
        per_reaction = [t4_buff_vol, mast_water, bsai_vol, t4_lig_vol]
        multiplier = float(clip_number + dead_vol/clip_vol)
        master_mix_df[VOL_COLUMN_1X] = per_reaction
        master_mix_df[VOL_COLUMN_TOTAL] = [round(x * multiplier, 2) for x in per_reaction]
        return master_mix_df


def generate_sources_paths_df(paths, sources_dict):
    """Generates a dataframe detailing source plate information.

    Args:
        paths (list): list of strings specifying paths to source plates.
        sources_dict (dict): dictionary mapping components to their source data.

    Returns:
        pd.DataFrame: DataFrame with deck positions, source plate names, and paths
    """
    source_plates_dict = {'Deck position': [], 'Source plate': [], 'Path': []}
    
    # Extract unique deck positions from the sources_dict
    deck_positions_used = set()
    for source_data in sources_dict.values():
        if len(source_data) >= 3:
            deck_positions_used.add(source_data[2])  # deck_position is at index 2
    
    # Create entries for each deck position used
    for deck_position in sorted(deck_positions_used):
        # Find which file contains this deck position
        file_index = None
        for i, path in enumerate(paths):
            # Check if this file contains the deck position
            with open(path, 'r') as csvfile:
                csv_reader = csv.reader(csvfile)
                for row in csv_reader:
                    if row and len(row) >= 1 and row[0].strip() == deck_position:
                        file_index = i
                        break
                if file_index is not None:
                    break
        
        if file_index is not None:
            source_plates_dict['Deck position'].append(deck_position)
            source_plates_dict['Source plate'].append(os.path.basename(paths[file_index]))
            source_plates_dict['Path'].append(paths[file_index])
    
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
    
    This function handles the new source data format where data_tuple contains:
    (well, concentration, deck_position, additional_data...)
    
    Args:
        data_tuple: Tuple or list containing source data
        
    Returns:
        Normalized 3-element tuple: (well, concentration, deck_position)
        
    Examples:
        >>> normalize_source_data(('A8', '', '1'))
        ('A8', '', '1')
        >>> normalize_source_data(('A8', '100', '2', 'extra_data'))
        ('A8', '100', '2')
    """
    if isinstance(data_tuple, list):
        data_tuple = data_tuple[0]
    
    if len(data_tuple) >= 3:
        # New format: (well, concentration, deck_position, ...)
        return (data_tuple[0], data_tuple[1], data_tuple[2])
    else:
        raise ValueError(f"Expected at least 3 elements (well, concentration, deck_position), got {len(data_tuple)}")


def validate_construct_data(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], unique_clips_df: pd.DataFrame) -> None:
    """
    Validate that construct data is properly formatted and within protocol limits.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing construct information
        unique_clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If constructs are malformed or exceed limits
    """
    if not constructs_dict:
        raise ValueError("No constructs found in input file")
    
    # Use unique_clips_df for accurate clip counting
    total_clips = unique_clips_df['number'].sum()
    
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


def validate_clips_data(unique_clips_df: pd.DataFrame) -> None:
    """
    Validate that CLIP reaction data is properly formatted.
    
    Args:
        unique_clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If CLIP data is malformed
    """
    if unique_clips_df.empty:
        raise ValueError("No CLIP reactions found")
    
    required_columns = ['prefixes', 'parts', 'suffixes', 'number', 'Clip_Well', 'plate']
    missing_columns = [col for col in required_columns if col not in unique_clips_df.columns]
    if missing_columns:
        raise ValueError(f"CLIP DataFrame missing required columns: {missing_columns}")
    
    # Check for negative or zero reaction numbers
    invalid_numbers = unique_clips_df[unique_clips_df['number'] <= 0]
    if not invalid_numbers.empty:
        raise ValueError(f"Found {len(invalid_numbers)} CLIP reactions with invalid numbers (≤0): {invalid_numbers.index.tolist()}")
    
    # Check for empty values in required columns
    for col in ['prefixes', 'parts', 'suffixes']:
        empty_values = unique_clips_df[unique_clips_df[col].isna() | (unique_clips_df[col] == '')]
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


def validate_tip_usage(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], unique_clips_df: pd.DataFrame) -> None:
    """
    Validate that the experiment doesn't exceed available tip capacity.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames containing construct information
        unique_clips_df: DataFrame containing CLIP reaction information
        
    Raises:
        ValueError: If tip usage exceeds limits
    """
    # Calculate total tips needed
    total_clip_tips = unique_clips_df['number'].sum()
    
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
                          unique_clips_df: pd.DataFrame,
                          sources_dict: Dict[str, Tuple[str, ...]]) -> None:
    """
    Log a summary of the processing results for user verification.
    
    Args:
        constructs_dict: Dictionary mapping construct positions to DataFrames
        unique_clips_df: CLIP reactions DataFrame
        sources_dict: Source locations dictionary
    """
    print("\n" + "="*60)
    print("PROCESSING SUMMARY")
    print("="*60)
    print(f"Total constructs: {len(constructs_dict)}")
    print(f"Total CLIP reactions: {unique_clips_df['number'].sum()}")
    print(f"Unique CLIP reactions: {len(unique_clips_df)}")
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


def generate_clip_dicts_for_constructs(constructs_dict, sources_dict, all_default_conc=False, start_plate_idx=1):
    """
    Given a constructs_dict, generate the unique_clips_df DataFrame, split into sets of up to 96 (whole plates), and return a list of clip dicts, a list of plate numbers, and the combined DataFrame.
    
    Args:
        constructs_dict: Dictionary of constructs
        sources_dict: Dictionary of sources
        all_default_conc: Whether all concentrations are default
        start_plate_idx: Starting plate index for numbering (default 1)
    """
    unique_clips_df = generate_unique_clips_df(constructs_dict, start_plate_idx)
    total_clips = unique_clips_df['number'].sum()
    clip_dict_list = []
    plate_numbers = []
    # Split into sets of 96
    start = 0
    plate_idx = start_plate_idx
    long_clip_df = []
    for idx, row in unique_clips_df.iterrows():
        for _ in range(row['number']):
            long_clip_df.append(row)
    long_clip_df = pd.DataFrame(long_clip_df)
    while start < total_clips:
        end = min(start + 96, total_clips)
        chunk_df = long_clip_df.iloc[start:end, :]
        clip_dict = generate_clips_dict(chunk_df, sources_dict, all_default_conc)
        clip_dict_list.append(clip_dict)
        plate_numbers.append(plate_idx)
        start = end
        plate_idx += 1
    return clip_dict_list, plate_numbers, unique_clips_df


def generate_optimised_clips_dict_list(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame], 
                                      sources_dict: Dict[str, Tuple[str, ...]],
                                      all_default_conc: bool = False):
    """
    New staged logic for clip distribution:
    1. Calculate total clips, error if >= max allowed
    2. Optimise distribution (try to split by assembly plates, as per user logic)
    3. Generate clip dicts for each plate (up to 96 clips at a time)
    4. Validate assignments
    """
    def count_clips_for_constructs(constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame]) -> int:
        """Count total number of clip wells needed for a set of constructs.
        This uses the same logic as generate_unique_clips_df to calculate how many wells/copies
        of each clip are needed, accounting for ASSEMBLY_FINAL_ASSEMBLIES_PER_CLIP.
        """
        if not constructs_dict:
            return 0
        unique_clips_df = generate_unique_clips_df(constructs_dict)
        return int(unique_clips_df['number'].sum())

    # --- Stage 1: Initial check ---
    total_clips = count_clips_for_constructs(constructs_dict)
    max_clips_allowed = PROTOCOL_CONFIG.ASSEMBLY_MAX_CLIPS_TOTAL * 2
    if total_clips >= max_clips_allowed:
        raise ValueError(f"Total number of CLIP reactions ({total_clips}) is greater than or equal to the maximum allowed ({max_clips_allowed}). Reduce the number of constructs or simplify construct designs.")

        # Get all assembly plates
    assembly_plates = sorted({plate for (_, plate, _) in constructs_dict.keys()})

    # --- Stage 2: Clip distribution optimisation ---
    # Default: no split
    best_constructs_dicts = [constructs_dict]

    if total_clips < 96 or len(assembly_plates) == 1:
        # Simple case: all fits on one plate or only one assembly plate
        print(f"No clip optimisation required. All assemblies: {assembly_plates}, Total clips: {total_clips}")
    else:
        # Try to split by assembly plates
        print(f"Attempting to optimise clip distribution...")
        plates = assembly_plates.copy()
        while True:
            mid = len(plates) // 2
            first_half_plates = plates[:mid]
            second_half_plates = plates[mid:]
            if not first_half_plates or not second_half_plates:
                print(f"Cannot split further. Plates: {plates}")
                break  # Can't split further
            # Build construct dicts for each half
            first_half = {k: v for k, v in constructs_dict.items() if k[1] in first_half_plates}
            second_half = {k: v for k, v in constructs_dict.items() if k[1] in second_half_plates}
            first_clips = count_clips_for_constructs(first_half)
            second_clips = count_clips_for_constructs(second_half)
            print(f"Clip plate 1. Associated assembly plates: {first_half_plates}, Constructs: {len(first_half)}, Clips: {first_clips}")
            print(f"Clip plate 2. Associated assembly plates: {second_half_plates}, Constructs: {len(second_half)}, Clips: {second_clips}")
            if first_clips <= 96 and second_clips <= 96:
                best_constructs_dicts = [first_half, second_half]
                print(f"Successful split found for optimised clip distribution.")
                break
            elif first_clips > 96 and second_clips > 96:
                print(f"Both halves exceed 96 clips. Reverting to non-optimised clip distribution.")
                break
            else:
                # Move the middle plate to the other half and try again
                if first_clips > 96:
                    move_plate = first_half_plates[-1]
                    print(f"Redistributing assembly plate {move_plate} from first to second clip plate.")
                    plates.remove(move_plate)
                else:
                    move_plate = second_half_plates[0]
                    print(f"Redistributing assembly plate {move_plate} from second to first clip plate.")
                    plates.remove(move_plate)
        # If no good split found, best_constructs_dicts remains as [constructs_dict]

    # --- Stage 3: Clip dict generation ---
    clips_dict_list = []
    all_clips_df = []
    assembly_to_clip_mapping = {}
    current_plate_idx = 1
    for i, sub_constructs_dict in enumerate(best_constructs_dicts):
        sub_clip_dicts, plate_numbers, sub_unique_clips_df = generate_clip_dicts_for_constructs(sub_constructs_dict, sources_dict, all_default_conc, current_plate_idx)
        clips_dict_list.extend(sub_clip_dicts)
        all_clips_df.append(sub_unique_clips_df)
        # Map assembly plates in this subdict to the corresponding plate numbers
        plates_in_subdict = sorted({k[1] for k in sub_constructs_dict.keys()})
        for plate in plates_in_subdict:
            # Each assembly plate maps to the list of clip plate numbers used for this subdict
            assembly_to_clip_mapping[plate] = plate_numbers
        current_plate_idx += len(plate_numbers)
    # Combine all clips DataFrames
    if all_clips_df:
        optimised_clips_df = pd.concat(all_clips_df, ignore_index=True)
    else:
        optimised_clips_df = pd.DataFrame()

            # --- Stage 4: Validation ---
        validate_clip_assignments(constructs_dict, optimised_clips_df, sources_dict, clips_dict_list)

    print(f"\nOptimisation complete:")
    print(f"  Total clip plates: {len(clips_dict_list)}")
    print(f"  Assembly to clip mapping: {assembly_to_clip_mapping}")

    return clips_dict_list, assembly_to_clip_mapping, optimised_clips_df 


def validate_clip_assignments(constructs_dict, clips_df, sources_dict, clips_dict_list=None):
    """
    Validate that the generated clip assignments can reconstruct all constructs.
    Checks that for each construct, all required clips are present in the clips_df.
    Also checks protocol limits: total clips <= 192, total plates <= 2.
    Raises a ValueError if any construct cannot be reconstructed or if protocol limits are exceeded.
    """
    # Protocol limit checks
    if clips_df is not None:
        total_clips_final = int(clips_df['number'].sum()) if not clips_df.empty else 0
        if total_clips_final > 192:
            raise ValueError(f"Total number of CLIP reactions ({total_clips_final}) exceeds the maximum allowed (192). Reduce the number of constructs or simplify construct designs.")
    if clips_dict_list is not None:
        total_plates_final = len(clips_dict_list)
        if total_plates_final > 2:
            raise ValueError(f"Total number of CLIP plates ({total_plates_final}) exceeds the maximum allowed (2). Reduce the number of constructs or simplify construct designs.")

    for position, construct_df in constructs_dict.items():
        for _, clip in construct_df.iterrows():
            matches = clips_df[
                (clips_df['prefixes'] == clip['prefixes']) &
                (clips_df['parts'] == clip['parts']) &
                (clips_df['suffixes'] == clip['suffixes'])
            ]
            if matches.empty:
                raise ValueError(
                    f"Clip assignment validation failed: "
                    f"Clip (prefix: {clip['prefixes']}, part: {clip['parts']}, suffix: {clip['suffixes']}) "
                    f"for construct at {position} not found in generated clips."
                )
    print("✓ Clip assignment validation passed")


def validate_final_assembly_reconstruction(assembly_dict: Dict[int, Dict[str, List]], 
                                   clips_dict_list: List[Dict[str, dict]], 
                                   constructs_dict: Dict[Tuple[int, int, str], pd.DataFrame],
                                   sources_dict: Dict[str, Tuple[str, ...]],
                                   sample_wells: List[Tuple[int, str]] = None) -> None:
    """
    Validate that final assemblies can correctly reconstruct the original constructs.
    For each sample well, print only: plate/well, (clip well, plate) pairs, clip details, and sequences with match status.
    """
    import random

    def collect_clip_components(clip_wells, clip_plates, clips_dict_list, well_to_component):
        components = []
        for i, (clip_well, clip_plate) in enumerate(zip(clip_wells, clip_plates)):
            clip_dict = clips_dict_list[clip_plate - 1]
            clip_info = clip_dict[clip_well]
            prefix_linker = find_component_by_well(clip_info['prefix_source_plate'], clip_info['prefix_source_well'], well_to_component)
            part = find_component_by_well(clip_info['part_source_plate'], clip_info['part_source_well'], well_to_component)
            suffix_linker = find_component_by_well(clip_info['suffix_source_plate'], clip_info['suffix_source_well'], well_to_component)
            components.append({
                'clip_well': clip_well,
                'clip_plate': clip_plate,
                'prefix_linker': prefix_linker,
                'prefix_source': (clip_info['prefix_source_plate'], clip_info['prefix_source_well']),
                'part': part,
                'part_source': (clip_info['part_source_plate'], clip_info['part_source_well']),
                'suffix_linker': suffix_linker,
                'suffix_source': (clip_info['suffix_source_plate'], clip_info['suffix_source_well'])
            })
        return components

    well_to_component = build_well_to_component_lookup(sources_dict)

    # Select sample wells if not provided
    if sample_wells is None:
        all_wells = [(plate_num, well) for plate_num, plate_dict in assembly_dict.items() for well in plate_dict.keys()]
        sample_wells = all_wells if len(all_wells) < 5 else random.sample(all_wells, 5)

    for plate_num, well in sample_wells:
        print(f"\nPlate {plate_num}, Well {well}:")
        if plate_num not in assembly_dict or well not in assembly_dict[plate_num]:
            raise ValueError(f"Construct at Plate {plate_num}, Well {well} not found in assembly dict")
        clip_info = assembly_dict[plate_num][well]
        if not clip_info or len(clip_info) != 2:
            raise ValueError(f"Invalid clip info format: {clip_info}")
        clip_wells, clip_plates = clip_info
        clip_components = collect_clip_components(clip_wells, clip_plates, clips_dict_list, well_to_component)
        print("  Clips:")
        for comp in clip_components:
            print(f"    Well {comp['clip_well']} (Plate {comp['clip_plate']})")
            print(f"      Prefix: {comp['prefix_linker']} [Plate {comp['prefix_source'][0]}, Well {comp['prefix_source'][1]}]")
            print(f"      Part:   {comp['part']} [Plate {comp['part_source'][0]}, Well {comp['part_source'][1]}]")
            print(f"      Suffix: {comp['suffix_linker']} [Plate {comp['suffix_source'][0]}, Well {comp['suffix_source'][1]}]")
        reconstructed_sequence = reconstruct_construct_sequence([[c['prefix_linker'], c['part'], c['suffix_linker']] for c in clip_components])
        # Find the original construct definition
        original_construct = None
        for (order, orig_plate, orig_well), construct_df in constructs_dict.items():
            if orig_plate == plate_num and orig_well == well:
                original_construct = construct_df
                break
        if original_construct is None:
            raise ValueError(f"Original construct definition not found for Plate {plate_num}, Well {well}")
        original_sequence = construct_df_to_sequence(original_construct)
        print(f"  Reconstructed sequence: {reconstructed_sequence}")
        print(f"  Original sequence:     {original_sequence}")
        if reconstructed_sequence == original_sequence:
            print("  MATCH")
        else:
            print("  MISMATCH")
            raise ValueError(
                f"Sequence mismatch for construct at Plate {plate_num}, Well {well}:\n"
                f"  Original: {original_sequence}\n"
                f"  Reconstructed: {reconstructed_sequence}"
            )
    print(f"\nAll {len(sample_wells)} sample constructs validated successfully!\n")


def reconstruct_construct_sequence(clip_components: List[List[str]]) -> List[str]:
    """
    Reconstruct the full construct sequence from clip components.
    Args:
        clip_components: List of [prefix_linker, part, suffix_linker] for each clip
    Returns:
        List representing the unified construct sequence [part1, linker1, part2, linker2, ...]
    """
    # Removed all debug print statements
    if not clip_components:
        return []
    sequence = []
    for i, clip in enumerate(clip_components):
        prefix_linker = clip[0]
        part = clip[1]
        suffix_linker = clip[2]
        if i == 0:
            sequence.extend([part, suffix_linker])
        else:
            previous_suffix = sequence[-1]
            previous_suffix_base = previous_suffix.replace('-S', '').replace('-P', '')
            current_prefix_base = prefix_linker.replace('-S', '').replace('-P', '')
            if previous_suffix_base == current_prefix_base:
                sequence[-1] = previous_suffix_base
                sequence.append(part)
                sequence.append(suffix_linker)
            else:
                sequence.append(prefix_linker)
                sequence.append(part)
                sequence.append(suffix_linker)
    if len(clip_components) > 1:
        first_prefix = clip_components[0][0]
        last_suffix = sequence[-1]
        first_prefix_base = first_prefix.replace('-S', '').replace('-P', '')
        last_suffix_base = last_suffix.replace('-S', '').replace('-P', '')
        if first_prefix_base == last_suffix_base:
            sequence[-1] = last_suffix_base
    return sequence


def construct_df_to_sequence(construct_df: pd.DataFrame) -> List[str]:
    """
    Convert a construct DataFrame to a sequence list.
    Args:
        construct_df: DataFrame with columns ['prefixes', 'parts', 'suffixes']
    Returns:
        List representing the construct sequence [part1, linker1, part2, linker2, ...]
    """
    # Removed all debug print statements
    sequence = []
    for i, (_, row) in enumerate(construct_df.iterrows()):
        prefix = row['prefixes'].replace('-P', '')
        part = row['parts']
        suffix = row['suffixes'].replace('-S', '')
        if not sequence:
            sequence.extend([part, suffix])
        else:
            current_suffix = sequence[-1]
            if current_suffix == prefix:
                sequence[-1] = prefix
                sequence.append(part)
                sequence.append(suffix)
            else:
                sequence.append(prefix)
                sequence.append(part)
                sequence.append(suffix)
    return sequence


# 1. Add utility function to convert DataFrame to list of dicts

def clips_df_to_dict_list(clips_df):
    """
    Convert the canonical clips DataFrame to a list of per-plate, per-well dictionaries.
    Each dict in the list corresponds to a plate, mapping well names to clip info dicts.
    """
    dict_list = []
    # Get all unique plate numbers (flatten if tuple/list)
    all_plates = set()
    for plates in clips_df['plate']:
        if isinstance(plates, (tuple, list)):
            all_plates.update(plates)
        else:
            all_plates.add(plates)
    for plate_num in sorted(all_plates, key=lambda x: int(x)):
        # Select all rows for this plate (handle tuple/list in 'plate' column)
        plate_rows = clips_df[clips_df['plate'].apply(lambda x: plate_num in (x if isinstance(x, (tuple, list)) else [x]))]
        plate_dict = {}
        for _, row in plate_rows.iterrows():
            wells = row['Clip_Well'] if isinstance(row['Clip_Well'], (tuple, list)) else [row['Clip_Well']]
            for well in wells:
                # Store a dict of all relevant info for this well
                plate_dict[well] = row.to_dict()
        dict_list.append(plate_dict)
    return dict_list

# 2. Refactor _process_input_files to only keep the optimised DataFrame
# (Remove all_long_clip_dfs, clips_dict_list, etc. from returned data structures)
# 3. Update all downstream code to use clips_df_to_dict_list(optimised_clips_df) where needed
# 4. Remove any code that stores or passes around per-plate DataFrames or dicts
# 5. Add comments to clarify the new approach

# ... existing code ...

def build_well_to_component_lookup(sources_dict):
    """
    Build a lookup dict mapping (plate, well) → component_name for fast reverse lookup.
    """
    lookup = {}
    for component_name, source_info in sources_dict.items():
        well = source_info[0]
        plate = normalize_source_data(source_info)[2]
        lookup[(plate, well)] = component_name
    return lookup

def find_component_by_well(plate: str, well: str, well_to_component: dict) -> str:
    """
    Find a component name by (plate, well) using a precomputed lookup dict.
    """
    key = (plate, well)
    if key in well_to_component:
        return well_to_component[key]
    raise ValueError(f"No component found at plate {plate}, well {well}")

if __name__ == '__main__':
    main()
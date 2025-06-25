import csv
import json
import argparse
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
import sys
import re

# Volume constants
FINAL_PLATE_VOLUME = 100  # µL
MIX_PLATE_VOLUME = 200    # µL
SUBCULTURE_TRANSFER = 10  # µL

# Reservoir capacity constants
RESERVOIR_CAPACITIES = {
    "Media reservoir": 20000,    # µL - 4ti0131_12_reservoir_21000ul
    "Supplement reservoir": 2000  # µL - 4ti0136_96_wellplate_2200ul
}

def calculate_supplement_volume(reps):
    """
    Calculate the volume of each supplement based on total volume needed.
    
    Args:
        reps (int): Number of replicate plates
        
    Returns:
        float: Volume of each supplement in µL
    """
    total_volume = (reps * (FINAL_PLATE_VOLUME - SUBCULTURE_TRANSFER) + 
                   (MIX_PLATE_VOLUME - SUBCULTURE_TRANSFER))
    return total_volume / 10  # 1/10 of total volume

def validate_well_format(well):
    """
    Validate that a well identifier is in the correct format (e.g., 'A1' or 'H12').
    
    Args:
        well (str): Well identifier to validate
        
    Returns:
        bool: True if well is valid, False otherwise
        
    Raises:
        ValueError: If well format is invalid
    """
    if not isinstance(well, str):
        raise ValueError(f"Well identifier must be a string, got {type(well)}")
    
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

def validate_unique_wells(wells, context="wells"):
    """
    Validate that all wells are unique.
    
    Args:
        wells (list): List of well identifiers
        context (str): Context for error message (e.g., "supplement map wells")
        
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

def load_supplement_map(sup_map_path):
    """
    Load the supplement and media map from CSV file.
    
    Args:
        sup_map_path (str): Path to supplement map CSV file
        
    Returns:
        dict: Dictionary mapping supplement/media names to well positions and plates
        
    Raises:
        ValueError: If CSV format is incorrect or contains invalid plate names
    """
    print(f"\n... Loading supplement map from: {sup_map_path}")
    supplement_map = {}
    valid_plates = ["Media reservoir", "Supplement reservoir"]
    
    try:
        with open(sup_map_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check required columns
            if not all(col in reader.fieldnames for col in ["Well", "Supplement", "Plate"]):
                raise ValueError("Supplement map CSV must contain 'Well', 'Supplement', and 'Plate' columns")
            
            # Track wells and plate combinations for uniqueness
            well_plate_combos = []
            
            for row in reader:
                plate_name = row["Plate"]
                if plate_name not in valid_plates:
                    raise ValueError(f"Invalid plate name: '{plate_name}'. Only 'Media reservoir' and 'Supplement reservoir' are allowed. Please check for extra spaces in your plate names.")
                
                # Validate well format
                validate_well_format(row["Well"])
                
                # Track well/plate combination
                well_plate_combos.append((row["Well"], plate_name))
                
                # Strip whitespace from supplement name
                supplement_name = row["Supplement"].strip()
                if supplement_name not in supplement_map:
                    supplement_map[supplement_name] = {
                        "wells": [row["Well"]],  # Store as list of wells
                        "plate": plate_name
                    }
                else:
                    # Add to existing wells list if same plate
                    if supplement_map[supplement_name]["plate"] == plate_name:
                        supplement_map[supplement_name]["wells"].append(row["Well"])
                    else:
                        raise ValueError(f"Supplement '{supplement_name}' found in multiple plates: {supplement_map[supplement_name]['plate']} and {plate_name}")
        
        # Validate unique well/plate combinations
        validate_unique_wells(well_plate_combos, "well/plate combinations in supplement map")
        
        print(f"... Successfully loaded supplement map with {len(supplement_map)} items")
        return supplement_map
        
    except Exception as e:
        print(f"[ERROR] Failed to load supplement map: {str(e)}")
        raise

def load_strain_map(strain_map_path):
    """
    Load the strain map from CSV file.
    
    Args:
        strain_map_path (str): Path to strain map CSV file
        
    Returns:
        dict: Dictionary mapping strain names to well positions
        
    Raises:
        ValueError: If CSV format is incorrect
    """
    print(f"\n... Loading strain map from: {strain_map_path}")
    strain_map = {}
    
    try:
        with open(strain_map_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check required columns
            if not all(col in reader.fieldnames for col in ["Well", "Strain"]):
                raise ValueError("Strain map CSV must contain 'Well' and 'Strain' columns")
            
            wells = []
            for row in reader:
                # Validate well format
                validate_well_format(row["Well"])
                wells.append(row["Well"])
                
                # Strip whitespace from strain name
                strain_name = row["Strain"].strip()
                strain_map[strain_name] = row["Well"]
        
        # Validate unique wells
        validate_unique_wells(wells, "wells in strain map")
        
        print(f"... Successfully loaded strain map with {len(strain_map)} items")
        return strain_map
        
    except Exception as e:
        print(f"[ERROR] Failed to load strain map: {str(e)}")
        raise

def calculate_total_tip_uses(experiment, strain_map):
    """
    Calculate the total number of tip uses for media component transfers and experiment entries.
    
    Args:
        experiment (list): List of experiment entries
        strain_map (dict): Mapping of strains to preculture wells
        
    Raises:
        ValueError: If total tip uses exceeds available tips (96)
    """
    # Count media component transfers
    media_components = set()
    for entry in experiment:
        # Add media type
        media_components.add(entry["media_type"])
        # Add supplements
        media_components.update(entry["supplements"])
        # Add water if needed
        if entry["media_type"][:2] == "2x" and entry["water_volume"] > 0:
            media_components.add("Water")
    
    # Total tip uses = number of media components + number of experiment entries
    total_tip_uses = len(media_components) + len(experiment)
    
    if total_tip_uses > 96:
        raise ValueError(f"Total tip uses ({total_tip_uses}) exceeds available tips (96). This includes media transfers and non-matching subculture transfers.")
    
    return total_tip_uses

def validate_components(experiment, strain_map, supplement_map):
    """
    Validate strains, supplements and media components.
    
    Args:
        experiment (list): List of experiment entries
        strain_map (dict): Mapping of strains to wells
        supplement_map (dict): Mapping of supplements to wells
        
    Raises:
        ValueError: If any invalid components are found
    """
    # Use sets to collect unique invalid items
    invalid_strains = set()
    invalid_supplements = set()
    
    # Get all valid supplements and media types once
    valid_supplements = set(supplement_map.keys())
    
    for entry in experiment:
        # Check strain
        strain = entry["strain"]
        if strain not in strain_map and strain not in ["WT", "Blank"]:
            invalid_strains.add(strain)
        
        # Check media type and supplements
        media_type = entry["media_type"]
        if media_type not in valid_supplements:
            invalid_supplements.add(media_type)
        
        # Check supplements
        invalid_supplements.update(sup for sup in entry["supplements"] if sup not in valid_supplements)
        
        # Check water if needed
        if media_type[:2] == "2x" and "Water" not in valid_supplements:
            invalid_supplements.add("Water")
    
    # Raise errors if any invalid items found
    if invalid_strains:
        raise ValueError(f"The following strains were not found in the strain map: {', '.join(sorted(invalid_strains))}")
    
    if invalid_supplements:
        raise ValueError(f"The following supplements and media were not found in the supplement map: {', '.join(sorted(invalid_supplements))}")

def validate_reservoir_volumes(experiment, supplement_map, supplement_volume):
    """
    Validate that reservoirs have sufficient capacity for all components.
    
    Args:
        experiment (list): List of experiment entries
        supplement_map (dict): Mapping of supplements to wells
        supplement_volume (float): Volume of each supplement in µL
        
    Raises:
        ValueError: If any reservoir has insufficient capacity
    """
    print("\nValidating reservoir volumes...")
    
    def calculate_required_volumes(experiment):
        """
        Calculate the total volume required for each media component and supplement.
        
        Args:
            experiment (list): List of experiment entries
            
        Returns:
            dict: Dictionary mapping component names to required volumes
        """
        required_volumes = {}
        
        for entry in experiment:
            # Add media volume
            media_type = entry["media_type"]
            if media_type not in required_volumes:
                required_volumes[media_type] = 0
            required_volumes[media_type] += entry["media_volume"]
            
            # Add water volume
            if entry["water_volume"] > 0:
                if "Water" not in required_volumes:
                    required_volumes["Water"] = 0
                required_volumes["Water"] += entry["water_volume"]
            
            # Add supplement volumes
            for sup in entry["supplements"]:
                if sup not in required_volumes:
                    required_volumes[sup] = 0
                required_volumes[sup] += supplement_volume  # Use the calculated supplement volume
        
        print("\nRequired volumes:")
        for component, volume in required_volumes.items():
            print(f"  {component}: {volume}µL")
        
        return required_volumes
    
    def check_well_capacity(wells, required_volume, well_capacity=20000):
        """
        Check if a set of wells has sufficient capacity for a required volume.
        
        Args:
            wells (list): List of wells in the reservoir
            required_volume (float): Volume required in µL
            well_capacity (float): Capacity of each well in µL (default 20mL)
            
        Returns:
            tuple: (bool, int) - (whether sufficient capacity, number of wells needed)
        """
        total_capacity = len(wells) * well_capacity
        wells_needed = (required_volume + well_capacity - 1) // well_capacity  # Ceiling division
        print(f"  Checking capacity:")
        print(f"    Required volume: {required_volume}µL")
        print(f"    Total capacity: {total_capacity}µL ({len(wells)} wells)")
        print(f"    Wells needed: {wells_needed}")
        
        # Warn if multiple wells are needed
        if wells_needed > 1:
            print(f"    [WARNING] {wells_needed} wells needed for this component")
        
        return total_capacity >= required_volume, wells_needed
    
    # Calculate required volumes
    required_volumes = calculate_required_volumes(experiment)
    
    # Only check volumes for components that are actually used
    insufficient_volumes = []
    multiple_wells_needed = []
    
    print("\nChecking reservoir capacities:")
    for component, volume in required_volumes.items():
        print(f"\nChecking {component}:")
        # Get reservoir info from supplement map
        if component not in supplement_map:
            raise ValueError(f"Component '{component}' not found in supplement map")
        
        # Get reservoir and well capacity from dictionary
        reservoir = supplement_map[component]["plate"]
        well_capacity = RESERVOIR_CAPACITIES[reservoir]
        wells = supplement_map[component]["wells"]
        print(f"  Reservoir: {reservoir}")
        print(f"  Well capacity: {well_capacity}µL")
        
        # Check if we have enough wells for the required volume
        has_capacity, wells_needed = check_well_capacity(wells, volume, well_capacity)
        if not has_capacity:
            insufficient_volumes.append(
                f"{component} in {reservoir}. Required: {volume}µL. "
                f"Current well number: {len(wells)}. "
                f"Required well number: {wells_needed}"
            )
        elif wells_needed > 1:
            multiple_wells_needed.append(
                f"{component} in {reservoir} needs {wells_needed} wells for {volume}µL"
            )
    
    if insufficient_volumes:
        print("\n[ERROR] Insufficient volume in reservoirs:")
        for msg in insufficient_volumes:
            print(f"  {msg}")
        raise ValueError("Insufficient volume in reservoirs:\n" + "\n".join(insufficient_volumes))
    
    if multiple_wells_needed:
        print("\n[WARNING] Multiple wells needed for some components:")
        for msg in multiple_wells_needed:
            print(f"  {msg}")
    
    print("\nAll reservoir volumes validated successfully!")

def validate_experiment(experiment, supplement_map, strain_map, reps=3):
    """
    Master function to validate all aspects of the experiment.
    
    Args:
        experiment (list): List of experiment entries
        supplement_map (dict): Mapping of supplements to wells
        strain_map (dict): Mapping of strains to wells
        reps (int): Number of replicate plates
        
    Raises:
        ValueError: If any validation check fails
    """
    # Validate wells
    wells = [entry["well"] for entry in experiment]
    validate_unique_wells(wells, "wells in experiment")
    
    # Check total number of wells
    if len(experiment) > 96:
        raise ValueError("Experiment cannot exceed 96 wells")
    
    # Validate components
    validate_components(experiment, strain_map, supplement_map)
    
    # Check tip usage
    calculate_total_tip_uses(experiment, strain_map)
    
    # Calculate supplement volume
    supplement_volume = calculate_supplement_volume(reps)
    
    # Validate reservoir volumes
    validate_reservoir_volumes(experiment, supplement_map, supplement_volume)

def validate_media_name(media_name):
    """
    Validate that a media name follows the correct format.
    
    Rules:
    1. Must start with either '1x' or '2x'
    2. Must contain the word 'media'
    """
    if not (media_name.startswith("1x") or media_name.startswith("2x")):
        raise ValueError(f"Invalid media name format: '{media_name}'. Media name must start with either '1x' or '2x'")
    
    return True

def parse_experiment_csv(csv_path, sup_map_path, strain_map_path, reps=3):
    """
    Parse the CSV file containing experiment parameters.
    
    Args:
        csv_path (str): Path to the CSV file
        sup_map_path (str): Path to supplement map CSV file
        strain_map_path (str): Path to strain map CSV file
        reps (int): Number of replicate plates (1-3)
        
    Returns:
        tuple: (experiment list, metadata DataFrame, list of invalid supplements/strains)
        
    Raises:
        ValueError: If CSV format is incorrect, exceeds 96 wells, contains invalid strains/supplements,
                   or if number of supplements exceeds available water volume
    """
    print(f"\nStarting CSV parsing from: {csv_path}")
    
    # Calculate total volume needed based on number of replicates
    TOTAL_VOLUME = (reps * (FINAL_PLATE_VOLUME - SUBCULTURE_TRANSFER) + 
                    (MIX_PLATE_VOLUME - SUBCULTURE_TRANSFER))
    
    # Calculate supplement volume
    VOL_PER_SUPPLEMENT = calculate_supplement_volume(reps)
    print(f"\nCalculated volumes:")
    print(f"Total volume per well: {TOTAL_VOLUME}µL")
    print(f"Volume per supplement: {VOL_PER_SUPPLEMENT}µL")
    
    # Load supplement and strain maps
    supplement_map = load_supplement_map(sup_map_path)
    strain_map = load_strain_map(strain_map_path)
    
    experiment = []
    metadata_rows = []
    invalid_supplements = []
    invalid_strains = []
    
    try:
        with open(csv_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check required columns
            required_columns = ["Well", "Strain", "Media Type", "Supplements"]
            if not all(col in reader.fieldnames for col in required_columns):
                raise ValueError(f"CSV must contain all required columns: {required_columns}")
            
            wells = []
            for row_num, row in enumerate(reader, start=2):  # Start from 2 to account for header row
                # Convert None values to empty strings
                row = {k: (v if v is not None else "") for k, v in row.items()}
                
                # Skip empty rows or rows with only well specified
                if not row["Well"].strip() or (row["Well"].strip() and not any(row[field].strip() for field in ["Strain", "Media Type", "Supplements"])):
                    continue
                
                # Validate well format
                validate_well_format(row["Well"])
                wells.append(row["Well"])
                
                # Validate strain and media type
                strain_name = row["Strain"].strip()
                media_type = row["Media Type"].strip()
                
                if not strain_name or not media_type:
                    raise ValueError(f"Row {row_num}: Missing required fields. For non-blank wells, both Strain and Media Type are required.")
                
                if strain_name not in strain_map:
                    invalid_strains.append(strain_name)
                
                try:
                    validate_media_name(media_type)
                except ValueError as e:
                    raise ValueError(f"Error in well {row['Well']}: {str(e)}")
                
                
                supplements = [s.strip() for s in row["Supplements"].split(';')] if row["Supplements"].strip() else []
                
                # Check for invalid supplements, media, and water
                for sup in supplements:
                    if sup not in supplement_map:
                        invalid_supplements.append(sup)
                
                # Check if media type is in supplement map
                if media_type not in supplement_map:
                    invalid_supplements.append(media_type)
                
                # For 2x media, check if water is in supplement map
                if media_type.startswith("2x_") and "Water" not in supplement_map:
                    invalid_supplements.append("Water")
                
                if media_type.startswith("1x"):
                    # Full composition media - no supplements or water needed
                    media_volume = TOTAL_VOLUME
                    water_volume = 0
                else:
                    # 2x minimal media - need supplements and water
                    media_volume = TOTAL_VOLUME // 2  # Half of total volume
                    supplement_volume = len(supplements) * VOL_PER_SUPPLEMENT
                    water_volume = (TOTAL_VOLUME // 2) - supplement_volume  # Remaining volume after supplements
                    
                    # Validate that supplements don't exceed available water volume
                    max_supplements = (TOTAL_VOLUME // 2) // VOL_PER_SUPPLEMENT
                    if len(supplements) > max_supplements:
                        raise ValueError(
                            f"Too many supplements in well {row['Well']}. "
                            f"Maximum of {max_supplements} supplements allowed out of {TOTAL_VOLUME//2}µL water volume. "
                            f"Found {len(supplements)} supplements."
                        )
                
                experiment.append({
                    "well": row["Well"],
                    "strain": strain_name,
                    "strain_well": strain_map[strain_name],
                    "media_type": media_type,
                    "supplements": supplements,
                    "media_volume": media_volume,
                    "water_volume": water_volume,
                    "supplement_volume": VOL_PER_SUPPLEMENT
                })
                
                # Create metadata row
                metadata_rows.append({
                    "Well": row["Well"],
                    "Strain": strain_name,
                    "Strain Well": strain_map[strain_name],
                    "Media Type": media_type,
                    "Supplements": row["Supplements"],
                    "Media Volume (µL)": media_volume,
                    "Water Volume (µL)": water_volume,
                    "Supplement Volume (µL)": len(supplements) * VOL_PER_SUPPLEMENT,
                    "Volume per Supplement (µL)": VOL_PER_SUPPLEMENT,
                    "Total Volume (µL)": TOTAL_VOLUME
                })
        
        # Run all validation checks
        validate_experiment(experiment, supplement_map, strain_map, reps)
        
        print(f"Successfully processed {len(experiment)} wells")
        
        # Create metadata DataFrame
        metadata_rows_df = pd.DataFrame(metadata_rows)
        
        # Create supplement map DataFrame with correct structure
        sup_map_rows = []
        for sup, info in supplement_map.items():
            sup_map_rows.append({
                "Well": info["wells"][0],
                "Supplement/Media": sup,
                "Plate": info["plate"]
            })
        metadata_sup_map_df = pd.DataFrame(sup_map_rows)
        
        # Create strain map DataFrame
        strain_map_rows = []
        for strain, well in strain_map.items():
            strain_map_rows.append({
                "Well": well,
                "Strain": strain
            })
        metadata_strain_map_df = pd.DataFrame(strain_map_rows)
        
        # Return both DataFrames
        return experiment, [metadata_rows_df, metadata_sup_map_df, metadata_strain_map_df], supplement_map, strain_map, invalid_supplements, invalid_strains
        
    except Exception as e:
        print(f"[ERROR] Failed to parse CSV: {str(e)}")
        raise

def generate_protocol(csv_path, template_path, output_dir, sup_map_path, strain_map_path, reps=3):
    """
    Generate Opentrons protocols from CSV input.
    
    Args:
        csv_path (str): Path to input CSV file
        template_path (str): Path to protocol template
        output_dir (str): Directory to save output files
        sup_map_path (str): Path to supplement map CSV file
        strain_map_path (str): Path to strain map CSV file
        reps (int): Number of replicate plates (1-3)
        
    Returns:
        list: Paths to generated protocol files
    """
    print(f"\nStarting protocol generation...")
    
    # Create timestamped output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    timestamped_dir = os.path.join(output_dir, f"protocol_{timestamp}")
    os.makedirs(timestamped_dir, exist_ok=True)
    
    # Parse the CSV
    experiment_plan, metadata_dfs, supplement_map, strain_map, invalid_supplements, invalid_strains = parse_experiment_csv(csv_path, sup_map_path, strain_map_path, reps)
    print("... Successfully parsed CSV!")
    
    # Get base directory of templates
    template_dir = Path(template_path).parent
    
    # Generate both Flex and OT2 protocols
    generated_paths = []
    for robot_type in ['Flex', 'OT2']:
        # Read the appropriate template
        template_file = f"Media_Bot_{robot_type}_sub_template.py"
        template_path = template_dir / template_file
        
        with open(template_path, 'r') as f:
            template = f.read()
        
        # Convert experiment plan and maps to JSON string
        experiment_json = json.dumps(experiment_plan, indent=4)
        supplement_map_json = json.dumps(supplement_map, indent=4)
        strain_map_json = json.dumps(strain_map, indent=4)
        
        # Replace the experiment plan and map loading in the template
        modified_protocol = template.replace(
            "with open('experiment_plan.json') as f:\n        experiment_plan = json.load(f)\n\n    with open('supplement_map.json') as f:\n        supplement_map = json.load(f)\n        \n    with open('strain_map.json') as f:\n        strain_map = json.load(f)",
            f"experiment_plan = {experiment_json}\n    supplement_map = {supplement_map_json}\n    strain_map = {strain_map_json}"
        )
        
        # Add reps value
        modified_protocol = modified_protocol.replace(
            "reps = 3",
            f"reps = {reps}"
        )
        
        # Generate output filename
        protocol_filename = f"Media_Bot_{robot_type}_protocol_{timestamp}.py"
        
        # Write the modified protocol
        protocol_path = os.path.join(timestamped_dir, protocol_filename)
        with open(protocol_path, 'w') as f:
            f.write(modified_protocol)
        
        generated_paths.append(protocol_path)
        print(f"Generated {robot_type} protocol: {protocol_path}")
    
    # Save metadata with titles
    metadata_filename = f"protocol_metadata_{timestamp}.csv"
    metadata_path = os.path.join(timestamped_dir, metadata_filename)
    with open(metadata_path, 'w', newline='') as f:
        # Write experiment data section
        f.write("=== Experiment Data ===\n")
        metadata_dfs[0].to_csv(f, index=False)
        
        # Write supplement map section
        f.write("\n=== Supplement and Media Map ===\n")
        metadata_dfs[1].to_csv(f, index=False)
        
        # Write strain map section
        f.write("\n=== Strain Map ===\n")
        metadata_dfs[2].to_csv(f, index=False)
    
    print(f"Metadata saved to {metadata_path}")
    
    return generated_paths

def main():
    print("\nStarting main function...")
    parser = argparse.ArgumentParser(
        description="Generate Opentrons protocols from CSV input",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python Media_Bot_parser.py --csv input.csv --template template.py --output protocols/ --sup_map supplement_map.csv --strain_map strain_map.csv --reps 2

The CSV file should contain the following columns:
    - Well: Target well location (e.g., A1)
    - Strain: Name of the strain (must match strain map)
    - Media Type: Either '2x_Media_minimal' or '1x_Media_full_comp'
    - Supplements: Semicolon-separated list of supplements (only used with 2x_Media_minimal)
        """
    )
    
    # Get the directory of the current script
    base_dir = Path(__file__).parent
    
    parser.add_argument(
        "--csv",
        type=str,
        default=str(base_dir / "Media_Bot_experiment_template.csv"),
        help="Path to input CSV file (default: Media_Bot_experiment_template.csv in script directory)"
    )
    
    parser.add_argument(
        "--template",
        type=str,
        default=str(base_dir / "Media_Bot_Flex_sub_template.py"),
        help="Path to protocol template (default: Media_Bot_Flex_sub_template.py in script directory)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        help="Directory to save output files (default: generated_protocols in the same directory as input CSV)"
    )
    
    parser.add_argument(
        "--sup_map",
        type=str,
        default=str(base_dir / "Media_Bot_Supplement_Map_template.csv"),
        help="Path to supplement map CSV file (default: Media_Bot_Supplement_Map_template.csv in script directory)"
    )
    
    parser.add_argument(
        "--strain_map",
        type=str,
        default=str(base_dir / "Media_Bot_strains_map_template.csv"),
        help="Path to strain map CSV file (default: Media_Bot_strains_map_template.csv in script directory)"
    )

    parser.add_argument(
        "--reps",
        type=int,
        default=3,
        help="Number of replicate plates (1-3, default: 3)"
    )
    
    args = parser.parse_args()
    
    # Set default output to be in the same directory as the input CSV if not specified
    if args.output is None:
        csv_path = Path(args.csv)
        args.output = str(csv_path.parent / "generated_protocols")
    
    print(f"\nParsed arguments:...")
    print(f"... CSV: {args.csv}")
    print(f"... Template: {args.template}")
    print(f"... Output: {args.output}")
    print(f"... Supplement Map: {args.sup_map}")
    print(f"... Strain Map: {args.strain_map}")
    print(f"... Replicates: {args.reps}")

    # Validate reps argument
    if not 1 <= args.reps <= 3:
        print("[ERROR] Number of replicates must be between 1 and 3")
        return 1
    
    try:
        # Generate both Flex and OT2 protocols
        generate_protocol(args.csv, args.template, args.output, args.sup_map, args.strain_map, args.reps)
    except Exception as e:
        print(f"[ERROR] Failed to generate protocol: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    print("Media_Bot script started...")
    exit(main())

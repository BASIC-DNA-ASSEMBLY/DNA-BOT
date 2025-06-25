#!/usr/bin/env python3
"""
Custom test script for DNA-BOT that uses multistage_builds/stage_2 as reference
and specific CSV files as inputs.
"""

import os
import sys
import shutil
import datetime
from pathlib import Path
import subprocess

# Add parent directory to path so we can import dnabot_app
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)  # Insert at start to override any installed packages

# Now import the local version of DNA-BOT
import dnabot.dnabot_app as dnabot

# Constants
TEST_RUNS_DIR = Path(__file__).parent / "test_runs"
STAGE2_INPUTS_DIR = Path(__file__).parent.parent / "648_constructs" / "multistage_builds" / "stage_2"

def create_test_directory() -> Path:
    """Creates a new timestamped directory for this test run.
    
    Returns:
        Path: Path to the new test directory
    """
    # Create test_runs directory if it doesn't exist
    TEST_RUNS_DIR.mkdir(exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    test_dir = TEST_RUNS_DIR / f"test_run_stage2_{timestamp}"
    test_dir.mkdir(exist_ok=True)
    return test_dir

def copy_stage2_inputs(test_dir: Path) -> None:
    """Copies stage2 input files to the test directory.
    
    Args:
        test_dir (Path): Path to the test directory
    """
    if not STAGE2_INPUTS_DIR.exists():
        raise FileNotFoundError(f"Stage2 inputs directory not found at {STAGE2_INPUTS_DIR}")
    
    # Copy the specific CSV files
    construct_file = STAGE2_INPUTS_DIR / "stage2_constructs.csv"
    parts_file = STAGE2_INPUTS_DIR / "stage2_parts.csv"
    
    if not construct_file.exists():
        raise FileNotFoundError(f"Construct file not found at {construct_file}")
    if not parts_file.exists():
        raise FileNotFoundError(f"Parts file not found at {parts_file}")
    
    shutil.copy2(construct_file, test_dir / "stage2_constructs.csv")
    shutil.copy2(parts_file, test_dir / "stage2_parts.csv")
    
    print(f"Copied {construct_file} to {test_dir / 'stage2_constructs.csv'}")
    print(f"Copied {parts_file} to {test_dir / 'stage2_parts.csv'}")

def run_dnabot(test_dir: Path) -> None:
    """Runs DNA-BOT on the test directory.
    
    Args:
        test_dir (Path): Path to the test directory containing input files
    """
    # Find construct and source files
    construct_file = test_dir / "stage2_constructs.csv"
    parts_file = test_dir / "stage2_parts.csv"
    
    if not construct_file.exists():
        raise FileNotFoundError(f"Construct file not found: {construct_file}")
    if not parts_file.exists():
        raise FileNotFoundError(f"Parts file not found: {parts_file}")
    
    print(f"\nDebug: Using construct file: {construct_file}")
    print(f"Debug: Using parts file: {parts_file}")
    
    # Change to test directory
    original_dir = os.getcwd()
    print(f"Debug: Original directory: {original_dir}")
    os.chdir(test_dir)
    print(f"Debug: Changed to directory: {os.getcwd()}")
    
    try:
        # Run DNA-BOT
        print("\nDebug: Setting up DNA-BOT arguments...")
        sys.argv = [
            "dnabot_app.py",
            "nogui",
            "--construct_path", str(construct_file),
            "--source_paths", str(parts_file),
            "--etoh_well", "A11",
            "--soc_column", "1"
        ]
        print(f"Debug: DNA-BOT arguments: {sys.argv}")
        
        print("\nDebug: Starting DNA-BOT main()...")
        dnabot.main()
        print("Debug: DNA-BOT main() completed")
    finally:
        # Return to original directory
        os.chdir(original_dir)
        print(f"Debug: Returned to directory: {os.getcwd()}")

def main():
    """Main function to run the stage2 test."""
    try:
        print("Starting DNA-BOT stage2 test...")
        print(f"Using reference folder: {STAGE2_INPUTS_DIR}")
        
        # Create test directory
        test_dir = create_test_directory()
        print(f"Created test directory: {test_dir}")
        
        # Copy stage2 inputs
        copy_stage2_inputs(test_dir)
        print("Copied stage2 inputs")
        
        # Run DNA-BOT
        print("\nRunning DNA-BOT...")
        run_dnabot(test_dir)
        print("DNA-BOT completed successfully")
        
        print(f"\nTest completed successfully!")
        print(f"Generated files are in: {test_dir}")
        
        # List generated files
        print("\nGenerated files:")
        for file in test_dir.glob("*"):
            if file.is_file():
                print(f"  - {file.name}")
        
        sys.exit(0)
            
    except Exception as e:
        print(f"\nError during test: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main() 
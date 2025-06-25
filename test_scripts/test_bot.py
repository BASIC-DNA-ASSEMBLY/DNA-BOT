#!/usr/bin/env python3
"""
Test script for DNA-BOT that:
1. Creates a timestamped test directory
2. Copies example input files
3. Runs DNA-BOT on the inputs
4. Compares generated scripts with previous versions
"""

import os
import sys
import shutil
import datetime
import filecmp
from pathlib import Path
import subprocess
import pandas as pd
import re

# Add parent directory to path so we can import dnabot_app
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)  # Insert at start to override any installed packages

# Now import the local version of DNA-BOT
import dnabot.dnabot_app as dnabot

# Constants
EXAMPLE_INPUTS_DIR = Path(__file__).parent / "example_inputs"
REFERENCE_OUTPUTS_DIR = Path(__file__).parent / "reference_outputs"
TEST_RUNS_DIR = Path(__file__).parent / "test_runs"
STANDARD_TEST_DIR = TEST_RUNS_DIR / "test_run_standard"
SCRIPT_EXTENSIONS = ['.py', '.csv', '.txt']

def create_test_directory() -> Path:
    """Creates a new timestamped directory for this test run.
    
    Returns:
        Path: Path to the new test directory
    """
    # Create test_runs directory if it doesn't exist
    TEST_RUNS_DIR.mkdir(exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    test_dir = TEST_RUNS_DIR / f"test_run_{timestamp}"
    test_dir.mkdir(exist_ok=True)
    return test_dir

def copy_example_inputs(test_dir: Path) -> None:
    """Copies example input files to the test directory.
    
    Args:
        test_dir (Path): Path to the test directory
    """
    if not EXAMPLE_INPUTS_DIR.exists():
        raise FileNotFoundError(f"Example inputs directory not found at {EXAMPLE_INPUTS_DIR}")
    
    for file in EXAMPLE_INPUTS_DIR.glob("*.csv"):
        if file.is_file():
            shutil.copy2(file, test_dir / file.name)

def run_dnabot(test_dir: Path) -> None:
    """Runs DNA-BOT on the test directory.
    
    Args:
        test_dir (Path): Path to the test directory containing input files
    """
    # Find construct and source files
    construct_file = next(test_dir.glob("*constructs*.csv"))
    source_files = list(test_dir.glob("*parts*.csv"))
    
    if not construct_file or not source_files:
        raise FileNotFoundError("Could not find required input files")
    
    print(f"\nDebug: Found construct file: {construct_file}")
    print(f"Debug: Found source files: {source_files}")
    
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
            "--source_paths", *[str(f) for f in source_files],
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

def compare_with_standard(test_dir: Path) -> bool:
    """Compares the test directory contents with the standard test directory.
    
    Args:
        test_dir (Path): Path to the test directory containing generated files
    
    Returns:
        bool: True if contents match, False otherwise
    """
    if not STANDARD_TEST_DIR.exists():
        print(f"Standard test directory not found at {STANDARD_TEST_DIR}")
        return False
    
    # Lists to store comparison results
    mismatched_files = []
    missing_files = []
    new_files = []
    
    def normalize_csv_content(file_path):
        """Normalize CSV content by replacing test run timestamps with a placeholder."""
        if file_path.suffix != '.csv':
            return file_path.read_text()
            
        content = file_path.read_text()
        # Replace test run timestamps with placeholder
        content = re.sub(r'test_run_\d{8}_\d{6}', 'test_run_TIMESTAMP', content)
        return content
    
    def compare_csv_files(standard_file, test_file):
        """Compare CSV files line by line and return differences."""
        standard_lines = normalize_csv_content(standard_file).splitlines()
        test_lines = normalize_csv_content(test_file).splitlines()
        
        differences = []
        for i, (std_line, test_line) in enumerate(zip(standard_lines, test_lines)):
            if std_line != test_line:
                differences.append(f"Line {i+1}:")
                differences.append(f"  Standard: {std_line}")
                differences.append(f"  Test:     {test_line}")
        
        # Check for different number of lines
        if len(standard_lines) != len(test_lines):
            differences.append(f"Files have different number of lines:")
            differences.append(f"  Standard: {len(standard_lines)} lines")
            differences.append(f"  Test:     {len(test_lines)} lines")
            
        return differences
    
    # Compare files
    for standard_file in STANDARD_TEST_DIR.rglob("*"):
        if standard_file.suffix not in SCRIPT_EXTENSIONS:
            continue
            
        rel_path = standard_file.relative_to(STANDARD_TEST_DIR)
        test_file = test_dir / rel_path
        
        if not test_file.exists():
            missing_files.append(str(rel_path))
        else:
            if standard_file.suffix == '.csv':
                differences = compare_csv_files(standard_file, test_file)
                if differences:
                    mismatched_files.append((str(rel_path), differences))
            else:
                # For non-CSV files, use direct comparison
                standard_content = normalize_csv_content(standard_file)
                test_content = normalize_csv_content(test_file)
                if standard_content != test_content:
                    mismatched_files.append((str(rel_path), ["Binary files differ"]))
    
    # Check for extra files in test directory
    for test_file in test_dir.rglob("*"):
        if test_file.suffix not in SCRIPT_EXTENSIONS:
            continue
            
        rel_path = test_file.relative_to(test_dir)
        standard_file = STANDARD_TEST_DIR / rel_path
        
        if not standard_file.exists():
            new_files.append(str(rel_path))
    
    # Print comparison results
    if mismatched_files:
        print("\nFiles that differ from standard:")
        for file_info in mismatched_files:
            file_path, differences = file_info
            print(f"\n{file_path}:")
            for diff in differences:
                print(f"  {diff}")
    
    if missing_files:
        print("\nFiles missing from test directory:")
        for file in missing_files:
            print(f"  - {file}")
    
    if new_files:
        print("\nNew files in test directory:")
        for file in new_files:
            print(f"  - {file}")
    
    # Return True only if all lists are empty
    return not (mismatched_files or missing_files or new_files)

def main():
    """Main function to run the test suite."""
    try:
        print("Starting DNA-BOT test...")
        
        # Create test directory
        test_dir = create_test_directory()
        print(f"Created test directory: {test_dir}")
        
        # Copy example inputs
        copy_example_inputs(test_dir)
        print("Copied example inputs")
        
        # Run DNA-BOT
        print("\nRunning DNA-BOT...")
        run_dnabot(test_dir)
        print("DNA-BOT completed successfully")
        
        # Compare with standard test directory
        print("\nComparing with standard test directory...")
        matches_standard = compare_with_standard(test_dir)
        print(f"Matches standard test directory: {matches_standard}")
        
        # Return success/failure
        sys.exit(0 if matches_standard else 1)
            
    except Exception as e:
        print(f"\nError during test: {str(e)}")
        print(f"Error type: {type(e)}")
        import traceback
        print("\nFull traceback:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main() 
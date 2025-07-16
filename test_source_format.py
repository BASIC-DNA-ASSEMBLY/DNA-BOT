#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dnabot'))

from dnabot.dnabot_app import generate_sources_dict

def test_new_source_format():
    """Test the new source file format with explicit deck positions."""
    try:
        # Test with the new format file
        sources = generate_sources_dict(['examples/part_linker_csvs/source_example_new_format.csv'])
        print(f"✓ Successfully loaded {len(sources)} components")
        
        # Show a sample component
        sample_component = list(sources.items())[0]
        print(f"Sample component: {sample_component[0]} -> {sample_component[1]}")
        
        # Verify the data structure
        well, concentration, deck_position = sample_component[1][:3]
        print(f"  Well: {well}")
        print(f"  Concentration: '{concentration}'")
        print(f"  Deck position: {deck_position}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

if __name__ == "__main__":
    success = test_new_source_format()
    sys.exit(0 if success else 1) 
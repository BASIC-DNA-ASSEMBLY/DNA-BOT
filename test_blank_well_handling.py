#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dnabot'))

from dnabot.dnabot_app import generate_sources_dict

def test_blank_well_handling():
    """Test that blank wells are properly handled in the new source format."""
    try:
        # Test with the file that has a blank well
        sources = generate_sources_dict(['test_scripts/clip_optimisation_test/example_parts.csv'])
        print(f"✓ Successfully loaded {len(sources)} components")
        
        # Check that the blank well (G1) was not included
        if 'G1' in [source_data[0] for source_data in sources.values()]:
            print("✗ ERROR: Blank well G1 was included when it should have been skipped")
            return False
        else:
            print("✓ Correctly skipped blank well G1")
        
        # Show a sample component
        sample_component = list(sources.items())[0]
        print(f"Sample component: {sample_component[0]} -> {sample_component[1]}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error: {e}")
        return False

if __name__ == "__main__":
    success = test_blank_well_handling()
    sys.exit(0 if success else 1) 
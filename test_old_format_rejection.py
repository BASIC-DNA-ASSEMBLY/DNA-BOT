#!/usr/bin/env python3

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dnabot'))

from dnabot.dnabot_app import generate_sources_dict

def test_old_format_rejection():
    """Test that the old source file format is properly rejected."""
    try:
        # Test with the old format file
        sources = generate_sources_dict(['examples/part_linker_csvs/BIOLEGIO_BASIC_STD_SET.csv'])
        print(f"✗ ERROR: Old format was accepted when it should have been rejected!")
        return False
        
    except ValueError as e:
        if "Deck position" in str(e) and "Component name" in str(e):
            print(f"✓ Correctly rejected old format: {e}")
            return True
        else:
            print(f"✗ Unexpected error: {e}")
            return False
    except Exception as e:
        print(f"✗ Unexpected error type: {type(e)} - {e}")
        return False

if __name__ == "__main__":
    success = test_old_format_rejection()
    sys.exit(0 if success else 1) 
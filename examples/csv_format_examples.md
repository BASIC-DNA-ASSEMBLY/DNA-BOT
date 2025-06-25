# DNA-BOT CSV Format Examples

This document provides examples of the expected CSV formats for DNA-BOT input files.

## Constructs CSV Format

The constructs CSV file defines the DNA constructs to be assembled. Each row represents one construct.

### Expected Format
```csv
Construct_ID,Component1,Component2,Component3,Component4,Component5,...
```

### Example
```csv
Construct_ID,Linker1,Part1,Linker2,Part2,Linker3
Construct1,UTR1-P,GFP,UTR1-S,,
Construct2,UTR2-P,LacZ,UTR2-S,,
Construct3,UTR3-P,mCherry,UTR3-S,,
```

### Notes
- **First column**: Construct identifier (can be any name)
- **Subsequent columns**: Components in order (linker, part, linker, part, ...)
- **Empty cells**: Use empty cells for unused component slots
- **Headers required**: The first row must contain column headers
- **Flexible column count**: You can have as many component columns as needed

## Source Parts CSV Format

The source parts CSV file defines the location of parts and linkers in source plates.

### Expected Format
```csv
Part_Name,Well,Concentration,Additional_Info
```

### Example
```csv
Part_Name,Well,Concentration,Notes
UTR1-P,A1,10.5,Universal terminator region 1 prefix
UTR1-S,A2,10.5,Universal terminator region 1 suffix
GFP,B1,5.0,Green fluorescent protein
LacZ,B2,5.0,Beta-galactosidase
mCherry,B3,5.0,Red fluorescent protein
UTR2-P,C1,10.5,Universal terminator region 2 prefix
UTR2-S,C2,10.5,Universal terminator region 2 suffix
```

### Notes
- **Part_Name**: Must match exactly with the names used in constructs CSV
- **Well**: Must be in A1-H12 format (96-well plate)
- **Concentration**: Optional, used for volume calculations
- **Additional_Info**: Optional, for notes or descriptions
- **Headers required**: The first row must contain column headers
- **Multiple source files**: You can use multiple CSV files for different source plates

## Validation Rules

### Well Format Validation
- Must be 2-3 characters long (e.g., "A1", "H12")
- First character must be uppercase letter A-H
- Remaining characters must be digits 1-12
- Examples: "A1", "B5", "H12" ✓
- Examples: "a1", "I1", "A13", "A" ✗

### Unique Well Validation
- Each well can only be used once per source plate
- Different source plates can use the same well coordinates

### Component Availability
- All parts and linkers referenced in constructs must exist in source files
- Case-sensitive matching (e.g., "GFP" ≠ "gfp")

## Error Messages

The improved validation provides clear error messages:

```
Error: Well identifier must be an uppercase letter, got 'a' in 'a1'
Error: Duplicate source plate wells found: A1, B2
Error: The following parts/linkers are required but not found in the source data: GFP, LacZ
Error: CSV file 'constructs.csv' must have at least one component column after the construct ID column
```

## Migration from Old Format

If you have existing CSV files without headers, add a header row:

### Before (Old Format)
```csv
Construct1,UTR1-P,GFP,UTR1-S
Construct2,UTR2-P,LacZ,UTR2-S
```

### After (New Format)
```csv
Construct_ID,Component1,Component2,Component3
Construct1,UTR1-P,GFP,UTR1-S
Construct2,UTR2-P,LacZ,UTR2-S
```

## Benefits of Header-Based Access

1. **Column Order Independence**: You can reorder columns without breaking the code
2. **Additional Columns**: You can add extra columns for notes or metadata
3. **Self-Documenting**: Column names make the data structure clear
4. **Error Prevention**: Missing or malformed headers are caught early
5. **Flexibility**: Different users can use different column names (as long as they're consistent)

## Tips for Creating CSV Files

1. **Use a spreadsheet program**: Excel, Google Sheets, or LibreOffice Calc
2. **Save as CSV**: Make sure to save in CSV format, not Excel format
3. **Check headers**: Ensure the first row contains descriptive column names
4. **Validate wells**: Use the A1-H12 format for all well coordinates
5. **Check spelling**: Ensure part names match exactly between files
6. **Test with small files**: Start with a few constructs to test the format 
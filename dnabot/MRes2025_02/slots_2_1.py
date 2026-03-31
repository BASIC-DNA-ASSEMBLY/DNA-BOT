import ast
import tabulate
import pandas as pd
from pathlib import Path
"""deck is a dictionary that stores the assigned locations of various labware on the Opentrons robot's workspace. It is built dynamically by extracting values from a provided script and ensuring that required slots (such as the trash and thermocycler on the Flex) are correctly assigned.

How deck Works:
deck is initialized as an empty dictionary.
The function get_positions_from_clip(fpath, robot_type) reads and parses a script to extract the assigned positions for:
linker_plates: linker plates
parts_plates: Parts plates
tip_racks: Tip racks
tube_rack: Tube rack
destination_plate: Destination plate
If a clip plate position isn’t explicitly assigned, a default position is selected (7 for OT-2, B3 for Flex), ensuring no conflicts.
For the Flex robot, two additional fixed positions are assigned:
Trash is always in slot "A3"
Thermocycler must be in "B1"
The deck dictionary is then returned, mapping labware names to their assigned deck positions.
Example Output of deck for Flex
python
Copy
Edit
{
    "linker_plates": ["A1", "A2"],
    "parts_plates": ["C1", "C2"],
    "tip_racks": ["D1", "D2"],
    "tube_rack": "B2",
    "destination_plate": "C3",
    "clip_plate": "B3",
    "trash": "A3",
    "thermocycler": "B1"
}
This structure ensures that all labware is positioned correctly for the robot type in use."""

# Preserve white space for markdown output
#   used within to_markdown method from pandas df
tabulate.PRESERVE_WHITESPACE = True
MAXLEN_PLATE_NAME = 25

# Define deck layouts
OT2_DECK_SLOTS = [str(i) for i in range(1, 12)]
FLEX_DECK_SLOTS = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3"]
FLEX_TRASH_SLOT = "A3"  # Fixed trash position
FLEX_THERMO_SLOT = "B1"  # Fixed thermocycler position

def get_positions_from_clip(fpath: Path, robot_type) -> dict:
    """Get labware slots from a clip reaction script"""
    DEFAULT_CLIP_PLATE_SLOT = "7" if robot_type == "OT-2" else "B3"
    deck = {}

    with open(fpath) as ifh:
        code = ast.parse(ifh.read())

    for node in ast.walk(code):
        if isinstance(node, ast.Assign) and "id" in node.targets[0]._fields:
            try:
                name = node.targets[0].id
                literal_value = ast.unparse(node.value)
                value = ast.literal_eval(literal_value)
                
                if name == "clips_dict":
                    deck["linker_plates"] = sorted(list(set(value["prefixes_plates"])))
                    deck["parts_plates"] = sorted(list(set(value["parts_plates"])))
                elif name == "CANDIDATE_TIPRACK_SLOTS":
                    deck["tip_racks"] = value
                elif name == "TUBE_RACK_POSITION":
                    deck["tube_rack"] = value
                elif name == "DESTINATION_PLATE_POSITION":
                    deck["destination_plate"] = value
            except Exception as e:
                ast.dump(node)
                raise e

    # Ensure clip plate is assigned
    if "clip_plate" not in deck:
        used_slots = set(deck.get(slot, []) for slot in deck)
        if DEFAULT_CLIP_PLATE_SLOT in used_slots:
            raise AssertionError(f"Slot {DEFAULT_CLIP_PLATE_SLOT} already used for clip plate.")
        deck["clip_plate"] = DEFAULT_CLIP_PLATE_SLOT

    # Assign fixed slots for Flex
    if robot_type == "Flex":
        deck["trash"] = FLEX_TRASH_SLOT
        deck["thermocycler"] = FLEX_THERMO_SLOT

    return deck

# Similar modifications would be applied to other functions

def format_deck_info(deck: dict, robot_type, section="Deck info") -> str:
    """Format deck info"""
    data = {"Plate": [], "Positions": []}
    for plate, position in deck.items():
        if isinstance(position, list):
            data["Plate"].append(plate.replace("_", " "))
            data["Positions"].append(",".join(position))
        elif isinstance(position, str):
            data["Plate"].append(plate.replace("_", " "))
            data["Positions"].append(position)
        else:
            raise NotImplementedError()
    df = pd.DataFrame(data=data)
    df = df.sort_values(by="Positions")
    plate_table = df.to_markdown(tablefmt="grid", index=False)

    # Deck representation
    deck_slots = OT2_DECK_SLOTS if robot_type == "OT-2" else FLEX_DECK_SLOTS
    data = ["-" for _ in range(len(deck_slots))] + ["bin"]
    for plate, position in deck.items():
        if isinstance(position, str):
            data[deck_slots.index(position)] = plate
        elif isinstance(position, list):
            for pos in position:
                data[deck_slots.index(pos)] = plate
    
    for i in range(len(data)):
        data[i] = data[i].replace("_", " ")
        data[i] = f"{deck_slots[i]:2} | {data[i]:^{MAXLEN_PLATE_NAME}}"
    df = pd.DataFrame([data[i:i+3] for i in range(0, len(data), 3)])
    df = df.iloc[::-1]
    deck_table = df.to_markdown(index=False, headers="", tablefmt="grid", stralign="center")

    sout = f"""## {section}

### Plate table

{plate_table}

### Deck representation

{deck_table}

"""
    return sout

import ast
import tabulate
import pandas as pd
from pathlib import Path

# Preserve white space for markdown output
#   used within to_markdown method from pandas df
tabulate.PRESERVE_WHITESPACE = True
MAXLEN_PLATE_NAME = 25

def get_positions_from_clip(fpath: Path) -> dict:
    """Get labware slots from a clip reaction script

    Parameters
    ----------
    fpath : Path
        script file to be parsed

    Returns
    -------
    dict
        deck positions
        {
            'biolegio_plate': list,
            'parts_plates': list,
            'tip_racks': list,
            'tube_rack': str,
            'destination_plate': str
        }

    Raises
    ------
    e
        parsing error
    """
    DEFAULT_CLIP_PLATE_SLOT = "7"
    deck = {}

    with open(fpath) as ifh:
        code = ast.parse(ifh.read())

    for node in ast.walk(code):
        if isinstance(node, ast.Assign) and "id" in node.targets[0]._fields:
            try:
                name = node.targets[0].id
                if name == "clips_dict":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["biolegio_plate"] = sorted(list(set(value["prefixes_plates"])))
                    deck["parts_plates"] = sorted(list(set(value["parts_plates"])))
                elif name == "CANDIDATE_TIPRACK_SLOTS":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tip_racks"] = value
                elif name == "TUBE_RACK_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tube_rack"] = value
                elif name == "DESTINATION_PLATE_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["destination_plate"] = value
                else:
                    pass
            except Exception as e:
                ast.dump(node)
                raise e
    # Destination plate is only explicitely defined in the 
    #   no thermocycler script, but not in the with thermo
    #   script. The thermocycler is always at position
    #   position 7
    if "clip_plate" not in deck:
        used_slots = []
        for item in deck:
            if isinstance(item, list):
                used_slots += item
            if isinstance(item, str):
                used_slots.append(item)
        if DEFAULT_CLIP_PLATE_SLOT in used_slots:
            raise AssertionError(f"Slot {DEFAULT_CLIP_PLATE_SLOT} already used for clip plate.")
        deck["clip_plate"] = DEFAULT_CLIP_PLATE_SLOT

    return deck


def get_positions_from_purif(fpath: Path) -> dict:
    """Get labware positions from a purification / magbead script

    Parameters
    ----------
    fpath : Path
        script file to be parsed

    Returns
    -------
    dict
        deck positions
        {
            'tip_racks': list,
            'magdeck_plate': str,
            'mix_plate': str,
            'reagant_container_plate': str,
            'bead_container_plate': str
        }

    Raises
    ------
    e
        parsing error
    """
    deck = {}
    with open(fpath) as ifh:
        code = ast.parse(ifh.read())
    for node in ast.walk(code):
        if isinstance(node, ast.Assign) and "id" in node.targets[0]._fields:
            try:
                name = node.targets[0].id
                if name == "CANDIDATE_TIPRACK_SLOTS":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tip_racks"] = value
                if name == "MAGDECK_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["clip_plate"] = value
                elif name == "MIX_PLATE_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["mix_plate"] = value
                elif name == "REAGENT_CONTAINER_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["reagant_plate"] = value
                elif name == "BEAD_CONTAINER_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["bead_plate"] = value
                else:
                    pass
            except Exception as e:
                ast.dump(node)
                raise e
    return deck


def get_positions_from_assembly(fpath: Path) -> dict:
    """Get labware positions from an assembly script

    Parameters
    ----------
    fpath : Path
        script file to be parsed

    Returns
    -------
    dict
        deck positions
        {
            'purified_sample_plate': str,
            'tip_racks': list,
            'destination_plate': str,
            'tube_rack': str
        }

    Raises
    ------
    e
        parsing error
    """
    DEFAULT_DESTINATION_PLATE_SLOT = "7"
    deck = {}
    with open(fpath) as ifh:
        code = ast.parse(ifh.read())
    for node in ast.walk(code):
        if isinstance(node, ast.Assign) and "id" in node.targets[0]._fields:
            try:
                name = node.targets[0].id
                if name == "CANDIDATE_TIPRACK_SLOTS":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tip_racks"] = value
                if name == "MAG_PLATE_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["purified_sample_plate"] = value
                elif name == "TUBE_RACK_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tube_rack"] = value
                elif name == "REAGENT_CONTAINER_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["reagant_plate"] = value
                elif name == "BEAD_CONTAINER_POSITION":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["bead_plate"] = value
                elif name == "TEMPDECK_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["assembly_plate"] = value
                else:
                    pass
            except Exception as e:
                ast.dump(node)
                raise e

    # Destination plate is only explicitely defined in the 
    #   no thermocycler script, but not in the with thermo
    #   script. The thermocycler is always at position
    #   position 7
    if "destination_plate" not in deck:
        used_slots = []
        for item in deck:
            if isinstance(item, list):
                used_slots += item
            if isinstance(item, str):
                used_slots.append(item)
        if DEFAULT_DESTINATION_PLATE_SLOT in used_slots:
            raise AssertionError(f"Destination slot {DEFAULT_DESTINATION_PLATE_SLOT} already used.")
        deck["destination_plate"] = DEFAULT_DESTINATION_PLATE_SLOT

    return deck


def get_positions_from_transfo(fpath: Path) -> dict:
    """Get labware positions from aa transformation script

    Parameters
    ----------
    fpath : Path
        script file to be parsed

    Returns
    -------
    dict
        deck positions
        {
            "agar_plate": str,
            "tip_racks_p20": list,
            "tip_racks_p300": list,
            "tube_rack": str,
            "soc_plate": str,
            "assembly_plate": str,
            "transformation_plate": str
        }

    Raises
    ------
    e
        parsing error
    """
    DEFAULT_TRANSFORMATION_PLATE_SLOT = "7"
    deck = {}
    with open(fpath) as ifh:
        code = ast.parse(ifh.read())
    for node in ast.walk(code):
        if isinstance(node, ast.Assign) and "id" in node.targets[0]._fields:
            try:
                name = node.targets[0].id
                if name == "CANDIDATE_p20_SLOTS":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tip_racks_p20"] = value
                if name == "CANDIDATE_P300_SLOTS":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tip_racks_p300"] = value
                if name == "TUBE_RACK_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["tube_rack"] = value
                if name == "ASSEMBLY_PLATE_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["assembly_plate"] = value
                elif name == "SOC_PLATE_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["SOC_plate"] = value
                elif name == "AGAR_PLATE_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["agar_plate"] = value
                elif name == "TEMPDECK_SLOT":
                    literal_value = ast.unparse(node.value)
                    value = ast.literal_eval(literal_value)
                    deck["transformation_plate"] = value
                else:
                    pass
            except Exception as e:
                ast.dump(node)
                raise e

    # Transformation plate is only explicitely defined in the 
    #   no thermocycler script, but not in the with thermo
    #   script. The thermocycler is always at position
    #   position 7
    if "transformation_plate" not in deck:
        used_slots = []
        for item in deck:
            if isinstance(item, list):
                used_slots += item
            if isinstance(item, str):
                used_slots.append(item)
        if DEFAULT_TRANSFORMATION_PLATE_SLOT in used_slots:
            raise AssertionError(f"Transformation slot {DEFAULT_TRANSFORMATION_PLATE_SLOT} already used.")
        deck["transformation_plate"] = DEFAULT_TRANSFORMATION_PLATE_SLOT

    return deck


def format_deck_info(deck: dict, section="Deck info") -> str:
    """Format deck info

    Parameters
    ----------
    deck : dict

    Returns
    -------
    str
        Formated string
    """
    # Table of position per plate type
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
    data = ["-" for i in range(11)] + ["bin"]
    for plate, position in deck.items():
        if isinstance(position, str):
            data[int(position)-1] = plate
        elif isinstance(position, list):
            for pos in position:
                data[int(pos)-1] = plate
        else:
            raise NotImplementedError
    for i in range(len(data)):  # Prettify
        data[i] = data[i].replace("_", " ")
        data[i] = f"{i+1:2d} | {data[i]:^{MAXLEN_PLATE_NAME}}"
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
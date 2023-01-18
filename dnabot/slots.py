import ast
import pandas as pd
from pathlib import Path


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
    slots = {"Plate": [], "Positions": []}
    for key, value in deck.items():
        if isinstance(value, list):
            slots["Plate"].append(key.replace("_", " "))
            slots["Positions"].append(",".join(value))
        elif isinstance(value, str):
            slots["Plate"].append(key.replace("_", " "))
            slots["Positions"].append(value)
        else:
            raise NotImplementedError()
    df = pd.DataFrame(data=slots)
    df = df.sort_values(by="Positions")
    deck_table = df.to_markdown(tablefmt="grid", index=False)
    s =f"""# {section}

{deck_table}

"""
    return s
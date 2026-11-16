"""Read literal deck assignments from this platform's generated protocols."""
import ast
from pathlib import Path
import pandas as pd

DECK_SLOTS = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3"]


def get_positions(fpath):
    tree = ast.parse(Path(fpath).read_text(encoding="utf-8"))
    deck = {"thermocycler": "B1", "trash": "A3"}
    if Path(fpath).name.startswith("2_"):
        deck.update(
            magnetic_block="D1",
            candidate_small_tipracks=["D3", "C3", "B3"],
            candidate_large_tipracks=["C2"],
        )
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not isinstance(
            node.targets[0], ast.Name
        ):
            continue
        name = node.targets[0].id
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, TypeError):
            continue
        if name == "clips_dict":
            for kind in ("prefixes", "suffixes", "parts"):
                deck[kind + "_plates"] = sorted(set(value[kind + "_plates"]))
        elif name.endswith(("_POSITION", "_SLOT", "_SLOTS")):
            deck[name.lower()] = value
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in ("load_labware", "load_module")
            and len(node.args) > 1
        ):
            location = node.args[1]
            if (
                isinstance(location, ast.Constant)
                and isinstance(location.value, str)
                and location.value in DECK_SLOTS
            ):
                deck["loaded_at_" + location.value] = location.value
    return deck


def get_positions_from_clip(fpath):
    return get_positions(fpath)


def get_positions_from_purif(fpath):
    return get_positions(fpath)


def get_positions_from_assembly(fpath):
    return get_positions(fpath)


def get_positions_from_transformation(fpath):
    return get_positions(fpath)


def format_deck_info(deck, section="Deck info"):
    rows = []
    for name, positions in deck.items():
        positions = positions if isinstance(positions, list) else [positions]
        if any(str(slot) not in DECK_SLOTS for slot in positions):
            raise ValueError("Invalid deck position: " + repr(positions))
        rows.append((name.replace("_", " "), ", ".join(map(str, positions))))
    return (
        "## "
        + section
        + "\n\n"
        + pd.DataFrame(
            rows, columns=["Labware / candidate positions", "Slots"]
        ).to_markdown(index=False)
        + "\n"
    )

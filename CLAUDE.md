# DNA-BOT

DNA-BOT generates Opentrons protocols for BASIC DNA assembly: clip reactions, purification,
assembly and transformation. A person always reviews the generated protocols, analyses them
in the Opentrons App and runs them on the robot. Claude never connects to a robot.

## Layout

| Folder | What it is |
|---|---|
| `DNABOT-OT2/` | **Current OT-2 app** (v2.3): `dnabot_app_2_3_OT2.py` (CLI + GUI launcher), `dnabot_gui_2_3_OT2.py`, `slots_2_3_OT2.py`, `mplates_OT2.py`, four protocol templates in `template_opentrons_scripts/`, custom labware in `labware/`, `validate.py` |
| `DNABOT-Flex/` | **Current Flex app** (v2.3): same structure with `_FLEX` names. Independent of the OT-2 app – change one without assuming the other matches |
| `dnabot/` | Legacy shared app (older versions). Kept for reference; don't develop here |
| `tests/` | `test_platform_apps.py`, `test_platform_gui.py` cover the v2.3 apps. `test_main.py` tests the legacy `dnabot/` app and currently fails because its templates were removed – not a regression |
| `examples/`, `docs/` | Example protocols and the older user instructions |

Each app's `README.md` describes its hardware, deck layout and workflow changes in detail – read it
before changing that app.

## Python environments – two separate ones

| App | Environment | Python | opentrons |
|---|---|---|---|
| OT-2 (and repo root) | `.venv` at the repo root | 3.12 | **9.0.0** |
| Flex | `DNABOT-Flex/.venv` | 3.12 | **10.1.0** |

- **Never upgrade the OT-2 environment past `opentrons==9.0.0`.** From 9.1 the PyPI package is
  Flex-only and rejects OT-2 protocols ("This protocol is designed for an OT-2 robot"). 9.0.0
  supports Protocol API up to 2.28, the same as OT-2 robot software 26.6.0.
- The Flex environment follows the Flex robots: update `DNABOT-Flex/requirements.txt` when they
  are updated (and later Python too). It is separate so this never affects the OT-2.
- `opentrons` is part of each app's `requirements.txt` on purpose: a new user installs that list
  and is ready to generate **and** simulate.
- Both `.venv` folders are git-ignored and were created from Python 3.12. Don't install packages
  into Anaconda or any global Python, and don't create other environments without asking.
- The bundled Python 3.10 (`.python31011`) was removed in Nov 2026; don't reintroduce it.

## Common commands

OT-2 (from `DNABOT-OT2/`, using the root `.venv`):

```
..\.venv\Scripts\python.exe dnabot_app_2_3_OT2.py
..\.venv\Scripts\python.exe validate.py --generate-only
..\.venv\Scripts\python.exe validate.py
```

Flex (from `DNABOT-Flex/`, using its own `.venv`):

```
.venv\Scripts\python.exe dnabot_app_2_3_FLEX.py
.venv\Scripts\python.exe validate.py
```

Tests (from the repo root): `.venv\Scripts\python.exe -m pytest tests`

`validate.py` generates the bundled example and simulates all four protocols. Run it after any
change to an app, its templates, settings or labware, and report the result. A pass means the code
is valid – not that deck geometry, liquid handling or agar depth are right on the real robot.

## Protocol templates

- All eight templates (four OT-2, four Flex) currently declare `apiLevel` **2.21**. Lab robots
  support 2.26, but don't raise the level unless asked; if you do, change all four templates of
  that app together, keep metadata and filenames consistent, and re-run `validate.py`.
- Custom labware definitions live in each app's `labware/` folder. Users must also add them to the
  Opentrons App before analysing protocols.

## Working here

- Explain changes in lab terms (volumes, wells, deck slots, pauses), not just code terms.
- This repo is shared with others outside the lab (BASIC-DNA-ASSEMBLY on GitHub). Keep changes
  focused; don't commit or push unless asked.

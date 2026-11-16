"""Generate bundled examples in isolation and optionally simulate all four."""
import argparse
import ast
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    app = next(root.glob("dnabot_app_2_3_*.py"))
    with tempfile.TemporaryDirectory() as output:
        subprocess.run(
            [
                sys.executable,
                str(app),
                "nogui",
                "--construct_path",
                str(root / "examples/constructs_temp.csv"),
                "--source_paths",
                str(root / "examples/BIOLEGIO_BASIC_RBS_EXT_SET.csv"),
                str(root / "examples/part_plate_2_230419.csv"),
                "--output_dir",
                output,
            ],
            check=True,
            cwd=output,
        )
        paths = sorted(Path(output).glob("*.py"))
        if len(paths) != 4:
            raise RuntimeError("Expected four generated protocols.")
        for path in paths:
            ast.parse(path.read_text(encoding="utf-8"))
        print("PASS: generated and parsed all four protocols.")
        if args.generate_only:
            return
        from opentrons import simulate

        for path in paths:
            with path.open(encoding="utf-8") as protocol, contextlib.redirect_stdout(
                io.StringIO()
            ):
                commands, _ = simulate.simulate(
                    protocol, custom_labware_paths=[str(root / "labware")]
                )
            print(f"PASS: {path.name} ({len(commands)} commands)")


if __name__ == "__main__":
    main()

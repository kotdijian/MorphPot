"""Thin CLI entry point; inherited algorithms remain independently callable."""
from __future__ import annotations
import argparse
import importlib
import sys
from pathlib import Path
from morphpot import __version__
from morphpot.asset_metadata import resolve_units

COMMANDS = {
    "volume": "vessel_voxel_volume",
    "sections": "pottery_radial_sections",
    "surface-trace": "pottery_surface_trace",
    "fragment-boundary": "pottery_fragment_boundary",
}


def main(argv=None):
    parser = argparse.ArgumentParser(description="MorphPot: analysis of normalized pottery PLY meshes")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("command", choices=COMMANDS)
    parser.add_argument("ply", type=Path)
    parser.add_argument("--unit", choices=["auto", "mm", "cm", "m"], default="auto")
    args, remaining = parser.parse_known_args(argv)
    if args.ply.suffix.lower() != ".ply":
        parser.error("The MorphPot entry point accepts normalized PLY meshes.")
    try:
        unit, _, _ = resolve_units(args.ply, args.unit)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    module = importlib.import_module(COMMANDS[args.command])
    old_argv = sys.argv
    try:
        sys.argv = [module.__name__, str(args.ply), "--unit", unit, *remaining]
        module.main()
    finally:
        sys.argv = old_argv
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

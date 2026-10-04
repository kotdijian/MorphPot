"""Read coordinate units; never apply a pose matrix to an already normalized PLY."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

UNIT_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _validate_unit(unit, factor):
    if unit not in UNIT_TO_MM:
        raise ValueError("Coordinate unit must be mm, cm or m.")
    if isinstance(factor, bool) or not isinstance(factor, (int, float)):
        raise ValueError("unit_to_mm must be a number.")
    if not math.isfinite(factor) or not math.isclose(factor, UNIT_TO_MM[unit], rel_tol=1e-12):
        raise ValueError("Coordinate unit and unit_to_mm disagree.")
    return unit, float(factor)


def resolve_units(ply_path, input_unit="auto"):
    """Prefer <stem>.asset.json, then legacy transform.json, then explicit unit.

    Modern metadata is bound to the PLY's bytes, including with an explicit unit.
    Legacy source.sha256 identifies the RAW input, so it is never compared to PLY.
    """
    path = Path(ply_path)
    explicit = input_unit not in (None, "auto")
    if explicit:
        _validate_unit(input_unit, UNIT_TO_MM.get(input_unit))
    sidecar = path.with_suffix(".asset.json")
    if sidecar.exists():
        data = json.loads(sidecar.read_text(encoding="utf-8"))
        if data.get("schema_version") != "1.0":
            raise ValueError("Unsupported asset metadata schema_version.")
        file = data.get("file", {})
        if file.get("geometry_type") != "triangle_mesh":
            raise ValueError("Asset metadata must identify a triangle_mesh.")
        if file.get("sha256") != sha256_file(path):
            raise ValueError("Asset metadata SHA-256 does not match the PLY.")
        coords = data.get("coordinates", {})
        unit, factor = _validate_unit(coords.get("unit"), coords.get("unit_to_mm"))
        if explicit and input_unit != unit:
            raise ValueError("Explicit unit conflicts with asset metadata.")
        return unit, factor, data
    if explicit:
        return input_unit, UNIT_TO_MM[input_unit], None
    legacy = path.parent / "transform.json"
    if legacy.exists():
        data = json.loads(legacy.read_text(encoding="utf-8"))
        source = data.get("source", {})
        if source.get("coordinate_values_rescaled", False) is not False:
            raise ValueError("Legacy rescaled coordinates need an explicit unit or asset metadata.")
        unit, factor = _validate_unit(source.get("input_unit"), source.get("unit_to_mm"))
        return unit, factor, data
    raise ValueError("Unit is unknown. Supply --unit mm|cm|m or matching metadata.")

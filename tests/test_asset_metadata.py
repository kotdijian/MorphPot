import hashlib
import json
import pytest
from morphpot.asset_metadata import resolve_units


def test_legacy_raw_hash_and_matrix_are_not_applied(tmp_path):
    ply = tmp_path / "specimen.ply"
    ply.write_bytes(b"normalized model")
    (tmp_path / "transform.json").write_text(json.dumps({
        "source": {"input_unit": "m", "unit_to_mm": 1000, "sha256": "raw-hash",
                   "coordinate_values_rescaled": False},
        "transform": {"matrix_4x4": [[0]]}
    }))
    assert resolve_units(ply)[:2] == ("m", 1000)
    assert ply.read_bytes() == b"normalized model"


def test_asset_mismatch_and_unit_conflict(tmp_path):
    ply = tmp_path / "specimen.ply"
    ply.write_bytes(b"normalized")
    data = {"schema_version": "1.0", "file": {
        "geometry_type": "triangle_mesh", "sha256": hashlib.sha256(ply.read_bytes()).hexdigest()},
        "coordinates": {"unit": "cm", "unit_to_mm": 10}}
    sidecar = ply.with_suffix(".asset.json")
    sidecar.write_text(json.dumps(data))
    assert resolve_units(ply)[:2] == ("cm", 10)
    with pytest.raises(ValueError, match="conflicts"):
        resolve_units(ply, "mm")
    ply.write_bytes(b"changed")
    with pytest.raises(ValueError, match="SHA-256"):
        resolve_units(ply, "cm")


def test_missing_or_inconsistent_units(tmp_path):
    ply = tmp_path / "specimen.ply"
    ply.write_bytes(b"model")
    with pytest.raises(ValueError, match="unknown"):
        resolve_units(ply)
    assert resolve_units(ply, "mm")[:2] == ("mm", 1)
    (tmp_path / "transform.json").write_text(json.dumps({
        "source": {"input_unit": "m", "unit_to_mm": 1}}))
    with pytest.raises(ValueError, match="disagree"):
        resolve_units(ply)

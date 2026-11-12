import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import trimesh

from morphpot.asset_metadata import sha256_file
from pose_comparison import load_pose, run_pose_comparison


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def exports(tmp_path):
    # Solid wall/bottom with a central open cavity. Coordinates are metres.
    profile = np.array([[0, 0], [24, 0], [24, 42], [20, 42], [20, 4], [0, 4]]) / 1000
    original = trimesh.creation.revolve(profile, sections=48)
    folders = {}
    for mode, angle in (("rim", 0), ("base", 12), ("slice", 6)):
        folder = tmp_path / mode
        folder.mkdir()
        mesh = original.copy()
        matrix = trimesh.transformations.rotation_matrix(math.radians(angle), [1, 0, 0])
        mesh.apply_transform(matrix)
        model = folder / "test_rev.ply"
        mesh.export(model)
        transform = folder / "transform.json"
        transform.write_text(json.dumps({
            "application": "ArtefactsOrthoMaker", "artifact_type": "pottery",
            "coordinate_system": {"z_axis": "up"},
            "source": {"input_unit": "m", "unit_to_mm": 1000,
                       "coordinate_values_rescaled": False, "sha256": "a" * 64},
            "posture": {"orientation_method": mode, "pose_quality": {
                "slice": {"status": "reference", "reason": "fixture"}}},
            "transform": {"transform_direction": "raw_to_normalized", "matrix_4x4": matrix.tolist()},
        }))
        (folder / "item_manifest.json").write_text(json.dumps({
            "schema_version": "1.0.0", "application": {"name": "ArtefactsOrthoMaker"},
            "artifact_type": "pottery", "item_id": "test", "revision_id": mode,
            "input_unit": "m", "inputs": {"source_sha256": "a" * 64},
            "outputs": [
                {"role": role, "path": path.name, "size_bytes": path.stat().st_size,
                 "sha256": sha256_file(path)}
                for role, path in (("normalized_model", model), ("transform", transform))
            ],
        }))
        folders[mode] = folder
    return folders


def run_cli(entry, *args):
    return subprocess.run([sys.executable, str(ROOT / entry), *map(str, args)],
                          cwd=ROOT, capture_output=True, text=True)


def test_compare_cli_and_direct_cli_keep_normalized_units_and_pose(tmp_path, exports):
    target = tmp_path / "comparison"
    proc = run_cli("MorphPot.py", "volume", "--compare-poses", *exports.values(),
                   "--unit", "auto", "--pitch", 2, "--output-dir", target)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    result = json.loads((target / "pose_comparison.json").read_text())
    rim, base, sliced = result["poses"]
    assert result["unit"] == "m"
    assert result["baseline_mode"] == "rim"
    assert rim["status"] == base["status"] == "ok"
    assert sliced["status"] == "skipped"
    assert base["angle_from_baseline_deg"] == pytest.approx(12)
    assert 30 < rim["volume_ml"] < 60
    assert base["volume_ml"] < rim["volume_ml"]
    assert (target / "pose_comparison.csv").is_file()

    direct = tmp_path / "direct"
    proc = run_cli("vessel_voxel_volume.py", exports["base"] / "test_rev.ply",
                   "--pitch", 2, "--output-dir", direct)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    computed = json.loads(next(direct.rglob("result.json")).read_text())
    assert computed["volume_ml"] == pytest.approx(base["volume_ml"])


def test_compare_rejects_inconsistent_inputs_before_output(tmp_path, exports):
    target = tmp_path / "comparison"
    with pytest.raises(ValueError, match="disagrees"):
        run_pose_comparison(list(exports.values()), target, 2, "mm", None, False, None)
    assert not target.exists()

    mesh = exports["base"] / "test_rev.ply"
    mesh.write_bytes(b"X" + mesh.read_bytes()[1:])
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_pose(exports["base"])

"""Compare spill volumes of the same mesh exported in different AOM postures."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

from morphpot.asset_metadata import UNIT_TO_MM, resolve_units, sha256_file


POSE_MODES = ("rim", "base", "slice")


@dataclass(frozen=True)
class PoseInput:
    folder: Path
    mesh: Path
    item_id: str
    revision_id: str
    unit: str
    source_sha256: str
    bundle_sha256: str | None
    mode: str
    up: tuple[float, float, float]
    slice_status: str | None
    slice_reason: str | None
    manual_flip_z: bool


def _check_digest(path: Path, record: dict) -> None:
    if path.stat().st_size != record["size_bytes"]:
        raise ValueError(f"AOM manifest size mismatch: {path}")
    if sha256_file(path) != record["sha256"]:
        raise ValueError(f"AOM manifest SHA-256 mismatch: {path}")


def load_pose(folder: Path) -> PoseInput:
    folder = folder.resolve()
    manifest = json.loads((folder / "item_manifest.json").read_text(encoding="utf-8"))
    if (manifest["schema_version"] != "1.0.0" or manifest["artifact_type"] != "pottery"
            or manifest["application"]["name"] != "ArtefactsOrthoMaker"):
        raise ValueError(f"Expected an AOM pottery export with manifest schema 1.0.0: {folder}")

    outputs = manifest["outputs"]
    models = [r for r in outputs if r["role"] == "normalized_model"]
    transforms = [r for r in outputs if r["role"] == "transform" and r["path"] == "transform.json"]
    if len(models) != 1 or len(transforms) != 1:
        raise ValueError(f"Expected one normalized model and transform.json: {folder}")
    name = models[0]["path"]
    if Path(name).name != name or Path(name).suffix.lower() != ".ply":
        raise ValueError(f"Expected a normalized PLY filename: {name}")
    mesh = folder / name
    transform_path = folder / "transform.json"
    _check_digest(mesh, models[0])
    _check_digest(transform_path, transforms[0])
    metadata = json.loads(transform_path.read_text(encoding="utf-8"))
    source = metadata["source"]
    inputs = manifest["inputs"]
    unit = manifest["input_unit"]
    if (metadata["application"] != "ArtefactsOrthoMaker"
            or metadata["artifact_type"] != "pottery"
            or metadata["coordinate_system"]["z_axis"] != "up"
            or unit not in {"mm", "cm", "m"}
            or unit != source["input_unit"]
            or source["coordinate_values_rescaled"] is not False
            or not inputs["source_sha256"]
            or inputs["source_sha256"] != source["sha256"]):
        raise ValueError(f"Inconsistent AOM source, coordinate unit or axis metadata: {folder}")

    if (source["unit_to_mm"] != UNIT_TO_MM[unit]
            or isinstance(source["unit_to_mm"], bool)):
        raise ValueError(f"Inconsistent AOM unit conversion factor: {folder}")
    resolved_unit, _, _ = resolve_units(mesh)
    if resolved_unit != unit:
        raise ValueError(f"Asset metadata unit disagrees with the AOM manifest: {folder}")

    posture = metadata["posture"]
    mode = posture["orientation_method"]
    if mode not in POSE_MODES:
        raise ValueError(f"Unsupported comparison mode {mode!r}: {folder}")
    transform = metadata["transform"]
    matrix = transform["matrix_4x4"]
    if (transform["transform_direction"] != "raw_to_normalized" or len(matrix) != 4
            or any(len(row) != 4 or not all(math.isfinite(x) for x in row) for row in matrix)):
        raise ValueError(f"Invalid AOM transform: {folder}")
    # The normalized PLY already includes this transform. Read the rotation only
    # to recover normalized +Z in source coordinates for the angle comparison.
    up = tuple(matrix[2][:3])
    if not math.isclose(math.hypot(*up), 1, abs_tol=1e-6):
        raise ValueError(f"Invalid AOM vertical direction: {folder}")

    slice_quality = posture.get("pose_quality", {}).get("slice", {})
    slice_status = slice_quality.get("status")
    if slice_status not in {None, "accepted", "reference", "held"}:
        raise ValueError(f"Unknown Slice assessment {slice_status!r}: {folder}")
    return PoseInput(folder, mesh, manifest["item_id"], manifest["revision_id"],
                     unit, source["sha256"], inputs.get("source_bundle_sha256"),
                     mode, up, slice_status, slice_quality.get("reason"),
                     posture.get("manual_flip_z", False))


def run_pose_comparison(folders, output_dir, pitch, unit, baseline, include_reference,
                        estimate_volume, **settings):
    """Run the existing voxel calculation at one pitch for each eligible export."""
    if len(folders) < 2:
        raise ValueError("Specify at least two AOM export folders.")
    if not math.isfinite(pitch) or pitch <= 0:
        raise ValueError("Comparison pitch must be finite and positive.")
    poses = [load_pose(folder) for folder in folders]
    first = poses[0]
    if len({p.mode for p in poses}) != len(poses):
        raise ValueError("Specify at most one AOM export folder for each posture mode.")
    if any((p.item_id, p.source_sha256, p.unit) != (first.item_id, first.source_sha256, first.unit)
           for p in poses) or len({p.bundle_sha256 for p in poses if p.bundle_sha256}) > 1:
        raise ValueError("Use exports of the same AOM item, source model and coordinate unit.")
    if unit is not None and unit != first.unit:
        raise ValueError(f"--unit disagrees with the AOM export ({first.unit}).")

    rows = {}
    for pose in poses:
        rows[pose.mode] = {
            "mode": pose.mode,
            "revision_id": pose.revision_id,
            "aom_export": str(pose.folder),
            "normalized_model": str(pose.mesh),
            "up_source_xyz": pose.up,
            "manual_flip_z": pose.manual_flip_z,
            "slice_status": pose.slice_status if pose.mode == "slice" else None,
            "reason": pose.slice_reason if pose.mode == "slice" else None,
            "status": "pending",
        }
        if pose.mode == "slice":
            row = rows["slice"]
            if pose.slice_status is None:
                row.update(status="skipped", reason="The AOM export has no Slice assessment.")
            elif pose.slice_status == "held":
                row["status"] = "skipped"
            elif pose.slice_status == "reference" and not include_reference:
                row["status"] = "skipped"
                row["reason"] = "Reference Slice requires --include-reference-slice. " + (row["reason"] or "")

    if "slice" not in rows:
        held = next((p for p in poses if p.slice_status == "held"), None)
        if held is not None:
            rows["slice"] = {"mode": "slice", "status": "skipped", "slice_status": "held",
                             "reason": held.slice_reason, "aom_export": str(held.folder),
                             "revision_id": held.revision_id}

    eligible = {p.mode for p in poses if rows[p.mode]["status"] == "pending"}
    if baseline is not None and (baseline not in eligible or
            (baseline == "slice" and rows["slice"]["slice_status"] != "accepted")):
        raise ValueError("The requested baseline must be an eligible primary comparison mode.")

    # Reserve a fresh output root so previous or concurrent runs cannot mix results.
    root = (output_dir or first.folder.with_name(f"{first.folder.name}_PoseComparison")).resolve()
    root.mkdir(parents=True, exist_ok=False)
    for pose in poses:
        row = rows[pose.mode]
        if pose.mode not in eligible:
            continue
        try:
            result = estimate_volume(input_path=pose.mesh, pitch=pitch, unit=pose.unit,
                                     output_dir=root / pose.mode, **settings)
            row.update(
                status="ok",
                volume_ml=result["volume_ml"],
                spill_level_mm=result["voxel_qa"]["spill_level_z_mm"],
                spill_depth_above_mesh_base_mm=(
                    result["voxel_qa"]["spill_level_z_mm"] - result["mesh_z_min_mm"]
                ),
                result_json=str(Path(result["run_directory"]) / "result.json"),
            )
        except Exception as exc:
            row.update(status="failed", error=str(exc))

    successful = [rows[mode] for mode in POSE_MODES
                  if mode in rows and rows[mode]["status"] == "ok"]
    primary = [row for row in successful if row["slice_status"] != "reference"]
    if baseline is None and primary:
        baseline = primary[0]["mode"]
    if baseline is not None and rows[baseline]["status"] == "ok":
        base = rows[baseline]
        for row in successful:
            dot = sum(a * b for a, b in zip(row["up_source_xyz"], base["up_source_xyz"]))
            row["angle_from_baseline_deg"] = (0.0 if row is base else
                                              math.degrees(math.acos(max(-1.0, min(1.0, dot)))))
            row["volume_delta_ml"] = row["volume_ml"] - base["volume_ml"]
            row["volume_delta_pct"] = (100 * row["volume_delta_ml"] / base["volume_ml"]
                                       if base["volume_ml"] else None)

    ordered = [rows[mode] for mode in POSE_MODES if mode in rows]
    summary = {
        "item_id": first.item_id,
        "source_sha256": first.source_sha256,
        "source_bundle_sha256": next((p.bundle_sha256 for p in poses if p.bundle_sha256), None),
        "unit": first.unit,
        "pitch_mm": pitch,
        "calculation_settings": settings,
        "baseline_mode": baseline,
        "include_reference_slice": include_reference,
        "poses": ordered,
    }
    (root / "pose_comparison.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    fields = ["mode", "status", "slice_status", "reason", "revision_id", "aom_export",
              "normalized_model", "manual_flip_z", "angle_from_baseline_deg",
              "spill_level_mm", "spill_depth_above_mesh_base_mm", "volume_ml",
              "volume_delta_ml", "volume_delta_pct", "result_json", "error"]
    with (root / "pose_comparison.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ordered)
    print(f"Pose comparison: {root}")
    if len(successful) < 2 or not primary or any(row["status"] == "failed" for row in ordered):
        raise RuntimeError("Pose comparison incomplete. See pose_comparison.json for the results and reasons.")

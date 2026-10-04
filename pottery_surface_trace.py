#!/usr/bin/env python3
"""
PotterySurfaceTraceAnalyzer
Version 0.1.0

Standalone CLI for parallel testing with PotteryFragmentBoundaryExtractor.

It extracts a common multiscale geometry feature layer using requested physical
scales in millimetres, then assigns preliminary *morphometric* trace states.
These are not archaeological identifications.

Default scales: 0.5, 1.0, 2.0, 4.0 mm
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
import time
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np
import trimesh

from pottery_multiscale_features import (
    __version__ as feature_layer_version,
    UNIT_SCALE_TO_MM,
    SURFACE_IGNORED,
    SURFACE_OUTER,
    SURFACE_INNER,
    SURFACE_TRANSITION,
    SURFACE_NAMES,
    load_ply_preserve_topology,
    build_edge_topology,
    face_components,
    component_orientation_and_activity,
    estimate_axis_xy,
    classify_surfaces,
    build_same_surface_diffusion,
    compute_multiscale_features_mm,
    robust_unit_interval,
    safe_quantiles,
    write_binary_point_scalar_ply,
)

__version__ = "0.1.0"

# Morphometric state codes. These are geometry descriptions, not archaeological labels.
TRACE_SMOOTH = 0
TRACE_FINE_RIDGE = 1
TRACE_FINE_VALLEY = 2
TRACE_MEDIUM_RIDGE = 3
TRACE_MEDIUM_VALLEY = 4
TRACE_BROAD_RIDGE = 5
TRACE_BROAD_VALLEY = 6
TRACE_STEP_EDGE = 7
TRACE_MULTISCALE_ANOMALY = 8
TRACE_IRREGULAR = 9

TRACE_CLASS_NAMES = {
    TRACE_SMOOTH: "smooth_or_low_response",
    TRACE_FINE_RIDGE: "fine_ridge_like",
    TRACE_FINE_VALLEY: "fine_valley_like",
    TRACE_MEDIUM_RIDGE: "medium_ridge_like",
    TRACE_MEDIUM_VALLEY: "medium_valley_like",
    TRACE_BROAD_RIDGE: "broad_ridge_like",
    TRACE_BROAD_VALLEY: "broad_valley_like",
    TRACE_STEP_EDGE: "step_edge_like",
    TRACE_MULTISCALE_ANOMALY: "multiscale_persistent_anomaly",
    TRACE_IRREGULAR: "irregular_or_mixed_response",
}


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)


def installed_version(name: str):
    try:
        return importlib_metadata.version(name)
    except Exception:
        return None


def parse_scales(text: str) -> tuple[float, ...]:
    try:
        vals = tuple(float(x.strip()) for x in text.split(",") if x.strip())
    except Exception as exc:
        raise argparse.ArgumentTypeError(str(exc))
    if len(vals) < 2 or any(v <= 0 for v in vals):
        raise argparse.ArgumentTypeError("Specify at least two positive scales, e.g. 0.5,1,2,4")
    return vals


def compact_face_subset(vertices, faces, face_ids, vertex_colors=None):
    face_ids = np.asarray(face_ids, dtype=np.int64)
    if len(face_ids) == 0:
        return trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=np.float64),
            faces=np.empty((0, 3), dtype=np.int64),
            process=False,
            validate=False,
        )
    sf = faces[face_ids]
    vids, inverse = np.unique(sf.reshape(-1), return_inverse=True)
    mesh = trimesh.Trimesh(
        vertices=vertices[vids],
        faces=inverse.reshape(-1, 3),
        process=False,
        validate=False,
    )
    if vertex_colors is not None:
        mesh.visual.vertex_colors = vertex_colors[vids]
    return mesh


def export_surface_classes(directory, vertices, faces, vertex_colors, surface_class):
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for code in (SURFACE_OUTER, SURFACE_INNER, SURFACE_TRANSITION, SURFACE_IGNORED):
        ids = np.flatnonzero(surface_class == code)
        p = directory / f"{SURFACE_NAMES[code]}_surface.ply"
        compact_face_subset(vertices, faces, ids, vertex_colors).export(str(p))
        outputs[SURFACE_NAMES[code]] = str(p)
    return outputs


def score_colors(class_code: np.ndarray):
    """Fixed colors for class QA only; scalar fields remain primary."""
    palette = np.array([
        [150, 150, 150, 255],  # smooth
        [255, 180, 0, 255],    # fine ridge
        [0, 150, 255, 255],    # fine valley
        [255, 120, 0, 255],    # med ridge
        [0, 100, 255, 255],    # med valley
        [255, 50, 0, 255],     # broad ridge
        [0, 50, 200, 255],     # broad valley
        [220, 0, 220, 255],    # step
        [255, 255, 0, 255],    # persistent
        [255, 255, 255, 255],  # irregular
    ], dtype=np.uint8)
    return palette[np.clip(class_code.astype(int), 0, len(palette)-1)]


def export_class_mesh(path, vertices, faces, mask, class_code):
    ids = np.flatnonzero(mask)
    mesh = compact_face_subset(vertices, faces, ids)
    if len(ids):
        mesh.visual.face_colors = score_colors(class_code[ids])
    mesh.export(str(path))


def normalize_features(features, surface_class):
    """Normalize each descriptor independently on OUTER and INNER surfaces."""
    scales = features.requested_scales_mm
    norm = {"normal": {}, "residual": {}}
    calibration = {"outer": {}, "inner": {}}
    for side_code, side_name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        mask = surface_class == side_code
        calibration[side_name]["scales"] = {}
        for s in scales:
            n, nq = robust_unit_interval(features.normal_dev_deg[s], mask)
            r, rq = robust_unit_interval(features.abs_residual_mm[s], mask)
            if s not in norm["normal"]:
                norm["normal"][s] = np.zeros(len(surface_class), dtype=np.float32)
                norm["residual"][s] = np.zeros(len(surface_class), dtype=np.float32)
            norm["normal"][s][mask] = n[mask]
            norm["residual"][s][mask] = r[mask]
            calibration[side_name]["scales"][str(s)] = {
                "normal_dev_deg": nq,
                "abs_residual_mm": rq,
            }
    return norm, calibration


def classify_morphometric_states(features, surface_class):
    """
    Rule-based morphology classification.

    The classifier describes scale/sign/persistence of geometric response.
    It deliberately does NOT emit archaeological labels such as 'join' or 'hake'.
    """
    scales = list(features.requested_scales_mm)
    n_faces = len(surface_class)
    valid = (surface_class == SURFACE_OUTER) | (surface_class == SURFACE_INNER)

    norm, calibration = normalize_features(features, surface_class)

    response = np.zeros((len(scales), n_faces), dtype=np.float32)
    for i, s in enumerate(scales):
        response[i] = 0.55 * norm["normal"][s] + 0.45 * norm["residual"][s]

    dominant_idx = np.argmax(response, axis=0)
    dominant_scale = np.asarray([scales[i] for i in dominant_idx], dtype=np.float32)
    dominant_response = response[dominant_idx, np.arange(n_faces)]

    # Persistence: how many physical scales show a substantial normalized response.
    persistence_count = np.sum(response >= 0.55, axis=0).astype(np.float32)
    persistence_fraction = persistence_count / float(len(scales))

    fine_response = response[0]
    broad_response = response[-1]
    if len(scales) >= 3:
        medium_response = np.max(response[1:-1], axis=0)
    else:
        medium_response = response[-1]

    # Ridge/valley sign at dominant scale from signed residual.
    signed_stack = np.vstack([features.signed_residual_mm[s] for s in scales])
    signed_dom = signed_stack[dominant_idx, np.arange(n_faces)]

    # Dihedral is intentionally secondary; robustly normalize per side.
    dih_norm = np.zeros(n_faces, dtype=np.float32)
    dihedral_cal = {}
    for side_code, side_name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        d, q = robust_unit_interval(features.dihedral_max_deg, surface_class == side_code)
        dih_norm[surface_class == side_code] = d[surface_class == side_code]
        dihedral_cal[side_name] = q

    # Continuous morphology scores useful for later task-specific classifiers.
    fine_trace_score = np.clip(fine_response * (1.0 - 0.55 * broad_response), 0.0, 1.0)
    medium_trace_score = np.clip(0.7 * medium_response + 0.3 * np.minimum(medium_response, broad_response), 0.0, 1.0)
    broad_trace_score = np.clip(0.75 * broad_response + 0.25 * dih_norm, 0.0, 1.0)
    persistent_trace_score = np.clip(
        0.65 * np.mean(response, axis=0) + 0.35 * persistence_fraction,
        0.0, 1.0,
    )
    step_edge_score = np.clip(
        0.60 * dih_norm + 0.25 * medium_response + 0.15 * broad_response,
        0.0, 1.0,
    )

    cls = np.full(n_faces, TRACE_SMOOTH, dtype=np.int16)
    strong = valid & (dominant_response >= 0.48)
    persistent = strong & (persistence_count >= max(2, int(math.ceil(len(scales) / 2))))
    step = strong & (step_edge_score >= 0.62) & (dih_norm >= 0.55)

    cls[persistent] = TRACE_MULTISCALE_ANOMALY
    cls[step] = TRACE_STEP_EDGE

    # For non-persistent/non-step responses, classify by dominant physical scale and sign.
    unresolved = strong & ~(persistent | step)
    # Scale bands are defined by index, not numeric midpoint, so arbitrary
    # requested scale sequences retain an explicit fine / intermediate / broad split.
    fine_mask = unresolved & (dominant_idx == 0)
    broad_mask = unresolved & (dominant_idx == (len(scales) - 1))
    med_mask = unresolved & ~(fine_mask | broad_mask)

    ridge = signed_dom >= 0.0
    cls[fine_mask & ridge] = TRACE_FINE_RIDGE
    cls[fine_mask & ~ridge] = TRACE_FINE_VALLEY
    cls[med_mask & ridge] = TRACE_MEDIUM_RIDGE
    cls[med_mask & ~ridge] = TRACE_MEDIUM_VALLEY
    cls[broad_mask & ridge] = TRACE_BROAD_RIDGE
    cls[broad_mask & ~ridge] = TRACE_BROAD_VALLEY

    # Strong but contradictory responses become irregular/mixed.
    contradictory = strong & (
        (fine_response >= 0.75) & (broad_response >= 0.75) & (persistence_count < 2)
    )
    cls[contradictory] = TRACE_IRREGULAR
    cls[~valid] = TRACE_SMOOTH

    fields = {
        "trace_class": cls.astype(np.float32),
        "dominant_scale_mm": dominant_scale,
        "dominant_response": dominant_response.astype(np.float32),
        "persistence_fraction": persistence_fraction.astype(np.float32),
        "fine_trace_score": fine_trace_score.astype(np.float32),
        "medium_trace_score": medium_trace_score.astype(np.float32),
        "broad_trace_score": broad_trace_score.astype(np.float32),
        "persistent_trace_score": persistent_trace_score.astype(np.float32),
        "step_edge_score": step_edge_score.astype(np.float32),
        "dihedral_max_deg": features.dihedral_max_deg.astype(np.float32),
        "dihedral_mean_deg": features.dihedral_mean_deg.astype(np.float32),
    }
    for s in scales:
        label = str(s).replace(".", "p")
        fields[f"normal_dev_{label}mm"] = features.normal_dev_deg[s]
        fields[f"residual_signed_{label}mm"] = features.signed_residual_mm[s]
        fields[f"residual_abs_{label}mm"] = features.abs_residual_mm[s]
        fields[f"response_{label}mm"] = response[scales.index(s)]

    calibration["dihedral_max_deg"] = dihedral_cal
    return cls, fields, calibration


def analyze(
    input_path: Path,
    unit: str,
    scales_mm: tuple[float, ...],
    max_iterations: int,
    output_dir: Path | None,
    surface_normal_threshold: float,
    min_component_area_mm2: float,
):
    started = time.perf_counter()
    input_path = Path(input_path)
    if unit not in UNIT_SCALE_TO_MM:
        raise ValueError("unit must be mm, cm, or m")
    scale_to_mm = UNIT_SCALE_TO_MM[unit]

    base = Path(output_dir) if output_dir else input_path.parent / f"{input_path.stem}_SurfaceTrace_v0"
    qa_dir = base / "qa"
    surf_dir = base / "surfaces"
    feat_dir = base / "features"
    class_dir = base / "classification"
    for p in (base, qa_dir, surf_dir, feat_dir, class_dir):
        p.mkdir(parents=True, exist_ok=True)

    print(f"=== PotterySurfaceTraceAnalyzer v{__version__} ===")
    print(f"input      : {input_path}")
    print(f"unit       : {unit}")
    print(f"scales mm  : {', '.join(f'{x:g}' for x in scales_mm)}")
    print(f"output     : {base}")

    # A. Load/topology
    t = time.perf_counter()
    data = load_ply_preserve_topology(input_path)
    mesh = trimesh.Trimesh(data.vertices, data.faces, process=False, validate=False)
    centers = np.asarray(mesh.triangles_center, dtype=np.float64)
    face_normals = np.asarray(mesh.face_normals, dtype=np.float64)
    edge_topology = build_edge_topology(data.faces)
    ncomp, comp_ids = face_components(len(data.faces), edge_topology)
    orient_sign, active_mask, comp_area, comp_face_count = component_orientation_and_activity(
        data.vertices, data.faces, comp_ids, min_component_area_mm2, scale_to_mm
    )
    topo_time = time.perf_counter() - t

    topo = {
        "vertices": int(len(data.vertices)),
        "faces": int(len(data.faces)),
        "unique_geometry_edges": int(len(edge_topology.edges)),
        "boundary_geometry_edges": int(np.count_nonzero(edge_topology.incidence == 1)),
        "nonmanifold_geometry_edges": int(np.count_nonzero(edge_topology.incidence > 2)),
        "geometry_components": int(ncomp),
        "active_faces": int(np.count_nonzero(active_mask)),
        "stored_vertex_normals_available": bool(data.stored_vertex_normals is not None),
        "reference_normal_source": "triangle winding / Trimesh face normals",
    }
    write_json(qa_dir / "input_topology.json", topo)

    # B. Outer/inner
    t = time.perf_counter()
    axis_xy = estimate_axis_xy(data.vertices, data.faces, active_mask)
    surface_class, radial_score, corrected_normals, face_radius = classify_surfaces(
        centers, face_normals, comp_ids, orient_sign, active_mask, axis_xy,
        surface_normal_threshold,
    )
    surface_time = time.perf_counter() - t
    surface_summary = {
        "axis_xy_native": [float(axis_xy[0]), float(axis_xy[1])],
        "outer_faces": int(np.count_nonzero(surface_class == SURFACE_OUTER)),
        "inner_faces": int(np.count_nonzero(surface_class == SURFACE_INNER)),
        "transition_faces": int(np.count_nonzero(surface_class == SURFACE_TRANSITION)),
        "ignored_faces": int(np.count_nonzero(surface_class == SURFACE_IGNORED)),
    }
    write_json(qa_dir / "surface_classification.json", surface_summary)
    surface_outputs = export_surface_classes(
        surf_dir, data.vertices, data.faces, data.vertex_colors, surface_class
    )

    # C. Physical-scale feature layer
    t = time.perf_counter()
    diffusion, same_a, same_b = build_same_surface_diffusion(
        len(data.faces), edge_topology, surface_class
    )
    features = compute_multiscale_features_mm(
        centers, corrected_normals, same_a, same_b, diffusion, scale_to_mm,
        scales_mm=scales_mm, max_iterations=max_iterations,
    )
    feature_time = time.perf_counter() - t

    resolution_report = {
        "requested_scales_mm": [float(x) for x in features.requested_scales_mm],
        "median_same_surface_face_step_mm": features.median_same_surface_face_step_mm,
        "iterations_by_requested_scale": {str(k): int(v) for k, v in features.iterations_by_scale.items()},
        "realized_diffusion_scale_mm": {str(k): float(v) for k, v in features.realized_scales_mm.items()},
        "clipped_by_max_iterations": {str(k): bool(v) for k, v in features.clipped_scales.items()},
        "resolution_ratio_scale_over_face_step": {
            str(k): float(k / features.median_same_surface_face_step_mm)
            for k in features.requested_scales_mm
        },
        "interpretation": (
            "Physical scale is an approximate graph-diffusion radius calibrated by median same-surface "
            "face-center spacing. Scales near one face step are diagnostic but not expected to support "
            "precise width measurement."
        ),
    }
    write_json(qa_dir / "physical_scale_calibration.json", resolution_report)

    # D. Preliminary morphology states
    t = time.perf_counter()
    trace_class, fields, calibration = classify_morphometric_states(features, surface_class)
    classify_time = time.perf_counter() - t
    write_json(qa_dir / "normalization_calibration.json", calibration)

    class_legend = [
        {"code": int(code), "name": name, "archaeological_identity": "not_assigned"}
        for code, name in TRACE_CLASS_NAMES.items()
    ]
    write_csv(class_dir / "trace_class_legend.csv", class_legend)

    side_outputs = {}
    for side_code, side_name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        ids = np.flatnonzero(surface_class == side_code)
        side_fields = {k: np.asarray(v)[ids] for k, v in fields.items()}
        out = feat_dir / f"{side_name}_multiscale_features_face_centers.ply"
        write_binary_point_scalar_ply(
            out, centers[ids], side_fields, normals=corrected_normals[ids]
        )
        class_mesh = class_dir / f"{side_name}_trace_classes_colored.ply"
        export_class_mesh(class_mesh, data.vertices, data.faces, surface_class == side_code, trace_class)
        side_outputs[side_name] = {
            "feature_scalar_ply": str(out),
            "class_colored_ply": str(class_mesh),
        }

    # Summary statistics
    class_counts = {}
    for side_code, side_name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        mask = surface_class == side_code
        class_counts[side_name] = {
            TRACE_CLASS_NAMES[c]: int(np.count_nonzero(mask & (trace_class == c)))
            for c in TRACE_CLASS_NAMES
        }

    feature_quantiles = {"outer": {}, "inner": {}}
    for side_code, side_name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        mask = surface_class == side_code
        for s in features.requested_scales_mm:
            feature_quantiles[side_name][str(s)] = {
                "normal_dev_deg": safe_quantiles(features.normal_dev_deg[s][mask]),
                "signed_residual_mm": safe_quantiles(features.signed_residual_mm[s][mask]),
                "abs_residual_mm": safe_quantiles(features.abs_residual_mm[s][mask]),
            }

    total_time = time.perf_counter() - started
    result = {
        "program": "PotterySurfaceTraceAnalyzer",
        "program_version": __version__,
        "feature_layer_version": feature_layer_version,
        "status": "ok",
        "input": str(input_path),
        "input_unit": unit,
        "method": "physical-mm-calibrated multiscale geometry descriptors + preliminary morphometric state classification",
        "uv_information_used": False,
        "archaeological_labeling": False,
        "interpretation_warning": (
            "Trace classes describe geometry scale/sign/persistence only. They are not archaeological "
            "identifications such as sherd join, hake-me, crack, or repair."
        ),
        "topology": topo,
        "surface_classification": surface_summary,
        "physical_scale_calibration": resolution_report,
        "trace_class_counts": class_counts,
        "feature_quantiles": feature_quantiles,
        "outputs": {
            "surfaces": surface_outputs,
            "features": side_outputs,
            "trace_class_legend_csv": str(class_dir / "trace_class_legend.csv"),
        },
        "timing_seconds": {
            "topology": float(topo_time),
            "surface_classification": float(surface_time),
            "multiscale_feature_extraction": float(feature_time),
            "morphometric_classification": float(classify_time),
            "total": float(total_time),
        },
    }
    write_json(base / "result.json", result)

    print("\n=== Physical-scale calibration ===")
    print(f"median same-surface face step : {features.median_same_surface_face_step_mm:.4f} mm")
    for s in features.requested_scales_mm:
        print(
            f"requested {s:g} mm -> {features.iterations_by_scale[s]} iterations "
            f"-> realized ~{features.realized_scales_mm[s]:.3f} mm"
            + (" [CLIPPED]" if features.clipped_scales[s] else "")
        )

    print("\n=== Surface classes ===")
    print(f"outer      : {surface_summary['outer_faces']:,}")
    print(f"inner      : {surface_summary['inner_faces']:,}")
    print(f"transition : {surface_summary['transition_faces']:,}")
    print(f"ignored    : {surface_summary['ignored_faces']:,}")

    print("\n=== Output ===")
    print(f"outer features : {side_outputs['outer']['feature_scalar_ply']}")
    print(f"inner features : {side_outputs['inner']['feature_scalar_ply']}")
    print(f"result JSON    : {base / 'result.json'}")
    print(f"total time     : {total_time:.2f} s")
    print("NOTE: trace classes are morphometric descriptions, not archaeological labels.")
    return result


def print_environment():
    print(f"PotterySurfaceTraceAnalyzer {__version__}")
    print(f"Python  : {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    for name in ("numpy", "scipy", "trimesh"):
        print(f"{name:8s}: {installed_version(name) or 'unknown'}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract physical-mm-calibrated multiscale geometry features from pottery OUTER/INNER surfaces "
            "and assign preliminary morphometric trace states. UV seams are not used."
        )
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--diagnose-env", action="store_true")
    parser.add_argument("input", type=Path, nargs="?", help="Input triangular PLY")
    parser.add_argument("--unit", choices=["mm", "cm", "m"], default="mm")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--scales-mm", type=parse_scales, default=(0.5, 1.0, 2.0, 4.0),
        help="Comma-separated physical analysis scales in mm. Default: 0.5,1,2,4"
    )
    parser.add_argument(
        "--max-diffusion-iters", type=int, default=160,
        help="Safety cap for graph diffusion. Default: 160"
    )
    parser.add_argument("--surface-normal-threshold", type=float, default=0.35)
    parser.add_argument("--min-component-area-mm2", type=float, default=10.0)
    args = parser.parse_args()

    if args.diagnose_env:
        print_environment(); return
    if args.input is None:
        parser.error("input PLY is required unless --diagnose-env is used")
    if args.max_diffusion_iters < 1:
        parser.error("--max-diffusion-iters must be positive")
    if not (0 <= args.surface_normal_threshold < 1):
        parser.error("--surface-normal-threshold must be in [0,1)")

    try:
        analyze(
            input_path=args.input,
            unit=args.unit,
            scales_mm=args.scales_mm,
            max_iterations=args.max_diffusion_iters,
            output_dir=args.output_dir,
            surface_normal_threshold=args.surface_normal_threshold,
            min_component_area_mm2=args.min_component_area_mm2,
        )
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()


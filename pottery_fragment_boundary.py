#!/usr/bin/env python3
"""
PotteryFragmentBoundaryExtractor
Version 0.3.0

Experimental geometry-only joining-boundary detector for archaeological pottery
meshes.

v0.3.0 deliberately does NOT use texture UV seams as joining evidence.
It preserves the original PLY geometry topology, separates OUTER and INNER
vessel surfaces, derives multi-scale geometric anomaly scores from triangle
geometry, and then uses the opposite vessel surface only as supporting evidence.

The primary output is a continuous score field, not a categorical assertion
that a face or edge is an archaeological joining boundary.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import struct
import sys
import time
from dataclasses import dataclass
from importlib import metadata as importlib_metadata
from pathlib import Path

try:
    import numpy as np
    import trimesh
    from scipy.sparse import coo_matrix, diags
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
except ModuleNotFoundError as exc:
    missing = exc.name or "unknown"
    print(
        "\nERROR: required Python module is missing.\n"
        f"missing module: {missing}\n"
        "Install at least:\n"
        "  python3 -m pip install numpy scipy trimesh\n",
        file=sys.stderr,
    )
    raise SystemExit(2)

from trimesh.exchange.ply import load_ply

__version__ = "0.3.0"

UNIT_SCALE_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}

SURFACE_IGNORED = -2
SURFACE_TRANSITION = 0
SURFACE_OUTER = 1
SURFACE_INNER = -1


# ----------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------

def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(jsonable(data), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in fieldnames})


def installed_version(name: str):
    try:
        return importlib_metadata.version(name)
    except Exception:
        return None


def safe_quantiles(values: np.ndarray, probs=(0.05, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99)):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return None
    q = np.quantile(values, probs)
    return {f"q{int(round(p * 100)):02d}": float(v) for p, v in zip(probs, q)}


def robust_unit_interval(values: np.ndarray, mask: np.ndarray, qlo=0.50, qhi=0.98):
    """Robustly map one feature to [0,1] independently within one surface class."""
    out = np.zeros(len(values), dtype=np.float32)
    valid = np.asarray(mask, dtype=bool) & np.isfinite(values)
    sample = np.asarray(values[valid], dtype=np.float64)
    if len(sample) == 0:
        return out, {"qlo": None, "qhi": None}
    lo, hi = np.quantile(sample, [qlo, qhi])
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return out, {"qlo": float(lo), "qhi": float(hi)}
    out[valid] = np.clip((values[valid] - lo) / (hi - lo), 0.0, 1.0)
    return out, {"qlo": float(lo), "qhi": float(hi)}


def compact_face_subset(vertices, faces, face_ids, vertex_colors=None):
    face_ids = np.asarray(face_ids, dtype=np.int64)
    if len(face_ids) == 0:
        return trimesh.Trimesh(
            vertices=np.empty((0, 3), dtype=np.float64),
            faces=np.empty((0, 3), dtype=np.int64),
            process=False,
            validate=False,
        )
    sub_faces_global = np.asarray(faces, dtype=np.int64)[face_ids]
    vertex_ids, inverse = np.unique(sub_faces_global.reshape(-1), return_inverse=True)
    sub_vertices = np.asarray(vertices, dtype=np.float64)[vertex_ids]
    sub_faces = inverse.reshape(-1, 3)
    mesh = trimesh.Trimesh(
        vertices=sub_vertices,
        faces=sub_faces,
        process=False,
        validate=False,
    )
    if vertex_colors is not None:
        vc = np.asarray(vertex_colors)[vertex_ids]
        mesh.visual.vertex_colors = vc
    return mesh


# ----------------------------------------------------------------------
# Original-topology PLY loader
# ----------------------------------------------------------------------

def _has_property(data, name: str) -> bool:
    return bool(data.dtype.names and name in data.dtype.names)


def _property_array(data, name: str) -> np.ndarray:
    field = data[name]
    if field.dtype.names and "f1" in field.dtype.names:
        return np.asarray(field["f1"])
    return np.asarray(field)


def _extract_fixed_list_property(structured: np.ndarray, name: str):
    if not _has_property(structured, name):
        raise KeyError(f"PLY property not found: {name}")
    field = structured[name]
    if field.dtype.names and "f0" in field.dtype.names and "f1" in field.dtype.names:
        return np.asarray(field["f0"], dtype=np.int64), np.asarray(field["f1"])
    if field.dtype == object:
        seq = [np.asarray(x) for x in field]
        counts = np.asarray([len(x) for x in seq], dtype=np.int64)
        if len(set(counts.tolist())) != 1:
            raise ValueError(f"Variable-length PLY property '{name}' is unsupported.")
        return counts, np.vstack(seq)
    arr = np.asarray(field)
    if arr.ndim == 2:
        counts = np.full(len(arr), arr.shape[1], dtype=np.int64)
        return counts, arr
    raise TypeError(f"Unsupported raw PLY representation for '{name}': {field.dtype}")


@dataclass
class PlyMeshData:
    vertices: np.ndarray
    faces: np.ndarray
    vertex_colors: np.ndarray | None
    stored_vertex_normals: np.ndarray | None
    header_vertex_count: int
    header_face_count: int


def load_ply_geometry(path: Path) -> PlyMeshData:
    path = Path(path)
    if path.suffix.lower() != ".ply":
        raise ValueError("v0.3.0 currently supports triangular PLY input only.")
    with path.open("rb") as f:
        loaded = load_ply(
            f,
            resolver=None,
            fix_texture=False,
            skip_materials=True,
            prefer_color="vertex",
        )
    metadata = loaded.get("metadata") or {}
    raw = metadata.get("_ply_raw")
    if not raw or "vertex" not in raw or "face" not in raw:
        raise RuntimeError("Raw PLY elements required for topology-preserving analysis are unavailable.")
    v_element = raw["vertex"]
    f_element = raw["face"]
    vdata = v_element["data"]
    fdata = f_element["data"]
    for axis in ("x", "y", "z"):
        if not _has_property(vdata, axis):
            raise ValueError(f"PLY vertex property '{axis}' is missing.")
    vertices = np.column_stack([
        np.asarray(_property_array(vdata, "x"), dtype=np.float64),
        np.asarray(_property_array(vdata, "y"), dtype=np.float64),
        np.asarray(_property_array(vdata, "z"), dtype=np.float64),
    ])
    counts, vals = _extract_fixed_list_property(fdata, "vertex_indices")
    if not np.all(counts == 3):
        observed = sorted(set(int(x) for x in counts))
        raise ValueError(f"v0.3.0 requires triangular source faces; observed {observed}.")
    faces = np.asarray(vals[:, :3], dtype=np.int64)
    if faces.min() < 0 or faces.max() >= len(vertices):
        raise ValueError("Face index is outside the source vertex array.")
    if not np.isfinite(vertices).all():
        raise ValueError("Vertex coordinates contain NaN or Inf.")

    vertex_colors = None
    if all(_has_property(vdata, n) for n in ("red", "green", "blue")):
        vertex_colors = np.column_stack([
            np.asarray(_property_array(vdata, "red"), dtype=np.uint8),
            np.asarray(_property_array(vdata, "green"), dtype=np.uint8),
            np.asarray(_property_array(vdata, "blue"), dtype=np.uint8),
        ])

    stored_normals = None
    if all(_has_property(vdata, n) for n in ("nx", "ny", "nz")):
        stored_normals = np.column_stack([
            np.asarray(_property_array(vdata, "nx"), dtype=np.float64),
            np.asarray(_property_array(vdata, "ny"), dtype=np.float64),
            np.asarray(_property_array(vdata, "nz"), dtype=np.float64),
        ])

    high_faces = np.asarray(loaded.get("faces"), dtype=np.int64)
    if high_faces.shape != faces.shape or not np.array_equal(high_faces, faces):
        raise RuntimeError("High-level face topology differs despite fix_texture=False.")
    if len(np.asarray(loaded.get("vertices"))) != len(vertices):
        raise RuntimeError("Unexpected vertex re-indexing occurred despite fix_texture=False.")

    return PlyMeshData(
        vertices=vertices,
        faces=faces,
        vertex_colors=vertex_colors,
        stored_vertex_normals=stored_normals,
        header_vertex_count=int(v_element["length"]),
        header_face_count=int(f_element["length"]),
    )


# ----------------------------------------------------------------------
# Geometry topology / face data
# ----------------------------------------------------------------------

@dataclass
class EdgeTopology:
    unique_edges: np.ndarray
    incidence: np.ndarray
    face_a: np.ndarray
    face_b: np.ndarray
    manifold_edges: np.ndarray
    nonmanifold_face_groups: list[np.ndarray]


def build_edge_topology(faces: np.ndarray) -> EdgeTopology:
    n = len(faces)
    edge_pairs = ((0, 1), (1, 2), (2, 0))
    edges = np.vstack([faces[:, [a, b]] for a, b in edge_pairs])
    occ_face = np.tile(np.arange(n, dtype=np.int64), 3)
    lo = np.minimum(edges[:, 0], edges[:, 1])
    hi = np.maximum(edges[:, 0], edges[:, 1])
    order = np.lexsort((hi, lo))
    los = lo[order]
    his = hi[order]
    starts = np.r_[0, np.flatnonzero((los[1:] != los[:-1]) | (his[1:] != his[:-1])) + 1]
    ends = np.r_[starts[1:], len(order)]
    counts = (ends - starts).astype(np.int64)
    first = order[starts]
    unique_edges = np.column_stack([lo[first], hi[first]]).astype(np.int64)

    two = counts == 2
    two_starts = starts[two]
    i0 = order[two_starts]
    i1 = order[two_starts + 1]
    face_a = occ_face[i0]
    face_b = occ_face[i1]
    manifold_edges = unique_edges[two]
    nonmanifold_face_groups = []
    for start, count in zip(starts[counts > 2], counts[counts > 2]):
        occ = order[start:start + count]
        nonmanifold_face_groups.append(np.unique(occ_face[occ]).astype(np.int64))
    return EdgeTopology(
        unique_edges, counts, face_a, face_b, manifold_edges, nonmanifold_face_groups
    )


def face_geometry(vertices: np.ndarray, faces: np.ndarray):
    tri = vertices[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    normals = np.zeros_like(cross, dtype=np.float64)
    valid = twice_area > np.finfo(np.float64).eps
    normals[valid] = cross[valid] / twice_area[valid, None]
    areas = 0.5 * twice_area
    centers = tri.mean(axis=1)
    signed_volume = np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])) / 6.0
    return centers, normals, areas, signed_volume


def topology_components(n_faces: int, edge_topology: EdgeTopology):
    fa = edge_topology.face_a
    fb = edge_topology.face_b
    row_parts = [fa, fb]
    col_parts = [fb, fa]
    for faces_here in edge_topology.nonmanifold_face_groups:
        if len(faces_here) > 1:
            anchor = faces_here[0]
            other = faces_here[1:]
            row_parts.extend([np.full(len(other), anchor, dtype=np.int64), other])
            col_parts.extend([other, np.full(len(other), anchor, dtype=np.int64)])
    rows = np.concatenate(row_parts) if row_parts else np.empty(0, dtype=np.int64)
    cols = np.concatenate(col_parts) if col_parts else np.empty(0, dtype=np.int64)
    graph = coo_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, cols)),
        shape=(n_faces, n_faces),
    ).tocsr()
    count, labels = connected_components(graph, directed=False, return_labels=True)
    return int(count), labels.astype(np.int64)


def component_statistics(vertices, faces, labels_raw, areas_native, signed_volume, scale_to_mm, min_area_mm2):
    raw_count = int(labels_raw.max()) + 1 if len(labels_raw) else 0
    face_counts = np.bincount(labels_raw, minlength=raw_count)
    area_sums = np.bincount(labels_raw, weights=areas_native, minlength=raw_count)
    volume_sums = np.bincount(labels_raw, weights=signed_volume, minlength=raw_count)
    old_ids = np.arange(raw_count, dtype=np.int64)
    order = np.lexsort((old_ids, -face_counts))
    remap = np.empty(raw_count, dtype=np.int64)
    remap[order] = np.arange(1, raw_count + 1, dtype=np.int64)
    component_ids = remap[labels_raw]
    rows = []
    active_ids = []
    orientation_sign = np.ones(raw_count + 1, dtype=np.float64)
    for old_id in order:
        cid = int(remap[old_id])
        fids = np.flatnonzero(labels_raw == old_id)
        vids = np.unique(faces[fids].reshape(-1))
        pts = vertices[vids]
        area_native = float(area_sums[old_id])
        area_mm2 = area_native * scale_to_mm ** 2
        vol = float(volume_sums[old_id])
        sign = 1.0 if vol >= 0 else -1.0
        orientation_sign[cid] = sign
        active = bool(area_mm2 >= min_area_mm2)
        if active:
            active_ids.append(cid)
        rows.append({
            "component_id": cid,
            "faces": int(len(fids)),
            "vertices": int(len(vids)),
            "surface_area_mm2": float(area_mm2),
            "signed_volume_native3": vol,
            "normal_orientation_multiplier": sign,
            "active_for_surface_analysis": active,
            "bbox_min_x": float(pts[:, 0].min()),
            "bbox_min_y": float(pts[:, 1].min()),
            "bbox_min_z": float(pts[:, 2].min()),
            "bbox_max_x": float(pts[:, 0].max()),
            "bbox_max_y": float(pts[:, 1].max()),
            "bbox_max_z": float(pts[:, 2].max()),
        })
    active_mask = np.isin(component_ids, np.asarray(active_ids, dtype=np.int64))
    return component_ids, rows, orientation_sign, active_mask


# ----------------------------------------------------------------------
# Outer / inner surface classification
# ----------------------------------------------------------------------

def estimate_axis_xy(vertices, faces, active_face_mask, z_bins=32):
    active_faces = faces[np.asarray(active_face_mask, dtype=bool)]
    if len(active_faces) == 0:
        raise RuntimeError("No active faces for vessel-axis estimation.")
    vids = np.unique(active_faces.reshape(-1))
    pts = vertices[vids]
    z = pts[:, 2]
    zlo, zhi = np.quantile(z, [0.05, 0.95])
    bins = np.linspace(zlo, zhi, z_bins + 1)
    centers = []
    for i in range(z_bins):
        mask = (z >= bins[i]) & ((z <= bins[i + 1]) if i == z_bins - 1 else (z < bins[i + 1]))
        slab = pts[mask]
        if len(slab) < 100:
            continue
        qx = np.quantile(slab[:, 0], [0.02, 0.98])
        qy = np.quantile(slab[:, 1], [0.02, 0.98])
        centers.append([0.5 * (qx[0] + qx[1]), 0.5 * (qy[0] + qy[1])])
    if not centers:
        bmin = pts.min(axis=0); bmax = pts.max(axis=0)
        return np.array([0.5 * (bmin[0] + bmax[0]), 0.5 * (bmin[1] + bmax[1])]), {
            "method": "active-geometry bbox center fallback", "slice_centers_used": 0,
        }
    centers = np.asarray(centers)
    axis = np.median(centers, axis=0)
    return axis, {
        "method": "median of robust per-Z-slice centers",
        "z_bins_requested": int(z_bins),
        "slice_centers_used": int(len(centers)),
    }


def classify_surfaces(face_centers, face_normals, component_ids, orientation_sign, active_mask, axis_xy, normal_threshold):
    corrected = np.asarray(face_normals, dtype=np.float64).copy()
    corrected *= orientation_sign[component_ids][:, None]
    radial_xy = face_centers[:, :2] - axis_xy[None, :]
    radius = np.linalg.norm(radial_xy, axis=1)
    radial_unit = np.zeros_like(corrected)
    valid = radius > np.finfo(np.float64).eps
    radial_unit[valid, 0] = radial_xy[valid, 0] / radius[valid]
    radial_unit[valid, 1] = radial_xy[valid, 1] / radius[valid]
    radial_score = np.einsum("ij,ij->i", corrected, radial_unit)
    cls = np.full(len(face_centers), SURFACE_IGNORED, dtype=np.int8)
    active = np.asarray(active_mask, dtype=bool)
    cls[active] = SURFACE_TRANSITION
    cls[active & (radial_score >= normal_threshold)] = SURFACE_OUTER
    cls[active & (radial_score <= -normal_threshold)] = SURFACE_INNER
    return cls, radial_score, corrected, radius


def export_surface_products(directory, vertices, faces, vertex_colors, surface_class):
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for name, cls in (("outer", SURFACE_OUTER), ("inner", SURFACE_INNER), ("transition", SURFACE_TRANSITION), ("ignored", SURFACE_IGNORED)):
        ids = np.flatnonzero(surface_class == cls)
        mesh = compact_face_subset(vertices, faces, ids, vertex_colors)
        p = directory / f"{name}_surface.ply"
        mesh.export(str(p))
        outputs[f"{name}_surface_ply"] = str(p)
    return outputs


# ----------------------------------------------------------------------
# Same-surface diffusion graph and geometry features
# ----------------------------------------------------------------------

def build_same_surface_diffusion(n_faces, edge_topology, surface_class):
    fa = edge_topology.face_a
    fb = edge_topology.face_b
    valid_class = (
        ((surface_class[fa] == SURFACE_OUTER) & (surface_class[fb] == SURFACE_OUTER))
        | ((surface_class[fa] == SURFACE_INNER) & (surface_class[fb] == SURFACE_INNER))
    )
    a = fa[valid_class]
    b = fb[valid_class]
    self_ids = np.arange(n_faces, dtype=np.int64)
    row = np.r_[a, b, self_ids]
    col = np.r_[b, a, self_ids]
    vals = np.ones(len(row), dtype=np.float32)
    adjacency = coo_matrix((vals, (row, col)), shape=(n_faces, n_faces)).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    diffusion = diags(1.0 / np.maximum(degree, 1.0)).dot(adjacency).tocsr()
    return diffusion, a, b


def compute_multiscale_geometry_features(
    centers_native,
    corrected_normals,
    same_face_a,
    same_face_b,
    diffusion,
    scale_to_mm,
    iterations=(2, 8, 32),
):
    iterations = tuple(sorted(set(int(x) for x in iterations if int(x) > 0)))
    if len(iterations) < 2:
        raise ValueError("At least two positive smoothing iteration targets are required.")

    normals0 = np.asarray(corrected_normals, dtype=np.float32)
    centers0 = np.asarray(centers_native, dtype=np.float32)
    sm_n = normals0.copy()
    sm_c = centers0.copy()
    normal_dev = {}
    residual_abs_mm = {}
    current = 0
    for target in iterations:
        while current < target:
            sm_n = diffusion.dot(sm_n)
            nn = np.linalg.norm(sm_n, axis=1)
            sm_n /= np.maximum(nn[:, None], 1e-12)
            sm_c = diffusion.dot(sm_c)
            current += 1
        dot = np.einsum("ij,ij->i", normals0, sm_n)
        normal_dev[target] = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))).astype(np.float32)
        delta = centers0 - sm_c
        signed = np.einsum("ij,ij->i", delta, sm_n) * float(scale_to_mm)
        residual_abs_mm[target] = np.abs(signed).astype(np.float32)

    # Local edge dihedral response on OUTER/INNER same-surface adjacency.
    dot = np.einsum(
        "ij,ij->i",
        corrected_normals[same_face_a],
        corrected_normals[same_face_b],
    )
    edge_angle = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))).astype(np.float32)
    dihedral_max = np.zeros(len(centers_native), dtype=np.float32)
    np.maximum.at(dihedral_max, same_face_a, edge_angle)
    np.maximum.at(dihedral_max, same_face_b, edge_angle)

    center_step_mm = np.linalg.norm(
        centers_native[same_face_a] - centers_native[same_face_b], axis=1
    ) * scale_to_mm
    median_step = float(np.median(center_step_mm)) if len(center_step_mm) else None

    return {
        "normal_dev_deg": normal_dev,
        "residual_abs_mm": residual_abs_mm,
        "dihedral_max_deg": dihedral_max,
        "median_same_surface_face_step_mm": median_step,
        "iterations": iterations,
        "approx_diffusion_scale_mm": {
            str(k): (float(median_step * math.sqrt(k)) if median_step is not None else None)
            for k in iterations
        },
    }


def build_geometry_score(features, surface_class, diffusion):
    """
    Build one continuous geometry-only score per face.

    The smallest scale is exported for QA but deliberately excluded from the
    primary score so that very fine scanning/striated texture is down-weighted.
    """
    iterations = features["iterations"]
    medium = iterations[-2]
    large = iterations[-1]
    n_med = features["normal_dev_deg"][medium]
    n_large = features["normal_dev_deg"][large]
    r_med = features["residual_abs_mm"][medium]
    r_large = features["residual_abs_mm"][large]
    dih = features["dihedral_max_deg"]

    score = np.zeros(len(surface_class), dtype=np.float32)
    normalization = {}
    feature_norm = {}

    for cls, name in ((SURFACE_OUTER, "outer"), (SURFACE_INNER, "inner")):
        mask = surface_class == cls
        nm, q_nm = robust_unit_interval(n_med, mask)
        nl, q_nl = robust_unit_interval(n_large, mask)
        rm, q_rm = robust_unit_interval(r_med, mask)
        rl, q_rl = robust_unit_interval(r_large, mask)
        dh, q_dh = robust_unit_interval(dih, mask)
        raw = 0.30 * nm + 0.25 * nl + 0.20 * rm + 0.15 * rl + 0.10 * dh
        persistence = np.minimum(nm, nl)
        raw = np.clip(0.85 * raw + 0.15 * persistence, 0.0, 1.0)
        score[mask] = raw[mask]
        normalization[name] = {
            f"normal_dev_iter_{medium}": q_nm,
            f"normal_dev_iter_{large}": q_nl,
            f"residual_iter_{medium}": q_rm,
            f"residual_iter_{large}": q_rl,
            "dihedral_max_deg": q_dh,
        }
        feature_norm[name] = {
            "normal_medium": nm,
            "normal_large": nl,
            "residual_medium": rm,
            "residual_large": rl,
            "dihedral": dh,
        }

    # One diffusion step turns face-level noise into a slightly more coherent
    # band without attempting to trace a final line.
    score_smoothed = diffusion.dot(score.astype(np.float32)).astype(np.float32)
    valid = (surface_class == SURFACE_OUTER) | (surface_class == SURFACE_INNER)
    score_smoothed[~valid] = 0.0
    return score_smoothed, normalization, feature_norm, medium, large


# ----------------------------------------------------------------------
# Opposite-surface correspondence and thickness support
# ----------------------------------------------------------------------

def cylindrical_embedding(theta, z_native, r_ref_native):
    return np.column_stack([
        r_ref_native * np.cos(theta),
        r_ref_native * np.sin(theta),
        z_native,
    ])


@dataclass
class SurfaceMatch:
    source_faces: np.ndarray
    target_faces: np.ndarray
    param_distance_mm: np.ndarray
    thickness_mm: np.ndarray
    opposition_cos: np.ndarray
    quality: np.ndarray


def match_opposite_surface(
    face_centers,
    corrected_normals,
    face_radius,
    surface_class,
    axis_xy,
    scale_to_mm,
    source_class,
    target_class,
    max_param_mm,
    thickness_min_mm,
    thickness_max_mm,
    min_opposition_cos,
    k_neighbors,
):
    source_ids = np.flatnonzero(surface_class == source_class)
    target_ids = np.flatnonzero(surface_class == target_class)
    if len(source_ids) == 0 or len(target_ids) == 0:
        return SurfaceMatch(
            source_ids,
            np.full(len(source_ids), -1, dtype=np.int64),
            np.full(len(source_ids), np.nan),
            np.full(len(source_ids), np.nan),
            np.full(len(source_ids), np.nan),
            np.zeros(len(source_ids)),
        )

    def theta(ids):
        rel = face_centers[ids, :2] - axis_xy[None, :]
        return np.arctan2(rel[:, 1], rel[:, 0])

    r_ref = float(np.median(np.r_[face_radius[source_ids], face_radius[target_ids]]))
    src_embed = cylindrical_embedding(theta(source_ids), face_centers[source_ids, 2], r_ref) * scale_to_mm
    dst_embed = cylindrical_embedding(theta(target_ids), face_centers[target_ids, 2], r_ref) * scale_to_mm
    tree = cKDTree(dst_embed)
    k = min(max(1, int(k_neighbors)), len(target_ids))
    distances, indices = tree.query(src_embed, k=k)
    if distances.ndim == 1:
        distances = distances[:, None]
        indices = indices[:, None]

    # Vectorized candidate matrices.
    cand_target_faces = target_ids[indices]
    src_r = face_radius[source_ids][:, None]
    dst_r = face_radius[cand_target_faces]
    if source_class == SURFACE_OUTER:
        thick = (src_r - dst_r) * scale_to_mm
    else:
        thick = (dst_r - src_r) * scale_to_mm

    src_n = corrected_normals[source_ids][:, None, :]
    dst_n = corrected_normals[cand_target_faces]
    opp = -np.einsum("nki,nki->nk", np.broadcast_to(src_n, dst_n.shape), dst_n)

    valid = (
        (distances <= max_param_mm)
        & (thick >= thickness_min_mm)
        & (thick <= thickness_max_mm)
        & (opp >= min_opposition_cos)
    )
    objective = distances + 1.5 * (1.0 - opp)
    objective = np.where(valid, objective, np.inf)
    best_q = np.argmin(objective, axis=1)
    best_obj = objective[np.arange(len(source_ids)), best_q]
    ok = np.isfinite(best_obj)

    target_face = np.full(len(source_ids), -1, dtype=np.int64)
    param = np.full(len(source_ids), np.nan, dtype=np.float64)
    thickness = np.full(len(source_ids), np.nan, dtype=np.float64)
    opposition = np.full(len(source_ids), np.nan, dtype=np.float64)
    target_face[ok] = cand_target_faces[np.arange(len(source_ids))[ok], best_q[ok]]
    param[ok] = distances[np.arange(len(source_ids))[ok], best_q[ok]]
    thickness[ok] = thick[np.arange(len(source_ids))[ok], best_q[ok]]
    opposition[ok] = opp[np.arange(len(source_ids))[ok], best_q[ok]]

    # A transparent support-quality score. It is NOT a join score.
    quality = np.zeros(len(source_ids), dtype=np.float32)
    if np.any(ok):
        q_param = np.clip(1.0 - param[ok] / max_param_mm, 0.0, 1.0)
        q_opp = np.clip((opposition[ok] - min_opposition_cos) / max(1e-9, 1.0 - min_opposition_cos), 0.0, 1.0)
        quality[ok] = (0.65 * q_param + 0.35 * q_opp).astype(np.float32)

    return SurfaceMatch(source_ids, target_face, param, thickness, opposition, quality)


def thickness_anomaly_score(match: SurfaceMatch, surface_class, diffusion):
    n = len(surface_class)
    thickness = np.full(n, np.nan, dtype=np.float32)
    valid = match.target_faces >= 0
    thickness[match.source_faces[valid]] = match.thickness_mm[valid].astype(np.float32)

    values = np.where(np.isfinite(thickness), thickness, 0.0).astype(np.float32)
    weights = np.isfinite(thickness).astype(np.float32)
    sv = values.copy(); sw = weights.copy()
    # Four local averaging passes. The numerator/denominator formulation avoids
    # pulling missing values toward zero.
    for _ in range(4):
        sv = diffusion.dot(sv)
        sw = diffusion.dot(sw)
    smooth = np.full(n, np.nan, dtype=np.float32)
    ok = sw > 1e-6
    smooth[ok] = sv[ok] / sw[ok]
    anomaly = np.full(n, np.nan, dtype=np.float32)
    both = np.isfinite(thickness) & np.isfinite(smooth)
    anomaly[both] = np.abs(thickness[both] - smooth[both])
    score = np.zeros(n, dtype=np.float32)
    source_mask = np.zeros(n, dtype=bool); source_mask[match.source_faces] = True
    score, norm = robust_unit_interval(anomaly, source_mask & np.isfinite(anomaly), 0.50, 0.98)
    return thickness, smooth, anomaly, score, norm


def combine_opposite_support(primary_score, opposite_score, match, thickness_score):
    combined = np.asarray(primary_score, dtype=np.float32).copy()
    support = np.zeros(len(primary_score), dtype=np.float32)
    valid = match.target_faces >= 0
    if np.any(valid):
        src = match.source_faces[valid]
        dst = match.target_faces[valid]
        paired = np.sqrt(
            np.clip(primary_score[src], 0.0, 1.0)
            * np.clip(opposite_score[dst], 0.0, 1.0)
        )
        support[src] = (paired * match.quality[valid]).astype(np.float32)
        # Opposite-side evidence is a bonus, never a prerequisite.
        combined[src] = np.clip(
            primary_score[src]
            + 0.20 * support[src]
            + 0.10 * thickness_score[src],
            0.0,
            1.0,
        )
    return combined, support


# ----------------------------------------------------------------------
# PLY/mesh score exports
# ----------------------------------------------------------------------

def score_colors(score: np.ndarray):
    """Simple perceptual-ish dark-purple -> red -> yellow mapping without matplotlib."""
    s = np.clip(np.asarray(score, dtype=np.float64), 0.0, 1.0)
    r = np.clip(1.6 * s, 0.08, 1.0)
    g = np.clip(1.8 * (s - 0.45), 0.02, 1.0)
    b = np.clip(0.35 - 0.25 * s, 0.02, 0.35)
    a = np.ones_like(s)
    return (255.0 * np.column_stack([r, g, b, a])).astype(np.uint8)


def export_face_score_mesh(path, vertices, faces, face_score, valid_mask):
    face_ids = np.flatnonzero(valid_mask)
    mesh = compact_face_subset(vertices, faces, face_ids)
    if len(face_ids):
        mesh.visual.face_colors = score_colors(face_score[face_ids])
    mesh.export(str(path))


def export_threshold_face_mesh(path, vertices, faces, face_ids):
    mesh = compact_face_subset(vertices, faces, face_ids)
    mesh.export(str(path))


def write_binary_point_scalar_ply(
    path: Path, points: np.ndarray, fields: dict[str, np.ndarray], normals: np.ndarray | None = None
):
    """Write face-center points with normals and float scalar properties for CloudCompare."""
    path.parent.mkdir(parents=True, exist_ok=True)
    points = np.asarray(points, dtype=np.float32)
    names = list(fields.keys())
    arrays = [np.asarray(fields[name], dtype=np.float32) for name in names]
    n = len(points)
    if any(len(a) != n for a in arrays):
        raise ValueError("Scalar field length mismatch.")
    normal_array = None if normals is None else np.asarray(normals, dtype=np.float32)
    if normal_array is not None and normal_array.shape != points.shape:
        raise ValueError("Normal array shape mismatch.")
    dtype = [("x", "<f4"), ("y", "<f4"), ("z", "<f4")]
    if normal_array is not None:
        dtype += [("nx", "<f4"), ("ny", "<f4"), ("nz", "<f4")]
    dtype += [(name, "<f4") for name in names]
    rec = np.empty(n, dtype=dtype)
    rec["x"] = points[:, 0]; rec["y"] = points[:, 1]; rec["z"] = points[:, 2]
    if normal_array is not None:
        rec["nx"] = normal_array[:, 0]; rec["ny"] = normal_array[:, 1]; rec["nz"] = normal_array[:, 2]
    for name, arr in zip(names, arrays):
        rec[name] = arr
    header = [
        "ply", "format binary_little_endian 1.0",
        "comment PotteryFragmentBoundaryExtractor geometry score face centers",
        f"element vertex {n}",
        "property float x", "property float y", "property float z",
    ]
    if normal_array is not None:
        header += ["property float nx", "property float ny", "property float nz"]
    header += [f"property float {name}" for name in names]
    header += ["end_header"]
    with path.open("wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        rec.tofile(f)


def export_score_outputs(
    directory,
    side_name,
    side_class,
    vertices,
    faces,
    centers,
    surface_class,
    geometry_score,
    combined_score,
    support_score,
    thickness_mm,
    thickness_anomaly_mm,
    features,
    medium_iter,
    large_iter,
    candidate_percentiles,
    corrected_normals,
):
    directory.mkdir(parents=True, exist_ok=True)
    mask = surface_class == side_class
    ids = np.flatnonzero(mask)
    outputs = {}

    mesh_path = directory / f"{side_name}_geometry_score_colored.ply"
    export_face_score_mesh(mesh_path, vertices, faces, geometry_score, mask)
    outputs["geometry_score_colored_mesh"] = str(mesh_path)

    comb_path = directory / f"{side_name}_combined_score_colored.ply"
    export_face_score_mesh(comb_path, vertices, faces, combined_score, mask)
    outputs["combined_score_colored_mesh"] = str(comb_path)

    scalar_path = directory / f"{side_name}_geometry_score_face_centers.ply"
    write_binary_point_scalar_ply(
        scalar_path,
        centers[ids],
        {
            "geometry_score": geometry_score[ids],
            "combined_score": combined_score[ids],
            "opposite_support": support_score[ids],
            "normal_dev_medium_deg": features["normal_dev_deg"][medium_iter][ids],
            "normal_dev_large_deg": features["normal_dev_deg"][large_iter][ids],
            "residual_medium_mm": features["residual_abs_mm"][medium_iter][ids],
            "residual_large_mm": features["residual_abs_mm"][large_iter][ids],
            "dihedral_max_deg": features["dihedral_max_deg"][ids],
            "wall_thickness_mm": np.nan_to_num(thickness_mm[ids], nan=-1.0),
            "thickness_anomaly_mm": np.nan_to_num(thickness_anomaly_mm[ids], nan=-1.0),
        },
        normals=corrected_normals[ids],
    )
    outputs["face_center_scalar_ply"] = str(scalar_path)

    thresholds = {}
    vals = geometry_score[ids]
    for pct in candidate_percentiles:
        q = float(np.quantile(vals, pct / 100.0))
        selected = ids[geometry_score[ids] >= q]
        p = directory / f"{side_name}_geometry_candidates_p{int(pct):02d}.ply"
        export_threshold_face_mesh(p, vertices, faces, selected)
        outputs[f"candidate_p{int(pct):02d}_mesh"] = str(p)
        thresholds[f"p{int(pct):02d}"] = {
            "score_threshold": q,
            "faces": int(len(selected)),
            "fraction_of_side_faces": float(len(selected) / len(ids)) if len(ids) else None,
        }

    return outputs, thresholds


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------

def analyze(
    input_path: Path,
    unit: str,
    output_dir: Path | None,
    surface_normal_threshold: float,
    min_component_area_mm2: float,
    axis_x: float | None,
    axis_y: float | None,
    smooth_iterations: tuple[int, ...],
    candidate_percentiles: tuple[int, ...],
    pair_max_param_mm: float,
    pair_thickness_min_mm: float,
    pair_thickness_max_mm: float,
    pair_min_opposition_cos: float,
    pair_k_neighbors: int,
):
    started = time.perf_counter()
    input_path = Path(input_path)
    scale_to_mm = UNIT_SCALE_TO_MM[unit]
    base = Path(output_dir) if output_dir is not None else input_path.parent / f"{input_path.stem}_FragmentBoundary_v0_3"
    qa_dir = base / "qa"
    surface_dir = base / "surfaces"
    score_dir = base / "geometry_scores"
    for d in (base, qa_dir, surface_dir, score_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"=== PotteryFragmentBoundaryExtractor v{__version__} ===")
    print("method     : geometry-only multi-scale surface anomaly")
    print("UV evidence: NOT USED")
    print(f"input      : {input_path}")
    print(f"input unit : {unit}")
    print(f"output     : {base}")

    # A. source geometry
    t = time.perf_counter()
    data = load_ply_geometry(input_path)
    vertices = data.vertices; faces = data.faces
    load_time = time.perf_counter() - t
    print("\n=== Stage A: source PLY geometry ===")
    print("texture topology fix : DISABLED (fix_texture=False)")
    print(f"vertices             : {len(vertices):,}")
    print(f"faces                : {len(faces):,}")
    print(f"stored vertex normals: {data.stored_vertex_normals is not None}")
    print(f"time                 : {load_time:.2f} s")

    # B. topology / face geometry
    t = time.perf_counter()
    edge_topo = build_edge_topology(faces)
    centers, normals, areas, signed_volume = face_geometry(vertices, faces)
    ncomp, raw_labels = topology_components(len(faces), edge_topo)
    component_ids, component_rows, orientation_sign, active_mask = component_statistics(
        vertices, faces, raw_labels, areas, signed_volume, scale_to_mm, min_component_area_mm2
    )
    topology_time = time.perf_counter() - t
    topo = {
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "unique_geometry_edges": int(len(edge_topo.unique_edges)),
        "boundary_geometry_edges": int(np.count_nonzero(edge_topo.incidence == 1)),
        "two_face_manifold_edges": int(np.count_nonzero(edge_topo.incidence == 2)),
        "nonmanifold_geometry_edges": int(np.count_nonzero(edge_topo.incidence > 2)),
        "geometry_components": int(ncomp),
        "active_faces": int(np.count_nonzero(active_mask)),
        "stored_vertex_normals_available": bool(data.stored_vertex_normals is not None),
        "reference_normal_source": "triangle winding / computed face normals",
        "uv_information_used": False,
    }
    write_json(qa_dir / "input_topology.json", topo)
    write_csv(qa_dir / "geometry_components.csv", component_rows)
    print("\n=== Stage B: geometry topology QA ===")
    print(f"boundary edges       : {topo['boundary_geometry_edges']:,}")
    print(f"nonmanifold edges    : {topo['nonmanifold_geometry_edges']:,}")
    print(f"components           : {ncomp:,}")
    print(f"active faces         : {topo['active_faces']:,}")
    print(f"time                 : {topology_time:.2f} s")

    # C. classify outer/inner
    t = time.perf_counter()
    if (axis_x is None) ^ (axis_y is None):
        raise ValueError("--axis-x and --axis-y must be supplied together.")
    if axis_x is None:
        axis_xy, axis_diag = estimate_axis_xy(vertices, faces, active_mask)
    else:
        axis_xy = np.array([axis_x, axis_y], dtype=np.float64)
        axis_diag = {"method": "explicit CLI axis"}
    surface_class, radial_score, corrected_normals, face_radius = classify_surfaces(
        centers, normals, component_ids, orientation_sign, active_mask, axis_xy, surface_normal_threshold
    )
    surface_outputs = export_surface_products(surface_dir, vertices, faces, data.vertex_colors, surface_class)
    surface_summary = {
        "axis_xy_native": [float(x) for x in axis_xy],
        "axis_unit": unit,
        "axis_diagnostics": axis_diag,
        "surface_normal_threshold": float(surface_normal_threshold),
        "outer_faces": int(np.count_nonzero(surface_class == SURFACE_OUTER)),
        "inner_faces": int(np.count_nonzero(surface_class == SURFACE_INNER)),
        "transition_faces": int(np.count_nonzero(surface_class == SURFACE_TRANSITION)),
        "ignored_faces": int(np.count_nonzero(surface_class == SURFACE_IGNORED)),
        "outputs": surface_outputs,
    }
    write_json(surface_dir / "surface_classification_summary.json", surface_summary)
    surface_time = time.perf_counter() - t
    print("\n=== Stage C: OUTER / INNER surface classification ===")
    print(f"outer      : {surface_summary['outer_faces']:,}")
    print(f"inner      : {surface_summary['inner_faces']:,}")
    print(f"transition : {surface_summary['transition_faces']:,}")
    print(f"ignored    : {surface_summary['ignored_faces']:,}")
    print(f"time       : {surface_time:.2f} s")

    # D. geometry-only multiscale score
    t = time.perf_counter()
    diffusion, same_a, same_b = build_same_surface_diffusion(len(faces), edge_topo, surface_class)
    features = compute_multiscale_geometry_features(
        centers, corrected_normals, same_a, same_b, diffusion, scale_to_mm, smooth_iterations
    )
    geometry_score, normalization, feature_norm, medium_iter, large_iter = build_geometry_score(
        features, surface_class, diffusion
    )
    geometry_time = time.perf_counter() - t
    geometry_summary = {
        "primary_evidence": "geometry only",
        "uv_information_used": False,
        "smooth_iterations": [int(x) for x in features["iterations"]],
        "median_same_surface_face_step_mm": features["median_same_surface_face_step_mm"],
        "approx_diffusion_scale_mm": features["approx_diffusion_scale_mm"],
        "score_formula": (
            "0.85*(0.30*Nmed + 0.25*Nlarge + 0.20*Rmed + 0.15*Rlarge + 0.10*D) "
            "+ 0.15*min(Nmed,Nlarge), then one same-surface diffusion step; "
            "features robustly normalized from median to 98th percentile per surface"
        ),
        "normalization": normalization,
        "outer_score_quantiles": safe_quantiles(geometry_score[surface_class == SURFACE_OUTER]),
        "inner_score_quantiles": safe_quantiles(geometry_score[surface_class == SURFACE_INNER]),
    }
    print("\n=== Stage D: geometry-only multi-scale score ===")
    print(f"smoothing iterations     : {features['iterations']}")
    print(f"median face step         : {features['median_same_surface_face_step_mm']:.3f} mm")
    print(f"approx scales            : {features['approx_diffusion_scale_mm']}")
    print(f"time                     : {geometry_time:.2f} s")

    # E. outer-inner support; never required for a geometry candidate
    t = time.perf_counter()
    outer_match = match_opposite_surface(
        centers, corrected_normals, face_radius, surface_class, axis_xy, scale_to_mm,
        SURFACE_OUTER, SURFACE_INNER, pair_max_param_mm, pair_thickness_min_mm,
        pair_thickness_max_mm, pair_min_opposition_cos, pair_k_neighbors,
    )
    inner_match = match_opposite_surface(
        centers, corrected_normals, face_radius, surface_class, axis_xy, scale_to_mm,
        SURFACE_INNER, SURFACE_OUTER, pair_max_param_mm, pair_thickness_min_mm,
        pair_thickness_max_mm, pair_min_opposition_cos, pair_k_neighbors,
    )
    outer_thick, outer_thick_smooth, outer_thick_anom, outer_thick_score, outer_thick_norm = thickness_anomaly_score(
        outer_match, surface_class, diffusion
    )
    inner_thick, inner_thick_smooth, inner_thick_anom, inner_thick_score, inner_thick_norm = thickness_anomaly_score(
        inner_match, surface_class, diffusion
    )
    outer_combined, outer_support = combine_opposite_support(
        geometry_score, geometry_score, outer_match, outer_thick_score
    )
    inner_combined, inner_support = combine_opposite_support(
        geometry_score, geometry_score, inner_match, inner_thick_score
    )
    combined_score = geometry_score.copy()
    support_score = np.zeros(len(faces), dtype=np.float32)
    combined_score[surface_class == SURFACE_OUTER] = outer_combined[surface_class == SURFACE_OUTER]
    combined_score[surface_class == SURFACE_INNER] = inner_combined[surface_class == SURFACE_INNER]
    support_score += outer_support + inner_support
    thickness_mm = np.where(np.isfinite(outer_thick), outer_thick, inner_thick)
    thickness_anom_mm = np.where(np.isfinite(outer_thick_anom), outer_thick_anom, inner_thick_anom)

    def match_summary(match):
        valid = match.target_faces >= 0
        return {
            "source_faces": int(len(match.source_faces)),
            "matched_faces": int(np.count_nonzero(valid)),
            "matched_fraction": float(np.mean(valid)) if len(valid) else None,
            "param_distance_mm_quantiles": safe_quantiles(match.param_distance_mm[valid]),
            "wall_thickness_mm_quantiles": safe_quantiles(match.thickness_mm[valid]),
            "normal_opposition_cos_quantiles": safe_quantiles(match.opposition_cos[valid]),
            "support_quality_quantiles": safe_quantiles(match.quality[valid]),
        }
    correspondence_summary = {
        "role": "supporting evidence only; opposite-side match is not required for a boundary candidate",
        "outer_to_inner": match_summary(outer_match),
        "inner_to_outer": match_summary(inner_match),
        "outer_thickness_anomaly_normalization": outer_thick_norm,
        "inner_thickness_anomaly_normalization": inner_thick_norm,
    }
    correspondence_time = time.perf_counter() - t
    print("\n=== Stage E: opposite-surface support ===")
    print(f"outer->inner matched      : {correspondence_summary['outer_to_inner']['matched_faces']:,} / {correspondence_summary['outer_to_inner']['source_faces']:,}")
    print(f"inner->outer matched      : {correspondence_summary['inner_to_outer']['matched_faces']:,} / {correspondence_summary['inner_to_outer']['source_faces']:,}")
    print("NOTE: opposite-side support is a score bonus, not a prerequisite.")
    print(f"time                      : {correspondence_time:.2f} s")

    # F. score exports
    t = time.perf_counter()
    outer_outputs, outer_thresholds = export_score_outputs(
        score_dir, "outer", SURFACE_OUTER, vertices, faces, centers, surface_class,
        geometry_score, combined_score, support_score, thickness_mm, thickness_anom_mm,
        features, medium_iter, large_iter, candidate_percentiles, corrected_normals,
    )
    inner_outputs, inner_thresholds = export_score_outputs(
        score_dir, "inner", SURFACE_INNER, vertices, faces, centers, surface_class,
        geometry_score, combined_score, support_score, thickness_mm, thickness_anom_mm,
        features, medium_iter, large_iter, candidate_percentiles, corrected_normals,
    )
    export_time = time.perf_counter() - t
    geometry_summary["outer_candidate_thresholds"] = outer_thresholds
    geometry_summary["inner_candidate_thresholds"] = inner_thresholds
    geometry_summary["outputs"] = {"outer": outer_outputs, "inner": inner_outputs}
    write_json(score_dir / "geometry_score_summary.json", geometry_summary)
    write_json(qa_dir / "outer_inner_correspondence_summary.json", correspondence_summary)

    total_time = time.perf_counter() - started
    result = {
        "program": "PotteryFragmentBoundaryExtractor",
        "program_version": __version__,
        "status": "ok",
        "input": str(input_path),
        "input_unit": unit,
        "up_axis": "+Z",
        "method": "geometry-only multi-scale surface anomaly + optional opposite-surface support",
        "uv_information_used": False,
        "interpretation_warning": (
            "v0.3.0 outputs a continuous geometry-based candidate score. "
            "Percentile candidate subsets are visualization aids, not automatic archaeological sherd boundaries."
        ),
        "topology": topo,
        "surface_classification": surface_summary,
        "geometry_score": geometry_summary,
        "outer_inner_correspondence": correspondence_summary,
        "timing_seconds": {
            "load_source_ply": float(load_time),
            "topology": float(topology_time),
            "surface_classification": float(surface_time),
            "geometry_features_and_score": float(geometry_time),
            "opposite_surface_support": float(correspondence_time),
            "score_exports": float(export_time),
            "total": float(total_time),
        },
    }
    write_json(base / "result.json", result)
    print("\n=== Result ===")
    print(f"result JSON               : {base / 'result.json'}")
    print(f"total time                : {total_time:.2f} s")
    print("Primary validation target : outer_geometry_score_face_centers.ply")
    print("UV seams were not read or used as evidence.")
    return result


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def parse_int_tuple(text: str, min_count=2):
    vals = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    vals = tuple(sorted(set(x for x in vals if x > 0)))
    if len(vals) < min_count:
        raise argparse.ArgumentTypeError(f"need at least {min_count} positive comma-separated integers")
    return vals


def print_environment():
    print(f"PotteryFragmentBoundaryExtractor {__version__}")
    print(f"Python  : {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    print(f"numpy   : {installed_version('numpy') or 'unknown'}")
    print(f"scipy   : {installed_version('scipy') or 'unknown'}")
    print(f"trimesh : {installed_version('trimesh') or 'unknown'}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Geometry-only experimental pottery joining-boundary scoring: preserve original PLY topology, "
            "separate OUTER/INNER surfaces, compute multi-scale normal/residual/dihedral anomalies, "
            "and use the opposite surface only as supporting evidence. UV seams are not used."
        )
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--diagnose-env", action="store_true")
    parser.add_argument("input", type=Path, nargs="?", help="Input triangular PLY.")
    parser.add_argument("--unit", choices=["mm", "cm", "m"], default="mm")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--surface-normal-threshold", type=float, default=0.35)
    parser.add_argument("--min-component-area-mm2", type=float, default=10.0)
    parser.add_argument("--axis-x", type=float, default=None)
    parser.add_argument("--axis-y", type=float, default=None)
    parser.add_argument(
        "--smooth-iters",
        type=lambda s: parse_int_tuple(s, 2),
        default=(2, 8, 32),
        help="Comma-separated same-surface diffusion iteration targets. Default: 2,8,32.",
    )
    parser.add_argument(
        "--candidate-percentiles",
        type=lambda s: parse_int_tuple(s, 1),
        default=(90, 95, 98),
        help="Percentile subsets exported only for visual QA. Default: 90,95,98.",
    )
    parser.add_argument("--pair-max-param-mm", type=float, default=4.0)
    parser.add_argument("--pair-thickness-min-mm", type=float, default=0.5)
    parser.add_argument("--pair-thickness-max-mm", type=float, default=20.0)
    parser.add_argument("--pair-min-opposition-cos", type=float, default=0.5)
    parser.add_argument("--pair-k-neighbors", type=int, default=12)
    args = parser.parse_args()

    if args.diagnose_env:
        print_environment(); return
    if args.input is None:
        parser.error("input PLY is required unless --diagnose-env is used.")
    if not (0.0 <= args.surface_normal_threshold < 1.0):
        parser.error("--surface-normal-threshold must be in [0,1).")
    if any(p <= 0 or p >= 100 for p in args.candidate_percentiles):
        parser.error("candidate percentiles must be between 1 and 99.")

    try:
        analyze(
            input_path=args.input,
            unit=args.unit,
            output_dir=args.output_dir,
            surface_normal_threshold=args.surface_normal_threshold,
            min_component_area_mm2=args.min_component_area_mm2,
            axis_x=args.axis_x,
            axis_y=args.axis_y,
            smooth_iterations=args.smooth_iters,
            candidate_percentiles=args.candidate_percentiles,
            pair_max_param_mm=args.pair_max_param_mm,
            pair_thickness_min_mm=args.pair_thickness_min_mm,
            pair_thickness_max_mm=args.pair_thickness_max_mm,
            pair_min_opposition_cos=args.pair_min_opposition_cos,
            pair_k_neighbors=args.pair_k_neighbors,
        )
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()


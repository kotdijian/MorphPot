#!/usr/bin/env python3
"""
Pottery multiscale geometry feature layer
Version 0.1.0

This module extracts scale-aware geometric descriptors from pottery surface meshes.
It is intentionally archaeological-label agnostic: it does not decide that a trace
is a sherd join, hake-me brushing, crack, etc.  It measures geometry first.

Core design
-----------
- Preserve original PLY geometry topology (fix_texture=False).
- Separate OUTER / INNER / TRANSITION faces using the vessel axis and face normals.
- Convert requested physical analysis scales [mm] into graph-diffusion iterations
  from the observed same-surface face spacing.
- Compute per-face descriptors at each scale:
    * normal deviation [deg]
    * signed normal residual [mm]
    * absolute residual [mm]
- Compute mesh-scale local dihedral response.
- Return sufficient metadata to compare meshes at different face densities.

The physical scale mapping is approximate:
    scale_mm ~= median_same_surface_face_step_mm * sqrt(iterations)
This is a graph-diffusion calibration, not a geodesic-radius guarantee.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from trimesh.exchange.ply import load_ply

__version__ = "0.1.0"

UNIT_SCALE_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}

SURFACE_IGNORED = 0
SURFACE_OUTER = 1
SURFACE_INNER = 2
SURFACE_TRANSITION = 3

SURFACE_NAMES = {
    SURFACE_IGNORED: "ignored",
    SURFACE_OUTER: "outer",
    SURFACE_INNER: "inner",
    SURFACE_TRANSITION: "transition",
}


@dataclass
class MeshData:
    vertices: np.ndarray
    faces: np.ndarray
    vertex_colors: np.ndarray | None
    stored_vertex_normals: np.ndarray | None


@dataclass
class EdgeTopology:
    edges: np.ndarray
    incidence: np.ndarray
    face_a: np.ndarray
    face_b: np.ndarray


@dataclass
class MultiscaleFeatures:
    requested_scales_mm: tuple[float, ...]
    iterations_by_scale: dict[float, int]
    realized_scales_mm: dict[float, float]
    clipped_scales: dict[float, bool]
    median_same_surface_face_step_mm: float
    normal_dev_deg: dict[float, np.ndarray]
    signed_residual_mm: dict[float, np.ndarray]
    abs_residual_mm: dict[float, np.ndarray]
    dihedral_max_deg: np.ndarray
    dihedral_mean_deg: np.ndarray


def load_ply_preserve_topology(path: Path) -> MeshData:
    """Load triangular PLY without UV-driven geometry re-indexing."""
    path = Path(path)
    if path.suffix.lower() != ".ply":
        raise ValueError("v0.1.0 currently supports PLY input only.")

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
        raise RuntimeError("Trimesh did not expose raw PLY vertex/face elements.")

    vdata = raw["vertex"]["data"]
    fdata = raw["face"]["data"]

    for name in ("x", "y", "z"):
        if vdata.dtype.names is None or name not in vdata.dtype.names:
            raise ValueError(f"PLY vertex property '{name}' is missing.")

    vertices = np.column_stack([
        np.asarray(vdata["x"], dtype=np.float64),
        np.asarray(vdata["y"], dtype=np.float64),
        np.asarray(vdata["z"], dtype=np.float64),
    ])

    field = fdata["vertex_indices"]
    if field.dtype.names and "f0" in field.dtype.names and "f1" in field.dtype.names:
        counts = np.asarray(field["f0"], dtype=np.int64)
        values = np.asarray(field["f1"], dtype=np.int64)
    elif field.dtype == object:
        seq = [np.asarray(x, dtype=np.int64) for x in field]
        counts = np.asarray([len(x) for x in seq], dtype=np.int64)
        if len(set(counts.tolist())) != 1:
            raise ValueError("Variable-length polygons are unsupported in v0.1.0.")
        values = np.vstack(seq)
    else:
        raise TypeError(f"Unsupported PLY face representation: {field.dtype}")

    if not np.all(counts == 3):
        raise ValueError("Input PLY must already be triangulated.")
    faces = np.asarray(values[:, :3], dtype=np.int64)

    vertex_colors = None
    if vdata.dtype.names and all(x in vdata.dtype.names for x in ("red", "green", "blue")):
        vertex_colors = np.column_stack([
            np.asarray(vdata["red"], dtype=np.uint8),
            np.asarray(vdata["green"], dtype=np.uint8),
            np.asarray(vdata["blue"], dtype=np.uint8),
        ])

    stored_vertex_normals = None
    if vdata.dtype.names and all(x in vdata.dtype.names for x in ("nx", "ny", "nz")):
        stored_vertex_normals = np.column_stack([
            np.asarray(vdata["nx"], dtype=np.float64),
            np.asarray(vdata["ny"], dtype=np.float64),
            np.asarray(vdata["nz"], dtype=np.float64),
        ])

    if len(vertices) == 0 or len(faces) == 0:
        raise ValueError("PLY contains no vertices or faces.")
    if not np.isfinite(vertices).all():
        raise ValueError("Vertex coordinates contain NaN or Inf.")
    if faces.min() < 0 or faces.max() >= len(vertices):
        raise ValueError("Face index outside vertex array.")

    return MeshData(vertices, faces, vertex_colors, stored_vertex_normals)


def build_edge_topology(faces: np.ndarray) -> EdgeTopology:
    """Return unique geometry edges, incidence, and 2-face adjacency."""
    faces = np.asarray(faces, dtype=np.int64)
    n = len(faces)
    edge_occ = np.vstack([
        faces[:, [0, 1]],
        faces[:, [1, 2]],
        faces[:, [2, 0]],
    ])
    edge_occ = np.sort(edge_occ, axis=1)
    occ_face = np.tile(np.arange(n, dtype=np.int64), 3)

    order = np.lexsort((edge_occ[:, 1], edge_occ[:, 0]))
    es = edge_occ[order]
    fs = occ_face[order]
    starts = np.r_[0, np.flatnonzero(np.any(es[1:] != es[:-1], axis=1)) + 1]
    ends = np.r_[starts[1:], len(es)]
    incidence = (ends - starts).astype(np.int64)
    edges = es[starts]

    two = incidence == 2
    ts = starts[two]
    face_a = fs[ts]
    face_b = fs[ts + 1]
    return EdgeTopology(edges, incidence, face_a, face_b)


def face_components(n_faces: int, edge_topology: EdgeTopology) -> tuple[int, np.ndarray]:
    a = edge_topology.face_a
    b = edge_topology.face_b
    if len(a) == 0:
        return n_faces, np.arange(n_faces, dtype=np.int64)
    rows = np.r_[a, b]
    cols = np.r_[b, a]
    graph = coo_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, cols)),
        shape=(n_faces, n_faces),
    ).tocsr()
    count, labels = connected_components(graph, directed=False, return_labels=True)
    return int(count), np.asarray(labels, dtype=np.int64)


def component_orientation_and_activity(
    vertices: np.ndarray,
    faces: np.ndarray,
    component_ids: np.ndarray,
    min_component_area_mm2: float,
    scale_to_mm: float,
):
    """Orient closed components consistently and reject tiny disconnected debris."""
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False, validate=False)
    ncomp = int(component_ids.max()) + 1 if len(component_ids) else 0
    orientation_sign = np.ones(ncomp, dtype=np.float64)
    component_area_mm2 = np.zeros(ncomp, dtype=np.float64)
    component_face_count = np.zeros(ncomp, dtype=np.int64)

    tri = np.asarray(mesh.triangles, dtype=np.float64)
    face_area = np.asarray(mesh.area_faces, dtype=np.float64) * (scale_to_mm ** 2)
    signed_tetra = np.einsum(
        "ij,ij->i",
        tri[:, 0],
        np.cross(tri[:, 1], tri[:, 2]),
    ) / 6.0

    for cid in range(ncomp):
        ids = np.flatnonzero(component_ids == cid)
        component_face_count[cid] = len(ids)
        component_area_mm2[cid] = float(face_area[ids].sum())
        vol = float(signed_tetra[ids].sum())
        orientation_sign[cid] = 1.0 if vol >= 0.0 else -1.0

    active_components = component_area_mm2 >= float(min_component_area_mm2)
    active_face_mask = active_components[component_ids]
    return orientation_sign, active_face_mask, component_area_mm2, component_face_count


def estimate_axis_xy(vertices: np.ndarray, faces: np.ndarray, active_face_mask: np.ndarray, z_bins: int = 32):
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
        return np.array([0.5 * (bmin[0] + bmax[0]), 0.5 * (bmin[1] + bmax[1])])
    return np.median(np.asarray(centers), axis=0)


def classify_surfaces(
    face_centers: np.ndarray,
    face_normals: np.ndarray,
    component_ids: np.ndarray,
    orientation_sign: np.ndarray,
    active_mask: np.ndarray,
    axis_xy: np.ndarray,
    normal_threshold: float,
):
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


def build_same_surface_diffusion(
    n_faces: int,
    edge_topology: EdgeTopology,
    surface_class: np.ndarray,
):
    fa = edge_topology.face_a
    fb = edge_topology.face_b
    same = (
        ((surface_class[fa] == SURFACE_OUTER) & (surface_class[fb] == SURFACE_OUTER))
        | ((surface_class[fa] == SURFACE_INNER) & (surface_class[fb] == SURFACE_INNER))
    )
    a = fa[same]
    b = fb[same]
    self_ids = np.arange(n_faces, dtype=np.int64)
    rows = np.r_[a, b, self_ids]
    cols = np.r_[b, a, self_ids]
    vals = np.ones(len(rows), dtype=np.float32)
    adjacency = coo_matrix((vals, (rows, cols)), shape=(n_faces, n_faces)).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    diffusion = diags(1.0 / np.maximum(degree, 1.0)).dot(adjacency).tocsr()
    return diffusion, a, b


def median_same_surface_step_mm(
    centers_native: np.ndarray,
    same_face_a: np.ndarray,
    same_face_b: np.ndarray,
    scale_to_mm: float,
) -> float:
    if len(same_face_a) == 0:
        raise RuntimeError("No same-surface face adjacency available.")
    step = np.linalg.norm(
        centers_native[same_face_a] - centers_native[same_face_b], axis=1
    ) * float(scale_to_mm)
    return float(np.median(step))


def physical_scales_to_iterations(
    scales_mm: tuple[float, ...],
    median_step_mm: float,
    max_iterations: int = 160,
):
    if median_step_mm <= 0:
        raise ValueError("median_step_mm must be positive.")
    requested = tuple(sorted(set(float(s) for s in scales_mm if float(s) > 0)))
    if len(requested) < 2:
        raise ValueError("At least two positive physical scales are required.")

    iterations_by_scale = {}
    realized = {}
    clipped = {}
    used = set()
    for scale in requested:
        raw = max(1, int(round((scale / median_step_mm) ** 2)))
        n = min(raw, int(max_iterations))
        # Ensure distinct iteration targets so each requested scale produces a new state.
        while n in used and n < max_iterations:
            n += 1
        if n in used:
            # Multiple requested scales collapsed to max_iterations. Keep the duplicate
            # mapping but mark it clipped; caller can report limited scale separation.
            pass
        used.add(n)
        iterations_by_scale[scale] = int(n)
        realized[scale] = float(median_step_mm * math.sqrt(n))
        clipped[scale] = bool(raw > max_iterations)
    return requested, iterations_by_scale, realized, clipped


def compute_multiscale_features_mm(
    centers_native: np.ndarray,
    corrected_normals: np.ndarray,
    same_face_a: np.ndarray,
    same_face_b: np.ndarray,
    diffusion,
    scale_to_mm: float,
    scales_mm=(0.5, 1.0, 2.0, 4.0),
    max_iterations: int = 160,
) -> MultiscaleFeatures:
    median_step = median_same_surface_step_mm(
        centers_native, same_face_a, same_face_b, scale_to_mm
    )
    requested, iterations_by_scale, realized, clipped = physical_scales_to_iterations(
        tuple(scales_mm), median_step, max_iterations=max_iterations
    )

    normals0 = np.asarray(corrected_normals, dtype=np.float32)
    centers0 = np.asarray(centers_native, dtype=np.float32)
    sm_n = normals0.copy()
    sm_c = centers0.copy()

    normal_dev = {}
    signed_residual = {}
    abs_residual = {}

    targets = sorted(set(iterations_by_scale.values()))
    scales_at_target = {}
    for s, it in iterations_by_scale.items():
        scales_at_target.setdefault(it, []).append(s)

    current = 0
    for target in targets:
        while current < target:
            sm_n = diffusion.dot(sm_n)
            nn = np.linalg.norm(sm_n, axis=1)
            sm_n /= np.maximum(nn[:, None], 1e-12)
            sm_c = diffusion.dot(sm_c)
            current += 1

        dot = np.einsum("ij,ij->i", normals0, sm_n)
        ndev = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))).astype(np.float32)
        delta = centers0 - sm_c
        signed = (np.einsum("ij,ij->i", delta, sm_n) * float(scale_to_mm)).astype(np.float32)
        absolute = np.abs(signed).astype(np.float32)
        for scale in scales_at_target[target]:
            normal_dev[scale] = ndev.copy()
            signed_residual[scale] = signed.copy()
            abs_residual[scale] = absolute.copy()

    dot_edge = np.einsum(
        "ij,ij->i",
        corrected_normals[same_face_a],
        corrected_normals[same_face_b],
    )
    angle = np.degrees(np.arccos(np.clip(dot_edge, -1.0, 1.0))).astype(np.float32)
    dihedral_max = np.zeros(len(centers_native), dtype=np.float32)
    dihedral_sum = np.zeros(len(centers_native), dtype=np.float32)
    dihedral_count = np.zeros(len(centers_native), dtype=np.float32)
    np.maximum.at(dihedral_max, same_face_a, angle)
    np.maximum.at(dihedral_max, same_face_b, angle)
    np.add.at(dihedral_sum, same_face_a, angle)
    np.add.at(dihedral_sum, same_face_b, angle)
    np.add.at(dihedral_count, same_face_a, 1.0)
    np.add.at(dihedral_count, same_face_b, 1.0)
    dihedral_mean = dihedral_sum / np.maximum(dihedral_count, 1.0)

    return MultiscaleFeatures(
        requested_scales_mm=requested,
        iterations_by_scale=iterations_by_scale,
        realized_scales_mm=realized,
        clipped_scales=clipped,
        median_same_surface_face_step_mm=median_step,
        normal_dev_deg=normal_dev,
        signed_residual_mm=signed_residual,
        abs_residual_mm=abs_residual,
        dihedral_max_deg=dihedral_max,
        dihedral_mean_deg=dihedral_mean.astype(np.float32),
    )


def robust_unit_interval(values: np.ndarray, mask: np.ndarray, qlo=0.50, qhi=0.98):
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


def safe_quantiles(values: np.ndarray, probs=(0.05, 0.25, 0.5, 0.75, 0.90, 0.95, 0.98, 0.99)):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return None
    q = np.quantile(values, probs)
    return {f"q{int(round(p * 100)):02d}": float(v) for p, v in zip(probs, q)}


def write_binary_point_scalar_ply(
    path: Path,
    points: np.ndarray,
    fields: dict[str, np.ndarray],
    normals: np.ndarray | None = None,
):
    """Write CloudCompare-friendly face-center point PLY with scalar fields."""
    path.parent.mkdir(parents=True, exist_ok=True)
    points = np.asarray(points, dtype=np.float32)
    n = len(points)
    names = list(fields.keys())
    arrays = [np.asarray(fields[k], dtype=np.float32) for k in names]
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
        "ply",
        "format binary_little_endian 1.0",
        "comment PotterySurfaceTraceAnalyzer multiscale face-center features",
        f"element vertex {n}",
        "property float x",
        "property float y",
        "property float z",
    ]
    if normal_array is not None:
        header += ["property float nx", "property float ny", "property float nz"]
    header += [f"property float {name}" for name in names]
    header += ["end_header"]
    with path.open("wb") as f:
        f.write(("\n".join(header) + "\n").encode("ascii"))
        rec.tofile(f)


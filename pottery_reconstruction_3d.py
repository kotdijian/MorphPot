#!/usr/bin/env python3
"""
PotteryReconstruction3D
Version 0.2.4

Strict fixed-placement reconstruction from multiple pottery sherd meshes.

Core policy:
1. Sherds are spatially labelled, but their input XYZ coordinates and orientations
   are immutable observations.
2. One global axis of revolution is estimated from the sherd ensemble.
3. Inner / outer surface evidence and a longitudinal r(z) profile are estimated
   directly from the unchanged placement.
4. A lightweight axisymmetric shell is generated in the SAME coordinate system as
   the input placement.
5. Deviations are evaluated by comparing each unchanged sherd sample to the fitted
   vessel surface.
6. No per-sherd translation, rotation, Z shift, reversal, or pose refinement is
   performed anywhere in the reconstruction or export pipeline.

The GUI also preserves the current camera when switching from the initial sherd
view to the reconstruction/residual overlay, so a camera reset cannot be mistaken
for a pose change.

v0.2.4 adds rigid registration of externally reconstructed ring/arc PLY sets.
Only the external guide geometry is transformed: its fitted center axis is aligned
to the fixed sherd global axis, then its axial position is optimized against the
fixed-placement inner/outer r(z) profile. Sherd coordinates remain immutable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable

import numpy as np

__version__ = "0.2.5"

UNIT_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}
MM_TO_UNIT = {"mm": 1.0, "cm": 0.1, "m": 0.001}


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

@dataclass
class FragmentCluster:
    fragment_id: int
    sample_count: int
    centroid_mm: list[float]
    bounds_min_mm: list[float]
    bounds_max_mm: list[float]


@dataclass
class GlobalAxisFit:
    axis_point_mm: np.ndarray
    axis_direction: np.ndarray
    prior_center_mm: np.ndarray
    prior_direction: np.ndarray
    objective: float
    radial_dispersion_mm: float
    median_normal_residual: float
    usable_fraction: float
    axis_tilt_from_input_z_deg: float
    center_shift_from_prior_mm: float
    n_fit_samples: int


@dataclass
class ProfileModel:
    z_mm: np.ndarray
    inner_observed_mm: np.ndarray
    outer_observed_mm: np.ndarray
    inner_model_mm: np.ndarray
    outer_model_mm: np.ndarray
    thickness_model_mm: np.ndarray
    support_inner_points: np.ndarray
    support_outer_points: np.ndarray
    support_inner_fragments: np.ndarray
    support_outer_fragments: np.ndarray
    orientation_mapping: str


@dataclass
class ReconstructionResult:
    source_path: Path
    input_unit: str
    export_unit: str
    clusters: list[FragmentCluster]
    display_points_input: np.ndarray
    display_points_mm: np.ndarray
    display_labels: np.ndarray
    display_colors: np.ndarray | None
    face_centroids_mm: np.ndarray
    face_normals: np.ndarray
    face_labels: np.ndarray
    axis: GlobalAxisFit
    profile: ProfileModel
    mesh_vertices_mm: np.ndarray
    mesh_faces: np.ndarray
    point_residual_mm: np.ndarray
    point_nearest_surface: np.ndarray
    fragment_fit_quality: dict[int, dict]
    metadata: dict


@dataclass
class GuideRing:
    source_path: Path
    source_unit: str
    kind: str
    points_mm: np.ndarray
    ordered_points_mm: np.ndarray
    center_mm: np.ndarray
    normal: np.ndarray
    radius_mm: float
    plane_rms_mm: float
    circle_rms_mm: float
    arc_span_deg: float


@dataclass
class GuideRingSet:
    folder: Path
    rings: list[GuideRing]
    axis_point_mm: np.ndarray | None
    axis_direction: np.ndarray | None
    axis_line_rms_mm: float | None


@dataclass
class GuideRegistration:
    status: str
    axis_orientation_sign: int
    rotation_matrix: np.ndarray
    translation_mm: np.ndarray
    matrix_mm: np.ndarray
    angular_correction_deg: float
    lateral_shift_mm: float
    axial_shift_mm: float
    radius_rmse_mm: float
    median_abs_radius_mm: float
    p90_abs_radius_mm: float
    matched_rings: int
    total_rings: int
    overlap_z_mm: float
    guide_z_min_mm: float
    guide_z_max_mm: float
    target_z_min_mm: float
    target_z_max_mm: float
    bottom_mismatch_mm: float
    top_mismatch_mm: float


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


def normalize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if n < eps:
        raise ValueError("Cannot normalize near-zero vector")
    return v / n


def orthonormal_perp_basis(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a = normalize(a)
    ref = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.8 else np.array([0.0, 1.0, 0.0])
    e1 = normalize(np.cross(a, ref))
    e2 = normalize(np.cross(a, e1))
    return e1, e2


def sample_indices(n_total: int, n_target: int, seed: int = 0) -> np.ndarray:
    n_total = int(n_total)
    n_target = min(int(n_target), n_total)
    if n_target >= n_total:
        return np.arange(n_total, dtype=np.int64)
    rng = np.random.default_rng(seed)
    step = n_total / n_target
    idx = (np.arange(n_target, dtype=float) * step + rng.random(n_target) * step).astype(np.int64)
    return np.clip(idx, 0, n_total - 1)


def unit_scale(input_unit: str) -> float:
    if input_unit not in UNIT_TO_MM:
        raise ValueError(f"Unsupported unit: {input_unit}")
    return UNIT_TO_MM[input_unit]


def export_scale(export_unit: str) -> float:
    if export_unit not in MM_TO_UNIT:
        raise ValueError(f"Unsupported export unit: {export_unit}")
    return MM_TO_UNIT[export_unit]


def fragment_color(fid: int) -> np.ndarray:
    rng = np.random.default_rng(1000 + int(fid))
    return rng.integers(45, 230, size=3, dtype=np.uint8)


def _rotation_matrix_axis_angle(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    a = normalize(axis)
    x, y, z = a
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    C = 1.0 - c
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C],
    ], dtype=float)


# -----------------------------------------------------------------------------
# Fast PLY reader + generic trimesh fallback
# -----------------------------------------------------------------------------

PLY_DTYPE = {
    "char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1",
    "short": "<i2", "int16": "<i2", "ushort": "<u2", "uint16": "<u2",
    "int": "<i4", "int32": "<i4", "uint": "<u4", "uint32": "<u4",
    "float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8",
}


@dataclass
class PLYHeader:
    fmt: str
    header_bytes: int
    vertex_count: int
    vertex_properties: list[tuple[str, str]]
    face_count: int
    face_count_type: str | None
    face_index_type: str | None


def parse_ply_header(path: Path) -> PLYHeader:
    with Path(path).open("rb") as f:
        if f.readline().strip() != b"ply":
            raise ValueError("Not a PLY file")
        fmt = None
        vertex_count = face_count = 0
        vertex_properties: list[tuple[str, str]] = []
        current_element = None
        face_count_type = face_index_type = None
        while True:
            line_b = f.readline()
            if not line_b:
                raise ValueError("Unexpected EOF in PLY header")
            line = line_b.decode("ascii", "replace").strip()
            if line == "end_header":
                header_bytes = f.tell()
                break
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "format":
                fmt = parts[1]
            elif parts[0] == "element" and len(parts) >= 3:
                current_element = parts[1]
                if current_element == "vertex":
                    vertex_count = int(parts[2])
                elif current_element == "face":
                    face_count = int(parts[2])
            elif parts[0] == "property":
                if current_element == "vertex" and len(parts) == 3:
                    vertex_properties.append((parts[2], parts[1]))
                elif current_element == "face" and len(parts) >= 5 and parts[1] == "list":
                    face_count_type = parts[2]
                    face_index_type = parts[3]
        if fmt is None:
            raise ValueError("PLY format line missing")
        return PLYHeader(fmt, header_bytes, vertex_count, vertex_properties,
                         face_count, face_count_type, face_index_type)


class MeshSourceBase:
    path: Path
    n_vertices: int
    n_faces: int

    def sample_vertices(self, n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray | None]:
        raise NotImplementedError

    def sample_face_centroids_normals(self, n: int, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError


class FastBinaryPLY(MeshSourceBase):
    """Memory-mapped binary little-endian triangular PLY reader."""

    def __init__(self, path: Path):
        self.path = Path(path)
        h = parse_ply_header(self.path)
        self.header = h
        if h.fmt != "binary_little_endian":
            raise ValueError("Fast PLY path requires binary_little_endian")
        fields = []
        for name, typ in h.vertex_properties:
            if typ not in PLY_DTYPE:
                raise ValueError(f"Unsupported PLY vertex type: {typ}")
            fields.append((name, PLY_DTYPE[typ]))
        names = {x[0] for x in fields}
        if not {"x", "y", "z"}.issubset(names):
            raise ValueError("PLY vertex x/y/z properties are required")
        self.vertex_dtype = np.dtype(fields, align=False)
        self.vertices_mmapped = np.memmap(
            self.path, mode="r", dtype=self.vertex_dtype,
            offset=h.header_bytes, shape=(h.vertex_count,), order="C"
        )
        self.n_vertices = h.vertex_count
        self.n_faces = h.face_count
        self.face_offset = h.header_bytes + self.vertex_dtype.itemsize * h.vertex_count
        self._face_indices_view = None
        if h.face_count > 0 and h.face_count_type in {"uchar", "uint8"} and h.face_index_type in {"int", "int32"}:
            expected = h.face_count * 13
            actual = self.path.stat().st_size - self.face_offset
            if actual >= expected:
                raw = np.memmap(self.path, mode="r", dtype=np.uint8,
                                offset=self.face_offset, shape=(expected,))
                test_i = sample_indices(h.face_count, min(2048, h.face_count), seed=17)
                if np.all(raw[test_i * 13] == 3):
                    self._face_indices_view = np.ndarray(
                        shape=(h.face_count, 3), dtype="<i4", buffer=raw,
                        offset=1, strides=(13, 4)
                    )
        if self.n_faces and self._face_indices_view is None:
            raise ValueError("Fast PLY path currently requires fixed triangular face records")

    def _xyz(self, idx: np.ndarray) -> np.ndarray:
        v = self.vertices_mmapped[idx]
        return np.column_stack((v["x"], v["y"], v["z"])).astype(np.float64, copy=False)

    def sample_vertices(self, n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray | None]:
        idx = sample_indices(self.n_vertices, n, seed)
        pts = self._xyz(idx)
        names = self.vertices_mmapped.dtype.names or ()
        colors = None
        if {"red", "green", "blue"}.issubset(names):
            v = self.vertices_mmapped[idx]
            colors = np.column_stack((v["red"], v["green"], v["blue"])).astype(np.uint8)
        return pts, colors

    def sample_face_centroids_normals(self, n: int, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
        if self.n_faces <= 0 or self._face_indices_view is None:
            raise ValueError("No triangular faces available")
        fi = sample_indices(self.n_faces, n, seed)
        tri_idx = np.asarray(self._face_indices_view[fi], dtype=np.int64)
        p0 = self._xyz(tri_idx[:, 0])
        p1 = self._xyz(tri_idx[:, 1])
        p2 = self._xyz(tri_idx[:, 2])
        c = (p0 + p1 + p2) / 3.0
        nn = np.cross(p1 - p0, p2 - p0)
        lens = np.linalg.norm(nn, axis=1)
        good = lens > 1e-15
        return c[good], nn[good] / lens[good, None]


class TrimeshSource(MeshSourceBase):
    def __init__(self, path: Path):
        try:
            import trimesh
        except ModuleNotFoundError as exc:
            raise RuntimeError("Generic OBJ/PLY loading requires trimesh") from exc
        self.path = Path(path)
        loaded = trimesh.load(str(path), process=False)
        if isinstance(loaded, trimesh.Scene):
            if not loaded.geometry:
                raise ValueError("Mesh scene is empty")
            loaded = trimesh.util.concatenate(tuple(loaded.geometry.values()))
        if not isinstance(loaded, trimesh.Trimesh):
            raise ValueError("Input did not resolve to a triangle mesh")
        self.mesh = loaded
        self.n_vertices = len(loaded.vertices)
        self.n_faces = len(loaded.faces)

    def sample_vertices(self, n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray | None]:
        idx = sample_indices(self.n_vertices, n, seed)
        pts = np.asarray(self.mesh.vertices[idx], dtype=float)
        colors = None
        try:
            vc = np.asarray(self.mesh.visual.vertex_colors)
            if len(vc) == self.n_vertices and vc.ndim == 2 and vc.shape[1] >= 3:
                colors = vc[idx, :3].astype(np.uint8)
        except Exception:
            pass
        return pts, colors

    def sample_face_centroids_normals(self, n: int, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
        fi = sample_indices(self.n_faces, n, seed)
        f = np.asarray(self.mesh.faces[fi], dtype=np.int64)
        v = np.asarray(self.mesh.vertices, dtype=float)
        p0, p1, p2 = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
        c = (p0 + p1 + p2) / 3.0
        nn = np.cross(p1 - p0, p2 - p0)
        lens = np.linalg.norm(nn, axis=1)
        good = lens > 1e-15
        return c[good], nn[good] / lens[good, None]


def open_mesh_source(path: Path) -> MeshSourceBase:
    path = Path(path)
    if path.suffix.lower() == ".ply":
        try:
            return FastBinaryPLY(path)
        except Exception as fast_exc:
            print(f"[loader] fast PLY path unavailable: {fast_exc}", file=sys.stderr)
            print("[loader] falling back to trimesh", file=sys.stderr)
    return TrimeshSource(path)


# -----------------------------------------------------------------------------
# External reconstructed-ring PLY loading
# -----------------------------------------------------------------------------

def _read_ply_header_text(path: Path, max_bytes: int = 65536) -> str:
    """Read only the ASCII PLY header, regardless of payload encoding."""
    with Path(path).open("rb") as f:
        raw = f.read(max_bytes)
    pos = raw.find(b"end_header")
    if pos < 0:
        return raw.decode("ascii", "replace")
    end = raw.find(b"\n", pos)
    if end < 0:
        end = pos + len(b"end_header")
    return raw[:end + 1].decode("ascii", "replace")


def detect_ply_coordinate_unit(path: Path) -> str | None:
    """Best-effort unit detection from PLY comment lines."""
    header = _read_ply_header_text(path).lower()
    import re
    patterns = [
        r"coordinate_unit\s*=\s*(mm|cm|m)\b",
        r"coordinate[_ ]unit\s*[:=]\s*(mm|cm|m)\b",
        r"coordinates?\s+in\s+(millimetres?|millimeters?|centimetres?|centimeters?|metres?|meters?)",
    ]
    for pat in patterns:
        m = re.search(pat, header)
        if not m:
            continue
        v = m.group(1)
        if v in {"mm", "cm", "m"}:
            return v
        if v.startswith("milli"):
            return "mm"
        if v.startswith("centi"):
            return "cm"
        if v.startswith("met"):
            return "m"
    return None


def load_ply_vertices_all(path: Path) -> np.ndarray:
    """Load every PLY vertex as XYZ, supporting common ASCII/binary little-endian files."""
    path = Path(path)
    h = parse_ply_header(path)
    names = [name for name, _ in h.vertex_properties]
    if not {"x", "y", "z"}.issubset(names):
        raise ValueError(f"{path.name}: x/y/z vertex properties are required")
    ix, iy, iz = names.index("x"), names.index("y"), names.index("z")

    if h.fmt == "binary_little_endian":
        fields = []
        for name, typ in h.vertex_properties:
            if typ not in PLY_DTYPE:
                raise ValueError(f"{path.name}: unsupported PLY vertex type {typ}")
            fields.append((name, PLY_DTYPE[typ]))
        dt = np.dtype(fields, align=False)
        arr = np.memmap(path, mode="r", dtype=dt, offset=h.header_bytes,
                        shape=(h.vertex_count,), order="C")
        pts = np.column_stack((arr["x"], arr["y"], arr["z"])).astype(np.float64)
        return pts

    if h.fmt == "ascii":
        pts = np.empty((h.vertex_count, 3), dtype=np.float64)
        with path.open("rb") as f:
            f.seek(h.header_bytes)
            for i in range(h.vertex_count):
                line = f.readline()
                if not line:
                    raise ValueError(f"{path.name}: unexpected EOF in vertex list")
                vals = line.decode("ascii", "replace").strip().split()
                if len(vals) < len(h.vertex_properties):
                    raise ValueError(f"{path.name}: malformed ASCII vertex row {i}")
                pts[i] = [float(vals[ix]), float(vals[iy]), float(vals[iz])]
        return pts

    # Rare encodings fall back to trimesh if available.
    try:
        import trimesh
    except ModuleNotFoundError as exc:
        raise RuntimeError(f"{path.name}: unsupported PLY format {h.fmt}") from exc
    obj = trimesh.load(str(path), process=False)
    if isinstance(obj, trimesh.Scene):
        verts = [np.asarray(g.vertices, float) for g in obj.geometry.values() if hasattr(g, "vertices")]
        if not verts:
            raise ValueError(f"{path.name}: no vertices found")
        return np.vstack(verts)
    if hasattr(obj, "vertices"):
        return np.asarray(obj.vertices, dtype=float)
    raise ValueError(f"{path.name}: no vertices found")


def fit_circle_3d(points_mm: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float, float, np.ndarray]:
    """Fit a plane + circle to a full or partial 3-D ring/arc.

    Returns center, normal, radius, plane RMS, circle RMS, ordered points.
    """
    p = np.asarray(points_mm, dtype=float)
    p = p[np.all(np.isfinite(p), axis=1)]
    if len(p) < 6:
        raise ValueError("At least 6 finite points are required for a 3-D circle fit")

    g = np.mean(p, axis=0)
    q = p - g
    cov = np.cov(q.T)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)
    normal = normalize(evecs[:, order[0]])
    e1 = normalize(evecs[:, order[2]])
    e2 = normalize(np.cross(normal, e1))

    x = q @ e1
    y = q @ e2
    A = np.column_stack((2.0 * x, 2.0 * y, np.ones(len(p))))
    b = x * x + y * y
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, c0 = [float(v) for v in sol]
    r2 = c0 + cx * cx + cy * cy
    if not np.isfinite(r2) or r2 <= 0:
        raise ValueError("Degenerate circle fit")
    radius = math.sqrt(r2)
    center = g + cx * e1 + cy * e2

    plane_d = (p - g) @ normal
    plane_rms = float(np.sqrt(np.mean(plane_d * plane_d)))
    rr = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    circle_rms = float(np.sqrt(np.mean((rr - radius) ** 2)))

    angles = np.mod(np.arctan2(y - cy, x - cx), 2.0 * math.pi)
    idx = np.argsort(angles)
    sa = angles[idx]
    gaps = np.diff(np.r_[sa, sa[0] + 2.0 * math.pi])
    cut = int(np.argmax(gaps))
    # Start immediately after the largest empty angular gap. This prevents a
    # partial arc from being drawn with a spurious chord across the missing sector.
    idx_ordered = np.r_[idx[cut + 1:], idx[:cut + 1]]
    ordered = p[idx_ordered]
    return center, normal, float(radius), plane_rms, circle_rms, ordered


def _ring_arc_span_deg(points_mm: np.ndarray, center_mm: np.ndarray, normal: np.ndarray) -> float:
    p = np.asarray(points_mm, float)
    if len(p) < 3:
        return 0.0
    e1, e2 = orthonormal_perp_basis(normal)
    q = p - np.asarray(center_mm, float)
    ang = np.mod(np.arctan2(q @ e2, q @ e1), 2.0 * math.pi)
    sa = np.sort(ang)
    gaps = np.diff(np.r_[sa, sa[0] + 2.0 * math.pi])
    return float(math.degrees(max(0.0, 2.0 * math.pi - float(np.max(gaps)))))


def _guide_ring_kind(path: Path) -> str:
    s = path.stem.lower()
    if "outer" in s:
        return "outer"
    if "inner" in s:
        return "inner"
    return "unknown"


def load_guide_ring_folder(folder: Path, fallback_unit: str = "m") -> GuideRingSet:
    """Load a folder of reconstructed circle/arc PLY files without registration.

    Source coordinates are converted only to internal millimetres. No rigid transform
    is applied here; registration to the sherd coordinate system is a separate step.
    """
    folder = Path(folder)
    if fallback_unit not in UNIT_TO_MM:
        raise ValueError(f"Unsupported fallback unit: {fallback_unit}")
    all_ply = sorted(p for p in folder.rglob("*.ply") if p.is_file())
    if not all_ply:
        raise ValueError(f"No PLY files found under: {folder}")

    # Prefer individual ring files when a folder also contains combined clouds,
    # reconstructed shells, axes, etc.
    individual = [p for p in all_ply if p.stem.lower().startswith("ring_")]
    files = individual if individual else all_ply

    rings: list[GuideRing] = []
    errors: list[str] = []
    for path in files:
        try:
            pts = load_ply_vertices_all(path)
            if len(pts) < 6:
                raise ValueError("fewer than 6 vertices")
            unit = detect_ply_coordinate_unit(path) or fallback_unit
            pmm = pts * unit_scale(unit)
            center, normal, radius, plane_rms, circle_rms, ordered = fit_circle_3d(pmm)
            arc_span = _ring_arc_span_deg(pmm, center, normal)
            rings.append(GuideRing(
                source_path=path,
                source_unit=unit,
                kind=_guide_ring_kind(path),
                points_mm=pmm,
                ordered_points_mm=ordered,
                center_mm=center,
                normal=normal,
                radius_mm=radius,
                plane_rms_mm=plane_rms,
                circle_rms_mm=circle_rms,
                arc_span_deg=arc_span,
            ))
        except Exception as exc:
            errors.append(f"{path.name}: {exc}")

    if not rings:
        detail = "; ".join(errors[:5])
        raise ValueError(f"No usable guide rings could be loaded. {detail}")

    axis_point = axis_direction = None
    axis_rms = None
    if len(rings) >= 2:
        centers = np.vstack([r.center_mm for r in rings])
        g = centers.mean(axis=0)
        cov = np.cov((centers - g).T)
        _, ev = np.linalg.eigh(cov)
        a = normalize(ev[:, -1])
        # Make ring normals point roughly along the fitted ring-stack axis for
        # deterministic display/diagnostics; the points themselves are unchanged.
        mean_n = np.mean([r.normal for r in rings], axis=0)
        if np.dot(a, mean_n) < 0:
            a = -a
        d = centers - g
        radial = d - (d @ a)[:, None] * a[None, :]
        axis_point = g
        axis_direction = a
        axis_rms = float(np.sqrt(np.mean(np.sum(radial * radial, axis=1))))

    out = GuideRingSet(folder=folder, rings=rings,
                       axis_point_mm=axis_point, axis_direction=axis_direction,
                       axis_line_rms_mm=axis_rms)
    # Attach non-fatal parsing errors for callers that want to log them without
    # expanding the dataclass/API surface.
    setattr(out, "load_errors", errors)
    return out


# -----------------------------------------------------------------------------
# External guide-axis registration (GUIDE moves; sherds never move)
# -----------------------------------------------------------------------------

def rotation_matrix_from_vectors(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
    """Return a proper 3x3 rotation mapping v_from onto v_to."""
    u = normalize(np.asarray(v_from, float))
    v = normalize(np.asarray(v_to, float))
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    if c > 1.0 - 1e-12:
        return np.eye(3, dtype=float)
    if c < -1.0 + 1e-12:
        # 180-degree rotation around any stable unit vector perpendicular to u.
        e1, _ = orthonormal_perp_basis(u)
        # Rodrigues at pi: R = -I + 2 kk^T
        return -np.eye(3, dtype=float) + 2.0 * np.outer(e1, e1)
    k = np.cross(u, v)
    s = float(np.linalg.norm(k))
    k /= s
    K = np.array([[0.0, -k[2], k[1]],
                  [k[2], 0.0, -k[0]],
                  [-k[1], k[0], 0.0]], dtype=float)
    return np.eye(3) + s * K + (1.0 - c) * (K @ K)


def apply_homogeneous(points_mm: np.ndarray, matrix_mm: np.ndarray) -> np.ndarray:
    p = np.asarray(points_mm, float)
    M = np.asarray(matrix_mm, float)
    return p @ M[:3, :3].T + M[:3, 3]


def _guide_radius_residuals_at_shift(guide: GuideRingSet, target_profile: ProfileModel,
                                     R: np.ndarray, guide_axis_point_mm: np.ndarray,
                                     target_axis_point_mm: np.ndarray, target_axis_direction: np.ndarray,
                                     axial_shift_mm: float) -> tuple[np.ndarray, np.ndarray]:
    """Radius residuals and target-z values for guide rings at one axial shift."""
    a = normalize(target_axis_direction)
    residuals = []
    zvals = []
    for ring in guide.rings:
        rel = R @ (np.asarray(ring.center_mm, float) - guide_axis_point_mm)
        z = float(np.dot(rel, a) + axial_shift_mm)
        zz = np.array([z], dtype=float)
        ro = float(interp_profile(zz, target_profile.z_mm, target_profile.outer_model_mm)[0])
        ri = float(interp_profile(zz, target_profile.z_mm, target_profile.inner_model_mm)[0])
        if ring.kind == "outer":
            pred = ro
        elif ring.kind == "inner":
            pred = ri
        else:
            vals = [x for x in (ro, ri) if np.isfinite(x)]
            if not vals:
                pred = np.nan
            else:
                pred = min(vals, key=lambda x: abs(float(ring.radius_mm) - x))
        if np.isfinite(pred):
            residuals.append(float(ring.radius_mm) - float(pred))
            zvals.append(z)
    return np.asarray(residuals, float), np.asarray(zvals, float)


def _guide_shift_objective(guide: GuideRingSet, target_profile: ProfileModel, R: np.ndarray,
                           guide_axis_point_mm: np.ndarray, target_axis_point_mm: np.ndarray,
                           target_axis_direction: np.ndarray, axial_shift_mm: float) -> float:
    d, _ = _guide_radius_residuals_at_shift(guide, target_profile, R, guide_axis_point_mm,
                                            target_axis_point_mm, target_axis_direction, axial_shift_mm)
    n_total = max(1, len(guide.rings))
    if len(d) < max(2, min(4, n_total)):
        return 1e6 + 1000.0 * (max(2, min(4, n_total)) - len(d))
    # Robust pseudo-Huber loss with 2 mm transition, plus a mild overlap penalty.
    delta = 2.0
    loss = delta * delta * (np.sqrt(1.0 + (d / delta) ** 2) - 1.0)
    robust = float(math.sqrt(max(0.0, 2.0 * np.mean(loss))))
    overlap_penalty = 2.0 * (1.0 - len(d) / n_total)
    return robust + overlap_penalty


def register_guide_axis_to_reconstruction(guide: GuideRingSet, result: ReconstructionResult,
                                          shift_step_mm: float = 1.0) -> GuideRegistration:
    """Rigidly register external rings to the fixed sherd reconstruction.

    The sherd geometry is never transformed.  For both possible guide-axis signs,
    the guide axis is rotated onto the sherd global axis and translated laterally
    onto the same line.  The remaining one-dimensional translation along the axis
    is optimized against the fixed sherd inner/outer radius profile.
    """
    if guide.axis_point_mm is None or guide.axis_direction is None or len(guide.rings) < 2:
        raise ValueError("Guide registration requires at least two fitted rings and a ring-stack axis")
    if shift_step_mm <= 0:
        raise ValueError("shift_step_mm must be positive")

    a_t = normalize(result.axis.axis_direction)
    c_t = np.asarray(result.axis.axis_point_mm, float)
    c_g = np.asarray(guide.axis_point_mm, float)
    a_g0 = normalize(guide.axis_direction)

    finite_target = np.isfinite(result.profile.outer_model_mm) | np.isfinite(result.profile.inner_model_mm)
    if np.count_nonzero(finite_target) < 2:
        raise ValueError("Target sherd profile has too little finite inner/outer support")
    target_z = np.asarray(result.profile.z_mm[finite_target], float)
    tzmin, tzmax = float(np.min(target_z)), float(np.max(target_z))

    best = None
    for sign in (1, -1):
        a_g = a_g0 * float(sign)
        R = rotation_matrix_from_vectors(a_g, a_t)
        # Guide ring axial coordinates after direction alignment, before axial shift.
        gz = np.array([float(np.dot(R @ (r.center_mm - c_g), a_t)) for r in guide.rings], float)
        gzmin, gzmax = float(np.min(gz)), float(np.max(gz))

        # Search all shifts giving at least possible overlap, with a small margin.
        lo = tzmin - gzmax - 20.0
        hi = tzmax - gzmin + 20.0
        if hi <= lo:
            continue
        grid = np.arange(lo, hi + 0.5 * shift_step_mm, shift_step_mm)
        vals = np.array([_guide_shift_objective(guide, result.profile, R, c_g, c_t, a_t, float(sh))
                         for sh in grid], float)
        if not np.any(np.isfinite(vals)):
            continue
        ib = int(np.nanargmin(vals))
        sh_best = float(grid[ib])
        obj_best = float(vals[ib])
        # Refine locally without making scipy a hard dependency for this feature.
        try:
            from scipy.optimize import minimize_scalar
            rlo = max(lo, sh_best - 2.5 * shift_step_mm)
            rhi = min(hi, sh_best + 2.5 * shift_step_mm)
            opt = minimize_scalar(lambda x: _guide_shift_objective(guide, result.profile, R, c_g, c_t, a_t, float(x)),
                                  bounds=(rlo, rhi), method="bounded", options={"xatol": 1e-4})
            if opt.success and float(opt.fun) <= obj_best:
                sh_best = float(opt.x); obj_best = float(opt.fun)
        except Exception:
            pass
        d, zmatch = _guide_radius_residuals_at_shift(guide, result.profile, R, c_g, c_t, a_t, sh_best)
        if len(d) < 2:
            continue
        item = (obj_best, sign, R, sh_best, d, zmatch, gzmin, gzmax)
        if best is None or item[0] < best[0]:
            best = item

    if best is None:
        raise RuntimeError("Guide-axis alignment could not find a usable axial profile overlap")

    _, sign, R, sh, d, zmatch, gzmin, gzmax = best
    t = c_t + sh * a_t - R @ c_g
    M = np.eye(4, dtype=float)
    M[:3, :3] = R
    M[:3, 3] = t

    # Diagnostics.
    source_axis = a_g0 * float(sign)
    ang = math.degrees(math.acos(float(np.clip(np.dot(source_axis, a_t), -1.0, 1.0))))
    rg0 = R @ c_g
    dc = c_t - rg0
    lateral = float(np.linalg.norm(dc - np.dot(dc, a_t) * a_t))
    absd = np.abs(d)
    rmse = float(np.sqrt(np.mean(d * d)))
    med = float(np.median(absd))
    p90 = float(np.quantile(absd, 0.90))
    overlap_z = float(np.ptp(zmatch)) if len(zmatch) >= 2 else 0.0
    z_aligned = np.array([float(np.dot(R @ (r.center_mm - c_g), a_t) + sh) for r in guide.rings], float)
    gmin, gmax = float(np.min(z_aligned)), float(np.max(z_aligned))
    return GuideRegistration(
        status="ok", axis_orientation_sign=int(sign), rotation_matrix=R, translation_mm=t, matrix_mm=M,
        angular_correction_deg=float(ang), lateral_shift_mm=lateral, axial_shift_mm=float(sh),
        radius_rmse_mm=rmse, median_abs_radius_mm=med, p90_abs_radius_mm=p90,
        matched_rings=int(len(d)), total_rings=int(len(guide.rings)), overlap_z_mm=overlap_z,
        guide_z_min_mm=gmin, guide_z_max_mm=gmax, target_z_min_mm=tzmin, target_z_max_mm=tzmax,
        bottom_mismatch_mm=float(gmin - tzmin), top_mismatch_mm=float(gmax - tzmax),
    )


def export_registered_guide_rings(guide: GuideRingSet, reg: GuideRegistration, out: Path,
                                  export_unit_name: str = "m") -> None:
    """Export aligned guide rings. Only guide coordinates are transformed."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    ringdir = out / "aligned_rings"; ringdir.mkdir(exist_ok=True)
    s = export_scale(export_unit_name)
    allp, allc = [], []
    rows = []
    for i, ring in enumerate(guide.rings):
        p = apply_homogeneous(ring.points_mm, reg.matrix_mm)
        if ring.kind == "outer":
            col = np.array([30, 170, 240], np.uint8)
        elif ring.kind == "inner":
            col = np.array([245, 145, 20], np.uint8)
        else:
            col = np.array([190, 70, 200], np.uint8)
        cc = np.tile(col, (len(p), 1))
        allp.append(p); allc.append(cc)
        stem = f"ring_{i:04d}_{ring.kind}_{ring.source_path.stem}"
        write_ply_point_cloud(ringdir / f"{stem}.ply", p * s, cc, coordinate_unit=export_unit_name)
        center = apply_homogeneous(ring.center_mm[None, :], reg.matrix_mm)[0]
        rows.append({
            "ring_index": i, "source_file": ring.source_path.name, "kind": ring.kind,
            "radius_mm": ring.radius_mm, "arc_span_deg": ring.arc_span_deg,
            "circle_rms_mm": ring.circle_rms_mm, "center_x_mm": center[0],
            "center_y_mm": center[1], "center_z_mm": center[2],
        })
    if allp:
        write_ply_point_cloud(out / "aligned_guide_rings.ply", np.vstack(allp) * s, np.vstack(allc),
                              coordinate_unit=export_unit_name)
    with (out / "aligned_ring_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        fields = ["ring_index","source_file","kind","radius_mm","arc_span_deg","circle_rms_mm",
                  "center_x_mm","center_y_mm","center_z_mm"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    meta = {
        "program": "PotteryReconstruction3D", "version": __version__,
        "policy": "external guide transformed; sherd coordinates unchanged",
        "source_folder": str(guide.folder), "export_unit": export_unit_name,
        "registration": jsonable(asdict(reg)),
    }
    (out / "guide_registration.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    with (out / "guide_registration_transform.csv").open("w", newline="", encoding="utf-8-sig") as f:
        fields = ["status","axis_orientation_sign","angular_correction_deg","lateral_shift_mm","axial_shift_mm",
                  "radius_rmse_mm","median_abs_radius_mm","p90_abs_radius_mm","matched_rings","total_rings",
                  "overlap_z_mm","bottom_mismatch_mm","top_mismatch_mm"] + [f"m{r}{c}" for r in range(4) for c in range(4)]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        row = {k: getattr(reg, k) for k in fields if hasattr(reg, k)}
        for r in range(4):
            for c in range(4):
                row[f"m{r}{c}"] = float(reg.matrix_mm[r, c])
        w.writerow(row)


# -----------------------------------------------------------------------------
# Sherd clustering: labels only; coordinates are never transformed
# -----------------------------------------------------------------------------

def _neighbor_offsets(radius_vox: float) -> list[tuple[int, int, int]]:
    r = int(math.ceil(radius_vox))
    out = []
    for dx in range(-r, r + 1):
        for dy in range(-r, r + 1):
            for dz in range(-r, r + 1):
                if dx == dy == dz == 0:
                    continue
                if math.sqrt(dx*dx + dy*dy + dz*dz) <= radius_vox + 1e-9:
                    if (dx, dy, dz) > (0, 0, 0):
                        out.append((dx, dy, dz))
    return out


def cluster_sample_points(points_mm: np.ndarray, voxel_mm: float = 1.5,
                          link_radius_vox: float = 2.1,
                          min_sample_points: int = 300) -> tuple[np.ndarray, list[FragmentCluster]]:
    pts = np.asarray(points_mm, dtype=float)
    if len(pts) < 3:
        raise ValueError("Not enough points to cluster")
    if voxel_mm <= 0:
        raise ValueError("voxel_mm must be positive")
    pmin = pts.min(axis=0)
    q = np.floor((pts - pmin) / float(voxel_mm)).astype(np.int32)
    vox, inv = np.unique(q, axis=0, return_inverse=True)
    nv = len(vox)
    parent = np.arange(nv, dtype=np.int32)
    rank = np.zeros(nv, dtype=np.uint8)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    lookup = {tuple(v.tolist()): i for i, v in enumerate(vox)}
    offsets = _neighbor_offsets(link_radius_vox)
    for i, v in enumerate(vox):
        x, y, z = map(int, v)
        for dx, dy, dz in offsets:
            j = lookup.get((x + dx, y + dy, z + dz))
            if j is not None:
                union(i, j)

    roots = np.array([find(i) for i in range(nv)], dtype=np.int32)
    _, root_inv = np.unique(roots, return_inverse=True)
    point_comp0 = root_inv[inv]
    comp_counts = np.bincount(point_comp0)
    keep_old = np.where(comp_counts >= int(min_sample_points))[0]
    keep_old = keep_old[np.argsort(comp_counts[keep_old])[::-1]]
    old_to_new = {int(old): int(new) for new, old in enumerate(keep_old, start=1)}
    labels = np.array([old_to_new.get(int(c), 0) for c in point_comp0], dtype=np.int32)

    clusters = []
    for fid in range(1, len(keep_old) + 1):
        sel = labels == fid
        p = pts[sel]
        clusters.append(FragmentCluster(
            fragment_id=fid,
            sample_count=int(len(p)),
            centroid_mm=[float(x) for x in p.mean(axis=0)],
            bounds_min_mm=[float(x) for x in p.min(axis=0)],
            bounds_max_mm=[float(x) for x in p.max(axis=0)],
        ))
    return labels, clusters


def assign_face_labels(face_centroids_mm: np.ndarray,
                       display_points_mm: np.ndarray,
                       display_labels: np.ndarray,
                       max_distance_mm: float = 4.0) -> np.ndarray:
    """Assign sampled face centroids to the nearest labelled display point."""
    from scipy.spatial import cKDTree
    labelled = display_labels > 0
    if np.count_nonzero(labelled) == 0:
        return np.zeros(len(face_centroids_mm), dtype=np.int32)
    tree = cKDTree(display_points_mm[labelled])
    dist, idx = tree.query(face_centroids_mm, k=1, workers=-1)
    lab = display_labels[labelled][idx].astype(np.int32)
    lab[np.asarray(dist) > float(max_distance_mm)] = 0
    return lab


# -----------------------------------------------------------------------------
# Global axis geometry
# -----------------------------------------------------------------------------

def axis_coordinates(points_mm: np.ndarray, normals: np.ndarray | None,
                     axis_point_mm: np.ndarray, axis_direction: np.ndarray):
    p = np.asarray(points_mm, float)
    a = normalize(axis_direction)
    c = np.asarray(axis_point_mm, float)
    d = p - c
    z = d @ a
    radial = d - z[:, None] * a[None, :]
    rho = np.linalg.norm(radial, axis=1)
    er = np.zeros_like(radial)
    good = rho > 1e-9
    er[good] = radial[good] / rho[good, None]
    out = {"z": z, "rho": rho, "radial_unit": er}
    if normals is not None:
        n = np.asarray(normals, float)
        nl = np.linalg.norm(n, axis=1)
        ngood = nl > 1e-12
        nn = np.zeros_like(n)
        nn[ngood] = n[ngood] / nl[ngood, None]
        tang = np.cross(np.broadcast_to(a, er.shape), er)
        tl = np.linalg.norm(tang, axis=1)
        tg = tl > 1e-12
        tang[tg] /= tl[tg, None]
        compatibility = np.ones(len(p), float)
        compatibility[good & ngood & tg] = np.abs(np.sum(nn[good & ngood & tg] * tang[good & ngood & tg], axis=1))
        radial_sign = np.sum(nn * er, axis=1)
        out["compatibility"] = compatibility
        out["radial_normal_sign"] = radial_sign
    return out


def _trimmed_mad(x: np.ndarray, qlo: float = 0.10, qhi: float = 0.90) -> tuple[float, float] | None:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return None
    lo, hi = np.quantile(x, [qlo, qhi])
    q = x[(x >= lo) & (x <= hi)]
    if len(q) < 6:
        return None
    med = float(np.median(q))
    mad = float(np.median(np.abs(q - med)) * 1.4826)
    return med, mad


def _balanced_axis_fit_subset(points: np.ndarray, normals: np.ndarray, labels: np.ndarray,
                              max_per_fragment: int = 4500, max_total: int = 30000,
                              seed: int = 101) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    pieces = []
    for fid in sorted(int(x) for x in np.unique(labels) if x > 0):
        idx = np.where(labels == fid)[0]
        if len(idx) > max_per_fragment:
            idx = rng.choice(idx, size=max_per_fragment, replace=False)
        pieces.append(idx)
    if not pieces:
        idx = sample_indices(len(points), min(max_total, len(points)), seed)
    else:
        idx = np.concatenate(pieces)
        if len(idx) > max_total:
            idx = rng.choice(idx, size=max_total, replace=False)
    return points[idx], normals[idx], labels[idx]


def _axis_profile_dispersion(points: np.ndarray, normals: np.ndarray,
                             c: np.ndarray, a: np.ndarray,
                             objective_bin_mm: float = 3.0,
                             normal_max: float = 0.30,
                             radial_sign_min: float = 0.18,
                             min_bin_samples: int = 18) -> tuple[float, float, float, int]:
    g = axis_coordinates(points, normals, c, a)
    z, rho = g["z"], g["rho"]
    comp, sign = g["compatibility"], g["radial_normal_sign"]
    usable = np.isfinite(z) & np.isfinite(rho) & (rho > 1.0) & (comp <= normal_max) & (np.abs(sign) >= radial_sign_min)
    if np.count_nonzero(usable) < max(250, min_bin_samples * 5):
        return 1e3, 1.0, 0.0, 0
    zu, ru, su, cu = z[usable], rho[usable], sign[usable], comp[usable]
    z0 = float(np.min(zu))
    bins = np.floor((zu - z0) / float(objective_bin_mm)).astype(np.int32)
    dispersions = []
    n_bins = 0
    for b in np.unique(bins):
        sb = bins == b
        for positive in (True, False):
            q = ru[sb & ((su > 0) if positive else (su < 0))]
            if len(q) < min_bin_samples:
                continue
            fit = _trimmed_mad(q)
            if fit is not None:
                dispersions.append(fit[1])
                n_bins += 1
    if len(dispersions) < 5:
        return 1e3, float(np.median(cu)), float(np.mean(usable)), n_bins
    return float(np.median(dispersions)), float(np.median(cu)), float(np.mean(usable)), n_bins


def _solve_axis_center_from_normals(points: np.ndarray, normals: np.ndarray,
                                    a: np.ndarray) -> tuple[np.ndarray, float]:
    """Robustly intersect projected surface-normal lines to locate the axis line."""
    a = normalize(a)
    P = np.asarray(points, float); N = np.asarray(normals, float)
    g = np.mean(P, axis=0)
    e1, e2 = orthonormal_perp_basis(a)
    ndota = N @ a
    nperp = N - ndota[:, None] * a[None, :]
    nplen = np.linalg.norm(nperp, axis=1)
    valid = nplen > 0.15
    if np.count_nonzero(valid) < 30:
        return g, 1e6
    pv = P[valid]
    dv = nperp[valid] / nplen[valid, None]
    q = pv - g
    q2 = np.column_stack((q @ e1, q @ e2))
    d2 = np.column_stack((dv @ e1, dv @ e2))
    dl = np.linalg.norm(d2, axis=1)
    good = dl > 1e-9
    q2 = q2[good]; d2 = d2[good] / dl[good, None]
    if len(q2) < 30:
        return g, 1e6
    # Each projected normal gives a 2-D line q + t*d.  The axis centre is the
    # robust intersection of those lines.
    M = np.column_stack((-d2[:, 1], d2[:, 0]))
    b = np.sum(M * q2, axis=1)
    w = np.ones(len(M), float); c2 = np.zeros(2, float)
    for _ in range(7):
        Aw = M * np.sqrt(w[:, None]); bw = b * np.sqrt(w)
        c2, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
        rr = M @ c2 - b
        scale = np.median(np.abs(rr)) * 1.4826 + 1e-6
        u = rr / (2.5 * scale)
        w = 1.0 / (1.0 + u*u)
    rr = M @ c2 - b
    line_scale = float(np.median(np.abs(rr)) * 1.4826)
    c = g + c2[0]*e1 + c2[1]*e2
    return c, line_scale


def _surface_like_mask_by_fragment(normals: np.ndarray, labels: np.ndarray,
                                   keep_fraction: float = 0.70,
                                   min_abs_dot: float = 0.55) -> tuple[np.ndarray, dict[int, dict]]:
    """
    Select likely vessel-wall faces before global-axis fitting.

    A pottery sherd is a thin shell. Outer and inner vessel surfaces therefore
    dominate two approximately antipodal normal populations, whereas fracture
    faces tend to be more nearly orthogonal and more dispersed.

    This selection is used ONLY for axis estimation. No sherd coordinate or pose
    is changed.
    """
    N = np.asarray(normals, float)
    L = np.asarray(labels, int)
    keep_fraction = float(np.clip(keep_fraction, 0.35, 0.95))
    min_abs_dot = float(np.clip(min_abs_dot, 0.0, 0.95))
    mask = np.zeros(len(N), dtype=bool)
    diag: dict[int, dict] = {}
    for fid in sorted(int(x) for x in np.unique(L) if int(x) > 0):
        idx = np.where(L == fid)[0]
        if len(idx) < 30:
            continue
        ns = N[idx]
        nl = np.linalg.norm(ns, axis=1)
        good = nl > 1e-9
        if np.count_nonzero(good) < 20:
            continue
        nsu = ns[good] / nl[good, None]
        C = (nsu.T @ nsu) / max(len(nsu), 1)
        vals, vecs = np.linalg.eigh(C)
        n0 = normalize(vecs[:, int(np.argmax(vals))])
        dots = np.abs(nsu @ n0)
        q = max(0.0, min(1.0, 1.0 - keep_fraction))
        threshold = max(min_abs_dot, float(np.quantile(dots, q)))
        selected_local = dots >= threshold
        good_idx = idx[good]
        mask[good_idx[selected_local]] = True
        diag[fid] = {
            "n_total": int(len(idx)),
            "n_surface_like": int(np.count_nonzero(selected_local)),
            "surface_like_fraction": float(np.mean(selected_local)),
            "normal_threshold_abs_dot": float(threshold),
            "dominant_unsigned_normal": n0.tolist(),
            "normal_eigenvalues": vals.tolist(),
        }
    return mask, diag


def _fibonacci_axis_directions(n: int) -> list[np.ndarray]:
    """Approximately uniform undirected axis directions over one hemisphere."""
    n = max(48, int(n))
    golden = math.pi * (3.0 - math.sqrt(5.0))
    out: list[np.ndarray] = []
    for i in range(n):
        z = (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - z*z))
        th = golden * i
        out.append(np.array([r*math.cos(th), r*math.sin(th), z], float))
    out.extend([
        np.array([1.0, 0.0, 0.0], float),
        np.array([0.0, 1.0, 0.0], float),
        np.array([0.0, 0.0, 1.0], float),
    ])
    return out


def _trimmed_rms_fraction(values: np.ndarray, keep_fraction: float = 0.65) -> float:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return 1e6
    keep_fraction = float(np.clip(keep_fraction, 0.1, 1.0))
    k = max(8, int(round(len(x) * keep_fraction)))
    k = min(k, len(x))
    y = np.partition(x, k - 1)[:k]
    return float(np.sqrt(np.mean(y*y)))


def _unrestricted_axis_candidate_score(points: np.ndarray, normals: np.ndarray,
                                       labels: np.ndarray, direction: np.ndarray,
                                       normal_max: float = 0.30) -> tuple[float, np.ndarray, dict]:
    """
    Geometry-only score for a global axis candidate.

    No input +Z prior and no bounding-box-centre prior are used.
    """
    a = normalize(direction)
    if a[2] < -1e-12 or (abs(a[2]) <= 1e-12 and (a[1] < -1e-12 or (abs(a[1]) <= 1e-12 and a[0] < 0))):
        a = -a

    c, line_scale = _solve_axis_center_from_normals(points, normals, a)
    g = axis_coordinates(points, normals, c, a)
    comp = np.asarray(g["compatibility"], float)
    rho = np.asarray(g["rho"], float)

    global_trimmed = _trimmed_rms_fraction(comp, 0.65)
    per_fragment = []
    compatible_fractions = []
    for fid in sorted(int(x) for x in np.unique(labels) if int(x) > 0):
        sel = labels == fid
        if np.count_nonzero(sel) < 20:
            continue
        q = comp[sel]
        per_fragment.append(_trimmed_rms_fraction(q, 0.65))
        compatible_fractions.append(float(np.mean(q <= normal_max)))

    if not per_fragment:
        return 1e6, c, {
            "global_trimmed_normal_rms": 1e6,
            "fragment_median_trimmed_normal_rms": 1e6,
            "median_compatible_fraction": 0.0,
            "line_scale_mm": float(line_scale),
            "median_radius_mm": 0.0,
        }

    frag_med = float(np.median(per_fragment))
    frac_med = float(np.median(compatible_fractions))
    rmed = max(float(np.median(rho[np.isfinite(rho)])), 1.0)
    line_norm = min(float(line_scale) / rmed, 3.0)

    score = (
        global_trimmed
        + 0.90 * frag_med
        + 0.30 * line_norm
        + 1.50 * max(0.0, 0.30 - frac_med)
    )
    return float(score), c, {
        "global_trimmed_normal_rms": float(global_trimmed),
        "fragment_median_trimmed_normal_rms": float(frag_med),
        "median_compatible_fraction": float(frac_med),
        "line_scale_mm": float(line_scale),
        "median_radius_mm": float(rmed),
        "per_fragment_trimmed_normal_rms": [float(x) for x in per_fragment],
        "per_fragment_compatible_fraction": [float(x) for x in compatible_fractions],
    }


def estimate_global_axis(points_mm: np.ndarray, normals: np.ndarray, face_labels: np.ndarray,
                         placement_points_mm: np.ndarray,
                         max_center_shift_mm: float = 20.0,
                         max_tilt_deg: float = 12.0,
                         objective_bin_mm: float = 3.0,
                         normal_max: float = 0.30,
                         radial_sign_min: float = 0.18,
                         max_fit_samples_per_fragment: int = 4500,
                         max_fit_samples_total: int = 30000,
                         seed: int = 1234,
                         search_directions: int = 480,
                         surface_keep_fraction: float = 0.70,
                         progress: Callable[[str], None] | None = None) -> GlobalAxisFit:
    """
    Estimate ONE global rotation axis in the current input coordinate system.

    v0.2.2:
    - sherd XYZ and pose remain immutable;
    - no input +Z direction prior is imposed;
    - the axis may lie outside the sherd bounding box;
    - fracture influence is reduced by per-sherd vessel-surface-like normals;
    - unrestricted hemispherical search is followed by local refinement.

    max_center_shift_mm and max_tilt_deg are accepted only for v0.2.1 API
    compatibility and are intentionally ignored.
    """
    from scipy.optimize import minimize

    P0 = np.asarray(points_mm, float)
    N0 = np.asarray(normals, float)
    L0 = np.asarray(face_labels, int)

    surface_mask, surface_diag = _surface_like_mask_by_fragment(
        N0, L0, keep_fraction=surface_keep_fraction, min_abs_dot=0.55
    )
    if np.count_nonzero(surface_mask) < 300:
        raise ValueError("Too few surface-like face samples for unrestricted global-axis estimation")

    P, N, L = _balanced_axis_fit_subset(
        P0[surface_mask], N0[surface_mask], L0[surface_mask],
        max_per_fragment=max_fit_samples_per_fragment,
        max_total=max_fit_samples_total, seed=seed
    )
    if len(P) < 300:
        raise ValueError("Too few balanced surface-like samples for global-axis estimation")

    pp = np.asarray(placement_points_mm, float)
    if len(pp) < 10:
        raise ValueError("Not enough retained placement points for axis estimation")
    bbox_center = (np.min(pp, axis=0) + np.max(pp, axis=0)) / 2.0

    if progress:
        progress(
            f"Global axis v0.2.2: unrestricted search; "
            f"surface-like samples={np.count_nonzero(surface_mask):,}, balanced fit={len(P):,}, "
            f"directions={int(search_directions)}"
        )
        progress("Global axis v0.2.2: input +Z is diagnostic only; no tilt/centre prior is applied")

    coarse = []
    for a in _fibonacci_axis_directions(search_directions):
        sc, c, st = _unrestricted_axis_candidate_score(P, N, L, a, normal_max=normal_max)
        if np.isfinite(sc):
            coarse.append((float(sc), normalize(a), np.asarray(c, float), st))
    if not coarse:
        raise RuntimeError("Global axis search produced no finite candidate")
    coarse.sort(key=lambda x: x[0])

    seeds: list[np.ndarray] = []
    for _, a, _, _ in coarse:
        if all(abs(float(np.dot(a, b))) < 0.985 for b in seeds):
            seeds.append(a)
        if len(seeds) >= 6:
            break
    if not seeds:
        seeds = [coarse[0][1]]

    refined = []
    for a0 in seeds:
        theta0 = math.acos(float(np.clip(a0[2], -1.0, 1.0)))
        phi0 = math.atan2(float(a0[1]), float(a0[0]))

        def objective_ang(x):
            theta = float(np.clip(x[0], 1e-6, math.pi/2))
            phi = float(x[1])
            a = np.array([
                math.sin(theta) * math.cos(phi),
                math.sin(theta) * math.sin(phi),
                math.cos(theta),
            ], float)
            sc, _, _ = _unrestricted_axis_candidate_score(P, N, L, a, normal_max=normal_max)
            return sc

        opt = minimize(
            objective_ang, np.array([theta0, phi0], float), method="Nelder-Mead",
            options={"maxiter": 120, "xatol": 1e-5, "fatol": 1e-5, "disp": False}
        )
        theta = float(np.clip(opt.x[0], 1e-6, math.pi/2))
        phi = float(opt.x[1])
        a = normalize(np.array([
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ], float))
        sc, c, st = _unrestricted_axis_candidate_score(P, N, L, a, normal_max=normal_max)
        refined.append((float(sc), a, np.asarray(c, float), st))

    candidates = coarse[:8] + refined
    candidates.sort(key=lambda x: x[0])
    fbest, a, c, best_stats = candidates[0]
    if a[2] < 0:
        a = -a

    dispersion, normal_med, usable_frac, _ = _axis_profile_dispersion(
        P, N, c, a, objective_bin_mm=objective_bin_mm,
        normal_max=normal_max, radial_sign_min=radial_sign_min
    )
    tilt = math.degrees(math.acos(float(np.clip(abs(a[2]), -1.0, 1.0))))
    d = bbox_center - c
    centre_offset = float(np.linalg.norm(d - np.dot(d, a) * a))

    if progress:
        progress(
            "Global axis v0.2.2: "
            f"centre={c.round(3).tolist()} mm; dir={a.round(6).tolist()}; "
            f"angle from input Z={tilt:.3f} deg"
        )
        progress(
            "Global axis v0.2.2: "
            f"geometric score={fbest:.4f}; radial dispersion={dispersion:.3f} mm; "
            f"normal residual median={normal_med:.4f}; usable={usable_frac:.3f}; "
            f"normal-line scale={best_stats['line_scale_mm']:.3f} mm"
        )

    return GlobalAxisFit(
        axis_point_mm=c,
        axis_direction=a,
        prior_center_mm=bbox_center,
        prior_direction=np.array([np.nan, np.nan, np.nan], float),
        objective=float(fbest),
        radial_dispersion_mm=float(dispersion),
        median_normal_residual=float(normal_med),
        usable_fraction=float(usable_frac),
        axis_tilt_from_input_z_deg=float(tilt),
        center_shift_from_prior_mm=float(centre_offset),
        n_fit_samples=int(len(P)),
    )


# -----------------------------------------------------------------------------
# Global inner / outer profile from current placement
# -----------------------------------------------------------------------------

def _finite_runs(mask: np.ndarray) -> list[np.ndarray]:
    mask = np.asarray(mask, bool)
    runs = []
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        if start is not None and (not v or i == len(mask) - 1):
            end = i if (v and i == len(mask) - 1) else i - 1
            if end >= start:
                runs.append(np.arange(start, end + 1, dtype=int))
            start = None
    return runs


def _fill_short_gaps(z: np.ndarray, r: np.ndarray, max_gap_mm: float) -> np.ndarray:
    out = np.asarray(r, float).copy()
    finite = np.isfinite(out)
    idx = np.where(finite)[0]
    if len(idx) < 2:
        return out
    for i0, i1 in zip(idx[:-1], idx[1:]):
        if i1 <= i0 + 1:
            continue
        if z[i1] - z[i0] <= max_gap_mm + 1e-9:
            ii = np.arange(i0 + 1, i1)
            out[ii] = np.interp(z[ii], [z[i0], z[i1]], [out[i0], out[i1]])
    return out


def _smooth_finite_runs(r: np.ndarray, window_bins: int = 3) -> np.ndarray:
    out = np.asarray(r, float).copy()
    if window_bins <= 1:
        return out
    half = window_bins // 2
    finite = np.isfinite(out)
    for run in _finite_runs(finite):
        if len(run) < 3:
            continue
        src = out.copy()
        for i in run:
            lo = max(run[0], i - half)
            hi = min(run[-1], i + half)
            q = src[lo:hi+1]
            q = q[np.isfinite(q)]
            if len(q):
                out[i] = float(np.median(q))
    return out


def build_global_profile(face_centroids_mm: np.ndarray, face_normals: np.ndarray,
                         face_labels: np.ndarray, axis: GlobalAxisFit,
                         bin_mm: float = 2.0,
                         compatibility_max: float = 0.30,
                         radial_sign_min: float = 0.18,
                         min_bin_samples: int = 18,
                         thickness_min_mm: float = 1.0,
                         thickness_max_mm: float = 20.0,
                         max_interpolation_gap_mm: float = 8.0,
                         smooth_window_bins: int = 3) -> ProfileModel:
    g = axis_coordinates(face_centroids_mm, face_normals, axis.axis_point_mm, axis.axis_direction)
    z, rho = g["z"], g["rho"]
    comp, sign = g["compatibility"], g["radial_normal_sign"]
    good = (face_labels > 0) & np.isfinite(z) & np.isfinite(rho) & (rho > 0.5) & (comp <= compatibility_max) & (np.abs(sign) >= radial_sign_min)
    if np.count_nonzero(good) < 200:
        raise ValueError("Too few global revolution-compatible surface samples")

    zmin = math.floor(float(np.min(z[good])) / bin_mm) * bin_mm
    zmax = math.ceil(float(np.max(z[good])) / bin_mm) * bin_mm
    centers = np.arange(zmin + 0.5 * bin_mm, zmax, bin_mm)
    pos = np.full(len(centers), np.nan)
    neg = np.full(len(centers), np.nan)
    pos_n = np.zeros(len(centers), int)
    neg_n = np.zeros(len(centers), int)
    pos_nf = np.zeros(len(centers), int)
    neg_nf = np.zeros(len(centers), int)

    for i, zc in enumerate(centers):
        selz = good & (z >= zc - 0.5*bin_mm) & (z < zc + 0.5*bin_mm)
        for is_pos in (True, False):
            sel = selz & ((sign > 0) if is_pos else (sign < 0))
            q = rho[sel]
            if len(q) < min_bin_samples:
                continue
            fit = _trimmed_mad(q)
            if fit is None:
                continue
            med, _ = fit
            labs = np.unique(face_labels[sel])
            labs = labs[labs > 0]
            if is_pos:
                pos[i] = med; pos_n[i] = len(q); pos_nf[i] = len(labs)
            else:
                neg[i] = med; neg_n[i] = len(q); neg_nf[i] = len(labs)

    both = np.isfinite(pos) & np.isfinite(neg)
    if np.count_nonzero(both) >= 2:
        if float(np.nanmedian(pos[both] - neg[both])) >= 0:
            outer_obs, inner_obs = pos.copy(), neg.copy()
            outer_n, inner_n = pos_n.copy(), neg_n.copy()
            outer_nf, inner_nf = pos_nf.copy(), neg_nf.copy()
            mapping = "positive radial normal = outer; negative = inner"
        else:
            outer_obs, inner_obs = neg.copy(), pos.copy()
            outer_n, inner_n = neg_n.copy(), pos_n.copy()
            outer_nf, inner_nf = neg_nf.copy(), pos_nf.copy()
            mapping = "negative radial normal = outer; positive = inner (global orientation reversed)"
    else:
        # Fallback: choose the larger global median as outer.
        mp = np.nanmedian(pos) if np.any(np.isfinite(pos)) else np.nan
        mn = np.nanmedian(neg) if np.any(np.isfinite(neg)) else np.nan
        if np.isfinite(mp) and np.isfinite(mn) and mp >= mn:
            outer_obs, inner_obs = pos.copy(), neg.copy(); outer_n, inner_n = pos_n.copy(), neg_n.copy(); outer_nf, inner_nf = pos_nf.copy(), neg_nf.copy()
            mapping = "positive radial normal = outer (global-median fallback)"
        else:
            outer_obs, inner_obs = neg.copy(), pos.copy(); outer_n, inner_n = neg_n.copy(), pos_n.copy(); outer_nf, inner_nf = neg_nf.copy(), pos_nf.copy()
            mapping = "negative radial normal = outer (global-median fallback)"

    outer_model = _fill_short_gaps(centers, outer_obs, max_interpolation_gap_mm)
    inner_model = _fill_short_gaps(centers, inner_obs, max_interpolation_gap_mm)
    outer_model = _smooth_finite_runs(outer_model, smooth_window_bins)
    inner_model = _smooth_finite_runs(inner_model, smooth_window_bins)

    thickness = outer_model - inner_model
    bad = np.isfinite(thickness) & ((thickness < thickness_min_mm) | (thickness > thickness_max_mm))
    # Preserve the better-supported side and invalidate the less-supported one.
    for i in np.where(bad)[0]:
        if outer_n[i] >= inner_n[i]:
            inner_model[i] = np.nan
        else:
            outer_model[i] = np.nan
    thickness = outer_model - inner_model
    thickness[~(np.isfinite(outer_model) & np.isfinite(inner_model))] = np.nan

    if np.count_nonzero(np.isfinite(outer_model)) < 3 and np.count_nonzero(np.isfinite(inner_model)) < 3:
        raise ValueError("Global profile contains fewer than 3 usable bins")

    return ProfileModel(
        z_mm=centers,
        inner_observed_mm=inner_obs, outer_observed_mm=outer_obs,
        inner_model_mm=inner_model, outer_model_mm=outer_model,
        thickness_model_mm=thickness,
        support_inner_points=inner_n, support_outer_points=outer_n,
        support_inner_fragments=inner_nf, support_outer_fragments=outer_nf,
        orientation_mapping=mapping,
    )


# -----------------------------------------------------------------------------
# Profile interpolation, reconstruction mesh, residuals
# -----------------------------------------------------------------------------

def interp_profile(zq: np.ndarray, z: np.ndarray, r: np.ndarray,
                   max_gap_mm: float = 8.0) -> np.ndarray:
    zq = np.asarray(zq, float)
    out = np.full(len(zq), np.nan, float)
    z = np.asarray(z, float); r = np.asarray(r, float)
    good = np.isfinite(z) & np.isfinite(r)
    for run in _finite_runs(good):
        zz, rr = z[run], r[run]
        if len(zz) >= 2:
            sel = (zq >= zz[0]) & (zq <= zz[-1])
            out[sel] = np.interp(zq[sel], zz, rr)
        elif len(zz) == 1:
            sel = np.abs(zq - zz[0]) <= max_gap_mm * 0.25
            out[sel] = rr[0]
    return out


def build_reconstructed_mesh(profile: ProfileModel, axis: GlobalAxisFit,
                             angle_step_deg: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    a = normalize(axis.axis_direction)
    e1, e2 = orthonormal_perp_basis(a)
    c = axis.axis_point_mm
    ntheta = max(12, int(round(360.0 / float(angle_step_deg))))
    theta = np.linspace(0.0, 2.0*math.pi, ntheta, endpoint=False)

    all_v, all_f = [], []
    offset = 0
    for r, inward in ((profile.outer_model_mm, False), (profile.inner_model_mm, True)):
        good = np.isfinite(r) & np.isfinite(profile.z_mm)
        for run in _finite_runs(good):
            if len(run) < 2:
                continue
            V = np.empty((len(run)*ntheta, 3), float)
            for iz, idx in enumerate(run):
                rad = float(r[idx]); zz = float(profile.z_mm[idx])
                sl = slice(iz*ntheta, (iz+1)*ntheta)
                V[sl] = (c[None, :] + zz*a[None, :] +
                         rad*np.cos(theta)[:, None]*e1[None, :] +
                         rad*np.sin(theta)[:, None]*e2[None, :])
            F = []
            for iz in range(len(run)-1):
                a0, a1 = iz*ntheta, (iz+1)*ntheta
                for j in range(ntheta):
                    jn = (j+1) % ntheta
                    if not inward:
                        F.append((a0+j, a0+jn, a1+jn)); F.append((a0+j, a1+jn, a1+j))
                    else:
                        F.append((a0+j, a1+jn, a0+jn)); F.append((a0+j, a1+j, a1+jn))
            all_v.append(V); all_f.append(np.asarray(F, np.int32)+offset); offset += len(V)
    if not all_v:
        return np.empty((0,3), float), np.empty((0,3), np.int32)
    return np.vstack(all_v), np.vstack(all_f) if all_f else np.empty((0,3), np.int32)


def evaluate_point_residuals(points_mm: np.ndarray, axis: GlobalAxisFit,
                             profile: ProfileModel) -> tuple[np.ndarray, np.ndarray]:
    g = axis_coordinates(points_mm, None, axis.axis_point_mm, axis.axis_direction)
    z, rho = g["z"], g["rho"]
    ro = interp_profile(z, profile.z_mm, profile.outer_model_mm)
    ri = interp_profile(z, profile.z_mm, profile.inner_model_mm)
    do = np.abs(rho - ro); di = np.abs(rho - ri)
    residual = np.full(len(points_mm), np.nan, float)
    surface = np.zeros(len(points_mm), np.uint8)  # 1 outer, 2 inner
    only_o = np.isfinite(do) & ~np.isfinite(di)
    only_i = np.isfinite(di) & ~np.isfinite(do)
    both = np.isfinite(do) & np.isfinite(di)
    residual[only_o] = do[only_o]; surface[only_o] = 1
    residual[only_i] = di[only_i]; surface[only_i] = 2
    choose_o = both & (do <= di); choose_i = both & (di < do)
    residual[choose_o] = do[choose_o]; surface[choose_o] = 1
    residual[choose_i] = di[choose_i]; surface[choose_i] = 2
    return residual, surface


def fragment_fit_quality(points_mm: np.ndarray, labels: np.ndarray,
                         residual_mm: np.ndarray) -> dict[int, dict]:
    out = {}
    for fid in sorted(int(x) for x in np.unique(labels) if x > 0):
        d = residual_mm[labels == fid]
        d = d[np.isfinite(d)]
        if len(d) < 20:
            out[fid] = {"status": "insufficient_overlap", "n": int(len(d))}
            continue
        # Trim the largest 10% because display vertices also contain fracture surfaces.
        hi = np.quantile(d, 0.90)
        q = d[d <= hi]
        out[fid] = {
            "status": "ok", "n": int(len(d)), "n_trimmed": int(len(q)),
            "rmse_mm": float(np.sqrt(np.mean(q*q))),
            "median_abs_mm": float(np.median(q)),
            "p90_abs_mm": float(np.quantile(q, 0.90)),
            "mean_abs_mm": float(np.mean(q)),
        }
    return out


# -----------------------------------------------------------------------------
# PLY export
# -----------------------------------------------------------------------------

def write_ply_point_cloud(path: Path, points: np.ndarray, colors: np.ndarray | None = None,
                          coordinate_unit: str = "mm",
                          fragment_ids: np.ndarray | None = None,
                          residual_mm: np.ndarray | None = None,
                          nearest_surface: np.ndarray | None = None) -> None:
    pts = np.asarray(points, dtype=np.float32)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    n = len(pts)
    has_color = colors is not None and len(colors) == n
    has_fid = fragment_ids is not None and len(fragment_ids) == n
    has_res = residual_mm is not None and len(residual_mm) == n
    has_surf = nearest_surface is not None and len(nearest_surface) == n
    fields = [("x","<f4"),("y","<f4"),("z","<f4")]
    header = ["ply", "format binary_little_endian 1.0",
              f"comment PotteryReconstruction3D v{__version__}",
              f"comment coordinate_unit={coordinate_unit}",
              f"element vertex {n}", "property float x", "property float y", "property float z"]
    if has_color:
        fields += [("red","u1"),("green","u1"),("blue","u1")]
        header += ["property uchar red","property uchar green","property uchar blue"]
    if has_fid:
        fields += [("fragment_id","<i4")]; header += ["property int fragment_id"]
    if has_res:
        fields += [("residual_mm","<f4")]; header += ["property float residual_mm"]
    if has_surf:
        fields += [("nearest_surface","u1")]; header += ["property uchar nearest_surface"]
    header += ["end_header"]
    arr = np.empty(n, dtype=np.dtype(fields))
    arr["x"], arr["y"], arr["z"] = pts[:,0], pts[:,1], pts[:,2]
    if has_color:
        c = np.asarray(colors,np.uint8)[:,:3]; arr["red"],arr["green"],arr["blue"] = c[:,0],c[:,1],c[:,2]
    if has_fid: arr["fragment_id"] = np.asarray(fragment_ids,np.int32)
    if has_res: arr["residual_mm"] = np.asarray(residual_mm,np.float32)
    if has_surf: arr["nearest_surface"] = np.asarray(nearest_surface,np.uint8)
    with path.open("wb") as f:
        f.write(("\n".join(header)+"\n").encode("ascii")); arr.tofile(f)


def write_ply_mesh(path: Path, vertices: np.ndarray, faces: np.ndarray,
                   coordinate_unit: str = "mm") -> None:
    v = np.asarray(vertices, dtype=np.float32); faces = np.asarray(faces, dtype=np.int32)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as f:
        header = ["ply","format binary_little_endian 1.0",
                  f"comment PotteryReconstruction3D v{__version__}",
                  f"comment coordinate_unit={coordinate_unit}",
                  f"element vertex {len(v)}","property float x","property float y","property float z",
                  f"element face {len(faces)}","property list uchar int vertex_indices","end_header"]
        f.write(("\n".join(header)+"\n").encode("ascii")); v.astype("<f4",copy=False).tofile(f)
        if len(faces):
            dt=np.dtype([("count","u1"),("idx","<i4",(3,))]); arr=np.empty(len(faces),dtype=dt); arr["count"]=3;arr["idx"]=faces;arr.tofile(f)


def write_axis_ply(path: Path, axis: GlobalAxisFit, zmin: float, zmax: float,
                   coordinate_unit: str = "mm") -> None:
    a=normalize(axis.axis_direction); c=axis.axis_point_mm
    zz=np.linspace(zmin,zmax,200); p=c[None,:]+zz[:,None]*a[None,:]
    write_ply_point_cloud(path,p*export_scale(coordinate_unit),coordinate_unit=coordinate_unit)



# -----------------------------------------------------------------------------
# Main reconstruction pipeline
# -----------------------------------------------------------------------------

def run_reconstruction(input_path: Path,
                       input_unit: str = "m", export_unit_name: str = "m",
                       display_points: int = 220_000, face_samples: int = 320_000,
                       cluster_voxel_mm: float = 1.5, cluster_link_radius_vox: float = 2.1,
                       min_fragment_sample_points: int = 400,
                       profile_bin_mm: float = 2.0,
                       compatibility_max: float = 0.30,
                       radial_sign_min: float = 0.18,
                       thickness_min_mm: float = 1.0, thickness_max_mm: float = 20.0,
                       profile_max_gap_mm: float = 8.0,
                       profile_smooth_window_bins: int = 3,
                       axis_center_max_shift_mm: float = 20.0,
                       axis_max_tilt_deg: float = 12.0,
                       axis_objective_bin_mm: float = 3.0,
                       axis_search_directions: int = 480,
                       axis_surface_keep_fraction: float = 0.70,
                       mesh_angle_step_deg: float = 10.0,
                       progress: Callable[[str], None] | None = None) -> ReconstructionResult:
    t0=time.perf_counter(); input_path=Path(input_path); ins=unit_scale(input_unit)
    def note(s):
        print(s)
        if progress: progress(s)

    note(f"Loading: {input_path}")
    src=open_mesh_source(input_path)
    p_input, colors=src.sample_vertices(display_points,seed=0); p_mm=p_input*ins
    note(f"Display sample: {len(p_mm):,} vertices")
    labels, clusters=cluster_sample_points(p_mm,voxel_mm=cluster_voxel_mm,
                                           link_radius_vox=cluster_link_radius_vox,
                                           min_sample_points=min_fragment_sample_points)
    if not clusters:
        raise ValueError("No retained sherd clusters. Lower min-fragment threshold or check input unit.")
    note(f"Retained sherds: {len(clusters)} (coordinates preserved)")

    fc_input, fn=src.sample_face_centroids_normals(face_samples,seed=2); fc_mm=fc_input*ins
    face_labels=assign_face_labels(fc_mm,p_mm,labels,max_distance_mm=max(3.0,cluster_voxel_mm*2.8))
    keep=face_labels>0
    note(f"Face samples assigned to retained sherds: {np.count_nonzero(keep):,}/{len(face_labels):,}")
    if np.count_nonzero(keep)<500:
        raise ValueError("Too few face samples could be assigned to retained sherds")

    placement=p_mm[labels>0]
    axis=estimate_global_axis(fc_mm[keep],fn[keep],face_labels[keep],placement,
                              max_center_shift_mm=axis_center_max_shift_mm,
                              max_tilt_deg=axis_max_tilt_deg,
                              objective_bin_mm=axis_objective_bin_mm,
                              normal_max=compatibility_max,
                              radial_sign_min=radial_sign_min,
                              progress=note)

    profile=build_global_profile(fc_mm,fn,face_labels,axis,
                                 bin_mm=profile_bin_mm,
                                 compatibility_max=compatibility_max,
                                 radial_sign_min=radial_sign_min,
                                 thickness_min_mm=thickness_min_mm,
                                 thickness_max_mm=thickness_max_mm,
                                 max_interpolation_gap_mm=profile_max_gap_mm,
                                 smooth_window_bins=profile_smooth_window_bins)
    note(f"Global profile: bins={len(profile.z_mm)}, outer={np.count_nonzero(np.isfinite(profile.outer_model_mm))}, inner={np.count_nonzero(np.isfinite(profile.inner_model_mm))}")
    note(f"Surface orientation: {profile.orientation_mapping}")

    mv,mf=build_reconstructed_mesh(profile,axis,angle_step_deg=mesh_angle_step_deg)
    residual,surf=evaluate_point_residuals(p_mm,axis,profile)
    quality=fragment_fit_quality(p_mm,labels,residual)

    metadata={
        "program":"PotteryReconstruction3D","version":__version__,
        "source_file":str(input_path),"source_sha256":sha256_file(input_path),
        "input_unit":input_unit,"internal_unit":"mm","export_unit":export_unit_name,
        "method":{
            "placement_policy":"STRICT IMMUTABLE: input sherd XYZ coordinates and orientations are never translated, rotated, shifted, reversed, or refined",
            "fragment_clustering":"spatial labels only; no coordinate transform",
            "coordinate_transform_applied_to_sherds":False,
            "axis":"one global axis estimated by unrestricted 3-D consensus search from surface-like sherd faces; no input-axis or centre prior",
            "profile":"inner/outer inferred in global axis coordinates from revolution-compatible face normals and radial-normal orientation",
            "mesh":"axisymmetric shell generated directly in original/current placement coordinate system",
            "old_v0_1_alignment":"not used: no independent fragment axes, no arbitrary Z shift, no Z reversal",
            "pose_refinement":"not implemented in v0.2.2; sherd poses are immutable",
        },
        "parameters":{
            "display_points":int(display_points),"face_samples":int(face_samples),
            "cluster_voxel_mm":float(cluster_voxel_mm),"cluster_link_radius_vox":float(cluster_link_radius_vox),
            "min_fragment_sample_points":int(min_fragment_sample_points),"profile_bin_mm":float(profile_bin_mm),
            "compatibility_max":float(compatibility_max),"radial_sign_min":float(radial_sign_min),
            "thickness_min_mm":float(thickness_min_mm),"thickness_max_mm":float(thickness_max_mm),
            "profile_max_gap_mm":float(profile_max_gap_mm),"profile_smooth_window_bins":int(profile_smooth_window_bins),
            "axis_center_max_shift_mm_deprecated":float(axis_center_max_shift_mm),"axis_max_tilt_deg_deprecated":float(axis_max_tilt_deg),
            "axis_objective_bin_mm":float(axis_objective_bin_mm),"axis_search_directions":int(axis_search_directions),
            "axis_surface_keep_fraction":float(axis_surface_keep_fraction),"mesh_angle_step_deg":float(mesh_angle_step_deg),
        },
        "axis":jsonable(asdict(axis)),
        "clusters":[asdict(c) for c in clusters],
        "fragment_fit_quality":jsonable(quality),
        "elapsed_seconds":float(time.perf_counter()-t0),
    }
    note(f"Done: {metadata['elapsed_seconds']:.2f} s")
    return ReconstructionResult(
        source_path=input_path,input_unit=input_unit,export_unit=export_unit_name,
        clusters=clusters,display_points_input=p_input,display_points_mm=p_mm,
        display_labels=labels,display_colors=colors,
        face_centroids_mm=fc_mm,face_normals=fn,face_labels=face_labels,
        axis=axis,profile=profile,mesh_vertices_mm=mv,mesh_faces=mf,
        point_residual_mm=residual,point_nearest_surface=surf,
        fragment_fit_quality=quality,metadata=metadata
    )


# -----------------------------------------------------------------------------
# Export complete result
# -----------------------------------------------------------------------------

def export_result(result: ReconstructionResult, out: Path) -> None:
    out=Path(out);out.mkdir(parents=True,exist_ok=True); fragdir=out/"fragment_samples";fragdir.mkdir(exist_ok=True)
    s=export_scale(result.export_unit)

    write_ply_mesh(out/"reconstructed_shell.ply",result.mesh_vertices_mm*s,result.mesh_faces,coordinate_unit=result.export_unit)
    rgb=np.zeros((len(result.display_points_mm),3),np.uint8);rgb[result.display_labels==0]=[170,170,170]
    for cl in result.clusters: rgb[result.display_labels==cl.fragment_id]=fragment_color(cl.fragment_id)
    write_ply_point_cloud(out/"current_placement_sample.ply",result.display_points_mm*s,rgb,
                          coordinate_unit=result.export_unit,fragment_ids=result.display_labels,
                          residual_mm=result.point_residual_mm,nearest_surface=result.point_nearest_surface)
    for cl in result.clusters:
        fid=cl.fragment_id;sel=result.display_labels==fid
        write_ply_point_cloud(fragdir/f"fragment_{fid:03d}_current_sample.ply",result.display_points_mm[sel]*s,
                              np.tile(fragment_color(fid),(np.count_nonzero(sel),1)),coordinate_unit=result.export_unit,
                              fragment_ids=np.full(np.count_nonzero(sel),fid,np.int32),residual_mm=result.point_residual_mm[sel],
                              nearest_surface=result.point_nearest_surface[sel])

    zfinite=result.profile.z_mm[np.isfinite(result.profile.z_mm)]
    if len(zfinite):
        write_axis_ply(out/"global_axis.ply",result.axis,float(zfinite.min()-10),float(zfinite.max()+10),coordinate_unit=result.export_unit)

    p=result.profile
    with (out/"global_profile.csv").open("w",newline="",encoding="utf-8-sig") as f:
        fields=["z_mm","inner_observed_mm","outer_observed_mm","inner_model_mm","outer_model_mm","thickness_model_mm",
                "support_inner_points","support_outer_points","support_inner_fragments","support_outer_fragments"]
        w=csv.writer(f);w.writerow(fields)
        for row in zip(p.z_mm,p.inner_observed_mm,p.outer_observed_mm,p.inner_model_mm,p.outer_model_mm,p.thickness_model_mm,
                       p.support_inner_points,p.support_outer_points,p.support_inner_fragments,p.support_outer_fragments):
            w.writerow([jsonable(x) for x in row])

    with (out/"fragment_fit_quality.csv").open("w",newline="",encoding="utf-8-sig") as f:
        fields=["fragment_id","status","n","n_trimmed","rmse_mm","median_abs_mm","p90_abs_mm","mean_abs_mm"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for fid in sorted(result.fragment_fit_quality):
            q=result.fragment_fit_quality[fid];w.writerow({"fragment_id":fid,**{k:q.get(k) for k in fields if k!="fragment_id"}})

    # Explicit residual export uses the exact same coordinates as current_placement_sample.ply.
    write_ply_point_cloud(out/"current_placement_residuals.ply",result.display_points_mm*s,rgb,
                          coordinate_unit=result.export_unit,fragment_ids=result.display_labels,
                          residual_mm=result.point_residual_mm,nearest_surface=result.point_nearest_surface)
    invariance={"sherd_pose_policy":"immutable","coordinate_transform_applied":False,
                "max_abs_position_change_mm":0.0,"note":"Only unit conversion is applied during geometry export."}
    (out/"coordinate_invariance_check.json").write_text(json.dumps(invariance,ensure_ascii=False,indent=2),encoding="utf-8")

    (out/"reconstruction_metadata.json").write_text(json.dumps(jsonable(result.metadata),ensure_ascii=False,indent=2),encoding="utf-8")


# -----------------------------------------------------------------------------
# GUI -- inherited interaction style, strict fixed-placement policy
# -----------------------------------------------------------------------------

def launch_gui(initial_input: Path | None = None):
    try:
        from PySide6 import QtCore, QtWidgets
        import pyvista as pv
        from pyvistaqt import QtInteractor
    except ModuleNotFoundError as exc:
        raise RuntimeError("GUI mode requires PySide6, pyvista, and pyvistaqt") from exc

    class MainWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__(); self.setWindowTitle(f"PotteryReconstruction3D v{__version__}"); self.resize(1540,940)
            self.source_path:Path|None=None; self.result:ReconstructionResult|None=None
            self.preview_points_input=None;self.preview_points_mm=None;self.preview_labels=None;self.preview_colors=None
            self.guide_ring_set:GuideRingSet|None=None; self.guide_registration:GuideRegistration|None=None
            central=QtWidgets.QWidget();self.setCentralWidget(central);layout=QtWidgets.QHBoxLayout(central)
            controls=QtWidgets.QWidget();controls.setMaximumWidth(430);form=QtWidgets.QVBoxLayout(controls);self.form=form

            self.open_btn=QtWidgets.QPushButton("1  PLY / OBJを開く");self.open_btn.clicked.connect(self.open_dialog);form.addWidget(self.open_btn)
            ur=QtWidgets.QHBoxLayout();ur.addWidget(QtWidgets.QLabel("入力単位"));self.input_unit=QtWidgets.QComboBox();self.input_unit.addItems(["m","mm","cm"]);ur.addWidget(self.input_unit);ur.addWidget(QtWidgets.QLabel("出力単位"));self.output_unit=QtWidgets.QComboBox();self.output_unit.addItems(["m","mm","cm"]);ur.addWidget(self.output_unit);form.addLayout(ur)

            self.display_points=self._spin("表示点数",20000,800000,220000,20000)
            self.face_samples=self._spin("法線解析 face sample",20000,1000000,320000,20000)
            self.cluster_voxel=self._dspin("破片分離 voxel [mm]",0.2,10.0,1.5,0.1,2)
            self.min_frag=self._spin("最小 sample 点数",30,50000,400,50)
            form.addWidget(self._separator("Global axis / current placement"))
            self.axis_search_dirs=self._spin("軸探索方向数",120,1200,480,60)
            self.axis_surface_keep=self._dspin("軸推定 surface-like face 採用率",0.40,0.90,0.70,0.05,2)
            self.compatibility=self._dspin("器面法線残差上限",0.05,0.8,0.30,0.01,2)
            self.radial_sign=self._dspin("inner/outer 法線radial成分下限",0.02,0.8,0.18,0.02,2)
            self.profile_bin=self._dspin("縦断面 bin [mm]",0.5,10.0,2.0,0.5,2)
            self.thick_min=self._dspin("最小器厚 [mm]",0.2,30.0,1.0,0.2,2)
            self.thick_max=self._dspin("最大器厚 [mm]",1.0,80.0,20.0,1.0,1)
            self.mesh_ang=self._dspin("復元mesh 角度step [deg]",2.0,30.0,10.0,1.0,1)
            fixed=QtWidgets.QLabel("🔒 破片座標・姿勢: 固定（解析中に移動・回転しません）");fixed.setStyleSheet("font-weight: bold; padding: 4px;");form.addWidget(fixed)

            self.preview_btn=QtWidgets.QPushButton("2  現在配置を表示 / 破片ラベル化");self.preview_btn.clicked.connect(self.load_preview);form.addWidget(self.preview_btn)
            self.analyze_btn=QtWidgets.QPushButton("3  共通回転軸・inner/outer・縦断面を推定");self.analyze_btn.clicked.connect(self.analyze);form.addWidget(self.analyze_btn)
            self.show_recon_btn=QtWidgets.QPushButton("4  初期配置固定＋推定器形＋偏差を表示");self.show_recon_btn.clicked.connect(self.show_reconstruction);form.addWidget(self.show_recon_btn)
            self.export_btn=QtWidgets.QPushButton("5  PLY / CSV / JSONを書き出す");self.export_btn.clicked.connect(self.export_dialog);form.addWidget(self.export_btn)

            form.addWidget(self._separator("External reconstructed rings / arcs"))
            gr=QtWidgets.QHBoxLayout();gr.addWidget(QtWidgets.QLabel("円弧PLY単位（未記録時）"));self.guide_unit=QtWidgets.QComboBox();self.guide_unit.addItems(["m","mm","cm"]);gr.addWidget(self.guide_unit);form.addLayout(gr)
            self.guide_load_btn=QtWidgets.QPushButton("6  復元円 / 円弧 PLY群を読み込む");self.guide_load_btn.clicked.connect(self.load_guide_rings_dialog);form.addWidget(self.guide_load_btn)
            self.guide_align_btn=QtWidgets.QPushButton("7  中心軸でアラインメント＋軸方向 r(z) fit");self.guide_align_btn.clicked.connect(self.align_guide_rings);form.addWidget(self.guide_align_btn)
            self.guide_export_btn=QtWidgets.QPushButton("8  アライン済み円PLYを書き出す");self.guide_export_btn.clicked.connect(self.export_aligned_guide);form.addWidget(self.guide_export_btn)
            self.guide_info=QtWidgets.QLabel("guide rings: 未読込");self.guide_info.setWordWrap(True);form.addWidget(self.guide_info)

            form.addWidget(QtWidgets.QLabel("破片 / fixed-placement fit"));self.table=QtWidgets.QTableWidget(0,5);self.table.setHorizontalHeaderLabels(["ID","sample","RMSE [mm]","median [mm]","p90 [mm]"]);self.table.horizontalHeader().setStretchLastSection(True);self.table.setMinimumHeight(210);form.addWidget(self.table)
            self.status=QtWidgets.QPlainTextEdit();self.status.setReadOnly(True);self.status.setMaximumBlockCount(800);form.addWidget(self.status,1)

            layout.addWidget(controls);self.plotter=QtInteractor(central);layout.addWidget(self.plotter.interactor,1);self.plotter.set_background("white");self.plotter.add_axes()
            if initial_input:
                self.source_path=Path(initial_input);self._log(f"Input: {self.source_path}");QtCore.QTimer.singleShot(150,self.load_preview)

        def _separator(self,text):
            lab=QtWidgets.QLabel(text);font=lab.font();font.setBold(True);lab.setFont(font);return lab
        def _spin(self,label,lo,hi,value,step=1):
            row=QtWidgets.QHBoxLayout();row.addWidget(QtWidgets.QLabel(label));w=QtWidgets.QSpinBox();w.setRange(lo,hi);w.setValue(value);w.setSingleStep(step);row.addWidget(w);self.form.addLayout(row);return w
        def _dspin(self,label,lo,hi,value,step,decimals):
            row=QtWidgets.QHBoxLayout();row.addWidget(QtWidgets.QLabel(label));w=QtWidgets.QDoubleSpinBox();w.setRange(lo,hi);w.setValue(value);w.setSingleStep(step);w.setDecimals(decimals);row.addWidget(w);self.form.addLayout(row);return w
        def _log(self,text): self.status.appendPlainText(str(text));QtWidgets.QApplication.processEvents()
        def open_dialog(self):
            fn,_=QtWidgets.QFileDialog.getOpenFileName(self,"土器破片3Dモデル","","Mesh (*.ply *.obj);;All files (*)")
            if fn:self.source_path=Path(fn);self.result=None;self._log(f"Input: {self.source_path}")

        def _plot_polyline(self, points, color="black", width=3, closed=False):
            """Display one connected polyline in PyVista.

            Plotter.add_lines(points) interprets points as independent endpoint pairs and
            therefore requires an even point count.  Guide circles/arcs are ordered
            polylines, so they must be supplied to VTK as one polyline cell instead.
            """
            pts=np.asarray(points,dtype=float)
            if len(pts)<2:
                return
            ids=np.arange(len(pts),dtype=np.int64)
            if closed:
                ids=np.concatenate([ids,np.array([0],dtype=np.int64)])
            cell=np.concatenate([np.array([len(ids)],dtype=np.int64),ids])
            line=pv.PolyData(pts)
            line.lines=cell
            self.plotter.add_mesh(line,color=color,line_width=width)

        def _add_guide_rings(self, gui_sc: float):
            if self.guide_ring_set is None:
                return
            M = self.guide_registration.matrix_mm if self.guide_registration is not None else None
            for ring in self.guide_ring_set.rings:
                pmm = apply_homogeneous(ring.ordered_points_mm, M) if M is not None else ring.ordered_points_mm
                pts = pmm * gui_sc
                if len(pts) < 2:
                    continue
                close = ring.arc_span_deg >= 345.0
                col = "deepskyblue" if ring.kind == "outer" else ("orange" if ring.kind == "inner" else "magenta")
                self._plot_polyline(pts,color=col,width=3,closed=close)
            if self.guide_ring_set.axis_point_mm is not None and self.guide_ring_set.axis_direction is not None:
                centers=np.vstack([r.center_mm for r in self.guide_ring_set.rings])
                g=self.guide_ring_set.axis_point_mm;a=self.guide_ring_set.axis_direction
                zz=(centers-g)@a
                lo=float(np.min(zz)-10);hi=float(np.max(zz)+10)
                raw=np.vstack([g+lo*a,g+hi*a])
                if M is not None:
                    raw=apply_homogeneous(raw,M)
                    axis_col="darkorange"
                else:
                    axis_col="red"
                self._plot_polyline(raw*gui_sc,color=axis_col,width=3,closed=False)

        def _render_preview(self, reset_camera: bool = True):
            if self.preview_points_input is None or self.preview_labels is None:
                return
            camera = self.plotter.camera_position
            self.plotter.clear()
            p=self.preview_points_input;labels=self.preview_labels
            poly=pv.PolyData(p);rgb=np.zeros((len(p),3),np.uint8);rgb[labels==0]=[170,170,170]
            for fid in np.unique(labels):
                if fid>0:rgb[labels==fid]=fragment_color(int(fid))
            poly["fragment_rgb"]=rgb;self.plotter.add_points(poly,scalars="fragment_rgb",rgb=True,point_size=3.0,render_points_as_spheres=False)
            gui_sc=1.0/unit_scale(self.input_unit.currentText())
            self._add_guide_rings(gui_sc)
            self.plotter.add_axes()
            if reset_camera:
                self.plotter.reset_camera();self.plotter.view_isometric()
            elif camera is not None:
                self.plotter.camera_position=camera
            self.plotter.render()

        def load_guide_rings_dialog(self):
            start = str(self.source_path.parent) if self.source_path is not None else str(Path.cwd())
            d=QtWidgets.QFileDialog.getExistingDirectory(self,"復元円 / 円弧 PLY群のフォルダを選択",start)
            if not d:
                return
            try:
                gs=load_guide_ring_folder(Path(d),fallback_unit=self.guide_unit.currentText())
                self.guide_ring_set=gs; self.guide_registration=None
                radii=np.asarray([r.radius_mm for r in gs.rings],float)
                spans=np.asarray([r.arc_span_deg for r in gs.rings],float)
                rms=np.asarray([r.circle_rms_mm for r in gs.rings],float)
                txt=f"guide rings: {len(gs.rings)} / radius {np.min(radii):.2f}–{np.max(radii):.2f} mm / median arc {np.median(spans):.1f}°"
                self.guide_info.setText(txt)
                self._log(f"Guide-ring folder: {gs.folder}")
                self._log(f"Guide rings loaded: {len(gs.rings)} (RAW coordinates; registration NOT applied)")
                self._log(f"  radius range: {np.min(radii):.3f} .. {np.max(radii):.3f} mm")
                self._log(f"  arc span median: {np.median(spans):.3f} deg; circle RMS median: {np.median(rms):.4f} mm")
                if gs.axis_direction is not None:
                    self._log(f"  ring-stack axis (guide coordinates): {gs.axis_direction.tolist()} / center-line RMS={gs.axis_line_rms_mm:.4f} mm")
                errors=getattr(gs,"load_errors",[])
                if errors:
                    self._log(f"  skipped PLY files: {len(errors)}")
                    for e in errors[:5]:self._log(f"    {e}")
                if self.result is not None:
                    self.show_reconstruction()
                elif self.preview_points_input is not None:
                    self._render_preview(reset_camera=False)
                else:
                    self.plotter.clear();self._add_guide_rings(1.0/unit_scale(self.input_unit.currentText()));self.plotter.add_axes();self.plotter.reset_camera();self.plotter.view_isometric();self.plotter.render()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"円弧PLY読み込みエラー",str(exc));self._log(f"GUIDE ERROR: {exc}")

        def align_guide_rings(self):
            if self.result is None:
                QtWidgets.QMessageBox.information(self,"未解析","先に 3 共通回転軸・inner/outer・縦断面を推定 を実行してください。")
                return
            if self.guide_ring_set is None:
                QtWidgets.QMessageBox.information(self,"円弧未読込","先に 6 復元円 / 円弧 PLY群を読み込む を実行してください。")
                return
            try:
                reg=register_guide_axis_to_reconstruction(self.guide_ring_set,self.result)
                self.guide_registration=reg
                self.guide_info.setText(
                    f"guide rings: {len(self.guide_ring_set.rings)} / ALIGNED / "
                    f"radius RMSE {reg.radius_rmse_mm:.3f} mm / axial shift {reg.axial_shift_mm:.2f} mm"
                )
                self._log("--- guide center-axis registration ---")
                self._log("🔒 Sherds remain fixed; ONLY the external guide rings are transformed.")
                self._log(f"  axis orientation sign: {reg.axis_orientation_sign:+d}")
                self._log(f"  angular correction: {reg.angular_correction_deg:.3f} deg")
                self._log(f"  lateral axis shift: {reg.lateral_shift_mm:.3f} mm")
                self._log(f"  axial profile shift: {reg.axial_shift_mm:.3f} mm")
                self._log(f"  radius fit: RMSE={reg.radius_rmse_mm:.3f} mm; median={reg.median_abs_radius_mm:.3f} mm; p90={reg.p90_abs_radius_mm:.3f} mm")
                self._log(f"  matched rings: {reg.matched_rings}/{reg.total_rings}; overlap Z={reg.overlap_z_mm:.3f} mm")
                self._log(f"  bottom/top extent mismatch: {reg.bottom_mismatch_mm:+.3f} / {reg.top_mismatch_mm:+.3f} mm")
                self.show_reconstruction()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"中心軸アラインメントエラー",str(exc));self._log(f"GUIDE ALIGN ERROR: {exc}")

        def export_aligned_guide(self):
            if self.guide_ring_set is None or self.guide_registration is None:
                QtWidgets.QMessageBox.information(self,"未アライン","先に 7 中心軸でアラインメント を実行してください。")
                return
            start=str(self.source_path.parent) if self.source_path is not None else str(Path.cwd())
            d=QtWidgets.QFileDialog.getExistingDirectory(self,"アライン済み円PLYの出力先を選択",start)
            if not d:
                return
            name=(self.source_path.stem if self.source_path is not None else "guide")+"_AlignedGuideRings_v024"
            out=Path(d)/name
            try:
                export_registered_guide_rings(self.guide_ring_set,self.guide_registration,out,self.output_unit.currentText())
                self._log(f"Aligned guide exported: {out}")
                QtWidgets.QMessageBox.information(self,"アライン済み円PLY書き出し完了",str(out))
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"円PLY書き出しエラー",str(exc));self._log(f"GUIDE EXPORT ERROR: {exc}")

        def load_preview(self):
            if self.source_path is None:self.open_dialog()
            if self.source_path is None:return
            try:
                src=open_mesh_source(self.source_path);p,c=src.sample_vertices(self.display_points.value(),seed=0);scale=unit_scale(self.input_unit.currentText());pmm=p*scale
                labels,clusters=cluster_sample_points(pmm,voxel_mm=self.cluster_voxel.value(),min_sample_points=self.min_frag.value())
                self.preview_points_input=p;self.preview_points_mm=pmm;self.preview_labels=labels;self.preview_colors=c
                self._render_preview(reset_camera=True);self._log(f"現在配置: {len(p):,} points / retained fragments={len(clusters)}; coordinates unchanged")
                self._populate_clusters(clusters)
            except Exception as exc:QtWidgets.QMessageBox.critical(self,"読み込みエラー",str(exc));self._log(f"ERROR: {exc}")
        def _populate_clusters(self,clusters):
            self.table.setRowCount(len(clusters))
            for r,cl in enumerate(clusters):
                vals=[cl.fragment_id,cl.sample_count,"","",""]
                for cc,v in enumerate(vals):self.table.setItem(r,cc,QtWidgets.QTableWidgetItem(str(v)))
        def analyze(self):
            if self.source_path is None:self.open_dialog()
            if self.source_path is None:return
            try:
                self._log("--- v0.2.4 STRICT fixed-placement reconstruction start ---")
                self.result=run_reconstruction(self.source_path,input_unit=self.input_unit.currentText(),export_unit_name=self.output_unit.currentText(),display_points=self.display_points.value(),face_samples=self.face_samples.value(),cluster_voxel_mm=self.cluster_voxel.value(),min_fragment_sample_points=self.min_frag.value(),profile_bin_mm=self.profile_bin.value(),compatibility_max=self.compatibility.value(),radial_sign_min=self.radial_sign.value(),thickness_min_mm=self.thick_min.value(),thickness_max_mm=self.thick_max.value(),axis_search_directions=self.axis_search_dirs.value(),axis_surface_keep_fraction=self.axis_surface_keep.value(),mesh_angle_step_deg=self.mesh_ang.value(),progress=self._log)
                if self.preview_points_input is not None and len(self.preview_points_input)==len(self.result.display_points_input):
                    max_delta=float(np.max(np.abs(self.preview_points_input-self.result.display_points_input)))
                    self._log(f"LOCK CHECK: max |analysis XYZ - preview XYZ| = {max_delta:.12g} input units")
                    if max_delta>1e-12:
                        raise RuntimeError(f"Fixed-placement invariant failed: sample coordinates changed by {max_delta}")
                self._populate_result_table();self.show_reconstruction()
            except Exception as exc:QtWidgets.QMessageBox.critical(self,"解析エラー",str(exc));self._log(f"ERROR: {exc}")
        def _populate_result_table(self):
            if self.result is None:return
            self.table.setRowCount(len(self.result.clusters))
            for rr,cl in enumerate(self.result.clusters):
                q=self.result.fragment_fit_quality.get(cl.fragment_id,{})
                vals=[
                    cl.fragment_id,
                    cl.sample_count,
                    "" if q.get("rmse_mm") is None else f"{q['rmse_mm']:.3f}",
                    "" if q.get("median_abs_mm") is None else f"{q['median_abs_mm']:.3f}",
                    "" if q.get("p90_abs_mm") is None else f"{q['p90_abs_mm']:.3f}",
                ]
                for cc,v in enumerate(vals):self.table.setItem(rr,cc,QtWidgets.QTableWidgetItem(str(v)))
        def _add_axis(self,sc):
            if self.result is None:return
            z=self.result.profile.z_mm;a=self.result.axis.axis_direction;c=self.result.axis.axis_point_mm
            if len(z):
                zz=np.linspace(float(np.min(z)-10),float(np.max(z)+10),100);pts=(c[None,:]+zz[:,None]*a[None,:])*sc
                self._plot_polyline(pts,color="black",width=3,closed=False)
        def show_reconstruction(self):
            if self.result is None:
                QtWidgets.QMessageBox.information(self,"未解析","先に解析を実行してください。")
                return

            # Preserve the exact current camera to avoid making unchanged sherds look
            # as if they were reoriented. GUI geometry is always shown in INPUT units.
            camera = self.plotter.camera_position
            self.plotter.clear()
            gui_sc = 1.0 / unit_scale(self.input_unit.currentText())

            # Exact unchanged sample coordinates from the input placement.
            p = self.result.display_points_input
            valid = np.isfinite(self.result.point_residual_mm)

            if np.any(valid):
                pv_valid = pv.PolyData(p[valid])
                pv_valid["residual_mm"] = self.result.point_residual_mm[valid]
                clim_hi = max(3.0,float(np.nanquantile(self.result.point_residual_mm[valid],0.90)))
                self.plotter.add_points(pv_valid,scalars="residual_mm",cmap="viridis",clim=[0.0,clim_hi],point_size=3.0,render_points_as_spheres=False)
            if np.any(~valid):
                pv_bad = pv.PolyData(p[~valid])
                self.plotter.add_points(pv_bad,color="lightgray",point_size=2.5,render_points_as_spheres=False)

            if len(self.result.mesh_vertices_mm) and len(self.result.mesh_faces):
                v=self.result.mesh_vertices_mm*gui_sc
                f=np.hstack([np.full((len(self.result.mesh_faces),1),3,np.int64),self.result.mesh_faces.astype(np.int64)]).ravel()
                mesh=pv.PolyData(v,f)
                self.plotter.add_mesh(mesh,color="lightgray",opacity=0.25,show_edges=True,edge_color="gray")

            self._add_axis(gui_sc)
            self._add_guide_rings(gui_sc)
            self.plotter.add_axes()
            if camera is not None:
                self.plotter.camera_position=camera
            else:
                self.plotter.reset_camera();self.plotter.view_isometric()
            self.plotter.render()
            self._log("🔒 初期配置のXYZ・姿勢は完全固定。推定器形と偏差だけを重ねて表示しました。")

        def export_dialog(self):
            if self.result is None:QtWidgets.QMessageBox.information(self,"未解析","先に解析を実行してください。");return
            d=QtWidgets.QFileDialog.getExistingDirectory(self,"出力先フォルダを選択",str(self.source_path.parent))
            if not d:return
            out=Path(d)/f"{self.source_path.stem}_Reconstruction3D_v024"
            try:export_result(self.result,out);self._log(f"Exported: {out}");QtWidgets.QMessageBox.information(self,"書き出し完了",str(out))
            except Exception as exc:QtWidgets.QMessageBox.critical(self,"書き出しエラー",str(exc))

    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv);win=MainWindow();win.show();sys.exit(app.exec())


# -----------------------------------------------------------------------------
# Self test: meaningful-placement synthetic vessel fragments
# -----------------------------------------------------------------------------

def self_test():
    rng=np.random.default_rng(42)
    true_c=np.array([12.0,-7.0,4.0]);true_a=normalize(np.array([0.04,-0.03,1.0]));e1,e2=orthonormal_perp_basis(true_a)
    allp=[];alln=[];alll=[]
    for fid,(t0,t1,z0,z1) in enumerate([(-0.7,-0.1,-35,10),(0.1,0.75,-15,35),(1.0,1.55,-30,25),(-1.55,-1.0,0,38)],start=1):
        n=1800;th=rng.uniform(t0,t1,n);z=rng.uniform(z0,z1,n);outer=rng.random(n)>.5;rmid=60+0.045*z;r=rmid+np.where(outer,2.5,-2.5)
        p=true_c[None,:]+z[:,None]*true_a[None,:]+r[:,None]*(np.cos(th)[:,None]*e1[None,:]+np.sin(th)[:,None]*e2[None,:])
        # normal in meridional plane, with opposite orientation on inner surface
        er=np.cos(th)[:,None]*e1[None,:]+np.sin(th)[:,None]*e2[None,:];nn=er-0.045*true_a[None,:];nn/=np.linalg.norm(nn,axis=1)[:,None];nn[~outer]*=-1
        # fracture-like contamination
        m=120;nn[:m]=rng.normal(size=(m,3));nn[:m]/=np.linalg.norm(nn[:m],axis=1)[:,None]
        allp.append(p);alln.append(nn);alll.append(np.full(n,fid,np.int32))
    P=np.vstack(allp);N=np.vstack(alln);L=np.concatenate(alll)
    axis=estimate_global_axis(P,N,L,P,max_center_shift_mm=15,max_tilt_deg=8,objective_bin_mm=3,normal_max=.32,radial_sign_min=.15,max_fit_samples_per_fragment=1600,max_fit_samples_total=6400,seed=9)
    ang=math.degrees(math.acos(float(np.clip(abs(np.dot(axis.axis_direction,true_a)),-1,1))))
    line_offset=np.linalg.norm((axis.axis_point_mm-true_c)-np.dot(axis.axis_point_mm-true_c,true_a)*true_a)
    if ang>4.0 or line_offset>8.0:raise SystemExit(f"SELF TEST FAILED axis: angle={ang:.3f}, line offset={line_offset:.3f}")
    prof=build_global_profile(P,N,L,axis,bin_mm=3,compatibility_max=.35,radial_sign_min=.12,min_bin_samples=20,thickness_min_mm=2,thickness_max_mm=10,max_interpolation_gap_mm=8,smooth_window_bins=3)
    th=prof.thickness_model_mm[np.isfinite(prof.thickness_model_mm)]
    if len(th)<5 or abs(float(np.median(th))-5.0)>1.5:raise SystemExit(f"SELF TEST FAILED thickness={np.median(th) if len(th) else np.nan}")
    v,f=build_reconstructed_mesh(prof,axis,angle_step_deg=15)
    if len(v)==0 or len(f)==0:raise SystemExit("SELF TEST FAILED mesh")

    # External ring/arc fit test: a partial tilted circle must recover centre/radius.
    rg_c=np.array([21.0,-14.0,8.0]);rg_n=normalize(np.array([0.3,-0.4,0.866]));rg_e1,rg_e2=orthonormal_perp_basis(rg_n)
    rg_t=np.linspace(-0.8,1.25,140);rg_r=73.0
    rg_p=rg_c[None,:]+rg_r*(np.cos(rg_t)[:,None]*rg_e1[None,:]+np.sin(rg_t)[:,None]*rg_e2[None,:])
    rg_fit_c,rg_fit_n,rg_fit_r,rg_plane,rg_circle,rg_ord=fit_circle_3d(rg_p)
    rg_center_err=float(np.linalg.norm(rg_fit_c-rg_c));rg_ang=math.degrees(math.acos(float(np.clip(abs(np.dot(rg_fit_n,rg_n)),-1,1))))
    if abs(rg_fit_r-rg_r)>1e-6 or rg_center_err>1e-5 or rg_ang>1e-4:
        raise SystemExit(f"SELF TEST FAILED guide ring: dr={rg_fit_r-rg_r}, dc={rg_center_err}, angle={rg_ang}")

    print(f"PotteryReconstruction3D v{__version__} SELF TEST PASSED")
    print(f"axis angular error : {ang:.3f} deg")
    print(f"axis line offset   : {line_offset:.3f} mm")
    print(f"radial dispersion : {axis.radial_dispersion_mm:.3f} mm")
    print(f"median thickness  : {np.median(th):.3f} mm")
    print(f"mesh vertices/faces: {len(v)}/{len(f)}")
    # Guide-axis registration test: arbitrary guide coordinate system must align
    # to the fixed target axis and recover the axial profile shift without moving target points.
    ta=normalize(np.array([0.20,-0.60,0.7745966692],float));tc=np.array([12.0,-8.0,20.0])
    gz_target=np.array([20.0,40.0,60.0,80.0]); gr_target=70.0+0.08*gz_target
    tp=ProfileModel(
        z_mm=np.arange(0.0,101.0,2.0),
        inner_observed_mm=np.full(51,np.nan), outer_observed_mm=70.0+0.08*np.arange(0.0,101.0,2.0),
        inner_model_mm=np.full(51,np.nan), outer_model_mm=70.0+0.08*np.arange(0.0,101.0,2.0),
        thickness_model_mm=np.full(51,np.nan), support_inner_points=np.zeros(51,int),
        support_outer_points=np.ones(51,int), support_inner_fragments=np.zeros(51,int),
        support_outer_fragments=np.ones(51,int), orientation_mapping="synthetic"
    )
    tax=GlobalAxisFit(tc,ta,tc,ta,0,0,0,1,0,0,100)
    ag=normalize(np.array([0.7,0.2,0.68556546],float));cg=np.array([-50.0,30.0,12.0]);ge1,ge2=orthonormal_perp_basis(ag)
    rings=[]
    centered=gz_target-np.mean(gz_target)
    for i,(zz,rad) in enumerate(zip(centered,gr_target)):
        th=np.linspace(0,2*np.pi,90,endpoint=False);cen=cg+zz*ag
        pts=cen+rad*np.cos(th)[:,None]*ge1+rad*np.sin(th)[:,None]*ge2
        rings.append(GuideRing(Path(f"synthetic_{i}.ply"),"mm","outer",pts,pts,cen,ag,rad,0,0,360.0))
    gs=GuideRingSet(Path("."),rings,cg,ag,0.0)
    fake=type("FakeResult",(),{})();fake.axis=tax;fake.profile=tp
    reg=register_guide_axis_to_reconstruction(gs,fake,shift_step_mm=1.0)
    centers_al=np.vstack([apply_homogeneous(r.center_mm[None,:],reg.matrix_mm)[0] for r in rings])
    target_centers=np.vstack([tc+z*ta for z in gz_target])
    center_err=float(np.max(np.linalg.norm(centers_al-target_centers,axis=1)))
    if center_err>0.05 or reg.radius_rmse_mm>0.05:
        raise SystemExit(f"SELF TEST FAILED guide registration: center={center_err}, radius={reg.radius_rmse_mm}")
    print(f"guide registration  : center max error {center_err:.4f} mm / radius RMSE {reg.radius_rmse_mm:.4f} mm")
    print(f"guide arc fit       : radius {rg_fit_r:.3f} mm / center error {rg_center_err:.3g} mm")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    ap=argparse.ArgumentParser(description="Reconstruct an axisymmetric vessel from sherds while preserving their current meaningful placement.")
    ap.add_argument("input",nargs="?",type=Path,help="Input PLY/OBJ. Omit to open in GUI.")
    ap.add_argument("--input-unit",choices=["mm","cm","m"],default="m")
    ap.add_argument("--export-unit",choices=["mm","cm","m"],default="m")
    ap.add_argument("--headless",action="store_true")
    ap.add_argument("--output-dir",type=Path)
    ap.add_argument("--display-points",type=int,default=220000)
    ap.add_argument("--face-samples",type=int,default=320000)
    ap.add_argument("--cluster-voxel-mm",type=float,default=1.5)
    ap.add_argument("--cluster-link-radius-vox",type=float,default=2.1)
    ap.add_argument("--min-fragment-sample-points",type=int,default=400)
    ap.add_argument("--profile-bin-mm",type=float,default=2.0)
    ap.add_argument("--compatibility-max",type=float,default=0.30)
    ap.add_argument("--radial-sign-min",type=float,default=0.18)
    ap.add_argument("--thickness-min-mm",type=float,default=1.0)
    ap.add_argument("--thickness-max-mm",type=float,default=20.0)
    ap.add_argument("--profile-max-gap-mm",type=float,default=8.0)
    ap.add_argument("--profile-smooth-window-bins",type=int,default=3)
    ap.add_argument("--axis-search-directions",type=int,default=480)
    ap.add_argument("--axis-surface-keep-fraction",type=float,default=0.70)
    ap.add_argument("--axis-objective-bin-mm",type=float,default=3.0)
    ap.add_argument("--mesh-angle-step-deg",type=float,default=10.0)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--version",action="version",version=f"%(prog)s {__version__}")
    args=ap.parse_args()
    if args.self_test:self_test();return
    if args.headless:
        if args.input is None:ap.error("--headless requires input")
        out=args.output_dir or args.input.parent/f"{args.input.stem}_Reconstruction3D_v023"
        r=run_reconstruction(args.input,input_unit=args.input_unit,export_unit_name=args.export_unit,display_points=args.display_points,face_samples=args.face_samples,cluster_voxel_mm=args.cluster_voxel_mm,cluster_link_radius_vox=args.cluster_link_radius_vox,min_fragment_sample_points=args.min_fragment_sample_points,profile_bin_mm=args.profile_bin_mm,compatibility_max=args.compatibility_max,radial_sign_min=args.radial_sign_min,thickness_min_mm=args.thickness_min_mm,thickness_max_mm=args.thickness_max_mm,profile_max_gap_mm=args.profile_max_gap_mm,profile_smooth_window_bins=args.profile_smooth_window_bins,axis_search_directions=args.axis_search_directions,axis_surface_keep_fraction=args.axis_surface_keep_fraction,axis_objective_bin_mm=args.axis_objective_bin_mm,mesh_angle_step_deg=args.mesh_angle_step_deg)
        export_result(r,out);print(f"Exported: {out}");return
    launch_gui(args.input)


if __name__=="__main__":
    main()

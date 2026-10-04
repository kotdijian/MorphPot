from __future__ import annotations

import sys
import math
import hashlib
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Optional, List, Any, Tuple

import numpy as np
import trimesh
from PIL import Image, ImageDraw, ImageFont
from scipy import sparse
from scipy.ndimage import gaussian_filter, laplace, map_coordinates, distance_transform_edt
from numba import njit

from PySide6.QtCore import Qt
from PySide6.QtGui import QImage, QPixmap
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QMainWindow, QMessageBox, QPushButton,
    QProgressBar, QScrollArea, QSizePolicy, QVBoxLayout, QWidget,
    QTabWidget, QToolBar, QListWidget, QListWidgetItem
)

APP_NAME = "Surface Enhancement Lab v0.7.3"

PROTOCOL_BASE_SCALE_MM = 0.80
PROTOCOL_BROAD_BASE_SCALE_MM = 1.60
PROTOCOL_SMOOTHING_ALPHA = 0.50
PROTOCOL_FINE_SCALE_MM = 0.35
PROTOCOL_TOP_POST_SIGMA_PX = 0.60
PROTOCOL_HIGHLIGHT_Q_LOW = 0.900
PROTOCOL_HIGHLIGHT_Q_HIGH = 0.995
PROTOCOL_HIGHLIGHT_STRENGTH = 0.22
PROTOCOL_HIGHLIGHT_GAMMA = 1.70
PROTOCOL_HIGHLIGHT_CAP_RATIO = 0.45
PROTOCOL_ALPHA_BLEND = 0.50
PROTOCOL_IMAGE_HEIGHT = 1200
PROTOCOL_MARGIN_RATIO = 0.03
PROTOCOL_DEPTH_MED_PCT = 0.0090
PROTOCOL_DEPTH_BROAD_PCT = 0.0240
PROTOCOL_DEPTH_MED_PERCENT = 0.90
PROTOCOL_DEPTH_BROAD_RATIO = PROTOCOL_DEPTH_BROAD_PCT / PROTOCOL_DEPTH_MED_PCT
PROTOCOL_MESH_TOP_SCALE_MM = 1.20
PROTOCOL_CURV_SIGMA_PERCENT = 6.0
PROTOCOL_THICKNESS_CONTRAST = 0.30
PROTOCOL_ROUGHNESS_RADIUS_MM = 1.00
PROTOCOL_ROUGHNESS_CONTRAST = 0.25
PROTOCOL_OPENNESS_RADIUS_MM = 2.00
PROTOCOL_HESSIAN_SIGMA_MM = 1.00
PROTOCOL_NORMAL_RADIUS_MM = 1.00
PROTOCOL_PRINCIPAL_SCALE_MM = 0.80
PROTOCOL_LINEARITY_SIGMA_PERCENT = 0.60
PROTOCOL_CONTINUITY_SPAN_PERCENT = 1.50
PROTOCOL_LINEARITY_WEIGHT = 0.75
PROTOCOL_CONTINUITY_WEIGHT = 1.25
TOP_SCALE_PRESETS = {"M8": 0.80, "M12": 1.20, "M16": 1.60}

LIGHT_DIR = np.array([0.35, 0.30, 0.886], dtype=np.float64)
LIGHT_DIR /= np.linalg.norm(LIGHT_DIR)
VIEW_DIR = np.array([0.0, 0.0, 1.0], dtype=np.float64)
DEFAULT_BG = 1.0

@dataclass(frozen=True)
class RenderParams:
    base_scale_mm: float = PROTOCOL_BASE_SCALE_MM
    broad_base_scale_mm: float = PROTOCOL_BROAD_BASE_SCALE_MM
    smoothing_alpha: float = PROTOCOL_SMOOTHING_ALPHA
    fine_scale_mm: float = PROTOCOL_FINE_SCALE_MM
    mesh_top_scale_mm: float = PROTOCOL_MESH_TOP_SCALE_MM
    depth_med_percent: float = PROTOCOL_DEPTH_MED_PERCENT
    top_post_sigma_px: float = PROTOCOL_TOP_POST_SIGMA_PX
    highlight_q_low: float = PROTOCOL_HIGHLIGHT_Q_LOW
    highlight_q_high: float = PROTOCOL_HIGHLIGHT_Q_HIGH
    highlight_strength: float = PROTOCOL_HIGHLIGHT_STRENGTH
    highlight_gamma: float = PROTOCOL_HIGHLIGHT_GAMMA
    highlight_cap_ratio: float = PROTOCOL_HIGHLIGHT_CAP_RATIO
    thickness_contrast: float = PROTOCOL_THICKNESS_CONTRAST
    roughness_radius_mm: float = PROTOCOL_ROUGHNESS_RADIUS_MM
    roughness_contrast: float = PROTOCOL_ROUGHNESS_CONTRAST
    curvature_sigma_percent: float = PROTOCOL_CURV_SIGMA_PERCENT
    openness_radius_mm: float = PROTOCOL_OPENNESS_RADIUS_MM
    hessian_sigma_mm: float = PROTOCOL_HESSIAN_SIGMA_MM
    normal_radius_mm: float = PROTOCOL_NORMAL_RADIUS_MM
    principal_scale_mm: float = PROTOCOL_PRINCIPAL_SCALE_MM
    linearity_sigma_percent: float = PROTOCOL_LINEARITY_SIGMA_PERCENT
    continuity_span_percent: float = PROTOCOL_CONTINUITY_SPAN_PERCENT
    linearity_weight: float = PROTOCOL_LINEARITY_WEIGHT
    continuity_weight: float = PROTOCOL_CONTINUITY_WEIGHT
    blend_mode: str = "multiply"
    alpha_blend: float = PROTOCOL_ALPHA_BLEND
    projection: str = "front"
    image_height: int = PROTOCOL_IMAGE_HEIGHT
    margin_ratio: float = PROTOCOL_MARGIN_RATIO
    base_choice: str = "Base reflect-suppressed"
    top1_choice: str = "Mesh-native medium relief"
    top2_choice: str = "None"
    auto_add_candidate: bool = True

@dataclass
class MeshState:
    path: Path
    sha256: str
    vertices: np.ndarray
    faces: np.ndarray
    edges_unique: np.ndarray
    bounds: np.ndarray
    median_edge: float
    operator: sparse.csr_matrix
    source_vertex_colors: Optional[np.ndarray] = None
    source_vertex_normals: Optional[np.ndarray] = None

@dataclass
class RenderBundle:
    params: RenderParams
    base_candidates: Dict[str, np.ndarray]
    top_candidates: Dict[str, np.ndarray]
    top_diagnostics: Dict[str, Dict[str, np.ndarray]]
    selected_base_name: str
    selected_top_names: List[str]
    result: np.ndarray
    mask: np.ndarray
    meta: Dict[str, Any]
    audit_text: str

# ---------------- numerical core ----------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def mm_to_iterations(scale_mm: float, median_edge_mm: float) -> int:
    return max(1, int(round((scale_mm / max(median_edge_mm, 1e-9)) ** 2)))


def build_operator(vertices: np.ndarray, edges: np.ndarray) -> Tuple[sparse.csr_matrix, float]:
    vi = edges[:, 0]
    vj = edges[:, 1]
    lengths = np.linalg.norm(vertices[vi] - vertices[vj], axis=1)
    lengths = np.maximum(lengths, 1e-6)
    weights = 1.0 / lengths
    n = len(vertices)
    rows = np.concatenate([vi, vj])
    cols = np.concatenate([vj, vi])
    data = np.concatenate([weights, weights])
    A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    P = sparse.diags(1.0 / row_sum) @ A
    return P.tocsr(), float(np.median(lengths))


def smooth_vertices(vertices: np.ndarray, P: sparse.csr_matrix, steps: int, alpha: float) -> np.ndarray:
    X = vertices.astype(np.float64).copy()
    for _ in range(int(steps)):
        PX = np.column_stack([P.dot(X[:, 0]), P.dot(X[:, 1]), P.dot(X[:, 2])])
        X = (1.0 - alpha) * X + alpha * PX
    return X


def vertex_normals_from_vertices(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    return np.asarray(mesh.vertex_normals, dtype=np.float64)


def projection_depth(vertices: np.ndarray, projection: str) -> np.ndarray:
    if projection == "front":
        return vertices[:, 2]
    if projection == "right":
        return vertices[:, 0]
    if projection == "bottom":
        return vertices[:, 1]
    raise ValueError(projection)


def normals_to_view(normals: np.ndarray, projection: str) -> np.ndarray:
    if projection == "front":
        return normals.copy()
    if projection == "right":
        return normals[:, [1, 2, 0]].copy()
    if projection == "bottom":
        return normals[:, [0, 2, 1]].copy()
    raise ValueError(projection)


def project_setup(vertices: np.ndarray, projection: str, height: int, margin: float):
    if projection == "front":
        a, b = vertices[:, 0], vertices[:, 1]
    elif projection == "right":
        a, b = vertices[:, 1], vertices[:, 2]
    elif projection == "bottom":
        a, b = vertices[:, 0], vertices[:, 2]
    else:
        raise ValueError(projection)
    depth = projection_depth(vertices, projection).astype(np.float64)
    amin, amax = float(a.min()), float(a.max())
    bmin, bmax = float(b.min()), float(b.max())
    amid = 0.5 * (amin + amax)
    bmid = 0.5 * (bmin + bmax)
    ar = max((amax - amin) * (1.0 + 2.0 * margin), 1e-9)
    br = max((bmax - bmin) * (1.0 + 2.0 * margin), 1e-9)
    width = max(1, int(round(height * ar / br)))
    amin2 = amid - ar / 2.0
    bmax2 = bmid + br / 2.0
    u = (a - amin2) / ar * (width - 1)
    v = (bmax2 - b) / br * (height - 1)
    return np.column_stack([u, v]).astype(np.float64), depth, width, height


@njit(cache=True)
def rasterize_normals(verts2d, depths, normals, faces, H, W):
    zbuf = np.full((H, W), -1e30, dtype=np.float64)
    nxbuf = np.zeros((H, W), dtype=np.float64)
    nybuf = np.zeros((H, W), dtype=np.float64)
    nzbuf = np.zeros((H, W), dtype=np.float64)
    mask = np.zeros((H, W), dtype=np.uint8)
    for fi in range(faces.shape[0]):
        i0, i1, i2 = faces[fi, 0], faces[fi, 1], faces[fi, 2]
        x0, y0 = verts2d[i0, 0], verts2d[i0, 1]
        x1, y1 = verts2d[i1, 0], verts2d[i1, 1]
        x2, y2 = verts2d[i2, 0], verts2d[i2, 1]
        minx = int(max(0, np.floor(min(x0, x1, x2))))
        maxx = int(min(W - 1, np.ceil(max(x0, x1, x2))))
        miny = int(max(0, np.floor(min(y0, y1, y2))))
        maxy = int(min(H - 1, np.ceil(max(y0, y1, y2))))
        denom = ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2))
        if abs(denom) < 1e-12:
            continue
        z0, z1, z2 = depths[i0], depths[i1], depths[i2]
        n0x, n0y, n0z = normals[i0, 0], normals[i0, 1], normals[i0, 2]
        n1x, n1y, n1z = normals[i1, 0], normals[i1, 1], normals[i1, 2]
        n2x, n2y, n2z = normals[i2, 0], normals[i2, 1], normals[i2, 2]
        for py in range(miny, maxy + 1):
            yy = py + 0.5
            for px in range(minx, maxx + 1):
                xx = px + 0.5
                w0 = ((y1 - y2) * (xx - x2) + (x2 - x1) * (yy - y2)) / denom
                if w0 < 0.0:
                    continue
                w1 = ((y2 - y0) * (xx - x2) + (x0 - x2) * (yy - y2)) / denom
                if w1 < 0.0:
                    continue
                w2 = 1.0 - w0 - w1
                if w2 < 0.0:
                    continue
                z = w0 * z0 + w1 * z1 + w2 * z2
                if z > zbuf[py, px]:
                    zbuf[py, px] = z
                    nx = w0 * n0x + w1 * n1x + w2 * n2x
                    ny = w0 * n0y + w1 * n1y + w2 * n2y
                    nz = w0 * n0z + w1 * n1z + w2 * n2z
                    nlen = (nx * nx + ny * ny + nz * nz) ** 0.5
                    if nlen > 1e-12:
                        nx /= nlen
                        ny /= nlen
                        nz /= nlen
                    nxbuf[py, px] = nx
                    nybuf[py, px] = ny
                    nzbuf[py, px] = nz
                    mask[py, px] = 1
    return zbuf, nxbuf, nybuf, nzbuf, mask


@njit(cache=True)
def rasterize_scalar(verts2d, scalar, depth, faces, H, W):
    zbuf = np.full((H, W), -1e30, dtype=np.float64)
    sbuf = np.zeros((H, W), dtype=np.float64)
    mask = np.zeros((H, W), dtype=np.uint8)
    for fi in range(faces.shape[0]):
        i0, i1, i2 = faces[fi, 0], faces[fi, 1], faces[fi, 2]
        x0, y0 = verts2d[i0, 0], verts2d[i0, 1]
        x1, y1 = verts2d[i1, 0], verts2d[i1, 1]
        x2, y2 = verts2d[i2, 0], verts2d[i2, 1]
        minx = int(max(0, np.floor(min(x0, x1, x2))))
        maxx = int(min(W - 1, np.ceil(max(x0, x1, x2))))
        miny = int(max(0, np.floor(min(y0, y1, y2))))
        maxy = int(min(H - 1, np.ceil(max(y0, y1, y2))))
        denom = ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2))
        if abs(denom) < 1e-12:
            continue
        s0, s1, s2 = scalar[i0], scalar[i1], scalar[i2]
        z0, z1, z2 = depth[i0], depth[i1], depth[i2]
        for py in range(miny, maxy + 1):
            yy = py + 0.5
            for px in range(minx, maxx + 1):
                xx = px + 0.5
                w0 = ((y1 - y2) * (xx - x2) + (x2 - x1) * (yy - y2)) / denom
                if w0 < 0.0:
                    continue
                w1 = ((y2 - y0) * (xx - x2) + (x0 - x2) * (yy - y2)) / denom
                if w1 < 0.0:
                    continue
                w2 = 1.0 - w0 - w1
                if w2 < 0.0:
                    continue
                z = w0 * z0 + w1 * z1 + w2 * z2
                if z > zbuf[py, px]:
                    zbuf[py, px] = z
                    sbuf[py, px] = w0 * s0 + w1 * s1 + w2 * s2
                    mask[py, px] = 1
    return sbuf, mask


@njit(cache=True)
def rasterize_bake_correspondence(verts2d, depths, faces, H, W):
    """Rasterize the visible triangle id and barycentric weights for each pixel.

    The Z-buffer convention is identical to the rendering core: larger depth is
    closer to the observer.  This map is used only for vertex-color baking.
    """
    zbuf = np.full((H, W), -1e30, dtype=np.float64)
    facebuf = np.full((H, W), -1, dtype=np.int32)
    w0buf = np.zeros((H, W), dtype=np.float32)
    w1buf = np.zeros((H, W), dtype=np.float32)
    w2buf = np.zeros((H, W), dtype=np.float32)

    for fi in range(faces.shape[0]):
        i0, i1, i2 = faces[fi, 0], faces[fi, 1], faces[fi, 2]
        x0, y0 = verts2d[i0, 0], verts2d[i0, 1]
        x1, y1 = verts2d[i1, 0], verts2d[i1, 1]
        x2, y2 = verts2d[i2, 0], verts2d[i2, 1]
        minx = int(max(0, np.floor(min(x0, x1, x2))))
        maxx = int(min(W - 1, np.ceil(max(x0, x1, x2))))
        miny = int(max(0, np.floor(min(y0, y1, y2))))
        maxy = int(min(H - 1, np.ceil(max(y0, y1, y2))))
        denom = ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2))
        if abs(denom) < 1e-12:
            continue
        z0, z1, z2 = depths[i0], depths[i1], depths[i2]

        for py in range(miny, maxy + 1):
            yy = py + 0.5
            for px in range(minx, maxx + 1):
                xx = px + 0.5
                w0 = ((y1 - y2) * (xx - x2) + (x2 - x1) * (yy - y2)) / denom
                if w0 < 0.0:
                    continue
                w1 = ((y2 - y0) * (xx - x2) + (x0 - x2) * (yy - y2)) / denom
                if w1 < 0.0:
                    continue
                w2 = 1.0 - w0 - w1
                if w2 < 0.0:
                    continue
                z = w0 * z0 + w1 * z1 + w2 * z2
                if z > zbuf[py, px]:
                    zbuf[py, px] = z
                    facebuf[py, px] = fi
                    w0buf[py, px] = w0
                    w1buf[py, px] = w1
                    w2buf[py, px] = w2

    return facebuf, w0buf, w1buf, w2buf


@njit(cache=True)
def accumulate_bake_to_vertices(facebuf, w0buf, w1buf, w2buf, image, faces, n_vertices):
    """Back-project visible pixels to vertices by barycentric weighted averaging."""
    accum = np.zeros(n_vertices, dtype=np.float64)
    weights = np.zeros(n_vertices, dtype=np.float64)
    H, W = image.shape

    for py in range(H):
        for px in range(W):
            fi = facebuf[py, px]
            if fi < 0:
                continue
            i0, i1, i2 = faces[fi, 0], faces[fi, 1], faces[fi, 2]
            val = image[py, px]
            w0 = float(w0buf[py, px])
            w1 = float(w1buf[py, px])
            w2 = float(w2buf[py, px])
            accum[i0] += w0 * val
            accum[i1] += w1 * val
            accum[i2] += w2 * val
            weights[i0] += w0
            weights[i1] += w1
            weights[i2] += w2

    values = np.ones(n_vertices, dtype=np.float64)
    for i in range(n_vertices):
        if weights[i] > 1e-12:
            values[i] = accum[i] / weights[i]
    return values, weights


def bake_result_to_vertex_colors(mesh: MeshState, bundle: RenderBundle):
    """Bake the current 2D Result to visible vertex colors.

    The current orthographic projection is reconstructed exactly from the
    RenderParams.  A Z-buffer chooses visible triangles, and each winning pixel
    is distributed to the triangle's vertices by barycentric weights.

    Unobserved vertices preserve source vertex colors when available; if the
    input PLY had no vertex colors, they remain opaque white.
    """
    p = bundle.params
    proj_verts, depth, W, H = project_setup(mesh.vertices, p.projection, p.image_height, p.margin_ratio)
    image = np.asarray(bundle.result, dtype=np.float64)
    if image.shape != (H, W):
        raise ValueError(
            f"Render/image geometry mismatch: result={image.shape}, expected={(H, W)}. "
            "Re-render the current mesh before baking."
        )

    facebuf, w0buf, w1buf, w2buf = rasterize_bake_correspondence(
        proj_verts.astype(np.float64),
        depth.astype(np.float64),
        mesh.faces.astype(np.int64),
        H, W,
    )
    gray, weights = accumulate_bake_to_vertices(
        facebuf, w0buf, w1buf, w2buf,
        np.clip(image, 0.0, 1.0).astype(np.float64),
        mesh.faces.astype(np.int64),
        len(mesh.vertices),
    )
    observed = weights > 1e-12

    if mesh.source_vertex_colors is not None and len(mesh.source_vertex_colors) == len(mesh.vertices):
        rgba = np.asarray(mesh.source_vertex_colors, dtype=np.uint8).copy()
        if rgba.shape[1] == 3:
            rgba = np.hstack([rgba, np.full((len(rgba), 1), 255, dtype=np.uint8)])
    else:
        rgba = np.full((len(mesh.vertices), 4), 255, dtype=np.uint8)

    g8 = np.clip(gray * 255.0 + 0.5, 0, 255).astype(np.uint8)
    rgba[observed, 0] = g8[observed]
    rgba[observed, 1] = g8[observed]
    rgba[observed, 2] = g8[observed]
    rgba[observed, 3] = 255

    colored = int(np.count_nonzero(observed))
    total = int(len(mesh.vertices))
    stats = {
        "colored_vertices": colored,
        "total_vertices": total,
        "colored_ratio": (colored / total) if total else 0.0,
        "projection": p.projection,
        "image_width": W,
        "image_height": H,
        "source_image": "current composite Result",
        "hidden_vertex_policy": "preserve source vertex colors; otherwise white",
        "bake_method": "visible-pixel barycentric weighted vertex-color bake",
    }
    return rgba, stats


def export_vertex_color_ply(mesh: MeshState, rgba: np.ndarray, outpath: Path):
    rgba = np.asarray(rgba, dtype=np.uint8)
    if rgba.shape != (len(mesh.vertices), 4):
        raise ValueError(f"vertex color array must be (N,4); got {rgba.shape}")
    baked = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=np.float64).copy(),
        faces=np.asarray(mesh.faces, dtype=np.int64).copy(),
        process=False,
    )
    baked.visual.vertex_colors = rgba
    baked.export(str(outpath), file_type="ply")


def masked_gaussian(values: np.ndarray, mask: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return values.copy()
    w = mask.astype(np.float64)
    denom = gaussian_filter(w, sigma=sigma, mode="nearest")
    denom = np.maximum(denom, 1e-8)
    numer = gaussian_filter(np.where(mask, values, 0.0), sigma=sigma, mode="nearest")
    return numer / denom


def render_base(depths: np.ndarray, faces: np.ndarray, normals_view: np.ndarray, proj_verts: np.ndarray, W: int, H: int):
    _z, nx, ny, nz, masku8 = rasterize_normals(proj_verts, depths.astype(np.float64), normals_view.astype(np.float64), faces.astype(np.int64), H, W)
    mask = masku8 > 0
    dot = np.clip(nx * LIGHT_DIR[0] + ny * LIGHT_DIR[1] + nz * LIGHT_DIR[2], 0.0, 1.0)
    rx = 2 * dot * nx - LIGHT_DIR[0]
    ry = 2 * dot * ny - LIGHT_DIR[1]
    rz = 2 * dot * nz - LIGHT_DIR[2]
    spec = np.clip(rx * VIEW_DIR[0] + ry * VIEW_DIR[1] + rz * VIEW_DIR[2], 0.0, 1.0) ** 18
    img = np.full((H, W), DEFAULT_BG, dtype=np.float64)
    val = 0.20 + 0.70 * dot + 0.10 * spec
    img[mask] = np.clip(val[mask], 0.0, 1.0)
    return img, mask


def reflect_suppress(base: np.ndarray, mask: np.ndarray, p: RenderParams):
    out = base.copy()
    vals = out[mask]
    if vals.size == 0:
        return out
    lo = np.quantile(vals, p.highlight_q_low)
    hi = np.quantile(vals, p.highlight_q_high)
    if hi <= lo:
        return out
    v = out[mask]
    t = np.clip((v - lo) / (hi - lo), 0.0, 1.0)
    v2 = v - p.highlight_strength * (t ** p.highlight_gamma)
    v2 = np.minimum(v2, lo + p.highlight_cap_ratio * (v2 - lo))
    out[mask] = np.clip(v2, 0.0, 1.0)
    return out


def diffuse_values(values: np.ndarray, P: sparse.csr_matrix, steps: int, alpha: float) -> np.ndarray:
    X = np.asarray(values, dtype=np.float64).copy()
    for _ in range(int(steps)):
        X = (1.0 - alpha) * X + alpha * P.dot(X)
    return X


def rasterize_feature(faces, proj_verts, depth_orig, feature_vertex, W, H, sigma_px=0.0):
    field, masku8 = rasterize_scalar(
        proj_verts,
        np.asarray(feature_vertex, dtype=np.float64),
        np.asarray(depth_orig, dtype=np.float64),
        np.asarray(faces, dtype=np.int64),
        H, W,
    )
    mask = masku8 > 0
    field[~mask] = 0.0
    if sigma_px > 0:
        field = masked_gaussian(field, mask, sigma_px)
        field[~mask] = 0.0
    return field, mask


def render_mesh_native_feature(faces, proj_verts, depth_orig, W, H, v_fine, v_top, sigma_px):
    n_top = vertex_normals_from_vertices(v_top, faces)
    relief = np.einsum("ij,ij->i", v_fine - v_top, n_top)
    relief = relief - np.median(relief)
    return rasterize_feature(faces, proj_verts, depth_orig, relief, W, H, sigma_px)


def normalize_magnitude(arr: np.ndarray, mask: np.ndarray, q=97.0):
    out = np.zeros_like(arr, dtype=np.float64)
    vals = np.abs(arr[mask])
    if vals.size == 0:
        return out, 1.0
    qv = max(float(np.percentile(vals, q)), 1e-12)
    out[mask] = np.clip(np.abs(arr[mask]) / qv, 0.0, 1.0)
    return out, qv


def line_score_from_feature(feature: np.ndarray, mask: np.ndarray, p: RenderParams):
    """Convert any Top feature response to a common black-line representation.

    The signed response is retained for audit/analysis, while display polarity is automatic:
    large magnitude + high local coherence + high oriented continuity -> black.
    """
    H, W = feature.shape
    dim = max(H, W)
    feature = np.asarray(feature, dtype=np.float64).copy()
    feature[~mask] = 0.0
    if p.top_post_sigma_px > 0:
        feature = masked_gaussian(feature, mask, p.top_post_sigma_px)
        feature[~mask] = 0.0

    mag, qv = normalize_magnitude(feature, mask, q=97.0)
    signed_norm = np.zeros_like(feature)
    signed_norm[mask] = np.clip(feature[mask] / qv, -1.0, 1.0)

    sigma = max(0.5, (p.linearity_sigma_percent / 100.0) * dim)
    edge_margin = max(2.0, 2.5 * sigma)
    dist = distance_transform_edt(mask)
    safe = mask & (dist > edge_margin)
    if not np.any(safe):
        safe = mask.copy()

    gy, gx = np.gradient(signed_norm)
    jxx = gaussian_filter(gx * gx, sigma=sigma, mode="nearest")
    jyy = gaussian_filter(gy * gy, sigma=sigma, mode="nearest")
    jxy = gaussian_filter(gx * gy, sigma=sigma, mode="nearest")
    delta = np.sqrt(np.maximum((jxx - jyy) ** 2 + 4.0 * jxy * jxy, 0.0))
    l1 = 0.5 * (jxx + jyy + delta)
    l2 = 0.5 * (jxx + jyy - delta)
    coherence = np.zeros_like(feature)
    coherence[safe] = np.clip((l1[safe] - l2[safe]) / (l1[safe] + l2[safe] + 1e-12), 0.0, 1.0)

    # Gradient principal direction from the structure tensor; the line tangent is perpendicular.
    theta_g = 0.5 * np.arctan2(2.0 * jxy, jxx - jyy)
    tx = -np.sin(theta_g)
    ty = np.cos(theta_g)

    span = max(1.0, (p.continuity_span_percent / 100.0) * dim)
    yy, xx = np.indices(feature.shape, dtype=np.float32)
    continuity_sum = np.zeros_like(feature)
    continuity_w = np.zeros_like(feature)
    safe_f = safe.astype(np.float64)
    for offset in np.linspace(-span, span, 9):
        ys = yy + offset * ty
        xs = xx + offset * tx
        sm = map_coordinates(mag, [ys, xs], order=1, mode="constant", cval=0.0)
        sw = map_coordinates(safe_f, [ys, xs], order=0, mode="constant", cval=0.0)
        continuity_sum += sm * sw
        continuity_w += sw
    continuity = np.zeros_like(feature)
    ok = continuity_w > 0
    continuity[ok] = np.clip(continuity_sum[ok] / continuity_w[ok], 0.0, 1.0)
    continuity[~safe] = 0.0

    lin_factor = np.ones_like(feature) if p.linearity_weight == 0 else np.power(np.clip(coherence, 0.0, 1.0), p.linearity_weight)
    con_factor = np.ones_like(feature) if p.continuity_weight == 0 else np.power(np.clip(continuity, 0.0, 1.0), p.continuity_weight)
    line_score = np.clip(mag * lin_factor * con_factor, 0.0, 1.0)
    line_score[~safe] = 0.0

    top = np.ones_like(feature)
    top[mask] = 1.0 - line_score[mask]
    diagnostics = {
        "signed_feature": signed_norm,
        "magnitude": mag,
        "linearity": coherence,
        "continuity": continuity,
        "line_score": line_score,
    }
    meta = {
        "normalization_q97": qv,
        "linearity_sigma_px": sigma,
        "continuity_span_px": span,
        "linearity_weight": p.linearity_weight,
        "continuity_weight": p.continuity_weight,
    }
    return top, diagnostics, meta


def visible_depth_image(vertices: np.ndarray, faces: np.ndarray, proj_verts: np.ndarray, depth_orig: np.ndarray, W: int, H: int, projection: str):
    scalar = projection_depth(vertices, projection).astype(np.float64)
    img, masku8 = rasterize_scalar(proj_verts, scalar, depth_orig.astype(np.float64), faces.astype(np.int64), H, W)
    mask = masku8 > 0
    img[~mask] = 0.0
    return img, mask


def projection_pixel_size_mm(vertices: np.ndarray, projection: str, W: int, H: int, margin: float) -> float:
    if projection == "front":
        a, b = vertices[:, 0], vertices[:, 1]
    elif projection == "right":
        a, b = vertices[:, 1], vertices[:, 2]
    elif projection == "bottom":
        a, b = vertices[:, 0], vertices[:, 2]
    else:
        raise ValueError(projection)
    ar = max((float(a.max()) - float(a.min())) * (1.0 + 2.0 * margin), 1e-9)
    br = max((float(b.max()) - float(b.min())) * (1.0 + 2.0 * margin), 1e-9)
    sx = ar / max(W - 1, 1)
    sy = br / max(H - 1, 1)
    return float(0.5 * (sx + sy))


def render_visible_side_thickness(depth_visible: np.ndarray, mask: np.ndarray, depth_min: float, depth_max: float, contrast: float):
    z0 = 0.5 * (depth_min + depth_max)
    half = max(0.5 * (depth_max - depth_min), 1e-12)
    d = np.zeros_like(depth_visible)
    d[mask] = (depth_visible[mask] - z0) / half
    out = np.ones_like(depth_visible)
    # 0 at the virtual mid-plane -> near white; observer-side protrusion -> darker gray.
    out[mask] = np.clip(0.95 - contrast * d[mask], 0.55, 1.0)
    return out, d, {"reference_mid_depth": z0, "half_depth_range": half}


def local_surface_roughness_vertices(vertices: np.ndarray, faces: np.ndarray, P: sparse.csr_matrix, radius_mm: float, median_edge: float, alpha: float):
    """Graph-scale tangent-plane residual RMS approximation.

    A locally smoothed surface supplies the reference tangent surface. Signed normal residuals of
    the original mesh are squared and diffused over the same graph scale, then square-rooted.
    """
    steps = mm_to_iterations(radius_mm, median_edge)
    v_ref = smooth_vertices(vertices, P, steps, alpha)
    n_ref = vertex_normals_from_vertices(v_ref, faces)
    residual = np.einsum("ij,ij->i", vertices - v_ref, n_ref)
    rms2 = diffuse_values(residual * residual, P, steps, alpha)
    return np.sqrt(np.maximum(rms2, 0.0)), steps


def render_roughness_base(rough_vertex: np.ndarray, faces, proj_verts, depth_orig, W, H, contrast: float):
    field, mask = rasterize_feature(faces, proj_verts, depth_orig, rough_vertex, W, H, 0.0)
    norm, qv = normalize_magnitude(field, mask, q=97.0)
    out = np.ones_like(field)
    out[mask] = np.clip(0.95 - contrast * norm[mask], 0.55, 1.0)
    return out, norm, mask, qv


def render_depth_feature_fields(
    depth_img: np.ndarray,
    mask: np.ndarray,
    H: int,
    W: int,
    depth_med_percent: float,
    curvature_sigma_percent: float,
):
    dim = max(H, W)
    med_pct = max(float(depth_med_percent), 0.01)
    broad_pct = med_pct * PROTOCOL_DEPTH_BROAD_RATIO
    sig_med = max(0.5, (med_pct / 100.0) * dim)
    sig_broad = max(sig_med + 0.5, (broad_pct / 100.0) * dim)
    sig_curv = max(0.5, (curvature_sigma_percent / 100.0) * dim)
    d_med = masked_gaussian(depth_img, mask, sig_med)
    d_broad = masked_gaussian(depth_img, mask, sig_broad)
    med_relief = depth_img - d_med
    broad_relief = d_med - d_broad
    mb_relief = 0.65 * med_relief + 0.35 * broad_relief
    curv = laplace(masked_gaussian(depth_img, mask, sig_curv))
    mb_relief[~mask] = 0.0
    curv[~mask] = 0.0
    return mb_relief, curv, {"depth_med_percent": med_pct, "depth_broad_percent": broad_pct, "sig_med_px": sig_med, "sig_broad_px": sig_broad, "sig_curv_px": sig_curv}


def _neighbor_shift(arr: np.ndarray, dy: int, dx: int, fill_value: float):
    out = np.full_like(arr, fill_value)
    H, W = arr.shape
    y0 = max(0, -dy); y1 = min(H, H - dy)
    x0 = max(0, -dx); x1 = min(W, W - dx)
    if y1 <= y0 or x1 <= x0:
        return out
    out[y0:y1, x0:x1] = arr[y0 + dy:y1 + dy, x0 + dx:x1 + dx]
    return out


def openness_fields(depth_img: np.ndarray, mask: np.ndarray, radius_mm: float, pixel_size_mm: float, directions: int = 8):
    rpx = max(1, int(round(radius_mm / max(pixel_size_mm, 1e-9))))
    nsteps = min(12, rpx)
    radii = np.unique(np.maximum(1, np.rint(np.linspace(1, rpx, nsteps)).astype(int)))
    pos_sum = np.zeros_like(depth_img, dtype=np.float64)
    neg_sum = np.zeros_like(depth_img, dtype=np.float64)
    dir_count = np.zeros_like(depth_img, dtype=np.float64)

    for ang in np.linspace(0.0, 2.0 * np.pi, directions, endpoint=False):
        max_slope = np.full_like(depth_img, -np.inf, dtype=np.float64)
        min_slope = np.full_like(depth_img, np.inf, dtype=np.float64)
        have = np.zeros_like(mask, dtype=bool)
        ca, sa = math.cos(ang), math.sin(ang)
        seen_offsets = set()
        for r in radii:
            dx = int(round(ca * r)); dy = int(round(sa * r))
            if dx == 0 and dy == 0:
                continue
            if (dy, dx) in seen_offsets:
                continue
            seen_offsets.add((dy, dx))
            nb = _neighbor_shift(depth_img, dy, dx, np.nan)
            nm = _neighbor_shift(mask.astype(np.float64), dy, dx, 0.0) > 0.5
            valid = mask & nm & np.isfinite(nb)
            dist_mm = max(math.hypot(dx, dy) * pixel_size_mm, 1e-9)
            slope = np.zeros_like(depth_img)
            slope[valid] = (nb[valid] - depth_img[valid]) / dist_mm
            max_slope[valid] = np.maximum(max_slope[valid], slope[valid])
            min_slope[valid] = np.minimum(min_slope[valid], slope[valid])
            have |= valid
        if np.any(have):
            pos = np.zeros_like(depth_img)
            neg = np.zeros_like(depth_img)
            pos[have] = (np.pi / 2.0) - np.arctan(max_slope[have])
            neg[have] = (np.pi / 2.0) + np.arctan(min_slope[have])
            pos_sum[have] += pos[have]
            neg_sum[have] += neg[have]
            dir_count[have] += 1.0

    ok = mask & (dir_count > 0)
    positive = np.full_like(depth_img, np.pi / 2.0)
    negative = np.full_like(depth_img, np.pi / 2.0)
    positive[ok] = pos_sum[ok] / dir_count[ok]
    negative[ok] = neg_sum[ok] / dir_count[ok]
    neutral = np.pi / 2.0
    f_pos = np.zeros_like(depth_img)
    f_neg = np.zeros_like(depth_img)
    f_comb = np.zeros_like(depth_img)
    f_pos[ok] = np.maximum(positive[ok] - neutral, 0.0)
    f_neg[ok] = np.maximum(negative[ok] - neutral, 0.0)
    f_comb[ok] = 0.5 * (positive[ok] - negative[ok])
    meta = {"radius_px": rpx, "radius_mm": radius_mm, "directions": directions, "radial_samples": int(len(radii))}
    return f_pos, f_neg, f_comb, meta


def hessian_line_feature(depth_img: np.ndarray, mask: np.ndarray, sigma_mm: float, pixel_size_mm: float):
    sigma_px = max(0.5, sigma_mm / max(pixel_size_mm, 1e-9))
    z = masked_gaussian(depth_img, mask, sigma_px)
    gy, gx = np.gradient(z, pixel_size_mm, pixel_size_mm)
    gxy, gxx = np.gradient(gx, pixel_size_mm, pixel_size_mm)
    gyy, gyx = np.gradient(gy, pixel_size_mm, pixel_size_mm)
    hxy = 0.5 * (gxy + gyx)
    tr = gxx + gyy
    disc = np.sqrt(np.maximum((gxx - gyy) ** 2 + 4.0 * hxy * hxy, 0.0))
    la = 0.5 * (tr + disc)
    lb = 0.5 * (tr - disc)
    major = np.where(np.abs(la) >= np.abs(lb), la, lb)
    minor = np.where(np.abs(la) >= np.abs(lb), lb, la)
    ratio = np.abs(minor) / (np.abs(major) + 1e-12)
    line_shape = np.exp(-0.5 * (ratio / 0.5) ** 2)
    feature = major * line_shape
    safe = mask & (distance_transform_edt(mask) > max(2.0, 2.5 * sigma_px))
    feature[~safe] = 0.0
    return feature, {"sigma_px": sigma_px, "sigma_mm": sigma_mm}


def normal_variation_feature(vertices: np.ndarray, faces: np.ndarray, P: sparse.csr_matrix, radius_mm: float, median_edge: float, alpha: float):
    normals = vertex_normals_from_vertices(vertices, faces)
    steps = mm_to_iterations(radius_mm, median_edge)
    navg = diffuse_values(normals, P, steps, alpha)
    lengths = np.linalg.norm(navg, axis=1)
    lengths[lengths < 1e-12] = 1.0
    navg = navg / lengths[:, None]
    dots = np.clip(np.abs(np.einsum("ij,ij->i", normals, navg)), 0.0, 1.0)
    variation = np.arccos(dots)
    return variation, steps


def principal_curvature_feature(vertices: np.ndarray, faces: np.ndarray):
    """Cotangent-Laplacian principal-curvature estimate for experimental ridge/valley response."""
    V = np.asarray(vertices, dtype=np.float64)
    F = np.asarray(faces, dtype=np.int64)
    nV = len(V)
    i, j, k = F[:, 0], F[:, 1], F[:, 2]
    vi, vj, vk = V[i], V[j], V[k]
    eij = vj - vi; eik = vk - vi
    eji = vi - vj; ejk = vk - vj
    eki = vi - vk; ekj = vj - vk
    cross = np.cross(eij, eik)
    dblA = np.linalg.norm(cross, axis=1)
    dblA_safe = np.maximum(dblA, 1e-12)
    area = 0.5 * dblA
    cot_i = np.einsum("ij,ij->i", eij, eik) / dblA_safe
    cot_j = np.einsum("ij,ij->i", eji, ejk) / dblA_safe
    cot_k = np.einsum("ij,ij->i", eki, ekj) / dblA_safe

    rows = np.concatenate([j, k, i, k, i, j])
    cols = np.concatenate([k, j, k, i, j, i])
    data = 0.5 * np.concatenate([cot_i, cot_i, cot_j, cot_j, cot_k, cot_k])
    Wmat = sparse.csr_matrix((data, (rows, cols)), shape=(nV, nV))
    row_sum = np.asarray(Wmat.sum(axis=1)).ravel()
    area_v = np.zeros(nV, dtype=np.float64)
    np.add.at(area_v, i, area / 3.0); np.add.at(area_v, j, area / 3.0); np.add.at(area_v, k, area / 3.0)
    area_v = np.maximum(area_v, 1e-12)
    lap = (Wmat.dot(V) - row_sum[:, None] * V) / area_v[:, None]
    normals = vertex_normals_from_vertices(V, F)
    Hmean = 0.5 * np.einsum("ij,ij->i", lap, normals)

    def angle(a, b):
        den = np.maximum(np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-12)
        return np.arccos(np.clip(np.einsum("ij,ij->i", a, b) / den, -1.0, 1.0))
    ai = angle(eij, eik); aj = angle(eji, ejk); ak = angle(eki, ekj)
    asum = np.zeros(nV, dtype=np.float64)
    np.add.at(asum, i, ai); np.add.at(asum, j, aj); np.add.at(asum, k, ak)
    K = (2.0 * np.pi - asum) / area_v

    # Suppress open-boundary curvature where the angle-defect formula changes convention.
    tm = trimesh.Trimesh(vertices=V, faces=F, process=False)
    counts = np.bincount(tm.edges_unique_inverse, minlength=len(tm.edges_unique))
    bedges = tm.edges_unique[counts == 1]
    boundary = np.zeros(nV, dtype=bool)
    if len(bedges):
        boundary[np.unique(bedges)] = True

    root = np.sqrt(np.maximum(Hmean * Hmean - K, 0.0))
    k1 = Hmean + root
    k2 = Hmean - root
    major = np.where(np.abs(k1) >= np.abs(k2), k1, k2)
    minor = np.where(np.abs(k1) >= np.abs(k2), k2, k1)
    ratio = np.abs(minor) / (np.abs(major) + 1e-12)
    line_shape = np.exp(-0.5 * (ratio / 0.5) ** 2)
    feature = major * line_shape
    feature[boundary] = 0.0
    return feature, {"boundary_vertices_suppressed": int(boundary.sum())}


def blend_pair(base: np.ndarray, top: np.ndarray, mask: np.ndarray, mode: str, alpha: float):
    out = np.ones_like(base)
    if mode == "multiply":
        comp = base * top
    elif mode == "alpha":
        comp = (1.0 - alpha) * base + alpha * top
    elif mode == "overlay":
        comp = np.where(base < 0.5, 2 * base * top, 1 - 2 * (1 - base) * (1 - top))
    elif mode == "screen":
        comp = 1.0 - (1.0 - base) * (1.0 - top)
    elif mode == "add":
        comp = np.clip(base + alpha * top, 0.0, 1.0)
    else:
        raise ValueError(mode)
    out[mask] = np.clip(comp[mask], 0.0, 1.0)
    return out


def select_top_names(p: RenderParams) -> List[str]:
    names: List[str] = []
    for name in (p.top1_choice, p.top2_choice):
        if name == "None":
            continue
        if name not in names:
            names.append(name)
    return names[:2]


def compute_render_bundle(mesh: MeshState, p: RenderParams, progress_cb) -> RenderBundle:
    progress_cb(5, "Projection setup…")
    proj_verts, depth_orig, W, H = project_setup(mesh.vertices, p.projection, p.image_height, p.margin_ratio)
    pixel_size_mm = projection_pixel_size_mm(mesh.vertices, p.projection, W, H, p.margin_ratio)
    it_base = mm_to_iterations(p.base_scale_mm, mesh.median_edge)
    it_broad = mm_to_iterations(p.broad_base_scale_mm, mesh.median_edge)
    it_fine = mm_to_iterations(p.fine_scale_mm, mesh.median_edge)
    selected_top_names = select_top_names(p)

    progress_cb(15, "Computing base smoothing…")
    v_base = smooth_vertices(mesh.vertices, mesh.operator, it_base, p.smoothing_alpha)
    v_broad = smooth_vertices(mesh.vertices, mesh.operator, it_broad, p.smoothing_alpha)

    progress_cb(30, "Rendering Base methods…")
    n_base = normals_to_view(vertex_normals_from_vertices(v_base, mesh.faces), p.projection)
    base_raw, base_mask = render_base(projection_depth(v_base, p.projection), mesh.faces, n_base, proj_verts, W, H)
    base_ref = reflect_suppress(base_raw, base_mask, p)

    n_broad = normals_to_view(vertex_normals_from_vertices(v_broad, mesh.faces), p.projection)
    broad_raw, _broad_mask = render_base(projection_depth(v_broad, p.projection), mesh.faces, n_broad, proj_verts, W, H)
    broad_ref = reflect_suppress(broad_raw, base_mask, p)

    base_candidates: Dict[str, np.ndarray] = {
        "Medium-smoothed shade": base_raw,
        "Base reflect-suppressed": base_ref,
        "Broad-base reflect-suppressed": broad_ref,
    }
    base_meta: Dict[str, Any] = {}

    # Visible-side Thickness Proxy: visible surface only; no hidden/back surface is queried.
    depth_visible_orig, visible_mask = visible_depth_image(mesh.vertices, mesh.faces, proj_verts, depth_orig, W, H, p.projection)
    depth_axis = projection_depth(mesh.vertices, p.projection)
    if p.base_choice == "Visible-side Thickness Proxy":
        if p.projection != "front":
            raise ValueError(
                "Visible-side Thickness Proxy is defined on the normalized Z thickness axis; "
                "select projection='front' for this Base method."
            )
        # In the normalized archaeological model, Z is the thickness axis.
        thickness_img, thickness_field, thickness_meta = render_visible_side_thickness(
            depth_visible_orig, visible_mask, float(mesh.vertices[:, 2].min()), float(mesh.vertices[:, 2].max()), p.thickness_contrast
        )
        base_candidates["Visible-side Thickness Proxy"] = thickness_img
        base_meta["Visible-side Thickness Proxy"] = thickness_meta

    if p.base_choice == "Local Surface Roughness":
        progress_cb(38, "Computing local surface roughness…")
        rough_v, rough_steps = local_surface_roughness_vertices(
            mesh.vertices, mesh.faces, mesh.operator, p.roughness_radius_mm, mesh.median_edge, p.smoothing_alpha
        )
        rough_img, rough_norm, rough_mask, rough_q = render_roughness_base(
            rough_v, mesh.faces, proj_verts, depth_orig, W, H, p.roughness_contrast
        )
        base_candidates["Local Surface Roughness"] = rough_img
        base_meta["Local Surface Roughness"] = {
            "radius_mm": p.roughness_radius_mm,
            "graph_diffusion_steps": rough_steps,
            "normal_residual_rms_q97_mm": rough_q,
            "implementation": "graph-diffused tangent-plane normal-residual RMS approximation",
        }

    selected_base_name = p.base_choice
    selected_base = base_candidates[selected_base_name]

    progress_cb(45, "Preparing Top feature fields…")
    top_candidates: Dict[str, np.ndarray] = {}
    top_diagnostics: Dict[str, Dict[str, np.ndarray]] = {}
    top_meta: Dict[str, Any] = {}
    feature_fields: Dict[str, Tuple[np.ndarray, np.ndarray, Dict[str, Any]]] = {}

    # Mesh-native medium relief. In v0.7.3 its medium scale is a sweepable
    # parameter rather than three fixed M8/M12/M16 UI methods.
    mesh_scales = {
        "Mesh-native medium relief": p.mesh_top_scale_mm,
        # Legacy names remain readable for older saved parameter sets.
        "Mesh-native top M8": 0.80,
        "Mesh-native top M12": 1.20,
        "Mesh-native top M16": 1.60,
    }
    needed_mesh = [n for n in selected_top_names if n in mesh_scales]
    if needed_mesh:
        v_fine = smooth_vertices(mesh.vertices, mesh.operator, it_fine, p.smoothing_alpha)
        for name in needed_mesh:
            scale = float(mesh_scales[name])
            v_top = smooth_vertices(
                mesh.vertices, mesh.operator,
                mm_to_iterations(scale, mesh.median_edge), p.smoothing_alpha
            )
            field, fmask = render_mesh_native_feature(
                mesh.faces, proj_verts, depth_orig, W, H, v_fine, v_top, 0.0
            )
            feature_fields[name] = (
                field, fmask,
                {"fine_mm": p.fine_scale_mm, "medium_mm": scale}
            )

    depth_based_names = {
        "Depth multi-scale relief", "Curvature proxy", "Positive Openness", "Negative Openness",
        "Combined Openness", "Hessian Line Response"
    }
    need_depth = any(n in depth_based_names for n in selected_top_names)
    depth_base_img = None
    dmask = None
    depth_meta = {}
    if need_depth:
        depth_base_img, dmask = visible_depth_image(v_base, mesh.faces, proj_verts, depth_orig, W, H, p.projection)

    if need_depth and any(n in selected_top_names for n in ("Depth multi-scale relief", "Curvature proxy")):
        mb, curv, depth_meta = render_depth_feature_fields(depth_base_img, dmask, H, W, p.depth_med_percent, p.curvature_sigma_percent)
        if "Depth multi-scale relief" in selected_top_names:
            feature_fields["Depth multi-scale relief"] = (mb, dmask, {"type": "multi-scale depth", **depth_meta})
        if "Curvature proxy" in selected_top_names:
            feature_fields["Curvature proxy"] = (curv, dmask, {"type": "depth Laplacian", **depth_meta, "sigma_percent": p.curvature_sigma_percent})

    if need_depth and any(n in selected_top_names for n in ("Positive Openness", "Negative Openness", "Combined Openness")):
        fpos, fneg, fcomb, ometa = openness_fields(depth_base_img, dmask, p.openness_radius_mm, pixel_size_mm, directions=8)
        if "Positive Openness" in selected_top_names:
            feature_fields["Positive Openness"] = (fpos, dmask, {"type": "positive openness", **ometa})
        if "Negative Openness" in selected_top_names:
            feature_fields["Negative Openness"] = (fneg, dmask, {"type": "negative openness", **ometa})
        if "Combined Openness" in selected_top_names:
            feature_fields["Combined Openness"] = (fcomb, dmask, {"type": "signed combined openness", **ometa})

    if need_depth and "Hessian Line Response" in selected_top_names:
        hf, hmeta = hessian_line_feature(depth_base_img, dmask, p.hessian_sigma_mm, pixel_size_mm)
        feature_fields["Hessian Line Response"] = (hf, dmask, {"type": "depth Hessian ridge/valley", **hmeta})

    if "Normal Variation" in selected_top_names:
        progress_cb(55, "Computing mesh normal variation…")
        nv, nsteps = normal_variation_feature(mesh.vertices, mesh.faces, mesh.operator, p.normal_radius_mm, mesh.median_edge, p.smoothing_alpha)
        field, fmask = rasterize_feature(mesh.faces, proj_verts, depth_orig, nv, W, H, 0.0)
        feature_fields["Normal Variation"] = (field, fmask, {"radius_mm": p.normal_radius_mm, "graph_diffusion_steps": nsteps})

    if "Principal Curvature Ridge / Valley" in selected_top_names:
        progress_cb(60, "Computing experimental principal curvature…")
        pc_steps = mm_to_iterations(p.principal_scale_mm, mesh.median_edge)
        v_pc = smooth_vertices(mesh.vertices, mesh.operator, pc_steps, p.smoothing_alpha)
        pc, pcmeta = principal_curvature_feature(v_pc, mesh.faces)
        field, fmask = rasterize_feature(mesh.faces, proj_verts, depth_orig, pc, W, H, 0.0)
        feature_fields["Principal Curvature Ridge / Valley"] = (
            field, fmask, {"scale_mm": p.principal_scale_mm, "smoothing_steps": pc_steps, **pcmeta}
        )

    progress_cb(72, "Evaluating linearity and continuity…")
    for name in selected_top_names:
        field, fmask, fmeta = feature_fields[name]
        top_img, diag, eval_meta = line_score_from_feature(field, fmask, p)
        top_candidates[name] = top_img
        top_diagnostics[name] = diag
        top_meta[name] = {"feature": fmeta, "evaluation": eval_meta}

    progress_cb(88, "Compositing selected layers…")
    result = selected_base.copy()
    result_mask = base_mask.copy()
    for tname in selected_top_names:
        result = blend_pair(result, top_candidates[tname], result_mask, p.blend_mode, p.alpha_blend)

    audit_text = (
        f"{APP_NAME}\n\n"
        f"Input\n"
        f"- file: {mesh.path}\n"
        f"- sha256: {mesh.sha256}\n"
        f"- projection: {p.projection}\n"
        f"- image size: {W} x {H}\n"
        f"- projected pixel size approx mm: {pixel_size_mm:.6f}\n\n"
        f"Base\n"
        f"- method: {selected_base_name}\n"
        f"- base_scale_mm: {p.base_scale_mm}\n"
        f"- broad_base_scale_mm: {p.broad_base_scale_mm}\n"
        f"- visible-side thickness contrast: {p.thickness_contrast}\n"
        f"- roughness_radius_mm: {p.roughness_radius_mm}\n"
        f"- roughness_contrast: {p.roughness_contrast}\n"
        f"- thickness proxy reference axis: normalized Z (front projection only)\n"
        f"- hidden/back surface used for thickness proxy: False\n"
        f"- depth-proxy matte base included: False\n\n"
        f"Top\n"
        f"- top1 method: {p.top1_choice}\n"
        f"- top2 method: {p.top2_choice}\n"
        f"- selected tops used: {selected_top_names if selected_top_names else 'None'}\n"
        f"- manual Raw/Invert polarity: removed\n"
        f"- display polarity: automatic, larger Line Score -> darker/black\n"
        f"- fine_scale_mm: {p.fine_scale_mm}\n"
        f"- mesh-native medium scale mm: {p.mesh_top_scale_mm}\n"
        f"- depth multi-scale medium sigma percent: {p.depth_med_percent}\n"
        f"- top antinoise sigma px: {p.top_post_sigma_px}\n"
        f"- curvature sigma percent: {p.curvature_sigma_percent}  [baseline=6.0%]\n"
        f"- openness radius mm: {p.openness_radius_mm}\n"
        f"- Hessian sigma mm: {p.hessian_sigma_mm}\n"
        f"- normal variation radius mm: {p.normal_radius_mm}\n"
        f"- principal curvature scale mm: {p.principal_scale_mm}\n\n"
        f"Common line evaluation\n"
        f"- magnitude normalization: abs(feature) / Q97\n"
        f"- linearity sigma percent: {p.linearity_sigma_percent}\n"
        f"- continuity span percent: {p.continuity_span_percent}\n"
        f"- linearity weight beta: {p.linearity_weight}\n"
        f"- continuity weight gamma: {p.continuity_weight}\n"
        f"- line score: magnitude * linearity^beta * continuity^gamma\n\n"
        f"Reflect suppression\n"
        f"- q_low: {p.highlight_q_low}\n"
        f"- q_high: {p.highlight_q_high}\n"
        f"- strength: {p.highlight_strength}\n"
        f"- gamma: {p.highlight_gamma}\n"
        f"- cap_ratio: {p.highlight_cap_ratio}\n\n"
        f"Blend\n"
        f"- mode: {p.blend_mode}\n"
        f"- alpha_blend: {p.alpha_blend}\n"
        f"- max overlay layers: 3 (Base + Top1 + Top2)\n\n"
        f"Derived\n"
        f"- median_edge_mm: {mesh.median_edge:.6f}\n"
        f"- base_iterations: {it_base}\n"
        f"- broad_base_iterations: {it_broad}\n"
        f"- mesh-native medium scale mm: {p.mesh_top_scale_mm}\n"
        f"- base method metadata: {base_meta.get(selected_base_name, {})}\n"
        f"- top method metadata: {top_meta}\n"
        f"- all outputs recomputed from source PLY in this run: True\n"
    )

    progress_cb(100, "Done")
    return RenderBundle(
        params=p,
        base_candidates=base_candidates,
        top_candidates=top_candidates,
        top_diagnostics=top_diagnostics,
        selected_base_name=selected_base_name,
        selected_top_names=selected_top_names,
        result=result,
        mask=result_mask,
        meta={"width_px": W, "height_px": H, "pixel_size_mm": pixel_size_mm, "base_meta": base_meta, "top_meta": top_meta},
        audit_text=audit_text,
    )


def to_rgb(gray: np.ndarray) -> Image.Image:
    arr = np.clip(gray * 255.0 + 0.5, 0, 255).astype(np.uint8)
    return Image.fromarray(arr, mode="L").convert("RGB")


def array_to_qimage(gray: np.ndarray) -> QImage:
    arr = np.ascontiguousarray(np.clip(gray * 255.0 + 0.5, 0, 255).astype(np.uint8))
    h, w = arr.shape
    qimg = QImage(arr.data, w, h, arr.strides[0], QImage.Format.Format_Grayscale8)
    return qimg.copy()


def array_to_pixmap(gray: np.ndarray) -> QPixmap:
    return QPixmap.fromImage(array_to_qimage(gray))


def make_candidate_comparison(candidates: List[RenderBundle]) -> Image.Image:
    try:
        font_title = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 18)
        font_text = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 13)
    except Exception:
        font_title = ImageFont.load_default()
        font_text = ImageFont.load_default()

    if not candidates:
        return Image.new("RGB", (600, 200), "white")

    images = [to_rgb(c.result) for c in candidates]
    w = max(im.width for im in images)
    h = max(im.height for im in images)
    pad = 16
    text_h = 70
    title_h = 36
    total_w = len(images) * (w + pad) + pad
    total_h = title_h + text_h + h + 2 * pad
    canvas = Image.new("RGB", (total_w, total_h), "white")
    dr = ImageDraw.Draw(canvas)
    dr.text((pad, 8), "Accumulated candidates comparison", fill="black", font=font_title)
    for i, (im, cand) in enumerate(zip(images, candidates), 1):
        x = pad + (i - 1) * (w + pad)
        dr.text((x, 42), f"Candidate {i}", fill="black", font=font_title)
        dr.multiline_text((x, 64), cand.candidate_label, fill=(30, 30, 30), font=font_text, spacing=2)
        canvas.paste(im, (x, title_h + text_h))
    return canvas


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())



# ---------------- Lab UI / workflow v0.7.2 ----------------

METHOD_ABBR = {
    "Medium-smoothed shade": "BMS",
    "Base reflect-suppressed": "BRS",
    "Broad-base reflect-suppressed": "BBR",
    "Visible-side Thickness Proxy": "VTP",
    "Local Surface Roughness": "LSR",
    "Mesh-native medium relief": "MNR",
    "Mesh-native top M8": "M8",       # legacy
    "Mesh-native top M12": "M12",     # legacy
    "Mesh-native top M16": "M16",     # legacy
    "Depth multi-scale relief": "DMR",
    "Curvature proxy": "CPX",
    "Positive Openness": "POP",
    "Negative Openness": "NOP",
    "Combined Openness": "COP",
    "Hessian Line Response": "HLR",
    "Normal Variation": "NVR",
    "Principal Curvature Ridge / Valley": "PCR",
}

BASE_METHODS = [
    "Medium-smoothed shade",
    "Base reflect-suppressed",
    "Broad-base reflect-suppressed",
    "Visible-side Thickness Proxy",
    "Local Surface Roughness",
]

TOP_METHODS = [
    "Mesh-native medium relief",
    "Depth multi-scale relief",
    "Curvature proxy",
    "Positive Openness",
    "Negative Openness",
    "Combined Openness",
    "Hessian Line Response",
    "Normal Variation",
    "Principal Curvature Ridge / Valley",
]


def make_spin(mn, mx, step, value, decimals=3, suffix=""):
    w = QDoubleSpinBox()
    w.setRange(float(mn), float(mx))
    w.setSingleStep(float(step))
    w.setDecimals(decimals)
    w.setValue(float(value))
    w.setSuffix(suffix)
    return w


def short_value(v: float) -> str:
    return f"{float(v):g}".replace("-", "m").replace(".", "p")


def method_filename_params(method: str, p: RenderParams) -> List[Tuple[str, float]]:
    common_top = [
        ("ts", p.top_post_sigma_px),
        ("ls", p.linearity_sigma_percent),
        ("cs", p.continuity_span_percent),
        ("lw", p.linearity_weight),
        ("cw", p.continuity_weight),
    ]
    if method == "Medium-smoothed shade":
        return [("bs", p.base_scale_mm), ("a", p.smoothing_alpha)]
    if method == "Base reflect-suppressed":
        return [
            ("bs", p.base_scale_mm), ("a", p.smoothing_alpha),
            ("ql", p.highlight_q_low), ("qh", p.highlight_q_high),
            ("hs", p.highlight_strength), ("hg", p.highlight_gamma),
            ("hc", p.highlight_cap_ratio),
        ]
    if method == "Broad-base reflect-suppressed":
        return [
            ("bb", p.broad_base_scale_mm), ("a", p.smoothing_alpha),
            ("ql", p.highlight_q_low), ("qh", p.highlight_q_high),
            ("hs", p.highlight_strength), ("hg", p.highlight_gamma),
            ("hc", p.highlight_cap_ratio),
        ]
    if method == "Visible-side Thickness Proxy":
        return [("tc", p.thickness_contrast)]
    if method == "Local Surface Roughness":
        return [("rr", p.roughness_radius_mm), ("rc", p.roughness_contrast), ("a", p.smoothing_alpha)]
    if method == "Mesh-native medium relief":
        return [("ms", p.mesh_top_scale_mm), ("fs", p.fine_scale_mm)] + common_top
    if method.startswith("Mesh-native top"):
        return [("fs", p.fine_scale_mm)] + common_top
    if method == "Curvature proxy":
        return [("cv", p.curvature_sigma_percent)] + common_top
    if "Openness" in method:
        return [("or", p.openness_radius_mm)] + common_top
    if method == "Hessian Line Response":
        return [("hs", p.hessian_sigma_mm)] + common_top
    if method == "Normal Variation":
        return [("nr", p.normal_radius_mm), ("a", p.smoothing_alpha)] + common_top
    if method == "Principal Curvature Ridge / Valley":
        return [("ps", p.principal_scale_mm), ("a", p.smoothing_alpha)] + common_top
    if method == "Depth multi-scale relief":
        return [("dm", p.depth_med_percent)] + common_top
    return common_top


def method_file_token(method: str, p: RenderParams) -> str:
    abbr = METHOD_ABBR.get(method, "MTH")
    params = "_".join(f"{k}{short_value(v)}" for k, v in method_filename_params(method, p))
    return f"{abbr}_{params}" if params else abbr


def annotated_render_image(image: np.ndarray, title: str, lines: List[str]) -> Image.Image:
    base = to_rgb(image)
    pad = 14
    try:
        font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 17)
        font_text = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 13)
    except Exception:
        font_title = ImageFont.load_default()
        font_text = ImageFont.load_default()
    line_h = 18
    footer_h = pad * 2 + 24 + line_h * len(lines)
    out = Image.new("RGB", (base.width, base.height + footer_h), "white")
    out.paste(base, (0, 0))
    draw = ImageDraw.Draw(out)
    y = base.height + pad
    draw.text((pad, y), title, fill="black", font=font_title)
    y += 24
    for line in lines:
        draw.text((pad, y), line, fill="black", font=font_text)
        y += line_h
    return out


class PreviewLabel(QLabel):
    def __init__(self, text="", parent=None, minimum=(320, 360)):
        super().__init__(text, parent)
        self._source: Optional[QPixmap] = None
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setMinimumSize(int(minimum[0]), int(minimum[1]))
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def set_array(self, image: np.ndarray):
        self._source = array_to_pixmap(image)
        self._refit()

    def clear_preview(self, text=""):
        self._source = None
        self.clear()
        if text:
            self.setText(text)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._refit()

    def _refit(self):
        if self._source is None or self._source.isNull():
            return
        self.setPixmap(
            self._source.scaled(
                self.contentsRect().size(),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
        )


@dataclass
class LabRecord:
    lab: str
    method: str
    params: RenderParams
    image: np.ndarray
    mask: np.ndarray
    audit_text: str
    meta: Dict[str, Any]

    @property
    def label(self) -> str:
        return f"{METHOD_ABBR.get(self.method, 'MTH')} | " + ", ".join(
            f"{k}={v:g}" for k, v in method_filename_params(self.method, self.params)
        )


def load_mesh_state_for_lab(path: Path) -> MeshState:
    tm = trimesh.load(path, force="mesh", process=False)
    if not isinstance(tm, trimesh.Trimesh):
        raise ValueError("Loaded object is not a single triangular mesh.")
    vertices = np.asarray(tm.vertices, dtype=np.float64).copy()
    faces = np.asarray(tm.faces, dtype=np.int64).copy()
    if len(vertices) == 0 or len(faces) == 0:
        raise ValueError("PLY contains no usable mesh vertices/faces.")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("Triangular faces are required.")
    edges_unique = tm.edges_unique.astype(np.int64)
    operator, median_edge = build_operator(vertices, edges_unique)
    source_colors = None
    try:
        vc = np.asarray(tm.visual.vertex_colors)
        if vc.ndim == 2 and len(vc) == len(vertices) and vc.shape[1] in (3, 4):
            vc = vc.astype(np.uint8, copy=True)
            if vc.shape[1] == 3:
                vc = np.hstack([vc, np.full((len(vc), 1), 255, dtype=np.uint8)])
            source_colors = vc
    except Exception:
        source_colors = None
    try:
        source_normals = np.asarray(tm.vertex_normals, dtype=np.float64).copy()
        if source_normals.shape != vertices.shape:
            source_normals = vertex_normals_from_vertices(vertices, faces)
    except Exception:
        source_normals = vertex_normals_from_vertices(vertices, faces)
    return MeshState(
        path=path,
        sha256=sha256_file(path),
        vertices=vertices,
        faces=faces,
        edges_unique=edges_unique,
        bounds=np.asarray(tm.bounds, dtype=np.float64),
        median_edge=median_edge,
        operator=operator,
        source_vertex_colors=source_colors,
        source_vertex_normals=source_normals,
    )


def bake_array_to_vertex_colors(mesh: MeshState, image: np.ndarray, p: RenderParams):
    proj_verts, depth, W, H = project_setup(mesh.vertices, p.projection, p.image_height, p.margin_ratio)
    image = np.asarray(image, dtype=np.float64)
    if image.shape != (H, W):
        raise ValueError(f"Render/image geometry mismatch: image={image.shape}, expected={(H, W)}")
    facebuf, w0buf, w1buf, w2buf = rasterize_bake_correspondence(
        proj_verts.astype(np.float64), depth.astype(np.float64), mesh.faces.astype(np.int64), H, W
    )
    gray, weights = accumulate_bake_to_vertices(
        facebuf, w0buf, w1buf, w2buf, np.clip(image, 0.0, 1.0),
        mesh.faces.astype(np.int64), len(mesh.vertices)
    )
    observed = weights > 1e-12
    if mesh.source_vertex_colors is not None and len(mesh.source_vertex_colors) == len(mesh.vertices):
        rgba = np.asarray(mesh.source_vertex_colors, dtype=np.uint8).copy()
    else:
        rgba = np.full((len(mesh.vertices), 4), 255, dtype=np.uint8)
    g8 = np.clip(gray * 255.0 + 0.5, 0, 255).astype(np.uint8)
    rgba[observed, 0:3] = g8[observed, None]
    rgba[observed, 3] = 255
    return rgba, {
        "colored_vertices": int(np.count_nonzero(observed)),
        "total_vertices": int(len(mesh.vertices)),
        "colored_ratio": float(np.count_nonzero(observed) / max(len(mesh.vertices), 1)),
        "projection": p.projection,
        "bake_method": "visible-pixel barycentric weighted vertex-color bake",
    }


def export_baked_vertex_color_ply(mesh: MeshState, rgba: np.ndarray, outpath: Path):
    """Write XYZ + normals + RGBA, deliberately omitting active UV s/t.

    Standard s/t properties can make trimesh and some viewers prioritize
    TextureVisuals, hiding vertex colors. The baked PLY is therefore a true
    vertex-color product. The original source file is left unchanged and keeps
    its UV data.
    """
    v = np.asarray(mesh.vertices, dtype=np.float64)
    f = np.asarray(mesh.faces, dtype=np.int64)
    n = mesh.source_vertex_normals
    if n is None or np.asarray(n).shape != v.shape:
        n = vertex_normals_from_vertices(v, f)
    n = np.asarray(n, dtype=np.float64)
    rgba = np.asarray(rgba, dtype=np.uint8)
    if rgba.shape != (len(v), 4):
        raise ValueError(f"RGBA must be (N,4), got {rgba.shape}")

    dtype_v = np.dtype([
        ("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
        ("nx", "<f4"), ("ny", "<f4"), ("nz", "<f4"),
        ("red", "u1"), ("green", "u1"), ("blue", "u1"), ("alpha", "u1"),
    ], align=False)
    vd = np.empty(len(v), dtype=dtype_v)
    for key, values in (
        ("x", v[:, 0]), ("y", v[:, 1]), ("z", v[:, 2]),
        ("nx", n[:, 0]), ("ny", n[:, 1]), ("nz", n[:, 2]),
    ):
        vd[key] = values.astype(np.float32)
    vd["red"], vd["green"], vd["blue"], vd["alpha"] = rgba.T

    dtype_f = np.dtype([("count", "u1"), ("indices", "<i4", (3,))], align=False)
    fd = np.empty(len(f), dtype=dtype_f)
    fd["count"] = 3
    fd["indices"] = f.astype(np.int32)

    header = "\n".join([
        "ply",
        "format binary_little_endian 1.0",
        f"comment {APP_NAME} vertex-color bake",
        "comment source_uv_omitted_to_prioritize_vertex_color 1",
        f"element vertex {len(v)}",
        "property float x", "property float y", "property float z",
        "property float nx", "property float ny", "property float nz",
        "property uchar red", "property uchar green", "property uchar blue", "property uchar alpha",
        f"element face {len(f)}",
        "property list uchar int vertex_indices",
        "end_header", "",
    ]).encode("ascii")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "wb") as fh:
        fh.write(header)
        vd.tofile(fh)
        fd.tofile(fh)


def verify_baked_ply(path: Path, expected_vertices: int, expected_faces: int) -> str:
    with open(path, "rb") as fh:
        head = fh.read(65536)
    end = head.find(b"end_header\n")
    if end < 0:
        raise ValueError("Written PLY has no end_header")
    text = head[: end + len(b"end_header\n")].decode("ascii", errors="replace")
    for prop in (
        "property float nx", "property float ny", "property float nz",
        "property uchar red", "property uchar green", "property uchar blue", "property uchar alpha",
    ):
        if prop not in text:
            raise ValueError(f"Missing required baked PLY property: {prop}")
    if "property double s" in text or "property double t" in text:
        raise ValueError("Baked vertex-color PLY unexpectedly contains active UV s/t")
    check = trimesh.load(path, force="mesh", process=False)
    if len(check.vertices) != expected_vertices or len(check.faces) != expected_faces:
        raise ValueError("Written PLY vertex/face counts changed")
    vc = np.asarray(check.visual.vertex_colors)
    if vc.shape != (expected_vertices, 4):
        raise ValueError(f"RGBA did not reload as vertex colors: {vc.shape}")
    return "normals + RGBA verified; RGBA reloads as ColorVisuals"


def normalize_for_composite(image: np.ndarray, mask: np.ndarray, qlow: float, qhigh: float) -> np.ndarray:
    out = np.ones_like(image, dtype=np.float64)
    vals = np.asarray(image, dtype=np.float64)[mask]
    if vals.size == 0:
        return np.asarray(image, dtype=np.float64).copy()
    lo = float(np.percentile(vals, qlow))
    hi = float(np.percentile(vals, qhigh))
    if hi <= lo + 1e-12:
        out[mask] = np.clip(vals, 0.0, 1.0)
    else:
        out[mask] = np.clip((vals - lo) / (hi - lo), 0.0, 1.0)
    return out


def contrast_adjust(image: np.ndarray, mask: np.ndarray, contrast: float) -> np.ndarray:
    out = np.asarray(image, dtype=np.float64).copy()
    out[mask] = np.clip((out[mask] - 0.5) * contrast + 0.5, 0.0, 1.0)
    return out


def blend_with_strength(base: np.ndarray, top: np.ndarray, mask: np.ndarray, mode: str, alpha: float) -> np.ndarray:
    b = np.asarray(base, dtype=np.float64)
    t = np.asarray(top, dtype=np.float64)
    if mode == "multiply":
        raw = b * t
    elif mode == "alpha":
        raw = t
    elif mode == "overlay":
        raw = np.where(b < 0.5, 2 * b * t, 1 - 2 * (1 - b) * (1 - t))
    elif mode == "screen":
        raw = 1.0 - (1.0 - b) * (1.0 - t)
    elif mode == "add":
        raw = np.clip(b + t, 0.0, 1.0)
    else:
        raise ValueError(mode)
    out = b.copy()
    out[mask] = np.clip((1.0 - alpha) * b[mask] + alpha * raw[mask], 0.0, 1.0)
    return out


# ---------------------------------------------------------------------------
# v0.7.3 five-variation Lab workflow
#
# Selecting a method creates five variants:
#   center-2step, center-step, center, center+step, center+2step
#
# The user selects one preview, optionally retains it for Composite Lab,
# and may use Customize to alter center/step and method-specific auxiliary
# parameters before re-rendering another five-variation set.
# ---------------------------------------------------------------------------

BASE_SWEEP_CONFIG = {
    "Medium-smoothed shade": {
        "primary": "base_scale_mm", "label": "Smoothing scale",
        "baseline": 0.80, "step": 0.15, "min": 0.10, "max": 1.50,
        "decimals": 2, "suffix": " mm",
        "aux": ["smoothing_alpha"],
    },
    "Base reflect-suppressed": {
        "primary": "highlight_strength", "label": "Reflect suppression strength",
        "baseline": 0.22, "step": 0.05, "min": 0.00, "max": 0.44,
        "decimals": 2, "suffix": "",
        "aux": [
            "base_scale_mm", "smoothing_alpha", "highlight_q_low",
            "highlight_q_high", "highlight_gamma", "highlight_cap_ratio",
        ],
    },
    "Broad-base reflect-suppressed": {
        "primary": "broad_base_scale_mm", "label": "Broad smoothing scale",
        "baseline": 1.60, "step": 0.30, "min": 0.80, "max": 3.00,
        "decimals": 2, "suffix": " mm",
        "aux": [
            "smoothing_alpha", "highlight_q_low", "highlight_q_high",
            "highlight_strength", "highlight_gamma", "highlight_cap_ratio",
        ],
    },
    "Visible-side Thickness Proxy": {
        "primary": "thickness_contrast", "label": "Thickness contrast",
        "baseline": 0.30, "step": 0.05, "min": 0.10, "max": 0.45,
        "decimals": 2, "suffix": "",
        "aux": [],
    },
    "Local Surface Roughness": {
        "primary": "roughness_radius_mm", "label": "Roughness radius",
        "baseline": 1.00, "step": 0.25, "min": 0.25, "max": 5.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["roughness_contrast", "smoothing_alpha"],
    },
}

TOP_COMMON_AUX = [
    "top_post_sigma_px",
    "linearity_sigma_percent", "continuity_span_percent",
    "linearity_weight", "continuity_weight",
]

TOP_SWEEP_CONFIG = {
    "Mesh-native medium relief": {
        "primary": "mesh_top_scale_mm", "label": "Medium relief scale",
        "baseline": 1.20, "step": 0.20, "min": 0.40, "max": 3.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["fine_scale_mm", "smoothing_alpha"] + TOP_COMMON_AUX,
    },
    "Depth multi-scale relief": {
        "primary": "depth_med_percent", "label": "Medium depth sigma",
        "baseline": 0.90, "step": 0.20, "min": 0.20, "max": 3.00,
        "decimals": 2, "suffix": " %",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Curvature proxy": {
        "primary": "curvature_sigma_percent", "label": "Curvature sigma",
        "baseline": 6.00, "step": 1.50, "min": 0.50, "max": 12.00,
        "decimals": 1, "suffix": " %",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Positive Openness": {
        "primary": "openness_radius_mm", "label": "Openness radius",
        "baseline": 2.00, "step": 0.50, "min": 0.25, "max": 10.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Negative Openness": {
        "primary": "openness_radius_mm", "label": "Openness radius",
        "baseline": 2.00, "step": 0.50, "min": 0.25, "max": 10.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Combined Openness": {
        "primary": "openness_radius_mm", "label": "Openness radius",
        "baseline": 2.00, "step": 0.50, "min": 0.25, "max": 10.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Hessian Line Response": {
        "primary": "hessian_sigma_mm", "label": "Hessian sigma",
        "baseline": 1.00, "step": 0.25, "min": 0.10, "max": 5.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["base_scale_mm"] + TOP_COMMON_AUX,
    },
    "Normal Variation": {
        "primary": "normal_radius_mm", "label": "Normal neighborhood radius",
        "baseline": 1.00, "step": 0.25, "min": 0.10, "max": 5.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["smoothing_alpha"] + TOP_COMMON_AUX,
    },
    "Principal Curvature Ridge / Valley": {
        "primary": "principal_scale_mm", "label": "Principal curvature scale",
        "baseline": 0.80, "step": 0.15, "min": 0.10, "max": 3.00,
        "decimals": 2, "suffix": " mm",
        "aux": ["smoothing_alpha"] + TOP_COMMON_AUX,
    },
}

# Parameters exposed inside Customize. The primary parameter of the selected
# method is controlled by the Center / Step widgets and is hidden from the
# auxiliary list to avoid duplicate controls.
BASE_PARAM_SPECS = {
    "base_scale_mm": ("Base smoothing", 0.10, 1.50, 0.05, PROTOCOL_BASE_SCALE_MM, 2, " mm"),
    "broad_base_scale_mm": ("Broad-base smoothing", 0.80, 3.00, 0.05, PROTOCOL_BROAD_BASE_SCALE_MM, 2, " mm"),
    "smoothing_alpha": ("Mesh smoothing alpha", 0.10, 0.90, 0.05, PROTOCOL_SMOOTHING_ALPHA, 2, ""),
    "thickness_contrast": ("Thickness contrast", 0.10, 0.45, 0.01, PROTOCOL_THICKNESS_CONTRAST, 2, ""),
    "roughness_radius_mm": ("Roughness radius", 0.25, 5.00, 0.05, PROTOCOL_ROUGHNESS_RADIUS_MM, 2, " mm"),
    "roughness_contrast": ("Roughness contrast", 0.10, 0.45, 0.01, PROTOCOL_ROUGHNESS_CONTRAST, 2, ""),
    "highlight_q_low": ("Reflect q low", 0.805, 0.995, 0.005, PROTOCOL_HIGHLIGHT_Q_LOW, 3, ""),
    "highlight_q_high": ("Reflect q high", 0.990, 1.000, 0.001, PROTOCOL_HIGHLIGHT_Q_HIGH, 3, ""),
    "highlight_strength": ("Reflect strength", 0.00, 0.44, 0.01, PROTOCOL_HIGHLIGHT_STRENGTH, 2, ""),
    "highlight_gamma": ("Reflect gamma", 0.40, 3.00, 0.05, PROTOCOL_HIGHLIGHT_GAMMA, 2, ""),
    "highlight_cap_ratio": ("Reflect cap", 0.10, 0.80, 0.01, PROTOCOL_HIGHLIGHT_CAP_RATIO, 2, ""),
}

TOP_PARAM_SPECS = {
    "base_scale_mm": ("Depth/base smoothing", 0.10, 1.50, 0.05, PROTOCOL_BASE_SCALE_MM, 2, " mm"),
    "smoothing_alpha": ("Mesh smoothing alpha", 0.10, 0.90, 0.05, PROTOCOL_SMOOTHING_ALPHA, 2, ""),
    "fine_scale_mm": ("Fine smoothing", 0.10, 0.60, 0.01, PROTOCOL_FINE_SCALE_MM, 2, " mm"),
    "mesh_top_scale_mm": ("Medium relief scale", 0.40, 3.00, 0.05, PROTOCOL_MESH_TOP_SCALE_MM, 2, " mm"),
    "depth_med_percent": ("Depth medium sigma", 0.20, 3.00, 0.05, PROTOCOL_DEPTH_MED_PERCENT, 2, " %"),
    "top_post_sigma_px": ("Top post sigma", 0.00, 1.20, 0.05, PROTOCOL_TOP_POST_SIGMA_PX, 2, " px"),
    "curvature_sigma_percent": ("Curvature sigma", 0.50, 12.00, 0.50, PROTOCOL_CURV_SIGMA_PERCENT, 1, " %"),
    "openness_radius_mm": ("Openness radius", 0.25, 10.00, 0.25, PROTOCOL_OPENNESS_RADIUS_MM, 2, " mm"),
    "hessian_sigma_mm": ("Hessian sigma", 0.10, 5.00, 0.10, PROTOCOL_HESSIAN_SIGMA_MM, 2, " mm"),
    "normal_radius_mm": ("Normal radius", 0.10, 5.00, 0.10, PROTOCOL_NORMAL_RADIUS_MM, 2, " mm"),
    "principal_scale_mm": ("Principal scale", 0.10, 3.00, 0.10, PROTOCOL_PRINCIPAL_SCALE_MM, 2, " mm"),
    "linearity_sigma_percent": ("Linearity sigma", 0.10, 2.00, 0.05, PROTOCOL_LINEARITY_SIGMA_PERCENT, 2, " %"),
    "continuity_span_percent": ("Continuity span", 0.25, 5.00, 0.25, PROTOCOL_CONTINUITY_SPAN_PERCENT, 2, " %"),
    "linearity_weight": ("Linearity weight", 0.00, 2.00, 0.05, PROTOCOL_LINEARITY_WEIGHT, 2, ""),
    "continuity_weight": ("Continuity weight", 0.00, 3.00, 0.05, PROTOCOL_CONTINUITY_WEIGHT, 2, ""),
}


class VariationCard(QGroupBox):
    def __init__(self, index: int, on_select, parent=None):
        super().__init__(parent)
        self.index = index
        self.on_select = on_select
        self.record: Optional[LabRecord] = None
        self.setMinimumWidth(205)
        lay = QVBoxLayout(self)
        self.step_label = QLabel(["-2 step", "-1 step", "baseline", "+1 step", "+2 step"][index])
        self.step_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.value_label = QLabel("—")
        self.value_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.preview = PreviewLabel("Not rendered", minimum=(185, 230))
        self.select_btn = QPushButton("Select")
        self.select_btn.clicked.connect(lambda: self.on_select(self.index))
        lay.addWidget(self.step_label)
        lay.addWidget(self.value_label)
        lay.addWidget(self.preview, 1)
        lay.addWidget(self.select_btn)
        self.set_selected(False)

    def clear_record(self):
        self.record = None
        self.value_label.setText("—")
        self.preview.clear_preview("Not rendered")
        self.select_btn.setEnabled(False)
        self.set_selected(False)

    def set_record(self, record: LabRecord, parameter_label: str, value: float, suffix: str):
        self.record = record
        self.value_label.setText(f"{parameter_label}: {value:g}{suffix}")
        self.preview.set_array(record.image)
        self.select_btn.setEnabled(True)

    def set_selected(self, selected: bool):
        if selected:
            self.setStyleSheet(
                "QGroupBox { border: 3px solid #555; border-radius: 5px; margin-top: 4px; }"
            )
            self.select_btn.setText("Selected")
        else:
            self.setStyleSheet(
                "QGroupBox { border: 1px solid #aaa; border-radius: 5px; margin-top: 4px; }"
            )
            self.select_btn.setText("Select")


class SweepLabBase(QWidget):
    LAB_NAME = "lab"
    METHODS: List[str] = []
    SWEEP_CONFIG: Dict[str, Dict[str, Any]] = {}
    PARAM_SPECS: Dict[str, Tuple[Any, ...]] = {}

    def __init__(self, main: "LabMainWindow"):
        super().__init__()
        self.main = main
        # records are deliberately only RETAINED candidates. Composite Lab
        # therefore sees only candidates explicitly retained by the user.
        self.records: List[LabRecord] = []
        self.sweep_records: List[Optional[LabRecord]] = [None] * 5
        self.selected_sweep_index: Optional[int] = None
        self.active_method: Optional[str] = None
        self.method_states: Dict[str, Dict[str, Any]] = {}
        self.param_widgets: Dict[str, QDoubleSpinBox] = {}
        self.param_rows: Dict[str, Tuple[QLabel, QWidget]] = {}
        self._build()

    def _build(self):
        root = QHBoxLayout(self)

        # -------- compact left workflow panel --------
        controls_widget = QWidget()
        controls = QVBoxLayout(controls_widget)
        controls_widget.setMaximumWidth(390)

        workflow = QGroupBox(f"{self.LAB_NAME.title()} workflow")
        wf = QFormLayout(workflow)
        self.method = QComboBox()
        self.method.addItem("— Select method —")
        self.method.addItems(self.METHODS)
        self.primary_summary = QLabel("Select a method to generate five variations.")
        self.primary_summary.setWordWrap(True)
        wf.addRow("Method", self.method)
        wf.addRow("Sweep", self.primary_summary)
        controls.addWidget(workflow)

        row = QHBoxLayout()
        self.render_btn = QPushButton("Render 5 variations")
        self.customize_btn = QPushButton("Customize")
        self.customize_btn.setCheckable(True)
        row.addWidget(self.render_btn)
        row.addWidget(self.customize_btn)
        controls.addLayout(row)

        selected_box = QGroupBox("Selected variation")
        sl = QVBoxLayout(selected_box)
        self.selected_info = QLabel("None")
        self.selected_info.setWordWrap(True)
        sl.addWidget(self.selected_info)
        row2 = QHBoxLayout()
        self.retain_btn = QPushButton("Retain for Composite")
        self.save_btn = QPushButton("Save annotated image")
        self.bake_btn = QPushButton("Bake vertex-color PLY")
        row2.addWidget(self.retain_btn)
        row2.addWidget(self.save_btn)
        row2.addWidget(self.bake_btn)
        sl.addLayout(row2)
        controls.addWidget(selected_box)

        retained_box = QGroupBox("Retained candidates")
        rl = QVBoxLayout(retained_box)
        self.list = QListWidget()
        self.remove_retained_btn = QPushButton("Remove retained")
        rl.addWidget(self.list, 1)
        rl.addWidget(self.remove_retained_btn)
        controls.addWidget(retained_box, 1)

        # Customize panel is hidden by default.
        self.custom_panel = QGroupBox("Customize sweep and parameters")
        custom = QVBoxLayout(self.custom_panel)
        sweep_form = QFormLayout()
        self.primary_name = QLabel("Primary parameter")
        self.center = QDoubleSpinBox()
        self.step = QDoubleSpinBox()
        sweep_form.addRow(self.primary_name, self.center)
        sweep_form.addRow("Step", self.step)
        custom.addLayout(sweep_form)

        self.param_form = QFormLayout()
        defaults = RenderParams()
        for field, spec in self.PARAM_SPECS.items():
            label, mn, mx, st, default, dec, suffix = spec
            widget = make_spin(mn, mx, st, default, dec, suffix)
            lab = QLabel(label)
            self.param_widgets[field] = widget
            self.param_rows[field] = (lab, widget)
            self.param_form.addRow(lab, widget)

        self.projection = QComboBox()
        self.projection.addItems(["front", "right", "bottom"])
        self.height = make_spin(600, 1800, 50, PROTOCOL_IMAGE_HEIGHT, 0, " px")
        self.margin = make_spin(0.0, 0.06, 0.005, PROTOCOL_MARGIN_RATIO, 3)
        self.param_form.addRow("Projection", self.projection)
        self.param_form.addRow("Image height", self.height)
        self.param_form.addRow("Margin", self.margin)
        custom.addLayout(self.param_form)
        controls.addWidget(self.custom_panel)
        self.custom_panel.setVisible(False)

        # -------- five previews on the right --------
        right_widget = QWidget()
        right = QVBoxLayout(right_widget)
        title = QLabel("Five-variation comparison")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        right.addWidget(title)

        cards_widget = QWidget()
        cards_row = QHBoxLayout(cards_widget)
        cards_row.setContentsMargins(0, 0, 0, 0)
        self.cards: List[VariationCard] = []
        for i in range(5):
            card = VariationCard(i, self.select_sweep)
            self.cards.append(card)
            cards_row.addWidget(card, 1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(cards_widget)
        right.addWidget(scroll, 1)

        self.sweep_info = QLabel(
            "Choose a method. With a PLY loaded, five variations are rendered automatically."
        )
        self.sweep_info.setWordWrap(True)
        right.addWidget(self.sweep_info)

        root.addWidget(controls_widget)
        root.addWidget(right_widget, 1)

        # Connections
        self.method.currentTextChanged.connect(self.on_method_changed)
        self.render_btn.clicked.connect(self.render_sweep)
        self.customize_btn.toggled.connect(self.toggle_customize)
        self.retain_btn.clicked.connect(self.retain_selected)
        self.save_btn.clicked.connect(self.save_selected)
        self.bake_btn.clicked.connect(self.bake_selected)
        self.remove_retained_btn.clicked.connect(self.remove_retained)
        self.list.currentRowChanged.connect(lambda _row: self._update_action_state())
        self.center.valueChanged.connect(lambda _v: self._refresh_primary_summary())
        self.step.valueChanged.connect(lambda _v: self._refresh_primary_summary())

        self.set_mesh_ready(False)
        self._update_action_state()

    def toggle_customize(self, checked: bool):
        self.custom_panel.setVisible(bool(checked))
        self.customize_btn.setText("Close Customize" if checked else "Customize")

    def _capture_active_state(self):
        method = self.active_method
        if not method or method not in self.SWEEP_CONFIG:
            return
        self.method_states[method] = {
            "center": float(self.center.value()),
            "step": float(self.step.value()),
            "params": {k: float(w.value()) for k, w in self.param_widgets.items()},
            "projection": self.projection.currentText(),
            "image_height": int(self.height.value()),
            "margin_ratio": float(self.margin.value()),
        }

    def _restore_method_state(self, method: str):
        cfg = self.SWEEP_CONFIG[method]
        state = self.method_states.get(method)

        self.center.blockSignals(True)
        self.step.blockSignals(True)
        self.center.setDecimals(int(cfg["decimals"]))
        self.center.setRange(float(cfg["min"]), float(cfg["max"]))
        self.center.setSingleStep(float(cfg["step"]))
        self.center.setSuffix(str(cfg["suffix"]))
        self.step.setDecimals(int(cfg["decimals"]))
        self.step.setRange(0.000001, max(float(cfg["max"]) - float(cfg["min"]), 0.000001))
        self.step.setSingleStep(float(cfg["step"]))
        self.step.setSuffix(str(cfg["suffix"]))
        self.center.setValue(float(state["center"]) if state else float(cfg["baseline"]))
        self.step.setValue(float(state["step"]) if state else float(cfg["step"]))
        self.center.blockSignals(False)
        self.step.blockSignals(False)

        if state:
            for k, v in state["params"].items():
                if k in self.param_widgets:
                    self.param_widgets[k].setValue(float(v))
            idx = self.projection.findText(state["projection"])
            if idx >= 0:
                self.projection.setCurrentIndex(idx)
            self.height.setValue(float(state["image_height"]))
            self.margin.setValue(float(state["margin_ratio"]))

        self.primary_name.setText(cfg["label"])
        self._refresh_aux_visibility()
        self._refresh_primary_summary()

    def _refresh_aux_visibility(self):
        method = self.active_method
        aux = set(self.SWEEP_CONFIG.get(method, {}).get("aux", []))
        primary = self.SWEEP_CONFIG.get(method, {}).get("primary")
        for field, (lab, widget) in self.param_rows.items():
            visible = field in aux and field != primary
            lab.setVisible(visible)
            widget.setVisible(visible)

    def _refresh_primary_summary(self):
        method = self.active_method
        if not method or method not in self.SWEEP_CONFIG:
            self.primary_summary.setText("Select a method to generate five variations.")
            return
        cfg = self.SWEEP_CONFIG[method]
        c = float(self.center.value())
        s = float(self.step.value())
        vals = [c + d * s for d in (-2, -1, 0, 1, 2)]
        self.primary_summary.setText(
            f"{cfg['label']}: "
            + " / ".join(f"{v:g}{cfg['suffix']}" for v in vals)
        )

    def on_method_changed(self, method: str):
        # Save customized state for the method being left.
        if self.active_method and self.active_method in self.SWEEP_CONFIG:
            self._capture_active_state()

        if method not in self.SWEEP_CONFIG:
            self.active_method = None
            self.clear_sweep()
            self._refresh_primary_summary()
            self.render_btn.setEnabled(False)
            return

        self.active_method = method
        self._restore_method_state(method)
        self.clear_sweep()
        self.render_btn.setEnabled(self.main.mesh_state is not None)

        # Method selection is the normal trigger for the initial five renders.
        if self.main.mesh_state is not None:
            QApplication.processEvents()
            self.render_sweep()

    def set_mesh_ready(self, ready: bool):
        self.method.setEnabled(ready)
        self.customize_btn.setEnabled(ready)
        self.render_btn.setEnabled(
            ready and self.active_method is not None and self.active_method in self.SWEEP_CONFIG
        )
        if not ready:
            self.clear_all(reset_method=True)

    def clear_sweep(self):
        self.sweep_records = [None] * 5
        self.selected_sweep_index = None
        for card in self.cards:
            card.clear_record()
        self.selected_info.setText("None")
        self.sweep_info.setText(
            "Choose a method. With a PLY loaded, five variations are rendered automatically."
        )
        self._update_action_state()

    def clear_all(self, reset_method=False):
        self.records.clear()
        self.list.clear()
        self.clear_sweep()
        if reset_method:
            self.method.blockSignals(True)
            self.method.setCurrentIndex(0)
            self.method.blockSignals(False)
            self.active_method = None
        if hasattr(self.main, "composite_lab"):
            self.main.composite_lab.refresh_sources()

    # compatibility with v0.7.2 MainWindow calls
    def clear_records(self):
        self.clear_all(reset_method=False)

    def _base_params(self, method: str) -> RenderParams:
        p = RenderParams()
        kwargs = {}
        for field, widget in self.param_widgets.items():
            kwargs[field] = float(widget.value())
        kwargs.update({
            "projection": self.projection.currentText(),
            "image_height": int(self.height.value()),
            "margin_ratio": float(self.margin.value()),
        })
        p = replace(p, **kwargs)
        if self.LAB_NAME == "base":
            p = replace(p, base_choice=method, top1_choice="None", top2_choice="None")
        else:
            p = replace(p, base_choice="Base reflect-suppressed", top1_choice=method, top2_choice="None")
        return p

    def _five_values(self) -> Optional[List[float]]:
        method = self.active_method
        if not method:
            return None
        cfg = self.SWEEP_CONFIG[method]
        center = float(self.center.value())
        step = float(self.step.value())
        lo = center - 2.0 * step
        hi = center + 2.0 * step
        if step <= 0:
            QMessageBox.warning(self, "Invalid step", "Step must be greater than zero.")
            return None
        if lo < float(cfg["min"]) - 1e-12 or hi > float(cfg["max"]) + 1e-12:
            QMessageBox.warning(
                self,
                "Sweep outside allowed range",
                f"{cfg['label']} five-variation range would be "
                f"{lo:g} to {hi:g}{cfg['suffix']}.\n"
                f"Allowed range is {cfg['min']:g} to {cfg['max']:g}{cfg['suffix']}.\n\n"
                "Adjust Center or Step in Customize.",
            )
            return None
        return [center + d * step for d in (-2, -1, 0, 1, 2)]

    def render_sweep(self):
        if self.main.mesh_state is None or not self.active_method:
            return
        values = self._five_values()
        if values is None:
            return

        self._capture_active_state()
        method = self.active_method
        cfg = self.SWEEP_CONFIG[method]
        primary = cfg["primary"]
        base_p = self._base_params(method)

        self.clear_sweep()
        self.render_btn.setEnabled(False)
        self.method.setEnabled(False)
        try:
            for i, value in enumerate(values):
                p = replace(base_p, **{primary: float(value)})
                bundle = self.main.compute_bundle(
                    p,
                    f"{self.LAB_NAME.title()} {i+1}/5: {method} | {cfg['label']}={value:g}{cfg['suffix']}",
                )
                if self.LAB_NAME == "base":
                    image = bundle.base_candidates[method].copy()
                else:
                    image = bundle.top_candidates[method].copy()

                meta = dict(bundle.meta)
                meta["sweep"] = {
                    "primary": primary,
                    "label": cfg["label"],
                    "value": float(value),
                    "center": float(self.center.value()),
                    "step": float(self.step.value()),
                    "offset": i - 2,
                }
                rec = LabRecord(
                    self.LAB_NAME, method, p, image,
                    bundle.mask.copy(), bundle.audit_text, meta
                )
                self.sweep_records[i] = rec
                self.cards[i].set_record(rec, cfg["label"], value, cfg["suffix"])
                QApplication.processEvents()

            # Baseline is selected by default after a completed sweep.
            self.select_sweep(2)
            self.sweep_info.setText(
                f"{method}: five variations rendered. Select one, then Retain for Composite, "
                "Save, Bake, or use the selected value as the center of the next sweep."
            )
        except Exception as exc:
            QMessageBox.critical(
                self, f"{self.LAB_NAME.title()} render error",
                f"{type(exc).__name__}: {exc}"
            )
        finally:
            self.render_btn.setEnabled(True)
            self.method.setEnabled(True)
            self._refresh_primary_summary()

    def select_sweep(self, index: int):
        if index < 0 or index >= len(self.sweep_records):
            return
        rec = self.sweep_records[index]
        if rec is None:
            return
        self.selected_sweep_index = index
        for i, card in enumerate(self.cards):
            card.set_selected(i == index)

        sweep = rec.meta.get("sweep", {})
        value = float(sweep.get("value", self.center.value()))
        cfg = self.SWEEP_CONFIG[self.active_method]
        self.center.blockSignals(True)
        self.center.setValue(value)
        self.center.blockSignals(False)
        self._capture_active_state()
        self._refresh_primary_summary()

        self.selected_info.setText(
            f"{rec.method}\n{cfg['label']}={value:g}{cfg['suffix']}\n"
            f"Selected variation {index+1}/5; this value is now the center for the next sweep."
        )
        self._update_action_state()

    def selected_record(self) -> Optional[LabRecord]:
        if self.selected_sweep_index is None:
            return None
        if not (0 <= self.selected_sweep_index < len(self.sweep_records)):
            return None
        return self.sweep_records[self.selected_sweep_index]

    def retain_selected(self):
        rec = self.selected_record()
        if rec is None:
            return
        # Snapshot arrays/metadata so later sweeps never mutate the retained result.
        retained = LabRecord(
            rec.lab, rec.method, rec.params,
            rec.image.copy(), rec.mask.copy(),
            str(rec.audit_text), dict(rec.meta)
        )
        self.records.append(retained)
        self.list.addItem(QListWidgetItem(retained.label))
        self.list.setCurrentRow(len(self.records) - 1)
        self.main.composite_lab.refresh_sources()
        self.sweep_info.setText(
            f"Retained {METHOD_ABBR.get(retained.method, 'MTH')} candidate "
            f"#{len(self.records)} for Composite Lab."
        )

    def remove_retained(self):
        row = self.list.currentRow()
        if row < 0 or row >= len(self.records):
            return
        del self.records[row]
        self.list.takeItem(row)
        self.main.composite_lab.refresh_sources()

    def save_selected(self):
        rec = self.selected_record()
        if rec is not None:
            self.main.save_record_image(rec)

    def bake_selected(self):
        rec = self.selected_record()
        if rec is not None:
            self.main.bake_record(rec)

    def _update_action_state(self):
        selected = self.selected_record() is not None
        self.retain_btn.setEnabled(selected)
        self.save_btn.setEnabled(selected)
        self.bake_btn.setEnabled(selected)
        self.remove_retained_btn.setEnabled(self.list.currentRow() >= 0)

    def _list_selection_changed(self):
        self._update_action_state()


class BaseLab(SweepLabBase):
    LAB_NAME = "base"
    METHODS = BASE_METHODS
    SWEEP_CONFIG = BASE_SWEEP_CONFIG
    PARAM_SPECS = BASE_PARAM_SPECS


class TopLab(SweepLabBase):
    LAB_NAME = "top"
    METHODS = TOP_METHODS
    SWEEP_CONFIG = TOP_SWEEP_CONFIG
    PARAM_SPECS = TOP_PARAM_SPECS

class CompositeLab(QWidget):
    def __init__(self, main: "LabMainWindow"):
        super().__init__(); self.main=main; self.current_image: Optional[np.ndarray]=None; self.current_mask: Optional[np.ndarray]=None; self.current_base: Optional[LabRecord]=None; self.current_top: Optional[LabRecord]=None; self._build()

    def _build(self):
        root=QHBoxLayout(self); controls_widget=QWidget(); controls=QVBoxLayout(controls_widget); form=QFormLayout()
        self.base_source=QComboBox(); self.top_source=QComboBox(); self.mode=QComboBox(); self.mode.addItems(["multiply","alpha","overlay","screen","add"])
        self.base_qlo=make_spin(0,20,0.5,1.0,1," %"); self.base_qhi=make_spin(80,100,0.5,99.0,1," %")
        self.top_qlo=make_spin(0,20,0.5,1.0,1," %"); self.top_qhi=make_spin(80,100,0.5,99.0,1," %")
        self.base_contrast=make_spin(0.1,3.0,0.05,1.0,2); self.top_contrast=make_spin(0.1,3.0,0.05,1.0,2); self.alpha=make_spin(0,1,0.05,0.5,2)
        for label,widget in (("Base candidate",self.base_source),("Top candidate",self.top_source),("Normalization Base low",self.base_qlo),("Normalization Base high",self.base_qhi),("Normalization Top low",self.top_qlo),("Normalization Top high",self.top_qhi),("Base contrast",self.base_contrast),("Top contrast",self.top_contrast),("Blend mode",self.mode),("Alpha / strength",self.alpha)): form.addRow(label,widget)
        controls.addLayout(form)
        row=QHBoxLayout(); self.render_btn=QPushButton("Render Composite"); self.save_btn=QPushButton("Save annotated image"); self.bake_btn=QPushButton("Bake vertex-color PLY"); row.addWidget(self.render_btn); row.addWidget(self.save_btn); row.addWidget(self.bake_btn); controls.addLayout(row); controls.addStretch(1); controls_widget.setMaximumWidth(430)
        right=QVBoxLayout(); self.preview=PreviewLabel("Render Base and Top candidates first."); self.info=QLabel(""); self.info.setWordWrap(True); right.addWidget(self.preview,1); right.addWidget(self.info); rw=QWidget(); rw.setLayout(right); root.addWidget(controls_widget); root.addWidget(rw,1)
        self.render_btn.clicked.connect(self.render_composite); self.save_btn.clicked.connect(self.save_current); self.bake_btn.clicked.connect(self.bake_current); self.refresh_sources()

    def clear(self):
        self.current_image=None; self.current_mask=None; self.current_base=None; self.current_top=None; self.preview.clear_preview("Render Base and Top candidates first."); self.info.setText(""); self.refresh_sources()

    def refresh_sources(self):
        bi=self.base_source.currentIndex() if hasattr(self,'base_source') else -1; ti=self.top_source.currentIndex() if hasattr(self,'top_source') else -1
        if not hasattr(self,'base_source'): return
        self.base_source.blockSignals(True); self.top_source.blockSignals(True); self.base_source.clear(); self.top_source.clear()
        for i,r in enumerate(self.main.base_lab.records): self.base_source.addItem(f"B{i+1}: {r.label}")
        for i,r in enumerate(self.main.top_lab.records): self.top_source.addItem(f"T{i+1}: {r.label}")
        if self.base_source.count(): self.base_source.setCurrentIndex(min(max(bi,0),self.base_source.count()-1))
        if self.top_source.count(): self.top_source.setCurrentIndex(min(max(ti,0),self.top_source.count()-1))
        self.base_source.blockSignals(False); self.top_source.blockSignals(False)
        ready=self.base_source.count()>0 and self.top_source.count()>0 and self.main.mesh_state is not None
        self.render_btn.setEnabled(ready); self.save_btn.setEnabled(self.current_image is not None); self.bake_btn.setEnabled(self.current_image is not None)

    def render_composite(self):
        if self.base_source.count()==0 or self.top_source.count()==0: return
        b=self.main.base_lab.records[self.base_source.currentIndex()]; t=self.main.top_lab.records[self.top_source.currentIndex()]
        if b.image.shape!=t.image.shape or b.params.projection!=t.params.projection or b.params.image_height!=t.params.image_height or abs(b.params.margin_ratio-t.params.margin_ratio)>1e-12:
            QMessageBox.warning(self,"Composite geometry mismatch","Base and Top must use the same projection, image height and margin. Re-render one candidate with matching settings."); return
        mask=b.mask & t.mask
        bn=normalize_for_composite(b.image,mask,self.base_qlo.value(),self.base_qhi.value()); tn=normalize_for_composite(t.image,mask,self.top_qlo.value(),self.top_qhi.value())
        bc=contrast_adjust(bn,mask,self.base_contrast.value()); tc=contrast_adjust(tn,mask,self.top_contrast.value())
        out=blend_with_strength(bc,tc,mask,self.mode.currentText(),self.alpha.value())
        self.current_image=out; self.current_mask=mask; self.current_base=b; self.current_top=t; self.preview.set_array(out); self.info.setText(f"{b.method} + {t.method} | mode={self.mode.currentText()} alpha={self.alpha.value():.2f} | baseC={self.base_contrast.value():.2f} topC={self.top_contrast.value():.2f}"); self.save_btn.setEnabled(True); self.bake_btn.setEnabled(True)

    def save_current(self):
        if self.current_image is None or self.current_base is None or self.current_top is None: return
        outdir=QFileDialog.getExistingDirectory(self,"Select output folder");
        if not outdir:return
        btoken=METHOD_ABBR.get(self.current_base.method,"B"); ttoken=METHOD_ABBR.get(self.current_top.method,"T")
        stem=f"{self.main.mesh_state.path.stem}_CMP_{btoken}_{ttoken}_a{short_value(self.alpha.value())}_bc{short_value(self.base_contrast.value())}_tc{short_value(self.top_contrast.value())}"
        path=Path(outdir)/f"{stem}.png"
        ann=annotated_render_image(self.current_image,"Composite",[f"Base: {self.current_base.method}",f"Top: {self.current_top.method}",f"mode={self.mode.currentText()} alpha={self.alpha.value():g}",f"base contrast={self.base_contrast.value():g}, top contrast={self.top_contrast.value():g}",f"base normalization={self.base_qlo.value():g}-{self.base_qhi.value():g} percentile",f"top normalization={self.top_qlo.value():g}-{self.top_qhi.value():g} percentile"]); ann.save(path); QMessageBox.information(self,"Saved",f"Saved:\n{path}")

    def bake_current(self):
        if self.current_image is None or self.current_base is None or self.current_top is None:return
        self.main.bake_composite(self.current_image,self.current_base,self.current_top,self)


class LabMainWindow(QMainWindow):
    def __init__(self):
        super().__init__(); self.mesh_state: Optional[MeshState]=None; self.setWindowTitle(APP_NAME)
        tb=QToolBar("Main"); self.addToolBar(tb); self.open_btn=QPushButton("Open PLY"); self.status=QLabel("No PLY loaded"); self.progress=QProgressBar(); self.progress.setRange(0,100); self.progress.setMaximumWidth(180); tb.addWidget(self.open_btn); tb.addWidget(self.status); tb.addWidget(self.progress); self.open_btn.clicked.connect(self.open_ply)
        self.tabs=QTabWidget(); self.base_lab=BaseLab(self); self.top_lab=TopLab(self); self.composite_lab=CompositeLab(self); self.tabs.addTab(self.base_lab,"Base Lab"); self.tabs.addTab(self.top_lab,"Top Lab"); self.tabs.addTab(self.composite_lab,"Composite Lab"); self.setCentralWidget(self.tabs); self.resize(1720,960)

    def open_ply(self):
        path_str,_=QFileDialog.getOpenFileName(self,"Open PLY","","PLY files (*.ply)");
        if not path_str:return
        try:
            self.status.setText("Loading PLY…"); QApplication.processEvents(); self.mesh_state=load_mesh_state_for_lab(Path(path_str)); bbox=self.mesh_state.bounds[1]-self.mesh_state.bounds[0]
            self.base_lab.clear_all(reset_method=True); self.top_lab.clear_all(reset_method=True); self.base_lab.set_mesh_ready(True); self.top_lab.set_mesh_ready(True); self.composite_lab.clear(); self.status.setText(f"{self.mesh_state.path.name} | V={len(self.mesh_state.vertices):,} F={len(self.mesh_state.faces):,} | bbox={bbox[0]:.2f}×{bbox[1]:.2f}×{bbox[2]:.2f} | choose a Base/Top method to render five variations"); self.progress.setValue(0)
        except Exception as exc:
            self.mesh_state=None; self.base_lab.set_mesh_ready(False); self.top_lab.set_mesh_ready(False); self.composite_lab.clear(); self.status.setText("PLY load failed"); QMessageBox.critical(self,"PLY load error",f"{type(exc).__name__}: {exc}")

    def compute_bundle(self,p:RenderParams,status_text:str)->RenderBundle:
        if self.mesh_state is None: raise RuntimeError("No mesh loaded")
        self.status.setText(status_text); self.progress.setValue(0); QApplication.processEvents()
        def cb(value,text): self.progress.setValue(int(value)); self.status.setText(text); QApplication.processEvents()
        bundle=compute_render_bundle(self.mesh_state,p,cb); self.progress.setValue(100); self.status.setText("Render complete"); return bundle

    def save_record_image(self,rec:LabRecord):
        outdir=QFileDialog.getExistingDirectory(self,"Select output folder");
        if not outdir:return
        stem=f"{self.mesh_state.path.stem}_{method_file_token(rec.method,rec.params)}"; path=Path(outdir)/f"{stem}.png"; sweep=rec.meta.get("sweep",{}); lines=[f"Method: {rec.method}",f"Projection: {rec.params.projection}"]; 
        if sweep: lines += [f"Sweep: {sweep.get('label', sweep.get('primary','parameter'))}={sweep.get('value','')} | center={sweep.get('center','')} | step={sweep.get('step','')} | offset={sweep.get('offset','')}"]
        lines += [f"{k}={v:g}" for k,v in method_filename_params(rec.method,rec.params)]; annotated_render_image(rec.image, f"{rec.lab.title()} Lab",lines).save(path); QMessageBox.information(self,"Saved",f"Saved:\n{path}")

    def bake_record(self,rec:LabRecord):
        outdir=QFileDialog.getExistingDirectory(self,"Select output folder");
        if not outdir:return
        path=Path(outdir)/f"{self.mesh_state.path.stem}_{method_file_token(rec.method,rec.params)}.ply"
        try:
            self.status.setText("Baking vertex colors…"); self.progress.setValue(10); QApplication.processEvents(); rgba,stats=bake_array_to_vertex_colors(self.mesh_state,rec.image,rec.params); self.progress.setValue(75); export_baked_vertex_color_ply(self.mesh_state,rgba,path); verify=verify_baked_ply(path,len(self.mesh_state.vertices),len(self.mesh_state.faces)); self.progress.setValue(100); self.status.setText("Vertex Color Bake complete"); QMessageBox.information(self,"Vertex Color Bake",f"Saved:\n{path}\n\nColored vertices: {stats['colored_vertices']:,}/{stats['total_vertices']:,} ({stats['colored_ratio']*100:.1f}%)\nProjection: {stats['projection']}\nNormals: preserved\nUV: intentionally omitted in baked PLY so vertex color remains active\nVerification: {verify}")
        except Exception as exc:
            self.progress.setValue(0); self.status.setText("Vertex Color Bake failed"); QMessageBox.critical(self,"Bake error",f"{type(exc).__name__}: {exc}")

    def bake_composite(self,image:np.ndarray,b:LabRecord,t:LabRecord,parent):
        outdir=QFileDialog.getExistingDirectory(parent,"Select output folder");
        if not outdir:return
        btoken=METHOD_ABBR.get(b.method,"B"); ttoken=METHOD_ABBR.get(t.method,"T"); path=Path(outdir)/f"{self.mesh_state.path.stem}_CMP_{btoken}_{ttoken}.ply"
        try:
            rgba,stats=bake_array_to_vertex_colors(self.mesh_state,image,b.params); export_baked_vertex_color_ply(self.mesh_state,rgba,path); verify=verify_baked_ply(path,len(self.mesh_state.vertices),len(self.mesh_state.faces)); QMessageBox.information(parent,"Composite Vertex Color Bake",f"Saved:\n{path}\n\nColored vertices: {stats['colored_vertices']:,}/{stats['total_vertices']:,} ({stats['colored_ratio']*100:.1f}%)\nVerification: {verify}")
        except Exception as exc: QMessageBox.critical(parent,"Bake error",f"{type(exc).__name__}: {exc}")


if __name__ == "__main__":
    app=QApplication(sys.argv)
    window=LabMainWindow(); window.show()
    sys.exit(app.exec())

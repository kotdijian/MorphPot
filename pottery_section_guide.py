#!/usr/bin/env python3
"""
PotterySectionGuide
Version 0.1.2

Generate 3-D circular section guides from a 2-D pottery longitudinal-section drawing.

Primary use case
----------------
The generated geometry is intended as a placement guide for archaeological pottery
sherd 3-D data in applications such as CloudCompare.

Design principles
-----------------
- Scale is always calibrated from two user-selected points of known distance.
- Rotation axis can be detected automatically (OpenCV Hough-line candidate) or
  specified manually by two points.
- Every readable longitudinal profile should be digitized independently:
    outer_left, outer_right, inner_left, inner_right.
- Radius is NOT derived from vessel width and is NOT averaged between left/right.
  For each profile, radius is the perpendicular distance from the rotation axis to
  that profile at the requested axial Z position.
- Cross-sections are true 3-D circles. The drawing itself remains a 2-D longitudinal
  section in the X-Z plane (Y=0) of the output local coordinate system.
- Guide levels may be created at a regular Z interval, at user-clicked positions,
  or by combining both.
- Internal measurement and QA values are kept in millimetres.
- OBJ/PLY geometry can be exported in mm, cm, or m for direct use with an existing 3-D project.

Supported sources
-----------------
Raster: PNG / JPEG / TIFF / BMP / WEBP
Vector/doc: SVG / PDF (rendered to pixels with PyMuPDF for digitizing)

CloudCompare-oriented outputs
-----------------------------
- OBJ polylines: rings, profiles, and rotation axis
- PLY point clouds: combined guide cloud and one file per ring
- CSV / JSON QA metadata
- PNG overlay showing the digitized source geometry and requested section levels
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageTk

__version__ = "0.1.2"

RASTER_EXTENSIONS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
VECTOR_DOCUMENT_EXTENSIONS = {".svg", ".pdf"}
SUPPORTED_EXTENSIONS = RASTER_EXTENSIONS | VECTOR_DOCUMENT_EXTENSIONS
PROFILE_KEYS = ("outer_left", "outer_right", "inner_left", "inner_right")
EXPORT_UNIT_SCALE_FROM_MM = {"mm": 1.0, "cm": 0.1, "m": 0.001}


@dataclass
class Calibration:
    p1_px: tuple[float, float]
    p2_px: tuple[float, float]
    known_length_mm: float
    pixel_length: float
    mm_per_pixel: float


@dataclass
class Axis2D:
    """Rotation axis represented by two pixel points and an oriented unit vector."""

    p_bottom_px: tuple[float, float]
    p_top_px: tuple[float, float]
    source: str
    detection_score: float | None = None

    def basis(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (e_r, e_z) in source-image pixel coordinates.

        e_z points upward along the rotation axis. e_r is perpendicular to e_z and
        points approximately to image-right. The local 3-D drawing plane is X-Z,
        with X along e_r and Y=0.
        """
        p0 = np.asarray(self.p_bottom_px, dtype=float)
        p1 = np.asarray(self.p_top_px, dtype=float)
        v = p1 - p0
        n = float(np.linalg.norm(v))
        if n <= 0:
            raise ValueError("Rotation-axis points must be different.")
        e_z = v / n
        # Ensure +Z points generally upward in the image (negative image Y).
        if e_z[1] > 0:
            e_z = -e_z
        # Perpendicular, chosen to point generally to image-right.
        e_r = np.array([-e_z[1], e_z[0]], dtype=float)
        if e_r[0] < 0:
            e_r = -e_r
        return e_r, e_z

    def project_to_axis(self, point_px: tuple[float, float]) -> tuple[float, float]:
        """Return (u_px, z_px) relative to p_bottom_px.

        u_px is signed perpendicular distance to axis; z_px is position along axis.
        """
        e_r, e_z = self.basis()
        p = np.asarray(point_px, dtype=float)
        p0 = np.asarray(self.p_bottom_px, dtype=float)
        d = p - p0
        return float(np.dot(d, e_r)), float(np.dot(d, e_z))

    def point_on_axis_at_z_px(self, z_px: float) -> np.ndarray:
        _, e_z = self.basis()
        return np.asarray(self.p_bottom_px, dtype=float) + float(z_px) * e_z


@dataclass
class PreparedProfile:
    key: str
    points_px: list[tuple[float, float]]
    z_mm: np.ndarray
    signed_u_mm: np.ndarray
    radius_mm: np.ndarray

    @property
    def surface(self) -> str:
        return "outer" if self.key.startswith("outer_") else "inner"

    @property
    def side(self) -> str:
        return "left" if self.key.endswith("_left") else "right"


@dataclass
class RingRecord:
    ring_id: int
    z_mm: float
    profile_key: str
    surface: str
    side: str
    radius_mm: float
    section_source: str
    point_count: int


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


def export_scale_from_mm(unit: str) -> float:
    """Return the multiplicative scale applied to internal millimetre coordinates."""
    unit = str(unit).lower()
    if unit not in EXPORT_UNIT_SCALE_FROM_MM:
        raise ValueError(f"Unsupported export unit: {unit}. Choose mm, cm, or m.")
    return float(EXPORT_UNIT_SCALE_FROM_MM[unit])


def geometry_for_export(points_mm: np.ndarray, unit: str) -> np.ndarray:
    """Convert internal mm geometry to the numeric coordinate unit written to OBJ/PLY."""
    return np.asarray(points_mm, dtype=float) * export_scale_from_mm(unit)


def load_source_image(path: Path, page_number: int = 1, render_dpi: float = 180.0) -> tuple[Image.Image, dict]:
    """Load raster/SVG/PDF into the pixel image used by the GUI."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported drawing format: {ext}. Supported: {sorted(SUPPORTED_EXTENSIONS)}")

    if ext in RASTER_EXTENSIONS:
        img = Image.open(path).convert("RGB")
        return img, {
            "source_kind": "raster_image",
            "source_extension": ext,
            "page_number": None,
            "render_dpi": None,
            "vector_semantics_used": False,
        }

    try:
        import pymupdf
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "SVG/PDF input requires PyMuPDF. Install with: python3 -m pip install pymupdf"
        ) from exc

    doc = pymupdf.open(str(path))
    if doc.page_count < 1:
        raise ValueError("SVG/PDF contains no renderable page.")
    page_number = int(page_number)
    if not (1 <= page_number <= doc.page_count):
        raise ValueError(f"page_number must be 1..{doc.page_count}")
    page = doc[page_number - 1]
    zoom = float(render_dpi) / 72.0
    pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    kind = "svg_render" if ext == ".svg" else "pdf_page_render"
    return img, {
        "source_kind": kind,
        "source_extension": ext,
        "page_number": page_number,
        "document_page_count": int(doc.page_count),
        "render_dpi": float(render_dpi),
        "vector_semantics_used": False,
        "note": "Rendered pixels are digitized; physical scale comes only from GUI calibration.",
    }


# -----------------------------------------------------------------------------
# Axis detection and 2-D -> local coordinates
# -----------------------------------------------------------------------------

def _axis_from_two_points(p1: tuple[float, float], p2: tuple[float, float], source: str,
                          detection_score: float | None = None) -> Axis2D:
    a = np.asarray(p1, dtype=float)
    b = np.asarray(p2, dtype=float)
    if float(np.linalg.norm(a - b)) <= 1e-9:
        raise ValueError("Rotation-axis points must be different.")
    # Bottom is the endpoint with larger image Y. If close to horizontal, keep the
    # endpoint that yields an upward vector with the stronger negative-Y component.
    if a[1] >= b[1]:
        bottom, top = a, b
    else:
        bottom, top = b, a
    return Axis2D(tuple(map(float, bottom)), tuple(map(float, top)), source, detection_score)


def detect_rotation_axis(source_image: Image.Image, max_vertical_deviation_deg: float = 25.0) -> Axis2D:
    """Detect a likely central rotation axis using OpenCV Hough line segments.

    This is a candidate detector, not an authoritative interpretation. The GUI shows
    the detected line and the user can replace it with a manual two-point axis.
    """
    try:
        import cv2
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Automatic axis detection requires opencv-python. "
            "Install with: python3 -m pip install opencv-python\n"
            "Manual two-point axis specification works without OpenCV."
        ) from exc

    rgb = np.asarray(source_image.convert("RGB"))
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    h, w = gray.shape[:2]
    blur = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = cv2.Canny(blur, 50, 150, apertureSize=3)

    min_len = max(30, int(round(h * 0.08)))
    max_gap = max(5, int(round(h * 0.03)))
    threshold = max(30, int(round(h * 0.035)))
    lines = cv2.HoughLinesP(
        edges,
        rho=1,
        theta=np.pi / 180.0,
        threshold=threshold,
        minLineLength=min_len,
        maxLineGap=max_gap,
    )
    if lines is None or len(lines) == 0:
        raise RuntimeError("No line candidate was detected. Specify the rotation axis manually with two points.")

    max_dev = math.radians(float(max_vertical_deviation_deg))
    best = None
    for item in lines[:, 0, :]:
        x1, y1, x2, y2 = map(float, item)
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy)
        if length <= 0:
            continue
        # deviation from vertical; vertical vector has |dy| / length == 1
        vertical_dev = math.acos(min(1.0, abs(dy) / length))
        if vertical_dev > max_dev:
            continue
        mx = 0.5 * (x1 + x2)
        center_dist = abs(mx - 0.5 * w) / max(1.0, 0.5 * w)
        # Prefer long, central and near-vertical segments. A dashed centre line may
        # contribute only one segment; this intentionally remains conservative.
        vertical_factor = max(0.0, 1.0 - vertical_dev / max(max_dev, 1e-9))
        center_factor = math.exp(-2.0 * center_dist * center_dist)
        score = length * (0.55 + 0.45 * vertical_factor) * center_factor
        if best is None or score > best[0]:
            best = (score, x1, y1, x2, y2)

    if best is None:
        raise RuntimeError(
            f"No near-vertical axis candidate within {max_vertical_deviation_deg:g}°. "
            "Specify the rotation axis manually with two points."
        )

    score, x1, y1, x2, y2 = best
    dx, dy = x2 - x1, y2 - y1
    if abs(dy) < 1e-9:
        raise RuntimeError("Detected line is numerically horizontal; use manual axis specification.")

    # Extend candidate to top and bottom of the rendered image for stable display.
    def x_at_y(y: float) -> float:
        return x1 + (y - y1) * dx / dy

    p_top = (float(x_at_y(0.0)), 0.0)
    p_bottom = (float(x_at_y(float(h - 1))), float(h - 1))
    return _axis_from_two_points(p_bottom, p_top, source="opencv_hough_candidate", detection_score=float(score))


def project_point_local_mm(point_px: tuple[float, float], axis: Axis2D,
                           origin_px: tuple[float, float], mm_per_pixel: float) -> tuple[float, float]:
    """Project image point to local drawing-plane coordinates (signed X, Z) in mm."""
    e_r, e_z = axis.basis()
    p = np.asarray(point_px, dtype=float)
    o = np.asarray(origin_px, dtype=float)
    d = p - o
    u_mm = float(np.dot(d, e_r) * mm_per_pixel)
    z_mm = float(np.dot(d, e_z) * mm_per_pixel)
    return u_mm, z_mm


def project_click_to_axis(point_px: tuple[float, float], axis: Axis2D) -> tuple[float, float]:
    """Orthogonally project an image click onto the axis line."""
    _, e_z = axis.basis()
    p = np.asarray(point_px, dtype=float)
    p0 = np.asarray(axis.p_bottom_px, dtype=float)
    z = float(np.dot(p - p0, e_z))
    q = p0 + z * e_z
    return float(q[0]), float(q[1])


def automatic_origin_from_profiles(axis: Axis2D, profiles_px: dict[str, list[tuple[float, float]]]) -> tuple[float, float]:
    """Set Z=0 to the lowest axial position represented by any digitized profile point."""
    _, e_z = axis.basis()
    p0 = np.asarray(axis.p_bottom_px, dtype=float)
    zvals: list[float] = []
    for pts in profiles_px.values():
        for p in pts:
            zvals.append(float(np.dot(np.asarray(p, dtype=float) - p0, e_z)))
    if not zvals:
        raise ValueError("No digitized profile points are available for automatic origin definition.")
    zmin = min(zvals)
    q = p0 + zmin * e_z
    return float(q[0]), float(q[1])


# -----------------------------------------------------------------------------
# Profile preparation and interpolation
# -----------------------------------------------------------------------------

def _dedupe_profile_z(z_mm: np.ndarray, u_mm: np.ndarray, tolerance_mm: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(z_mm)
    z = np.asarray(z_mm, dtype=float)[order]
    u = np.asarray(u_mm, dtype=float)[order]
    out_z: list[float] = []
    out_u: list[float] = []
    i = 0
    while i < len(z):
        j = i + 1
        while j < len(z) and abs(float(z[j] - z[i])) <= tolerance_mm:
            j += 1
        out_z.append(float(np.mean(z[i:j])))
        out_u.append(float(np.median(u[i:j])))
        i = j
    return np.asarray(out_z, dtype=float), np.asarray(out_u, dtype=float)


def prepare_profile(key: str, points_px: Iterable[tuple[float, float]], axis: Axis2D,
                    origin_px: tuple[float, float], mm_per_pixel: float) -> PreparedProfile:
    if key not in PROFILE_KEYS:
        raise ValueError(f"Unknown profile key: {key}")
    pts = list(points_px)
    if len(pts) < 2:
        raise ValueError(f"{key} requires at least two points.")

    z_vals: list[float] = []
    u_vals: list[float] = []
    for p in pts:
        u, z = project_point_local_mm(p, axis, origin_px, mm_per_pixel)
        u_vals.append(u)
        z_vals.append(z)
    z, u = _dedupe_profile_z(np.asarray(z_vals), np.asarray(u_vals))
    if len(z) < 2 or float(np.ptp(z)) <= 0:
        raise ValueError(f"{key} does not span a measurable axial Z range.")

    # Side label is a QA expectation. We keep the measured signed geometry rather
    # than mirroring or averaging it.
    med = float(np.median(u))
    if key.endswith("_left") and med > 0:
        print(f"WARNING: {key} lies mostly on positive local X; check axis/profile labeling.", file=sys.stderr)
    if key.endswith("_right") and med < 0:
        print(f"WARNING: {key} lies mostly on negative local X; check axis/profile labeling.", file=sys.stderr)

    return PreparedProfile(
        key=key,
        points_px=[(float(x), float(y)) for x, y in pts],
        z_mm=z,
        signed_u_mm=u,
        radius_mm=np.abs(u),
    )


def interpolate_profile_radius(profile: PreparedProfile, z_mm: float, tolerance_mm: float = 1e-9) -> float | None:
    z = profile.z_mm
    if z_mm < float(z.min()) - tolerance_mm or z_mm > float(z.max()) + tolerance_mm:
        return None
    u = float(np.interp(float(z_mm), z, profile.signed_u_mm))
    return abs(u)


def interpolate_profile_signed_u(profile: PreparedProfile, z_mm: float, tolerance_mm: float = 1e-9) -> float | None:
    z = profile.z_mm
    if z_mm < float(z.min()) - tolerance_mm or z_mm > float(z.max()) + tolerance_mm:
        return None
    return float(np.interp(float(z_mm), z, profile.signed_u_mm))


def profile_z_range(profiles: dict[str, PreparedProfile]) -> tuple[float, float]:
    if not profiles:
        raise ValueError("No prepared profiles.")
    zmin = min(float(p.z_mm.min()) for p in profiles.values())
    zmax = max(float(p.z_mm.max()) for p in profiles.values())
    return zmin, zmax


# -----------------------------------------------------------------------------
# Guide-level selection and 3-D geometry
# -----------------------------------------------------------------------------

def regular_section_levels(zmin_mm: float, zmax_mm: float, interval_mm: float,
                           start_mm: float | None = None, end_mm: float | None = None) -> list[float]:
    interval = float(interval_mm)
    if interval <= 0:
        raise ValueError("Section interval must be > 0 mm.")
    lo = float(zmin_mm if start_mm is None else start_mm)
    hi = float(zmax_mm if end_mm is None else end_mm)
    if hi < lo:
        lo, hi = hi, lo
    # If no explicit start is supplied, anchor the regular series at zero whenever
    # zero lies in/near the profile range; otherwise start at the first interval
    # multiple inside the range.
    if start_mm is None:
        k = math.ceil((lo - 1e-9) / interval)
        first = k * interval
    else:
        first = lo
    vals: list[float] = []
    z = first
    # Numerically bounded loop.
    for _ in range(1_000_000):
        if z > hi + 1e-9:
            break
        if z >= lo - 1e-9:
            vals.append(float(z))
        z += interval
    return vals


def merge_section_levels(regular: Iterable[float], custom: Iterable[float], tolerance_mm: float = 1e-6) -> list[tuple[float, str]]:
    tagged = [(float(z), "regular") for z in regular] + [(float(z), "custom") for z in custom]
    tagged.sort(key=lambda x: x[0])
    out: list[tuple[float, str]] = []
    for z, src in tagged:
        if out and abs(z - out[-1][0]) <= tolerance_mm:
            # Preserve information that a level was explicitly requested.
            if src == "custom" and out[-1][1] != "custom":
                out[-1] = (out[-1][0], "regular+custom")
            elif src == "regular" and out[-1][1] == "custom":
                out[-1] = (out[-1][0], "regular+custom")
            continue
        out.append((z, src))
    return out


def circle_points(radius_mm: float, z_mm: float, point_count: int = 180) -> np.ndarray:
    r = float(radius_mm)
    n = int(point_count)
    if r < 0:
        raise ValueError("radius_mm must be non-negative")
    if n < 8:
        raise ValueError("point_count must be >= 8")
    theta = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    return np.column_stack((r * np.cos(theta), r * np.sin(theta), np.full(n, float(z_mm))))


def build_ring_records(profiles: dict[str, PreparedProfile], section_levels: list[tuple[float, str]],
                       point_count: int) -> tuple[list[RingRecord], dict[int, np.ndarray]]:
    records: list[RingRecord] = []
    geometries: dict[int, np.ndarray] = {}
    rid = 0
    for z, source in section_levels:
        for key in PROFILE_KEYS:
            profile = profiles.get(key)
            if profile is None:
                continue
            radius = interpolate_profile_radius(profile, z)
            if radius is None:
                continue
            record = RingRecord(
                ring_id=rid,
                z_mm=float(z),
                profile_key=key,
                surface=profile.surface,
                side=profile.side,
                radius_mm=float(radius),
                section_source=source,
                point_count=int(point_count),
            )
            records.append(record)
            geometries[rid] = circle_points(radius, z, point_count)
            rid += 1
    return records, geometries


def section_summary_rows(profiles: dict[str, PreparedProfile], section_levels: list[tuple[float, str]]) -> list[dict]:
    rows: list[dict] = []
    for z, source in section_levels:
        radii: dict[str, float | None] = {}
        for key in PROFILE_KEYS:
            profile = profiles.get(key)
            radii[key] = None if profile is None else interpolate_profile_radius(profile, z)
        row = {
            "z_mm": float(z),
            "section_source": source,
            "outer_left_radius_mm": radii["outer_left"],
            "outer_right_radius_mm": radii["outer_right"],
            "inner_left_radius_mm": radii["inner_left"],
            "inner_right_radius_mm": radii["inner_right"],
        }
        for side in ("left", "right"):
            ro = radii[f"outer_{side}"]
            ri = radii[f"inner_{side}"]
            row[f"wall_thickness_{side}_radial_mm"] = None if ro is None or ri is None else float(ro - ri)
            row[f"outer_inner_order_{side}_ok"] = None if ro is None or ri is None else bool(ro >= ri)
        rows.append(row)
    return rows


# -----------------------------------------------------------------------------
# Geometry writers
# -----------------------------------------------------------------------------

def _safe_name(text: str) -> str:
    return "".join(c if c.isalnum() or c in "_-" else "_" for c in text)


def write_obj_polylines(path: Path, named_polylines: list[tuple[str, np.ndarray, bool]], export_unit: str = "mm") -> None:
    """Write OBJ vertices and line elements. Input points are mm; output coordinates use export_unit."""
    scale = export_scale_from_mm(export_unit)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "# PotterySectionGuide OBJ polyline export",
        f"# version {__version__}",
        f"# coordinate_unit: {export_unit}",
        f"# scale_from_internal_mm: {scale:.12g}",
    ]
    v_index = 1
    for name, pts, closed in named_polylines:
        pts = geometry_for_export(pts, export_unit)
        if pts.ndim != 2 or pts.shape[1] != 3 or len(pts) < 2:
            continue
        lines.append(f"g {_safe_name(name)}")
        for x, y, z in pts:
            lines.append(f"v {x:.9f} {y:.9f} {z:.9f}")
        ids = list(range(v_index, v_index + len(pts)))
        if closed:
            ids.append(ids[0])
        lines.append("l " + " ".join(map(str, ids)))
        v_index += len(pts)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_ply_xyz(path: Path, points: np.ndarray, export_unit: str = "mm") -> None:
    pts = geometry_for_export(points, export_unit)
    scale = export_scale_from_mm(export_unit)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write(f"comment PotterySectionGuide point cloud; coordinate_unit={export_unit}\n")
        f.write(f"comment scale_from_internal_mm={scale:.12g}\n")
        f.write(f"element vertex {len(pts)}\n")
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write("end_header\n")
        for x, y, z in pts:
            f.write(f"{x:.9f} {y:.9f} {z:.9f}\n")


def write_combined_guide_ply(path: Path, records: list[RingRecord], geometries: dict[int, np.ndarray], export_unit: str = "mm") -> None:
    total = sum(len(geometries[r.ring_id]) for r in records)
    scale = export_scale_from_mm(export_unit)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write(f"comment PotterySectionGuide combined ring cloud; coordinate_unit={export_unit}\n")
        f.write(f"comment scale_from_internal_mm={scale:.12g}\n")
        f.write(f"element vertex {total}\n")
        f.write("property float x\nproperty float y\nproperty float z\n")
        f.write("property int ring_id\nproperty uchar surface_id\nproperty uchar side_id\n")
        f.write("end_header\n")
        for rec in records:
            sid = 1 if rec.surface == "outer" else 2
            sideid = 1 if rec.side == "left" else 2
            for x, y, z in geometry_for_export(geometries[rec.ring_id], export_unit):
                f.write(f"{x:.9f} {y:.9f} {z:.9f} {rec.ring_id:d} {sid:d} {sideid:d}\n")


# -----------------------------------------------------------------------------
# CSV / QA outputs
# -----------------------------------------------------------------------------

def write_profiles_csv(path: Path, profiles: dict[str, PreparedProfile]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["profile_key", "point_order", "x_source_px", "y_source_px", "signed_x_from_axis_mm", "radius_mm", "z_mm"])
        for key in PROFILE_KEYS:
            p = profiles.get(key)
            if p is None:
                continue
            # Preserve raw click order, recomputing projected coordinates by nearest
            # prepared samples is inappropriate. The export caller writes a separate
            # raw point table from source points; here prepared sorted samples are used.
            for i, (z, u, r) in enumerate(zip(p.z_mm, p.signed_u_mm, p.radius_mm), start=1):
                w.writerow([key, i, "", "", float(u), float(r), float(z)])


def write_raw_profile_csv(path: Path, profiles_px: dict[str, list[tuple[float, float]]], axis: Axis2D,
                          origin_px: tuple[float, float], mm_per_pixel: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["profile_key", "point_order", "x_source_px", "y_source_px", "signed_x_from_axis_mm", "radius_mm", "z_mm"])
        for key in PROFILE_KEYS:
            for i, p in enumerate(profiles_px.get(key, []), start=1):
                u, z = project_point_local_mm(p, axis, origin_px, mm_per_pixel)
                w.writerow([key, i, float(p[0]), float(p[1]), float(u), abs(float(u)), float(z)])


def write_sections_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "z_mm", "section_source",
        "outer_left_radius_mm", "outer_right_radius_mm",
        "inner_left_radius_mm", "inner_right_radius_mm",
        "wall_thickness_left_radial_mm", "outer_inner_order_left_ok",
        "wall_thickness_right_radial_mm", "outer_inner_order_right_ok",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_rings_csv(path: Path, records: list[RingRecord]) -> None:
    fields = ["ring_id", "z_mm", "profile_key", "surface", "side", "radius_mm", "section_source", "point_count"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow(asdict(r))


def axis_line_points_3d(zmin_mm: float, zmax_mm: float) -> np.ndarray:
    pad = max(5.0, 0.03 * max(1.0, zmax_mm - zmin_mm))
    return np.array([[0.0, 0.0, zmin_mm - pad], [0.0, 0.0, zmax_mm + pad]], dtype=float)


def prepared_profile_polyline_3d(profile: PreparedProfile) -> np.ndarray:
    return np.column_stack((profile.signed_u_mm, np.zeros(len(profile.z_mm)), profile.z_mm))


def _line_segment_for_z_on_image(axis: Axis2D, origin_px: tuple[float, float], z_mm: float,
                                 mm_per_pixel: float, image_size: tuple[int, int]) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return a long line perpendicular to the axis at local z, for QA overlay."""
    e_r, e_z = axis.basis()
    center = np.asarray(origin_px, dtype=float) + e_z * (float(z_mm) / mm_per_pixel)
    half = 1.5 * math.hypot(*image_size)
    p1 = center - e_r * half
    p2 = center + e_r * half
    return (float(p1[0]), float(p1[1])), (float(p2[0]), float(p2[1]))


def save_qc_overlay(source_image: Image.Image, path: Path, calibration: Calibration, axis: Axis2D,
                    origin_px: tuple[float, float], profiles_px: dict[str, list[tuple[float, float]]],
                    section_levels: list[tuple[float, str]]) -> None:
    im = source_image.copy().convert("RGB")
    draw = ImageDraw.Draw(im)
    width = max(2, int(round(max(im.size) / 900)))
    rr = max(3, 2 * width)

    # Calibration
    draw.line([calibration.p1_px, calibration.p2_px], fill=(180, 40, 40), width=width)
    for x, y in (calibration.p1_px, calibration.p2_px):
        draw.ellipse((x-rr, y-rr, x+rr, y+rr), outline=(180, 40, 40), width=width)

    # Axis
    p0 = np.asarray(axis.p_bottom_px, float)
    p1 = np.asarray(axis.p_top_px, float)
    v = p1 - p0
    n = np.linalg.norm(v)
    if n > 0:
        v = v / n
        extension = 2.0 * math.hypot(*im.size)
        a = p0 - v * extension
        b = p0 + v * extension
        draw.line([(float(a[0]), float(a[1])), (float(b[0]), float(b[1]))], fill=(30, 100, 210), width=width)

    # Origin
    ox, oy = origin_px
    draw.ellipse((ox-rr*1.5, oy-rr*1.5, ox+rr*1.5, oy+rr*1.5), outline=(20, 20, 20), width=width)

    # Profiles
    colors = {
        "outer_left": (10, 110, 50),
        "outer_right": (20, 140, 80),
        "inner_left": (190, 100, 20),
        "inner_right": (210, 140, 30),
    }
    for key in PROFILE_KEYS:
        pts = profiles_px.get(key, [])
        col = colors[key]
        if len(pts) >= 2:
            draw.line(pts, fill=col, width=max(2, width))
        for x, y in pts:
            draw.ellipse((x-rr, y-rr, x+rr, y+rr), fill=col)

    # Requested section levels
    for z, source in section_levels:
        a, b = _line_segment_for_z_on_image(axis, origin_px, z, calibration.mm_per_pixel, im.size)
        fill = (120, 50, 160) if source in {"custom", "regular+custom"} else (100, 100, 100)
        draw.line([a, b], fill=fill, width=max(1, width // 2))

    try:
        font = ImageFont.load_default()
        draw.text((10, 10), f"scale: {calibration.mm_per_pixel:.6g} mm/px", fill=(0, 0, 0), font=font)
        draw.text((10, 28), f"axis: {axis.source}", fill=(0, 0, 0), font=font)
        draw.text((10, 46), f"sections: {len(section_levels)}", fill=(0, 0, 0), font=font)
    except Exception:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    im.save(path)


def save_3d_preview(path: Path, records: list[RingRecord], geometries: dict[int, np.ndarray],
                    profiles: dict[str, PreparedProfile], zmin_mm: float, zmax_mm: float) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    for rec in records:
        pts = geometries[rec.ring_id]
        closed = np.vstack([pts, pts[0]])
        ax.plot(closed[:, 0], closed[:, 1], closed[:, 2], linewidth=0.7)
    for p in profiles.values():
        pts = prepared_profile_polyline_3d(p)
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], linewidth=1.0)
    axispts = axis_line_points_3d(zmin_mm, zmax_mm)
    ax.plot(axispts[:, 0], axispts[:, 1], axispts[:, 2], linestyle="--", linewidth=1.0)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm]")
    ax.set_title("PotterySectionGuide 3-D guide preview")
    # Equal-ish data scaling.
    all_r = [r.radius_mm for r in records] or [1.0]
    rmax = max(all_r)
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(-rmax, rmax)
    ax.set_zlim(zmin_mm, zmax_mm)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main export workflow
# -----------------------------------------------------------------------------

def export_guide(output_dir: Path, source_path: Path, source_image: Image.Image, source_meta: dict,
                 calibration: Calibration, axis: Axis2D, profiles_px: dict[str, list[tuple[float, float]]],
                 custom_section_z_mm: list[float], regular_interval_mm: float | None,
                 regular_start_mm: float | None, regular_end_mm: float | None,
                 origin_px: tuple[float, float] | None, ring_points_count: int, export_unit: str = "mm") -> dict:
    export_unit = str(export_unit).lower()
    export_scale = export_scale_from_mm(export_unit)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rings_dir = output_dir / "rings"
    rings_dir.mkdir(exist_ok=True)

    available = {k: list(v) for k, v in profiles_px.items() if len(v) >= 2}
    if not available:
        raise ValueError("Digitize at least one outer/inner profile with two or more points.")

    if origin_px is None:
        origin_px = automatic_origin_from_profiles(axis, available)
        origin_source = "automatic_lowest_digitized_profile_projection"
    else:
        origin_px = project_click_to_axis(origin_px, axis)
        origin_source = "manual_axis_projection"

    profiles: dict[str, PreparedProfile] = {
        key: prepare_profile(key, pts, axis, origin_px, calibration.mm_per_pixel)
        for key, pts in available.items()
    }
    zmin, zmax = profile_z_range(profiles)

    regular: list[float] = []
    if regular_interval_mm is not None:
        regular = regular_section_levels(zmin, zmax, regular_interval_mm, regular_start_mm, regular_end_mm)
    levels = merge_section_levels(regular, custom_section_z_mm)
    if not levels:
        raise ValueError("No section levels are defined. Add regular sections and/or custom section positions.")

    records, geometries = build_ring_records(profiles, levels, ring_points_count)
    if not records:
        raise ValueError("No rings could be generated because no requested level intersects a digitized profile range.")

    section_rows = section_summary_rows(profiles, levels)

    source_image.save(output_dir / "source_render.png")
    save_qc_overlay(source_image, output_dir / "guide_qc_overlay.png", calibration, axis, origin_px, available, levels)
    write_raw_profile_csv(output_dir / "profiles_raw.csv", available, axis, origin_px, calibration.mm_per_pixel)
    write_profiles_csv(output_dir / "profiles_prepared.csv", profiles)
    write_sections_csv(output_dir / "sections_summary.csv", section_rows)
    write_rings_csv(output_dir / "rings_summary.csv", records)

    # Combined OBJ: axis + profiles + every ring.
    obj_items: list[tuple[str, np.ndarray, bool]] = []
    obj_items.append(("rotation_axis", axis_line_points_3d(zmin, zmax), False))
    for key in PROFILE_KEYS:
        if key in profiles:
            obj_items.append((f"profile_{key}", prepared_profile_polyline_3d(profiles[key]), False))
    for rec in records:
        obj_items.append((f"ring_{rec.ring_id:04d}_{rec.profile_key}_z{rec.z_mm:+010.3f}mm", geometries[rec.ring_id], True))
    write_obj_polylines(output_dir / "guide_polylines.obj", obj_items, export_unit)

    # Rings-only OBJ is convenient when profiles are not wanted in CloudCompare.
    ring_items = [
        (f"ring_{r.ring_id:04d}_{r.profile_key}_z{r.z_mm:+010.3f}mm", geometries[r.ring_id], True)
        for r in records
    ]
    write_obj_polylines(output_dir / "guide_rings.obj", ring_items, export_unit)

    # Profiles and axis as a separate OBJ.
    profile_items = [("rotation_axis", axis_line_points_3d(zmin, zmax), False)]
    for key in PROFILE_KEYS:
        if key in profiles:
            profile_items.append((f"profile_{key}", prepared_profile_polyline_3d(profiles[key]), False))
    write_obj_polylines(output_dir / "guide_profiles_axis.obj", profile_items, export_unit)

    write_combined_guide_ply(output_dir / "guide_rings_cloud.ply", records, geometries, export_unit)

    for rec in records:
        base = f"ring_{rec.ring_id:04d}_{rec.profile_key}_z{rec.z_mm:+010.3f}mm_r{rec.radius_mm:09.3f}mm"
        pts = geometries[rec.ring_id]
        write_ply_xyz(rings_dir / f"{base}.ply", pts, export_unit)
        write_obj_polylines(rings_dir / f"{base}.obj", [(base, pts, True)], export_unit)

    save_3d_preview(output_dir / "guide_3d_preview.png", records, geometries, profiles, zmin, zmax)

    # QA warnings, especially outer/inner order when both are available.
    warnings: list[str] = []
    for row in section_rows:
        for side in ("left", "right"):
            ok = row[f"outer_inner_order_{side}_ok"]
            if ok is False:
                warnings.append(
                    f"At Z={row['z_mm']:.3f} mm, {side} outer radius is smaller than inner radius; "
                    "check digitizing, axis, or drawing interpretation."
                )

    metadata = {
        "program": "PotterySectionGuide",
        "version": __version__,
        "purpose": "3-D circular placement guides reconstructed from a 2-D pottery longitudinal-section drawing",
        "source_file": str(source_path),
        "source_sha256": sha256_file(source_path),
        "source": source_meta,
        "rendered_image_size_px": [int(source_image.width), int(source_image.height)],
        "coordinate_system": {
            "internal_measurement_unit": "mm",
            "export_geometry_unit": export_unit,
            "export_scale_from_mm": export_scale,
            "X": "signed perpendicular direction from rotation axis within the source drawing plane",
            "Y": "perpendicular to the source drawing plane",
            "Z": "upward along the rotation axis",
            "drawing_plane": "Y=0",
            "origin_px": [float(origin_px[0]), float(origin_px[1])],
            "origin_source": origin_source,
        },
        "calibration": asdict(calibration),
        "rotation_axis": asdict(axis),
        "radius_definition": (
            "For each digitized profile independently, radius at Z is the absolute perpendicular distance "
            "from the rotation axis to the interpolated profile line. Left/right radii are not averaged."
        ),
        "profiles_available": list(profiles.keys()),
        "profile_model": "piecewise-linear interpolation in local axial Z",
        "regular_sections": {
            "enabled": regular_interval_mm is not None,
            "interval_mm": regular_interval_mm,
            "start_mm": regular_start_mm,
            "end_mm": regular_end_mm,
        },
        "custom_section_z_mm": [float(v) for v in custom_section_z_mm],
        "section_levels": [{"z_mm": float(z), "source": src} for z, src in levels],
        "ring_sampling": {
            "point_count_per_circle": int(ring_points_count),
            "angular_step_deg": float(360.0 / ring_points_count),
            "geometry": "true circle in plane Z=constant",
        },
        "ring_count": len(records),
        "geometry_export": {
            "unit": export_unit,
            "scale_from_internal_mm": export_scale,
            "note": "OBJ/PLY formats do not enforce physical units; numeric coordinates are pre-scaled to the selected export unit for CloudCompare use.",
        },
        "warnings": warnings,
        "outputs": {
            "combined_obj": "guide_polylines.obj",
            "rings_obj": "guide_rings.obj",
            "profiles_axis_obj": "guide_profiles_axis.obj",
            "combined_ring_cloud": "guide_rings_cloud.ply",
            "individual_rings": "rings/",
            "sections_csv": "sections_summary.csv",
            "rings_csv": "rings_summary.csv",
            "qa_overlay": "guide_qc_overlay.png",
            "preview": "guide_3d_preview.png",
        },
    }
    (output_dir / "guide_metadata.json").write_text(
        json.dumps(jsonable(metadata), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return {
        "output_dir": output_dir,
        "metadata": metadata,
        "profiles": profiles,
        "section_levels": levels,
        "records": records,
        "geometries": geometries,
    }


# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

PROFILE_LABELS_JA = {
    "outer_left": "外面 左",
    "outer_right": "外面 右",
    "inner_left": "内面 左",
    "inner_right": "内面 右",
}
PROFILE_COLORS = {
    "outer_left": "#087334",
    "outer_right": "#15945a",
    "inner_left": "#bf6414",
    "inner_right": "#d88d1c",
}


class SectionGuideGUI:
    def __init__(self, root, source_path: Path | None, page_number: int, render_dpi: float,
                 output_dir: Path | None, ring_points_count: int, export_unit: str = "mm"):
        import tkinter as tk
        from tkinter import filedialog, messagebox, simpledialog, ttk

        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.simpledialog = simpledialog
        self.ttk = ttk
        self.root = root
        self.root.title(f"PotterySectionGuide v{__version__}")
        self.root.geometry("1380x900")

        self.source_path = Path(source_path) if source_path else None
        self.page_number = int(page_number)
        self.render_dpi = float(render_dpi)
        self.output_dir_arg = Path(output_dir) if output_dir else None
        self.ring_points_count = int(ring_points_count)
        self.initial_export_unit = str(export_unit).lower()
        export_scale_from_mm(self.initial_export_unit)

        self.source_image: Image.Image | None = None
        self.source_meta: dict = {}
        self.tk_image = None
        self.zoom = 1.0
        self.mode: str | None = None

        self.scale_points: list[tuple[float, float]] = []
        self.calibration: Calibration | None = None
        self.axis_points: list[tuple[float, float]] = []
        self.axis: Axis2D | None = None
        self.origin_click_px: tuple[float, float] | None = None
        self.profiles_px: dict[str, list[tuple[float, float]]] = {k: [] for k in PROFILE_KEYS}
        self.custom_section_clicks_px: list[tuple[float, float]] = []
        self.custom_section_z_mm: list[float] = []
        self.regular_interval_mm: float | None = 10.0
        self.regular_start_mm: float | None = None
        self.regular_end_mm: float | None = None

        self.status_var = tk.StringVar(value="断面実測図を開いてください。")
        self.scale_var = tk.StringVar(value="scale: 未設定")
        self.axis_var = tk.StringVar(value="axis: 未設定")
        self.origin_var = tk.StringVar(value="origin: 自動")
        self.interval_var = tk.StringVar(value="10")
        self.start_var = tk.StringVar(value="")
        self.end_var = tk.StringVar(value="")
        self.zoom_var = tk.StringVar(value="zoom: 100%")
        self.export_unit_var = tk.StringVar(value=self.initial_export_unit)

        # Toolbar row 1
        top1 = ttk.Frame(root)
        top1.pack(side="top", fill="x")
        buttons1 = [
            ("図面を開く", self.open_source_dialog),
            ("1 縮尺 2点", self.start_scale),
            ("2 軸 自動検出", self.auto_detect_axis),
            ("3 軸 2点指定", self.start_axis_manual),
            ("4 原点 Z=0 (任意)", self.start_origin),
            ("1点戻す", self.undo),
        ]
        for text, cmd in buttons1:
            ttk.Button(top1, text=text, command=cmd).pack(side="left", padx=3, pady=4)
        ttk.Button(top1, text="−", width=3, command=lambda: self.change_zoom(0.8)).pack(side="left", padx=2)
        ttk.Button(top1, text="＋", width=3, command=lambda: self.change_zoom(1.25)).pack(side="left", padx=2)
        ttk.Button(top1, text="100%", width=6, command=self.zoom_100).pack(side="left", padx=2)
        ttk.Button(top1, text="Fit", width=5, command=self.fit_to_window).pack(side="left", padx=2)
        ttk.Label(top1, textvariable=self.zoom_var).pack(side="left", padx=6)
        ttk.Label(top1, textvariable=self.scale_var).pack(side="left", padx=8)
        ttk.Label(top1, textvariable=self.axis_var).pack(side="left", padx=8)
        ttk.Label(top1, textvariable=self.origin_var).pack(side="left", padx=8)

        # Toolbar row 2: profile digitizing
        top2 = ttk.Frame(root)
        top2.pack(side="top", fill="x")
        ttk.Label(top2, text="プロファイル:").pack(side="left", padx=(5, 2))
        for key in PROFILE_KEYS:
            ttk.Button(top2, text=PROFILE_LABELS_JA[key], command=lambda k=key: self.start_profile(k)).pack(side="left", padx=3, pady=4)
        ttk.Button(top2, text="現在のプロファイルを消去", command=self.clear_current_profile).pack(side="left", padx=8)

        # Toolbar row 3: section definitions and export
        top3 = ttk.Frame(root)
        top3.pack(side="top", fill="x")
        ttk.Label(top3, text="等間隔 [mm]").pack(side="left", padx=(5, 2))
        ttk.Entry(top3, width=7, textvariable=self.interval_var).pack(side="left")
        ttk.Label(top3, text="開始Z").pack(side="left", padx=(8, 2))
        ttk.Entry(top3, width=8, textvariable=self.start_var).pack(side="left")
        ttk.Label(top3, text="終了Z").pack(side="left", padx=(8, 2))
        ttk.Entry(top3, width=8, textvariable=self.end_var).pack(side="left")
        ttk.Button(top3, text="5 等間隔を設定", command=self.set_regular_sections).pack(side="left", padx=5)
        ttk.Button(top3, text="6 指定位置を追加", command=self.start_custom_sections).pack(side="left", padx=5)
        ttk.Button(top3, text="指定位置を消去", command=self.clear_custom_sections).pack(side="left", padx=5)
        ttk.Label(top3, text="OBJ/PLY出力単位").pack(side="left", padx=(10, 2))
        ttk.Combobox(top3, width=5, textvariable=self.export_unit_var, values=("mm", "cm", "m"), state="readonly").pack(side="left")
        ttk.Button(top3, text="7 3Dガイド生成", command=self.generate_guide).pack(side="left", padx=12)

        body = ttk.Frame(root)
        body.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(body, background="white")
        xbar = ttk.Scrollbar(body, orient="horizontal", command=self.canvas.xview)
        ybar = ttk.Scrollbar(body, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=xbar.set, yscrollcommand=ybar.set)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        ybar.grid(row=0, column=1, sticky="ns")
        xbar.grid(row=1, column=0, sticky="ew")
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)

        ttk.Label(root, textvariable=self.status_var, anchor="w").pack(side="bottom", fill="x", padx=6, pady=3)
        self.canvas.bind("<Button-1>", self.on_click)
        self.canvas.bind("<Button-3>", lambda e: self.finish_current())
        self.root.bind("<Return>", lambda e: self.finish_current())
        self.root.bind("<Control-z>", lambda e: self.undo())
        # Zoom keeps the source point under the cursor fixed.  Middle-button drag
        # pans the current view; scrollbars remain available on trackpads/mice
        # without a middle button.
        self.canvas.bind("<MouseWheel>", self.on_mousewheel)
        self.canvas.bind("<Button-4>", lambda e: self.zoom_at_event(e, 1.12))
        self.canvas.bind("<Button-5>", lambda e: self.zoom_at_event(e, 0.89))
        self.canvas.bind("<ButtonPress-2>", self.start_pan)
        self.canvas.bind("<B2-Motion>", self.drag_pan)

        if self.source_path:
            self.load_source(self.source_path)

    # ----- source and coordinate helpers -----
    def open_source_dialog(self):
        filename = self.filedialog.askopenfilename(
            title="土器断面実測図を選択",
            filetypes=[
                ("対応図面", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.webp *.svg *.pdf"),
                ("画像", "*.png *.jpg *.jpeg *.tif *.tiff *.bmp *.webp"),
                ("SVG", "*.svg"), ("PDF", "*.pdf"), ("すべて", "*.*")
            ],
        )
        if filename:
            self.load_source(Path(filename))

    def load_source(self, path: Path):
        page = self.page_number
        if path.suffix.lower() == ".pdf":
            try:
                import pymupdf
                d = pymupdf.open(str(path))
                if d.page_count > 1:
                    val = self.simpledialog.askinteger(
                        "PDFページ", f"ページ番号 (1-{d.page_count})",
                        initialvalue=min(page, d.page_count), minvalue=1, maxvalue=d.page_count,
                    )
                    if val is None:
                        return
                    page = int(val)
            except Exception:
                pass
        try:
            image, meta = load_source_image(path, page_number=page, render_dpi=self.render_dpi)
        except Exception as exc:
            self.messagebox.showerror("読み込みエラー", str(exc))
            return
        self.source_path = path
        self.page_number = page
        self.source_image = image
        self.source_meta = meta
        self.scale_points = []
        self.calibration = None
        self.axis_points = []
        self.axis = None
        self.origin_click_px = None
        self.profiles_px = {k: [] for k in PROFILE_KEYS}
        self.custom_section_clicks_px = []
        self.custom_section_z_mm = []
        self.scale_var.set("scale: 未設定")
        self.axis_var.set("axis: 未設定")
        self.origin_var.set("origin: 自動")
        self.zoom = min(1.0, 1150 / max(1, image.width), 680 / max(1, image.height))
        self.zoom = max(self.zoom, 0.05)
        self.zoom_var.set(f"zoom: {self.zoom * 100:.0f}%")
        self.redraw()
        self.canvas.xview_moveto(0.0)
        self.canvas.yview_moveto(0.0)
        self.status_var.set("1) 縮尺を2点で指定してください。次に回転軸を自動検出または2点指定します。")

    def display_xy_to_source(self, event) -> tuple[float, float]:
        x = self.canvas.canvasx(event.x) / self.zoom
        y = self.canvas.canvasy(event.y) / self.zoom
        return float(x), float(y)

    def source_to_display(self, p):
        return float(p[0] * self.zoom), float(p[1] * self.zoom)

    # ----- modes -----
    def start_scale(self):
        if self.source_image is None:
            return
        self.mode = "scale"
        self.scale_points = []
        self.status_var.set("縮尺: 既知距離の1点目・2点目を順にクリックしてください。")
        self.redraw()

    def auto_detect_axis(self):
        if self.source_image is None:
            return
        try:
            self.axis = detect_rotation_axis(self.source_image)
            self.axis_points = []
            self.axis_var.set("axis: 自動候補")
            self.status_var.set("自動検出した回転軸候補を青線で表示しました。誤っている場合は『軸 2点指定』で置き換えてください。")
            self.redraw()
        except Exception as exc:
            self.messagebox.showwarning("回転軸の自動検出", str(exc))

    def start_axis_manual(self):
        if self.source_image is None:
            return
        self.mode = "axis_manual"
        self.axis_points = []
        self.status_var.set("回転軸: 軸上の離れた2点をクリックしてください。順序は問いません。")
        self.redraw()

    def start_origin(self):
        if self.source_image is None or self.axis is None:
            self.messagebox.showwarning("回転軸未設定", "先に回転軸を設定してください。")
            return
        self.mode = "origin"
        self.status_var.set("Z=0にしたい位置をクリックしてください。クリック点を回転軸へ直交投影して原点にします。")

    def start_profile(self, key: str):
        if self.source_image is None:
            return
        if self.axis is None:
            self.messagebox.showwarning("回転軸未設定", "先に回転軸を設定してください。")
            return
        self.mode = key
        self.status_var.set(
            f"{PROFILE_LABELS_JA[key]}: 線に沿って順にクリック。Enter/右クリックで終了。"
            "判読可能なouter/innerはそれぞれ独立に取得してください。"
        )

    def start_custom_sections(self):
        if self.source_image is None or self.axis is None or self.calibration is None:
            self.messagebox.showwarning("未設定", "縮尺と回転軸を先に設定してください。")
            return
        self.mode = "custom_section"
        self.status_var.set("指定横断面: 必要な高さをクリックしてください。クリック位置の軸方向Zを使用します。Enter/右クリックで終了。")

    def finish_current(self):
        if self.mode in PROFILE_KEYS:
            n = len(self.profiles_px[self.mode])
            self.status_var.set(f"{PROFILE_LABELS_JA[self.mode]}を{n}点で確定。")
        elif self.mode == "custom_section":
            self.status_var.set(f"指定位置を{len(self.custom_section_z_mm)}個保持しています。")
        self.mode = None

    def clear_current_profile(self):
        if self.mode in PROFILE_KEYS:
            key = self.mode
        else:
            # Ask user which profile to clear.
            val = self.simpledialog.askstring(
                "プロファイル消去",
                "消去するキーを入力:\nouter_left / outer_right / inner_left / inner_right",
            )
            key = val.strip() if val else ""
        if key in PROFILE_KEYS:
            self.profiles_px[key] = []
            self.redraw()
            self.status_var.set(f"{PROFILE_LABELS_JA[key]}を消去しました。")

    def clear_custom_sections(self):
        self.custom_section_clicks_px = []
        self.custom_section_z_mm = []
        self.redraw()
        self.status_var.set("指定位置の横断面を消去しました。")

    def set_regular_sections(self):
        try:
            interval = float(self.interval_var.get())
            if interval <= 0:
                raise ValueError
        except Exception:
            self.messagebox.showwarning("入力エラー", "等間隔 [mm] は正の数で指定してください。")
            return
        self.regular_interval_mm = interval
        try:
            self.regular_start_mm = float(self.start_var.get()) if self.start_var.get().strip() else None
            self.regular_end_mm = float(self.end_var.get()) if self.end_var.get().strip() else None
        except Exception:
            self.messagebox.showwarning("入力エラー", "開始Z / 終了Zは空欄または数値 [mm] にしてください。")
            return
        self.status_var.set(
            f"等間隔横断面: {interval:g} mm"
            + (f", start={self.regular_start_mm:g}" if self.regular_start_mm is not None else "")
            + (f", end={self.regular_end_mm:g}" if self.regular_end_mm is not None else "")
        )
        self.redraw()

    # ----- click handling -----
    def on_click(self, event):
        if self.source_image is None:
            return
        x, y = self.display_xy_to_source(event)
        if not (0 <= x < self.source_image.width and 0 <= y < self.source_image.height):
            return

        if self.mode == "scale":
            self.scale_points.append((x, y))
            if len(self.scale_points) == 2:
                p1, p2 = self.scale_points
                pix = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
                known = self.simpledialog.askfloat("実寸", "選択した2点間の実寸 [mm]", minvalue=1e-9)
                if known is None or pix <= 0:
                    self.scale_points = []
                else:
                    self.calibration = Calibration(p1, p2, float(known), float(pix), float(known / pix))
                    self.scale_var.set(f"scale: {self.calibration.mm_per_pixel:.6g} mm/px")
                    self.status_var.set("2) 回転軸を自動検出するか、2点で手動指定してください。")
                self.mode = None

        elif self.mode == "axis_manual":
            self.axis_points.append((x, y))
            if len(self.axis_points) == 2:
                try:
                    self.axis = _axis_from_two_points(self.axis_points[0], self.axis_points[1], source="manual_two_point")
                    self.axis_var.set("axis: 手動2点")
                    self.status_var.set("回転軸を確定しました。outer / innerの判読可能な線をデジタイズしてください。")
                except Exception as exc:
                    self.messagebox.showerror("回転軸", str(exc))
                    self.axis_points = []
                self.mode = None

        elif self.mode == "origin":
            if self.axis is not None:
                self.origin_click_px = project_click_to_axis((x, y), self.axis)
                self.origin_var.set("origin: 手動")
                self.status_var.set("Z=0原点を設定しました。")
            self.mode = None

        elif self.mode in PROFILE_KEYS:
            self.profiles_px[self.mode].append((x, y))

        elif self.mode == "custom_section":
            if self.axis is None or self.calibration is None:
                return
            # Origin for interactive custom Z: use manual origin if available;
            # otherwise use provisional automatic origin from already digitized profiles.
            try:
                available = {k: v for k, v in self.profiles_px.items() if len(v) >= 2}
                if self.origin_click_px is not None:
                    origin = project_click_to_axis(self.origin_click_px, self.axis)
                else:
                    origin = automatic_origin_from_profiles(self.axis, available)
                _, z = project_point_local_mm((x, y), self.axis, origin, self.calibration.mm_per_pixel)
                self.custom_section_clicks_px.append((x, y))
                self.custom_section_z_mm.append(float(z))
                self.status_var.set(f"指定横断面 Z={z:.3f} mm を追加。現在 {len(self.custom_section_z_mm)} 個。")
            except Exception as exc:
                self.messagebox.showwarning("指定位置", f"先に少なくとも1本のプロファイルを2点以上デジタイズしてください。\n{exc}")

        self.redraw()

    def undo(self):
        if self.mode in PROFILE_KEYS and self.profiles_px[self.mode]:
            self.profiles_px[self.mode].pop()
        elif self.mode == "axis_manual" and self.axis_points:
            self.axis_points.pop()
        elif self.mode == "scale" and self.scale_points:
            self.scale_points.pop()
        elif self.mode == "custom_section" and self.custom_section_z_mm:
            self.custom_section_z_mm.pop()
            if self.custom_section_clicks_px:
                self.custom_section_clicks_px.pop()
        else:
            # Most recent non-empty profile in fixed priority order.
            for key in reversed(PROFILE_KEYS):
                if self.profiles_px[key]:
                    self.profiles_px[key].pop()
                    break
        self.redraw()

    # ----- view -----
    def _viewport_center(self) -> tuple[float, float]:
        """Return viewport-center coordinates in canvas-window pixels."""
        self.root.update_idletasks()
        return (max(1.0, self.canvas.winfo_width()) / 2.0,
                max(1.0, self.canvas.winfo_height()) / 2.0)

    def _set_zoom_preserve_anchor(self, new_zoom: float, anchor_view_xy: tuple[float, float]):
        """Set display zoom while keeping one source-image point under the anchor.

        `anchor_view_xy` is expressed in visible canvas-widget coordinates (for
        example the mouse cursor).  Digitized/calibration coordinates remain in
        original source-image pixels and are therefore numerically unchanged by
        zooming.
        """
        if self.source_image is None:
            return
        old_zoom = float(self.zoom)
        new_zoom = max(0.03, min(8.0, float(new_zoom)))
        if abs(new_zoom - old_zoom) <= 1e-12:
            return

        ax, ay = map(float, anchor_view_xy)
        # Source pixel located under the chosen visible anchor before zoom.
        source_x = float(self.canvas.canvasx(ax)) / old_zoom
        source_y = float(self.canvas.canvasy(ay)) / old_zoom

        self.zoom = new_zoom
        self.zoom_var.set(f"zoom: {self.zoom * 100:.0f}%")
        self.redraw()
        self.root.update_idletasks()

        content_w = max(1.0, float(self.source_image.width) * self.zoom)
        content_h = max(1.0, float(self.source_image.height) * self.zoom)
        left = source_x * self.zoom - ax
        top = source_y * self.zoom - ay
        self.canvas.xview_moveto(max(0.0, min(1.0, left / content_w)))
        self.canvas.yview_moveto(max(0.0, min(1.0, top / content_h)))

    def change_zoom(self, factor):
        """Zoom from toolbar buttons, anchored at the visible viewport centre."""
        if self.source_image is None:
            return
        self._set_zoom_preserve_anchor(self.zoom * float(factor), self._viewport_center())

    def zoom_at_event(self, event, factor: float):
        """Cursor-centred zoom used by mouse-wheel events."""
        if self.source_image is None:
            return "break"
        self._set_zoom_preserve_anchor(self.zoom * float(factor), (event.x, event.y))
        return "break"

    def on_mousewheel(self, event):
        if event.delta == 0:
            return "break"
        return self.zoom_at_event(event, 1.12 if event.delta > 0 else 0.89)

    def zoom_100(self):
        """Set 1 source pixel = 1 display pixel, preserving the current view centre."""
        if self.source_image is None:
            return
        self._set_zoom_preserve_anchor(1.0, self._viewport_center())

    def fit_to_window(self):
        """Fit the entire source image into the current canvas viewport."""
        if self.source_image is None:
            return
        self.root.update_idletasks()
        vw = max(20.0, float(self.canvas.winfo_width()) - 4.0)
        vh = max(20.0, float(self.canvas.winfo_height()) - 4.0)
        fit = min(vw / max(1.0, float(self.source_image.width)),
                  vh / max(1.0, float(self.source_image.height)))
        self.zoom = max(0.03, min(8.0, fit))
        self.zoom_var.set(f"zoom: {self.zoom * 100:.0f}%")
        self.redraw()
        self.root.update_idletasks()
        self.canvas.xview_moveto(0.0)
        self.canvas.yview_moveto(0.0)

    def start_pan(self, event):
        """Start middle-button canvas panning."""
        self.canvas.scan_mark(event.x, event.y)
        return "break"

    def drag_pan(self, event):
        """Pan using middle-button drag without altering any source coordinates."""
        self.canvas.scan_dragto(event.x, event.y, gain=1)
        return "break"

    def redraw(self):
        self.canvas.delete("all")
        if self.source_image is None:
            return
        self.zoom_var.set(f"zoom: {self.zoom * 100:.0f}%")
        w = max(1, int(round(self.source_image.width * self.zoom)))
        h = max(1, int(round(self.source_image.height * self.zoom)))
        display = self.source_image.resize((w, h), Image.Resampling.LANCZOS)
        self.tk_image = ImageTk.PhotoImage(display)
        self.canvas.create_image(0, 0, image=self.tk_image, anchor="nw")
        self.canvas.configure(scrollregion=(0, 0, w, h))

        # scale
        if self.scale_points:
            dp = [self.source_to_display(p) for p in self.scale_points]
            for x, y in dp:
                self.canvas.create_oval(x-4, y-4, x+4, y+4, outline="red", width=2)
            if len(dp) == 2:
                self.canvas.create_line(*dp[0], *dp[1], fill="red", width=2)

        # axis candidate / manual
        if self.axis is not None:
            p0 = np.asarray(self.axis.p_bottom_px, float)
            p1 = np.asarray(self.axis.p_top_px, float)
            v = p1 - p0
            n = np.linalg.norm(v)
            if n > 0:
                v = v / n
                ext = 2.0 * math.hypot(self.source_image.width, self.source_image.height)
                a = self.source_to_display(p0 - v * ext)
                b = self.source_to_display(p0 + v * ext)
                self.canvas.create_line(*a, *b, fill="blue", width=2, dash=(6, 4))
        elif self.axis_points:
            dp = [self.source_to_display(p) for p in self.axis_points]
            for x, y in dp:
                self.canvas.create_oval(x-4, y-4, x+4, y+4, outline="blue", width=2)
            if len(dp) == 2:
                self.canvas.create_line(*dp[0], *dp[1], fill="blue", width=2)

        # origin
        if self.origin_click_px is not None:
            x, y = self.source_to_display(self.origin_click_px)
            self.canvas.create_oval(x-6, y-6, x+6, y+6, outline="black", width=2)

        # profiles
        for key in PROFILE_KEYS:
            pts = self.profiles_px[key]
            col = PROFILE_COLORS[key]
            dp = [self.source_to_display(p) for p in pts]
            if len(dp) >= 2:
                flat = [v for p in dp for v in p]
                self.canvas.create_line(*flat, fill=col, width=2)
            for x, y in dp:
                self.canvas.create_oval(x-3, y-3, x+3, y+3, fill=col, outline=col)

        # custom section click markers and perpendicular lines, when coordinate system is usable
        if self.axis is not None and self.calibration is not None:
            try:
                available = {k: v for k, v in self.profiles_px.items() if len(v) >= 2}
                if self.origin_click_px is not None:
                    origin = project_click_to_axis(self.origin_click_px, self.axis)
                elif available:
                    origin = automatic_origin_from_profiles(self.axis, available)
                else:
                    origin = None
                if origin is not None:
                    # custom
                    for z in self.custom_section_z_mm:
                        a, b = _line_segment_for_z_on_image(self.axis, origin, z, self.calibration.mm_per_pixel, self.source_image.size)
                        a = self.source_to_display(a); b = self.source_to_display(b)
                        self.canvas.create_line(*a, *b, fill="#7832a0", width=1)
                    # regular preview only after at least one prepared profile is possible
                    if self.regular_interval_mm is not None and available:
                        prepped = {
                            k: prepare_profile(k, pts, self.axis, origin, self.calibration.mm_per_pixel)
                            for k, pts in available.items()
                        }
                        zmin, zmax = profile_z_range(prepped)
                        levels = regular_section_levels(zmin, zmax, self.regular_interval_mm,
                                                        self.regular_start_mm, self.regular_end_mm)
                        for z in levels:
                            a, b = _line_segment_for_z_on_image(self.axis, origin, z, self.calibration.mm_per_pixel, self.source_image.size)
                            a = self.source_to_display(a); b = self.source_to_display(b)
                            self.canvas.create_line(*a, *b, fill="#777777", width=1)
            except Exception:
                pass

    # ----- export -----
    def generate_guide(self):
        if self.source_path is None or self.source_image is None:
            self.messagebox.showwarning("未設定", "図面を開いてください。")
            return
        if self.calibration is None:
            self.messagebox.showwarning("未設定", "縮尺を指定してください。")
            return
        if self.axis is None:
            self.messagebox.showwarning("未設定", "回転軸を自動検出または2点指定してください。")
            return
        available = {k: v for k, v in self.profiles_px.items() if len(v) >= 2}
        if not available:
            self.messagebox.showwarning("未設定", "少なくとも1本のouter/innerプロファイルを2点以上デジタイズしてください。")
            return

        # Read current regular-section fields.
        try:
            interval_text = self.interval_var.get().strip()
            interval = float(interval_text) if interval_text else None
            if interval is not None and interval <= 0:
                raise ValueError
            start = float(self.start_var.get()) if self.start_var.get().strip() else None
            end = float(self.end_var.get()) if self.end_var.get().strip() else None
        except Exception:
            self.messagebox.showwarning("入力エラー", "等間隔/開始Z/終了Zの値を確認してください。")
            return

        out = self.output_dir_arg or (self.source_path.parent / f"{self.source_path.stem}_SectionGuide")
        try:
            result = export_guide(
                out,
                self.source_path,
                self.source_image,
                self.source_meta,
                self.calibration,
                self.axis,
                self.profiles_px,
                list(self.custom_section_z_mm),
                interval,
                start,
                end,
                self.origin_click_px,
                self.ring_points_count,
                self.export_unit_var.get(),
            )
        except Exception as exc:
            self.messagebox.showerror("ガイド生成エラー", str(exc))
            return

        n = len(result["records"])
        warnings = result["metadata"].get("warnings", [])
        export_unit = result["metadata"]["geometry_export"]["unit"]
        self.status_var.set(f"3Dガイド生成完了: {n} rings / export unit={export_unit} - {out}")
        msg = (
            f"3Dガイドを生成しました。\n\n"
            f"ring数: {n}\n"
            f"横断面レベル数: {len(result['section_levels'])}\n"
            f"OBJ/PLY出力単位: {export_unit}\n"
            f"出力先:\n{out}\n\n"
            f"CloudCompareでは guide_polylines.obj または guide_rings_cloud.ply を読み込めます。"
        )
        if warnings:
            msg += f"\n\nQA warning: {len(warnings)}件。guide_metadata.json / sections_summary.csvを確認してください。"
        self.messagebox.showinfo("生成完了", msg)


# -----------------------------------------------------------------------------
# Self-test
# -----------------------------------------------------------------------------

def self_test() -> None:
    """Headless geometry self-test.

    Synthetic drawing:
      scale = 0.5 mm/px
      vertical axis x=200 px
      z=0 at y=300 px
      outer radius = 50 mm (x=100 / 300)
      inner radius = 45 mm (x=110 / 290)
      height = 100 mm
    """
    mmpp = 0.5
    axis = _axis_from_two_points((200.0, 300.0), (200.0, 100.0), source="self_test")
    origin = (200.0, 300.0)

    profiles_px = {
        "outer_left": [(100.0, 300.0), (100.0, 200.0), (100.0, 100.0)],
        "outer_right": [(300.0, 300.0), (300.0, 200.0), (300.0, 100.0)],
        "inner_left": [(110.0, 300.0), (110.0, 200.0), (110.0, 100.0)],
        "inner_right": [(290.0, 300.0), (290.0, 200.0), (290.0, 100.0)],
    }
    profiles = {k: prepare_profile(k, v, axis, origin, mmpp) for k, v in profiles_px.items()}
    levels = merge_section_levels(regular_section_levels(0.0, 100.0, 25.0), [37.5])
    records, geoms = build_ring_records(profiles, levels, point_count=180)

    # Every requested level should yield four rings.
    expected_levels = 6  # 0,25,37.5,50,75,100
    if len(records) != expected_levels * 4:
        raise SystemExit(f"SELF TEST FAILED: expected {expected_levels*4} rings, got {len(records)}")

    for rec in records:
        expected_r = 50.0 if rec.surface == "outer" else 45.0
        if abs(rec.radius_mm - expected_r) > 1e-9:
            raise SystemExit(f"SELF TEST FAILED: {rec.profile_key} radius {rec.radius_mm} != {expected_r}")
        pts = geoms[rec.ring_id]
        radial = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
        if float(np.max(np.abs(radial - expected_r))) > 1e-9:
            raise SystemExit("SELF TEST FAILED: circle geometry radius mismatch")
        if float(np.max(np.abs(pts[:, 2] - rec.z_mm))) > 1e-9:
            raise SystemExit("SELF TEST FAILED: ring Z mismatch")

    rows = section_summary_rows(profiles, [(50.0, "regular")])
    if abs(float(rows[0]["wall_thickness_left_radial_mm"]) - 5.0) > 1e-9:
        raise SystemExit("SELF TEST FAILED: wall thickness mismatch")

    # Tilted-axis distance test: point 10 px perpendicular from axis must be 5 mm.
    tilted = _axis_from_two_points((200.0, 300.0), (220.0, 100.0), source="self_test_tilt")
    e_r, _ = tilted.basis()
    q = np.asarray(tilted.p_bottom_px) + e_r * 10.0
    u, _ = project_point_local_mm(tuple(q), tilted, tilted.p_bottom_px, mmpp)
    if abs(abs(u) - 5.0) > 1e-9:
        raise SystemExit("SELF TEST FAILED: tilted-axis perpendicular distance mismatch")

    # Export-unit conversion test: 50 mm radius must become 0.05 m.
    sample_mm = np.array([[50.0, 0.0, 100.0]], dtype=float)
    sample_m = geometry_for_export(sample_mm, "m")
    if float(np.max(np.abs(sample_m - np.array([[0.05, 0.0, 0.1]])))) > 1e-12:
        raise SystemExit("SELF TEST FAILED: mm-to-m export scaling mismatch")

    print(f"PotterySectionGuide v{__version__} SELF TEST PASSED")
    print(f"levels          : {expected_levels}")
    print(f"rings           : {len(records)}")
    print("outer radius    : 50.0 mm")
    print("inner radius    : 45.0 mm")
    print("wall thickness  : 5.0 mm")
    print("tilted-axis distance test: passed")
    print("export unit scaling test: passed (50 mm -> 0.05 m)")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate CloudCompare-oriented 3-D circular guide geometry from a 2-D pottery longitudinal-section drawing."
        )
    )
    parser.add_argument("input", nargs="?", type=Path,
                        help="PNG/JPEG/TIFF/BMP/WEBP/SVG/PDF drawing. Omit to choose in GUI.")
    parser.add_argument("--page", type=int, default=1,
                        help="PDF page number, 1-based (default: 1; GUI asks for multipage PDF)")
    parser.add_argument("--render-dpi", type=float, default=180.0,
                        help="SVG/PDF rendering DPI for digitizing (default: 180)")
    parser.add_argument("--ring-points", type=int, default=180,
                        help="Points per 3-D circle (default: 180 = 2 degree spacing)")
    parser.add_argument("--output-dir", type=Path,
                        help="Output directory (default: <source>_SectionGuide)")
    parser.add_argument("--export-unit", choices=["mm", "cm", "m"], default="mm",
                        help="Numeric coordinate unit written to OBJ/PLY (default: mm). Internal measurements and CSV remain mm.")
    parser.add_argument("--self-test", action="store_true", help="Run headless numerical/geometry self-test and exit")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return
    if args.render_dpi <= 0:
        parser.error("--render-dpi must be positive")
    if args.ring_points < 8:
        parser.error("--ring-points must be >= 8")

    try:
        import tkinter as tk
    except ModuleNotFoundError:
        parser.exit(2, "ERROR: tkinter is required for GUI mode.\n")

    root = tk.Tk()
    SectionGuideGUI(root, args.input, args.page, args.render_dpi, args.output_dir, args.ring_points, args.export_unit)
    root.mainloop()


if __name__ == "__main__":
    main()


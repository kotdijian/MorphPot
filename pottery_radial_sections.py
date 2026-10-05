#!/usr/bin/env python3
"""
PotteryRadialSections
Version 0.9.0

Estimate a pottery rotation axis from a series of horizontal XY sections,
then extract longitudinal sections through that estimated Z-parallel axis.

Default workflow:
1. Divide vessel height into 20 equal intervals and section at interval centers.
2. Stitch each horizontal mesh-plane intersection into closed contours.
3. Identify outer and, when present, inner contours by enclosed area.
4. Resample contours uniformly by arc length.
5. Robustly fit ellipses to the contours.
6. Use inner ellipse centers by default (outer/both modes are available).
7. Remove 2-D center outliers using a MAD rule.
8. Average retained centers to define the final XY position of the Z-parallel axis.
9. Extract 30-degree longitudinal sections through the final axis.
10. Optionally estimate capacity using one or more selectable methods:
    single selected radial profile, optimized axisymmetric model,
    angular integration without axial-symmetry assumption, and horizontal ellipses.
11. Export CSV/JSON/PLY/PNG QC products.
12. Numerical profile/section integration is delegated to pottery_volume_core.py.

Important:
- +Z is assumed to be the vessel vertical direction.
- Axis fitting estimates only the XY axis position; it does NOT rotate/tilt the model.
- A center-drift/tilt diagnostic is reported but not automatically corrected.
- In outer mode, wall thickness is NOT used for axis estimation. It can be supplied
  in a later wall-thickness reconstruction stage; current 3D capacity methods require
  an inner surface to be present in the input mesh.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import trimesh
import pottery_volume_core as volume_core
from morphpot.section_overlay import export_section_overlay
from morphpot.rim_standardization import export_rim_standardization
from scipy.optimize import least_squares

__version__ = "0.9.0"

UNIT_SCALE_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}


def load_geometry(path: Path) -> trimesh.Trimesh:
    """Load geometry without texture-UV re-indexing for PLY files."""
    suffix = path.suffix.lower()
    if suffix == ".ply":
        from trimesh.exchange.ply import load_ply
        from trimesh.resolvers import FilePathResolver
        with path.open("rb") as f:
            data = load_ply(
                f,
                resolver=FilePathResolver(str(path.parent)),
                fix_texture=False,
                skip_materials=True,
                prefer_color="vertex",
            )
        mesh = trimesh.Trimesh(
            vertices=np.asarray(data["vertices"], dtype=np.float64),
            faces=np.asarray(data["faces"], dtype=np.int64),
            process=False,
            validate=False,
        )
        vc = data.get("vertex_colors")
        if vc is not None and len(vc) == len(mesh.vertices):
            try:
                mesh.visual.vertex_colors = np.asarray(vc, dtype=np.uint8)
            except Exception:
                pass
    else:
        loaded = trimesh.load(str(path), process=False, force="mesh")
        if isinstance(loaded, trimesh.Scene):
            geoms = [g for g in loaded.geometry.values()
                     if isinstance(g, trimesh.Trimesh) and len(g.vertices) and len(g.faces)]
            if not geoms:
                raise ValueError("No triangular mesh geometry found.")
            mesh = trimesh.util.concatenate(geoms)
        elif isinstance(loaded, trimesh.Trimesh):
            mesh = loaded
        else:
            raise TypeError(f"Unsupported mesh type: {type(loaded)}")
    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        raise ValueError("Input mesh has no vertices or faces.")
    if not np.isfinite(mesh.vertices).all():
        raise ValueError("Input mesh contains NaN or Inf coordinates.")
    return mesh


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def polygon_area_xy(points: np.ndarray) -> float:
    p = np.asarray(points, dtype=np.float64)
    if len(p) < 3:
        return 0.0
    if np.linalg.norm(p[0] - p[-1]) < 1e-12:
        p = p[:-1]
    x, y = p[:, 0], p[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def point_in_polygon_xy(point: np.ndarray, polygon: np.ndarray) -> bool:
    """Ray-crossing point-in-polygon test. Suitable for simple section contours."""
    x, y = map(float, point[:2])
    p = np.asarray(polygon, dtype=np.float64)
    if np.linalg.norm(p[0] - p[-1]) < 1e-12:
        p = p[:-1]
    inside = False
    n = len(p)
    j = n - 1
    for i in range(n):
        xi, yi = p[i, 0], p[i, 1]
        xj, yj = p[j, 0], p[j, 1]
        hit = ((yi > y) != (yj > y)) and (
            x < (xj - xi) * (y - yi) / ((yj - yi) if abs(yj - yi) > 1e-30 else 1e-30) + xi
        )
        if hit:
            inside = not inside
        j = i
    return inside


def resample_closed_polyline(points_xy: np.ndarray, spacing_native: float, min_points: int = 48) -> np.ndarray:
    """Uniform arc-length resampling of a closed XY contour."""
    p = np.asarray(points_xy, dtype=np.float64)
    if len(p) < 3:
        return np.empty((0, 2), dtype=np.float64)
    if np.linalg.norm(p[0] - p[-1]) > 1e-12:
        p = np.vstack([p, p[0]])
    seg_len = np.linalg.norm(np.diff(p, axis=0), axis=1)
    keep = seg_len > 1e-15
    if not np.any(keep):
        return np.empty((0, 2), dtype=np.float64)
    # retain all vertices; zero-length segments are skipped during interpolation
    cumulative = np.r_[0.0, np.cumsum(seg_len)]
    total = float(cumulative[-1])
    if total <= 0:
        return np.empty((0, 2), dtype=np.float64)
    n = max(min_points, int(math.ceil(total / max(spacing_native, 1e-15))))
    s_values = np.linspace(0.0, total, n, endpoint=False)
    out = np.empty((n, 2), dtype=np.float64)
    j = 0
    for k, s in enumerate(s_values):
        while j + 1 < len(cumulative) - 1 and cumulative[j + 1] <= s:
            j += 1
        denom = cumulative[j + 1] - cumulative[j]
        if denom <= 1e-15:
            out[k] = p[j]
        else:
            t = (s - cumulative[j]) / denom
            out[k] = p[j] * (1.0 - t) + p[j + 1] * t
    return out


def horizontal_section_segments(mesh: trimesh.Trimesh, z_native: float) -> np.ndarray:
    seg = trimesh.intersections.mesh_plane(
        mesh=mesh,
        plane_normal=np.array([0.0, 0.0, 1.0]),
        plane_origin=np.array([0.0, 0.0, float(z_native)]),
        return_faces=False,
    )
    seg = np.asarray(seg, dtype=np.float64)
    if seg.size == 0:
        return np.empty((0, 2, 3), dtype=np.float64)
    return seg.reshape((-1, 2, 3))


def closed_horizontal_contours(segments: np.ndarray, closure_tol_native: float) -> list[dict]:
    """Stitch mesh-plane segments and return closed contours sorted by area."""
    if len(segments) == 0:
        return []
    path = trimesh.load_path(np.asarray(segments, dtype=np.float64))
    contours = []
    for discrete in path.discrete:
        pts = np.asarray(discrete, dtype=np.float64)
        if len(pts) < 4:
            continue
        gap = float(np.linalg.norm(pts[0] - pts[-1]))
        if gap > closure_tol_native:
            continue
        xy = pts[:, :2]
        area = polygon_area_xy(xy)
        if area <= 0:
            continue
        contours.append({"points_xyz": pts, "points_xy": xy, "area_native2": area, "closure_gap_native": gap})
    contours.sort(key=lambda r: r["area_native2"], reverse=True)
    return contours


def select_outer_inner(contours: list[dict], min_inner_area_ratio: float = 0.05):
    if not contours:
        return None, None
    outer = contours[0]
    inner = None
    outer_area = float(outer["area_native2"])
    for candidate in contours[1:]:
        ratio = float(candidate["area_native2"] / outer_area) if outer_area > 0 else 0.0
        if ratio < min_inner_area_ratio:
            continue
        test_point = np.asarray(candidate["points_xy"], dtype=float)[0]
        if point_in_polygon_xy(test_point, np.asarray(outer["points_xy"], dtype=float)):
            inner = candidate
            break
    return outer, inner


def ellipse_residual(params: np.ndarray, points_xy: np.ndarray) -> np.ndarray:
    """Approximate geometric radial residual, in native length units."""
    cx, cy, log_a, log_b, phi = params
    a = math.exp(float(log_a))
    b = math.exp(float(log_b))
    cp, sp = math.cos(float(phi)), math.sin(float(phi))
    q = np.asarray(points_xy, dtype=np.float64) - np.array([cx, cy])
    x = cp * q[:, 0] + sp * q[:, 1]
    y = -sp * q[:, 0] + cp * q[:, 1]
    rho = np.sqrt((x / a) ** 2 + (y / b) ** 2 + 1e-24)
    return (rho - 1.0) * math.sqrt(a * b)


def ellipse_initial_guess(points_xy: np.ndarray) -> np.ndarray:
    p = np.asarray(points_xy, dtype=np.float64)
    center = np.median(p, axis=0)
    q = p - center
    cov = np.cov(q.T)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vecs = vecs[:, order]
    proj = q @ vecs
    a = max(float(np.percentile(np.abs(proj[:, 0]), 95)), 1e-12)
    b = max(float(np.percentile(np.abs(proj[:, 1]), 95)), 1e-12)
    phi = math.atan2(float(vecs[1, 0]), float(vecs[0, 0]))
    return np.array([center[0], center[1], math.log(a), math.log(b), phi], dtype=np.float64)


def robust_ellipse_fit(points_xy: np.ndarray, scale_to_mm: float) -> dict:
    p = np.asarray(points_xy, dtype=np.float64)
    if len(p) < 12:
        return {"valid": False, "reason": "too_few_points"}
    guess = ellipse_initial_guess(p)
    a0, b0 = math.exp(guess[2]), math.exp(guess[3])
    f_scale = max(min(a0, b0) * 0.01, 1e-9)
    try:
        fit1 = least_squares(
            ellipse_residual, guess, args=(p,), loss="soft_l1", f_scale=f_scale,
            max_nfev=3000,
        )
    except Exception as exc:
        return {"valid": False, "reason": f"least_squares_failed:{exc}"}

    r1 = np.abs(ellipse_residual(fit1.x, p))
    med = float(np.median(r1))
    mad = float(np.median(np.abs(r1 - med)))
    sigma = 1.4826 * mad
    threshold = med + 3.0 * sigma
    if threshold <= 1e-12:
        threshold = max(float(np.percentile(r1, 90)), 1e-12)
    inliers = r1 <= threshold
    fit = fit1
    if int(np.count_nonzero(inliers)) >= 12 and int(np.count_nonzero(inliers)) < len(p):
        try:
            fit = least_squares(
                ellipse_residual, fit1.x, args=(p[inliers],), loss="soft_l1",
                f_scale=max(threshold, 1e-9), max_nfev=3000,
            )
        except Exception:
            fit = fit1

    cx, cy, log_a, log_b, phi = map(float, fit.x)
    a, b = math.exp(log_a), math.exp(log_b)
    if not all(np.isfinite([cx, cy, a, b, phi])) or a <= 0 or b <= 0:
        return {"valid": False, "reason": "nonfinite_or_nonpositive_ellipse"}
    if b > a:
        a, b = b, a
        phi += math.pi / 2.0
    # normalize angle to [0,180)
    angle_deg = (math.degrees(phi) % 180.0 + 180.0) % 180.0
    residual_all = np.abs(ellipse_residual(np.array([cx, cy, math.log(a), math.log(b), math.radians(angle_deg)]), p))
    rmse_native = float(np.sqrt(np.mean(residual_all ** 2)))
    mean_radius_native = 0.5 * (a + b)
    ecc = math.sqrt(max(0.0, 1.0 - (b * b) / (a * a)))
    return {
        "valid": True,
        "reason": "",
        "center_x_native": cx,
        "center_y_native": cy,
        "center_x_mm": cx * scale_to_mm,
        "center_y_mm": cy * scale_to_mm,
        "semi_major_native": a,
        "semi_minor_native": b,
        "semi_major_mm": a * scale_to_mm,
        "semi_minor_mm": b * scale_to_mm,
        "ellipse_angle_deg": angle_deg,
        "eccentricity": ecc,
        "fit_rmse_native": rmse_native,
        "fit_rmse_mm": rmse_native * scale_to_mm,
        "fit_rmse_relative": rmse_native / mean_radius_native if mean_radius_native > 0 else float("nan"),
        "fit_points": int(len(p)),
        "fit_inlier_fraction": float(np.count_nonzero(inliers) / len(p)),
    }


def ellipse_points(fit: dict, z_native: float, count: int = 240) -> np.ndarray:
    if not fit or not fit.get("valid"):
        return np.empty((0, 3), dtype=np.float64)
    t = np.linspace(0.0, 2.0 * math.pi, count, endpoint=False)
    a = float(fit["semi_major_native"])
    b = float(fit["semi_minor_native"])
    phi = math.radians(float(fit["ellipse_angle_deg"]))
    cp, sp = math.cos(phi), math.sin(phi)
    x0 = a * np.cos(t)
    y0 = b * np.sin(t)
    x = float(fit["center_x_native"]) + cp * x0 - sp * y0
    y = float(fit["center_y_native"]) + sp * x0 + cp * y0
    z = np.full_like(x, float(z_native))
    return np.column_stack([x, y, z])


def horizontal_z_positions(zmin_native: float, zmax_native: float, scale_to_mm: float,
                           z_sections: int | None, z_step_mm: float | None):
    height_native = zmax_native - zmin_native
    height_mm = height_native * scale_to_mm
    if height_native <= 0:
        raise ValueError("Mesh Z extent must be positive.")
    if z_step_mm is not None:
        if z_step_mm <= 0:
            raise ValueError("--z-step-mm must be positive.")
        step_native = z_step_mm / scale_to_mm
        vals = []
        z = zmin_native + 0.5 * step_native
        upper = zmax_native - 0.5 * step_native
        while z <= upper + 1e-12:
            vals.append(float(z))
            z += step_native
        if len(vals) < 3:
            raise ValueError("--z-step-mm produced fewer than 3 horizontal sections.")
        return vals, {"mode": "absolute_step", "z_step_mm": float(z_step_mm), "count": len(vals)}
    n = 20 if z_sections is None else int(z_sections)
    if n < 3:
        raise ValueError("--z-sections must be at least 3.")
    vals = [zmin_native + (i + 0.5) * height_native / n for i in range(n)]
    return vals, {"mode": "relative_equal_divisions", "z_sections": n, "count": n, "height_mm": height_mm}


def choose_axis_center(records: list[dict], axis_surface: str, outlier_mad_k: float):
    if axis_surface not in {"inner", "outer", "both"}:
        raise ValueError("axis_surface must be inner, outer, or both")
    source = "inner" if axis_surface in {"inner", "both"} else "outer"
    valid = []
    for r in records:
        fit = r.get(f"{source}_fit")
        if fit and fit.get("valid"):
            valid.append(r)
    if len(valid) < 3:
        raise RuntimeError(f"Fewer than 3 valid {source} ellipse centers; cannot estimate rotation axis.")

    centers = np.array([[r[f"{source}_fit"]["center_x_mm"], r[f"{source}_fit"]["center_y_mm"]] for r in valid], dtype=float)
    med_center = np.median(centers, axis=0)
    d = np.linalg.norm(centers - med_center[None, :], axis=1)
    med_d = float(np.median(d))
    mad_d = float(np.median(np.abs(d - med_d)))
    robust_sigma = 1.4826 * mad_d
    threshold = med_d + float(outlier_mad_k) * robust_sigma
    if robust_sigma <= 1e-12:
        threshold = max(med_d, 1e-9)
    inlier = d <= threshold + 1e-12
    if int(np.count_nonzero(inlier)) < 3:
        # Safety fallback: retain all valid sections rather than producing an unstable 2-point axis.
        inlier[:] = True
        threshold = float("inf")

    final_center = centers[inlier].mean(axis=0)
    for r in records:
        fit = r.get(f"{source}_fit")
        if not fit or not fit.get("valid"):
            r["center_outlier"] = False
            r["used_for_axis"] = False
            r["exclude_reason"] = f"no_valid_{source}_ellipse"
            continue
        idx = valid.index(r)
        r["center_outlier"] = bool(not inlier[idx])
        r["used_for_axis"] = bool(inlier[idx])
        r["exclude_reason"] = "center_outlier" if not inlier[idx] else ""
        r["selected_center_x_mm"] = float(fit["center_x_mm"])
        r["selected_center_y_mm"] = float(fit["center_y_mm"])
        r["delta_x_mm"] = float(fit["center_x_mm"] - final_center[0])
        r["delta_y_mm"] = float(fit["center_y_mm"] - final_center[1])
        r["radial_offset_mm"] = float(math.hypot(r["delta_x_mm"], r["delta_y_mm"]))

    used_records = [r for r in records if r.get("used_for_axis")]
    used_centers = np.array([[r["selected_center_x_mm"], r["selected_center_y_mm"]] for r in used_records], dtype=float)
    z_used = np.array([r["z_mm"] for r in used_records], dtype=float)
    dx = used_centers[:, 0] - final_center[0]
    dy = used_centers[:, 1] - final_center[1]
    radial = np.hypot(dx, dy)

    def mad1(a):
        a = np.asarray(a, dtype=float)
        med = np.median(a)
        return float(np.median(np.abs(a - med)))

    if len(used_records) >= 2 and np.ptp(z_used) > 0:
        slope_x, intercept_x = np.polyfit(z_used, used_centers[:, 0], 1)
        slope_y, intercept_y = np.polyfit(z_used, used_centers[:, 1], 1)
        tilt_deg = math.degrees(math.atan(math.hypot(slope_x, slope_y)))
    else:
        slope_x = slope_y = intercept_x = intercept_y = tilt_deg = float("nan")

    summary = {
        "axis_surface_mode": axis_surface,
        "axis_center_source": source,
        "attempted_sections": int(len(records)),
        "valid_source_sections": int(len(valid)),
        "used_for_axis": int(len(used_records)),
        "excluded_center_outliers": int(len(valid) - len(used_records)),
        "axis_center_x_mm": float(final_center[0]),
        "axis_center_y_mm": float(final_center[1]),
        "center_outlier_mad_k": float(outlier_mad_k),
        "center_outlier_threshold_mm": None if not np.isfinite(threshold) else float(threshold),
        "sd_x_mm": float(np.std(used_centers[:, 0], ddof=1)) if len(used_centers) > 1 else 0.0,
        "sd_y_mm": float(np.std(used_centers[:, 1], ddof=1)) if len(used_centers) > 1 else 0.0,
        "mad_x_mm": mad1(used_centers[:, 0]),
        "mad_y_mm": mad1(used_centers[:, 1]),
        "mean_radial_offset_mm": float(np.mean(radial)),
        "median_radial_offset_mm": float(np.median(radial)),
        "rms_radial_offset_mm": float(np.sqrt(np.mean(radial ** 2))),
        "p95_radial_offset_mm": float(np.percentile(radial, 95)),
        "max_radial_offset_mm": float(np.max(radial)),
        "axis_drift_x_mm_per_100mm": float(slope_x * 100.0),
        "axis_drift_y_mm_per_100mm": float(slope_y * 100.0),
        "estimated_centerline_tilt_deg": float(tilt_deg),
        "tilt_is_diagnostic_only": True,
        "tilt_is_not_auto_corrected": True,
    }
    return (float(final_center[0]), float(final_center[1])), summary


def angle_series(step_deg: float, start_deg: float, span_deg: float) -> list[float]:
    if step_deg <= 0:
        raise ValueError("angle step must be positive.")
    n_float = span_deg / step_deg
    n = int(round(n_float))
    if not math.isclose(n_float, n, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError(f"angle step must divide {span_deg:g} degrees exactly.")
    return [float((start_deg + i * step_deg) % 360.0) for i in range(n)]


def plane_basis(angle_deg: float) -> tuple[np.ndarray, np.ndarray]:
    a = math.radians(angle_deg)
    radial = np.array([math.cos(a), math.sin(a), 0.0], dtype=np.float64)
    normal = np.array([-math.sin(a), math.cos(a), 0.0], dtype=np.float64)
    return radial, normal


def section_segments(mesh: trimesh.Trimesh, center_xy_native: tuple[float, float], angle_deg: float) -> np.ndarray:
    _radial, normal = plane_basis(angle_deg)
    zmid = float(mesh.bounds[:, 2].mean())
    origin = np.array([center_xy_native[0], center_xy_native[1], zmid], dtype=np.float64)
    segments = trimesh.intersections.mesh_plane(mesh=mesh, plane_normal=normal, plane_origin=origin, return_faces=False)
    segments = np.asarray(segments, dtype=np.float64)
    if segments.size == 0:
        return np.empty((0, 2, 3), dtype=np.float64)
    return segments.reshape((-1, 2, 3))


def signed_r(points: np.ndarray, center_xy_native: tuple[float, float], angle_deg: float) -> np.ndarray:
    radial, _ = plane_basis(angle_deg)
    d = np.asarray(points, dtype=np.float64).copy()
    d[..., 0] -= center_xy_native[0]
    d[..., 1] -= center_xy_native[1]
    return np.tensordot(d, radial, axes=([-1], [0]))


def clip_segments_to_positive_half(segments: np.ndarray, center_xy_native: tuple[float, float], ray_angle_deg: float,
                                   eps: float = 1e-12) -> np.ndarray:
    if len(segments) == 0:
        return segments.copy()
    s = signed_r(segments, center_xy_native, ray_angle_deg)
    out = []
    for seg, rr in zip(segments, s):
        p0, p1 = seg
        r0, r1 = float(rr[0]), float(rr[1])
        in0, in1 = r0 >= -eps, r1 >= -eps
        if in0 and in1:
            out.append(np.array([p0, p1], dtype=np.float64)); continue
        if (not in0) and (not in1):
            continue
        denom = r0 - r1
        if abs(denom) < eps:
            continue
        t = min(1.0, max(0.0, r0 / denom))
        cross = p0 + t * (p1 - p0)
        out.append(np.array([p0, cross] if in0 else [cross, p1], dtype=np.float64))
    return np.stack(out, axis=0) if out else np.empty((0, 2, 3), dtype=np.float64)


def sample_segments(segments: np.ndarray, spacing_native: float) -> np.ndarray:
    if len(segments) == 0:
        return np.empty((0, 3), dtype=np.float64)
    pts = []
    for p0, p1 in segments:
        length = float(np.linalg.norm(p1 - p0))
        n = max(2, int(math.ceil(length / spacing_native)) + 1)
        t = np.linspace(0.0, 1.0, n)[:, None]
        pts.append(p0[None, :] * (1.0 - t) + p1[None, :] * t)
    return np.vstack(pts)


def hsv_color(index: int, total: int) -> np.ndarray:
    import colorsys
    hue = (index / max(1, total)) % 1.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.95)
    return np.array([round(255*r), round(255*g), round(255*b), 255], dtype=np.uint8)


def write_points_ply(path: Path, points: np.ndarray, color: np.ndarray | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pts = np.asarray(points, dtype=np.float64)
    if color is None:
        pc = trimesh.points.PointCloud(pts)
    else:
        c = np.asarray(color, dtype=np.uint8).reshape(1, 4)
        pc = trimesh.points.PointCloud(pts, colors=np.repeat(c, len(pts), axis=0))
    pc.export(str(path))


def write_combined_points(path: Path, groups: list[tuple[np.ndarray, np.ndarray]]) -> None:
    points, colors = [], []
    for pts, color in groups:
        if len(pts) == 0:
            continue
        points.append(np.asarray(pts, dtype=np.float64))
        colors.append(np.repeat(np.asarray(color, dtype=np.uint8).reshape(1,4), len(pts), axis=0))
    pc = trimesh.points.PointCloud(np.vstack(points), colors=np.vstack(colors)) if points else trimesh.points.PointCloud(np.empty((0,3)))
    path.parent.mkdir(parents=True, exist_ok=True)
    pc.export(str(path))


def write_edge_ply(path: Path, segments: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    points = segments.reshape((-1, 3)) if len(segments) else np.empty((0,3), dtype=float)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write(f"element vertex {len(points)}\n")
        f.write("property double x\nproperty double y\nproperty double z\n")
        f.write(f"element edge {len(segments)}\n")
        f.write("property int vertex1\nproperty int vertex2\nend_header\n")
        for p in points: f.write(f"{p[0]:.12g} {p[1]:.12g} {p[2]:.12g}\n")
        for i in range(len(segments)): f.write(f"{2*i} {2*i+1}\n")


def write_axis_edge_ply(path: Path, center_native: tuple[float,float], zmin: float, zmax: float) -> None:
    seg = np.array([[[center_native[0], center_native[1], zmin], [center_native[0], center_native[1], zmax]]], dtype=float)
    write_edge_ply(path, seg)


def write_segment_csv(path: Path, segments: np.ndarray, center_xy_native: tuple[float,float], angle_deg: float,
                      scale_to_mm: float, z_min_native: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rr = signed_r(segments, center_xy_native, angle_deg) if len(segments) else np.empty((0,2))
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["segment_id","endpoint","angle_deg","side","x_input","y_input","z_input",
                    "x_mm","y_mm","z_mm","signed_r_mm","radial_r_mm","z_from_bottom_mm"])
        for i, seg in enumerate(segments):
            for j in (0,1):
                r_native=float(rr[i,j]); p=seg[j]
                w.writerow([i,j,f"{angle_deg:.10g}","positive" if r_native>=0 else "negative",
                            f"{p[0]:.12g}",f"{p[1]:.12g}",f"{p[2]:.12g}",
                            f"{p[0]*scale_to_mm:.12g}",f"{p[1]*scale_to_mm:.12g}",f"{p[2]*scale_to_mm:.12g}",
                            f"{r_native*scale_to_mm:.12g}",f"{abs(r_native)*scale_to_mm:.12g}",
                            f"{(p[2]-z_min_native)*scale_to_mm:.12g}"])



# ----------------------------------------------------------------------
# Capacity calculation helpers (v0.4)
# ----------------------------------------------------------------------

def _dedupe_sorted(values: np.ndarray, tolerance: float) -> np.ndarray:
    v = np.sort(np.asarray(values, dtype=np.float64))
    if len(v) <= 1:
        return v
    out = [float(v[0])]
    for x in v[1:]:
        if float(x) - out[-1] > tolerance:
            out.append(float(x))
        else:
            # Multiple mesh-plane segments can meet at one geometric point.
            # Keep a stable representative rather than counting the vertex twice.
            out[-1] = 0.5 * (out[-1] + float(x))
    return np.asarray(out, dtype=np.float64)


def _fill_short_nan_gaps(values: np.ndarray, max_gap: int = 2) -> np.ndarray:
    a = np.asarray(values, dtype=np.float64).copy()
    good = np.isfinite(a)
    n = len(a)
    i = 0
    while i < n:
        if good[i]:
            i += 1
            continue
        j = i
        while j < n and not good[j]:
            j += 1
        gap = j - i
        if i > 0 and j < n and gap <= max_gap and good[i-1] and good[j]:
            a[i:j] = np.linspace(a[i-1], a[j], gap + 2)[1:-1]
            good[i:j] = True
        i = j
    return a


def _longest_finite_run(values: np.ndarray) -> tuple[int | None, int | None]:
    good = np.isfinite(np.asarray(values, dtype=np.float64))
    best = (None, None)
    best_len = 0
    i = 0
    while i < len(good):
        if not good[i]:
            i += 1
            continue
        j = i + 1
        while j < len(good) and good[j]:
            j += 1
        if j - i > best_len:
            best = (i, j)
            best_len = j - i
        i = j
    return best


def extract_ray_inner_profile(
    full_segments: np.ndarray,
    center_xy_native: tuple[float, float],
    ray_angle_deg: float,
    z_centers_native: np.ndarray,
    z_step_mm: float,
    scale_to_mm: float,
    dedupe_tolerance_mm: float = 0.05,
    max_gap_cells: int = 2,
) -> dict:
    """
    Sample a longitudinal mesh-plane intersection at regular Z positions.

    On a positive radial ray, a normal vessel-wall section has two crossings
    while the ray begins inside the cavity: the smaller positive radius is the
    inner surface and the larger is the outer surface. Below the inner bottom,
    normally only the outer crossing exists, so no inner radius is assigned.

    The longest continuous run of valid inner radii is retained. This rejects
    isolated intersections caused by local geometry or numerical degeneracy.
    """
    zc = np.asarray(z_centers_native, dtype=np.float64)
    inner = np.full(len(zc), np.nan, dtype=np.float64)
    outer = np.full(len(zc), np.nan, dtype=np.float64)
    counts = np.zeros(len(zc), dtype=np.int32)
    if len(full_segments) == 0:
        return {
            "angle_deg": float(ray_angle_deg), "z_native": zc,
            "z_mm": zc * scale_to_mm, "inner_radius_mm": inner,
            "outer_radius_mm": outer, "positive_intersections": counts,
            "run_start_index": None, "run_end_index": None,
            "bottom_z_mm": None, "rim_z_mm": None, "valid_cells": 0,
        }

    seg = np.asarray(full_segments, dtype=np.float64)
    z0, z1 = seg[:, 0, 2], seg[:, 1, 2]
    rr = signed_r(seg, center_xy_native, ray_angle_deg)
    r0, r1 = rr[:, 0], rr[:, 1]
    dz = z1 - z0
    eps_z = max(1e-12, (z_step_mm / scale_to_mm) * 1e-8)
    tol_native = dedupe_tolerance_mm / scale_to_mm

    for k, z in enumerate(zc):
        nonhorizontal = np.abs(dz) > eps_z
        lo = np.minimum(z0, z1) - eps_z
        hi = np.maximum(z0, z1) + eps_z
        mask = nonhorizontal & (z >= lo) & (z <= hi)
        if not np.any(mask):
            continue
        t = (z - z0[mask]) / dz[mask]
        ok = (t >= -1e-9) & (t <= 1.0 + 1e-9)
        if not np.any(ok):
            continue
        rad = r0[mask][ok] + t[ok] * (r1[mask][ok] - r0[mask][ok])
        rad = rad[rad >= -tol_native]
        if len(rad) == 0:
            continue
        rad = np.maximum(rad, 0.0)
        rad = _dedupe_sorted(rad, tol_native)
        counts[k] = len(rad)
        # Two or more positive crossings indicate cavity -> inner wall -> body -> outer wall.
        if len(rad) >= 2:
            inner[k] = float(rad[0] * scale_to_mm)
            outer[k] = float(rad[-1] * scale_to_mm)

    inner = _fill_short_nan_gaps(inner, max_gap=max_gap_cells)
    outer = _fill_short_nan_gaps(outer, max_gap=max_gap_cells)
    i0, i1 = _longest_finite_run(inner)
    if i0 is None:
        inner[:] = np.nan
        outer[:] = np.nan
        bottom = rim = None
        valid_cells = 0
    else:
        keep = np.zeros(len(inner), dtype=bool)
        keep[i0:i1] = True
        inner[~keep] = np.nan
        outer[~keep] = np.nan
        z_mm = zc * scale_to_mm
        bottom = float(z_mm[i0] - 0.5 * z_step_mm)
        rim = float(z_mm[i1 - 1] + 0.5 * z_step_mm)
        valid_cells = int(i1 - i0)

    return {
        "angle_deg": float(ray_angle_deg),
        "z_native": zc,
        "z_mm": zc * scale_to_mm,
        "inner_radius_mm": inner,
        "outer_radius_mm": outer,
        "positive_intersections": counts,
        "run_start_index": None if i0 is None else int(i0),
        "run_end_index": None if i1 is None else int(i1),
        "bottom_z_mm": bottom,
        "rim_z_mm": rim,
        "valid_cells": valid_cells,
    }


def write_dict_rows_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8-sig")
        return
    fields = []
    seen = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key); fields.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: (json.dumps(v, ensure_ascii=False) if isinstance(v, (dict,list,tuple)) else v) for k,v in row.items()})


def export_inner_profile_points(path: Path, ray_profiles: list[dict], center_native: tuple[float,float], scale_to_mm: float) -> None:
    groups = []
    for i, rp in enumerate(ray_profiles):
        rmm = np.asarray(rp["inner_radius_mm"], dtype=float)
        z_native = np.asarray(rp["z_native"], dtype=float)
        valid = np.isfinite(rmm)
        if not np.any(valid):
            continue
        a = math.radians(float(rp["angle_deg"]))
        rn = rmm[valid] / scale_to_mm
        pts = np.column_stack([
            center_native[0] + rn * math.cos(a),
            center_native[1] + rn * math.sin(a),
            z_native[valid],
        ])
        groups.append((pts, hsv_color(i, len(ray_profiles))))
    write_combined_points(path, groups)


def calculate_capacity_modes(
    mesh: trimesh.Trimesh,
    center_native: tuple[float,float],
    scale_to_mm: float,
    ray_angles: list[float],
    full_cache: dict[int, np.ndarray],
    full_angles: list[float],
    horizontal_records: list[dict],
    requested_modes: list[str],
    single_angle_deg: float,
    volume_z_step_mm: float,
    min_angular_valid_fraction: float,
    volume_dir: Path,
    make_visualization: bool,
) -> dict:
    """Acquire 3-D profiles, then delegate numerical capacity calculations to PotteryVolumeCore."""
    if volume_z_step_mm <= 0:
        raise ValueError("--volume-z-step-mm must be positive.")
    if not (0.0 < min_angular_valid_fraction <= 1.0):
        raise ValueError("--min-angular-valid-fraction must be in (0,1].")
    modes = set(requested_modes)
    if "none" in modes:
        return {
            "requested_modes": ["none"], "results": [], "implemented": True,
            "volume_core_version": volume_core.__version__,
        }
    if "all" in modes:
        modes = {"single", "optimized", "angular", "ellipse"}
    allowed = {"single", "optimized", "angular", "ellipse"}
    unknown = modes - allowed
    if unknown:
        raise ValueError(f"Unknown volume mode(s): {sorted(unknown)}")

    volume_dir.mkdir(parents=True, exist_ok=True)
    zmin_native, zmax_native = map(float, mesh.bounds[:,2])
    step_native = volume_z_step_mm / scale_to_mm
    z_centers_native = np.arange(zmin_native + 0.5*step_native, zmax_native, step_native, dtype=float)
    if len(z_centers_native) < 3:
        raise ValueError("Volume Z step produced fewer than 3 cells.")
    z_mm = z_centers_native * scale_to_mm

    # Geometry acquisition remains a responsibility of the radial front end.
    ray_profiles = []
    for j, angle in enumerate(ray_angles):
        base_index = j % len(full_angles)
        rp = extract_ray_inner_profile(
            full_cache[base_index], center_native, angle, z_centers_native,
            volume_z_step_mm, scale_to_mm,
        )
        ray_profiles.append(rp)

    spans = np.array([
        (rp["rim_z_mm"] - rp["bottom_z_mm"]) if rp["rim_z_mm"] is not None and rp["bottom_z_mm"] is not None else np.nan
        for rp in ray_profiles
    ], dtype=float)
    finite_spans = spans[np.isfinite(spans)]
    median_span = float(np.median(finite_spans)) if len(finite_spans) else float("nan")
    for rp, span in zip(ray_profiles, spans):
        rp["profile_span_mm"] = None if not np.isfinite(span) else float(span)
        rp["profile_reliable"] = bool(np.isfinite(span) and (not np.isfinite(median_span) or span >= 0.5 * median_span))

    profile_rows = []
    for rp in ray_profiles:
        for k in range(len(z_mm)):
            profile_rows.append({
                "angle_deg": rp["angle_deg"],
                "z_mm": float(z_mm[k]),
                "z_from_mesh_bottom_mm": float((z_centers_native[k]-zmin_native)*scale_to_mm),
                "inner_radius_mm": None if not np.isfinite(rp["inner_radius_mm"][k]) else float(rp["inner_radius_mm"][k]),
                "outer_radius_mm": None if not np.isfinite(rp["outer_radius_mm"][k]) else float(rp["outer_radius_mm"][k]),
                "positive_intersections": int(rp["positive_intersections"][k]),
                "profile_reliable": rp["profile_reliable"],
            })
    write_dict_rows_csv(volume_dir/"radial_inner_profiles.csv", profile_rows)
    export_inner_profile_points(volume_dir/"radial_inner_profile_points.ply", ray_profiles, center_native, scale_to_mm)

    reliable = [rp for rp in ray_profiles if rp.get("profile_reliable") and rp.get("rim_z_mm") is not None]
    physical_spill = min((rp["rim_z_mm"] for rp in reliable), default=None)
    physical_bottom = max((rp["bottom_z_mm"] for rp in reliable), default=None)
    min_valid_rays = max(3, int(math.ceil(len(ray_profiles) * min_angular_valid_fraction)))
    results = []

    def failed(method: str, note: str, **extra):
        row={"method":method,"status":"failed","core_version":volume_core.__version__,
             "volume_mm3":None,"volume_ml":None,"volume_l":None,"note":note}
        row.update(extra)
        return row

    if "single" in modes:
        angle = float(single_angle_deg % 360.0)
        opposite = float((angle + 180.0) % 360.0)
        matched = {}
        for target in (angle, opposite):
            for rp in ray_profiles:
                if abs(((rp["angle_deg"] - target + 180.0) % 360.0) - 180.0) < 1e-9:
                    matched[target] = rp
                    break
        if len(matched) < 2:
            seg = section_segments(mesh, center_native, angle)
            matched[angle] = extract_ray_inner_profile(seg, center_native, angle, z_centers_native, volume_z_step_mm, scale_to_mm)
            matched[opposite] = extract_ray_inner_profile(seg, center_native, opposite, z_centers_native, volume_z_step_mm, scale_to_mm)
        rp1, rp2 = matched[angle], matched[opposite]
        r1 = np.asarray(rp1["inner_radius_mm"], float)
        r2 = np.asarray(rp2["inner_radius_mm"], float)
        bottoms = [x for x in (rp1.get("bottom_z_mm"), rp2.get("bottom_z_mm")) if x is not None]
        rims = [x for x in (rp1.get("rim_z_mm"), rp2.get("rim_z_mm")) if x is not None]
        lo = max(bottoms) if len(bottoms) == 2 else None
        hi = min(rims) if len(rims) == 2 else None
        good1=np.isfinite(r1); good2=np.isfinite(r2)
        if lo is None or hi is None or np.count_nonzero(good1)<2 or np.count_nonzero(good2)<2:
            results.append(failed("single","Both opposing inner profiles are required for the selected longitudinal plane.",angle_deg=angle,opposite_angle_deg=opposite))
        else:
            try:
                calc=volume_core.calculate_single_two_sided(
                    z_mm[good1], r1[good1], z_mm[good2], r2[good2],
                    z_step_mm=volume_z_step_mm, bottom_z_mm=lo, spill_z_mm=hi, method="single",
                )
                summary=dict(calc["summary"])
                summary.update({
                    "angle_deg": angle, "opposite_angle_deg": opposite,
                    "side_a_volume_l": summary.get("left_only_volume_l"),
                    "side_b_volume_l": summary.get("right_only_volume_l"),
                    "note": "One full longitudinal plane. Numerical integration is performed by pottery_volume_core.calculate_single_two_sided.",
                })
                results.append(summary)
                pr=calc["profile"]
                rows=[]
                for z,a,b,r,aa in zip(pr["z_mm"],pr["left_radius_mm"],pr["right_radius_mm"],pr["equivalent_radius_mm"],pr["area_mm2"]):
                    rows.append({"z_mm":float(z),"radius_side_a_mm":float(a),"radius_side_b_mm":float(b),"equivalent_radius_mm":float(r),"area_mm2":float(aa)})
                stem=f"single_section_{int(round(angle))%360:03d}_{int(round(opposite))%360:03d}_profile.csv"
                write_dict_rows_csv(volume_dir/stem,rows)
            except Exception as exc:
                results.append(failed("single",f"PotteryVolumeCore failed: {exc}",angle_deg=angle,opposite_angle_deg=opposite))

    radii_matrix = np.vstack([np.asarray(rp["inner_radius_mm"], float) for rp in ray_profiles])

    if "optimized" in modes:
        calc=volume_core.calculate_optimized_axisymmetric(
            z_mm, radii_matrix, z_step_mm=volume_z_step_mm,
            bottom_z_mm=physical_bottom, spill_z_mm=physical_spill,
            min_valid_rays=min_valid_rays, mad_k=3.0, method="optimized",
        )
        results.append(calc["summary"])
        pr=calc.get("profile",{})
        if pr:
            rows=[]
            for z,rr,aa,nv,ni in zip(pr["z_mm"],pr["optimized_radius_mm"],pr["area_mm2"],pr["valid_rays"],pr["radius_inliers"]):
                rows.append({"z_mm":float(z),"optimized_radius_mm":None if not np.isfinite(rr) else float(rr),"area_mm2":None if not np.isfinite(aa) else float(aa),"valid_rays":int(nv),"radius_inliers":int(ni)})
            write_dict_rows_csv(volume_dir/"optimized_rotation_profile.csv",rows)

    if "angular" in modes:
        angles=np.asarray([rp["angle_deg"] for rp in ray_profiles],float)
        calc=volume_core.calculate_angular_integration(
            z_mm, angles, radii_matrix, z_step_mm=volume_z_step_mm,
            bottom_z_mm=physical_bottom, spill_z_mm=physical_spill,
            min_valid_rays=min_valid_rays, method="angular",
        )
        results.append(calc["summary"])
        pr=calc.get("profile",{})
        if pr:
            rows=[]
            for z,aa,nv,ni in zip(pr["z_mm"],pr["area_mm2"],pr["valid_rays"],pr["interpolated_rays"]):
                rows.append({"z_mm":float(z),"area_mm2":None if not np.isfinite(aa) else float(aa),"valid_rays":int(nv),"interpolated_rays":int(ni)})
            write_dict_rows_csv(volume_dir/"angular_area_profile.csv",rows)

    if "ellipse" in modes:
        ez=[]; ea=[]; eb=[]; eids=[]
        for rec in horizontal_records:
            fit=rec.get("inner_fit")
            if fit and fit.get("valid"):
                ez.append(float(rec["z_mm"])); ea.append(float(fit["semi_major_mm"])); eb.append(float(fit["semi_minor_mm"])); eids.append(rec.get("section_id"))
        if len(ez)<3:
            results.append(failed("ellipse","Fewer than 3 valid inner horizontal ellipses."))
        else:
            calc=volume_core.calculate_ellipse_integration(
                ez,ea,eb,z_grid_mm=z_mm,z_step_mm=volume_z_step_mm,
                bottom_z_mm=physical_bottom,spill_z_mm=physical_spill,
                section_ids=eids,method="ellipse",
            )
            results.append(calc["summary"])
            pr=calc.get("profile",{})
            if pr:
                ids=pr.get("section_id")
                rows=[]
                for i,(z,a) in enumerate(zip(pr["section_z_mm"],pr["ellipse_area_mm2"])):
                    sid=None if ids is None else ids[i]
                    rows.append({"section_id":sid,"z_mm":float(z),"ellipse_area_mm2":float(a)})
                write_dict_rows_csv(volume_dir/"horizontal_ellipse_area_samples.csv",rows)

    write_dict_rows_csv(volume_dir/"volume_summary.csv",results)
    write_json(volume_dir/"volume_summary.json",{
        "requested_modes": sorted(modes),
        "volume_core_version": volume_core.__version__,
        "volume_z_step_mm": float(volume_z_step_mm),
        "radial_angles_deg": [float(x) for x in ray_angles],
        "min_angular_valid_fraction": float(min_angular_valid_fraction),
        "min_valid_rays": int(min_valid_rays),
        "physical_spill_z_mm_estimate": physical_spill,
        "physical_common_bottom_z_mm_estimate": physical_bottom,
        "spill_definition": "lowest reliable radial inner-profile rim height; resolution approximately +/- half volume_z_step_mm",
        "results": results,
    })

    if make_visualization and results:
        try:
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            ok=[r for r in results if r.get("status")=="ok" and r.get("volume_l") is not None]
            if ok:
                fig,ax=plt.subplots(figsize=(7.2,4.8))
                ax.bar([r["method"] for r in ok],[r["volume_l"] for r in ok])
                ax.set_ylabel("Volume [L]"); ax.set_title("Capacity estimates by method")
                fig.tight_layout(); fig.savefig(volume_dir/"volume_method_comparison.png",dpi=220); plt.close(fig)
            fig,ax=plt.subplots(figsize=(7.0,7.4))
            for rp in ray_profiles:
                rr=np.asarray(rp["inner_radius_mm"],float); good=np.isfinite(rr)
                if np.any(good): ax.plot(rr[good],z_mm[good],linewidth=.7,alpha=.5)
            ax.set_xlabel("Inner radius [mm]"); ax.set_ylabel("Z [mm]"); ax.set_title("Inner radial profiles used for capacity calculation")
            fig.tight_layout(); fig.savefig(volume_dir/"radial_inner_profiles.png",dpi=220); plt.close(fig)
        except Exception as exc:
            write_json(volume_dir/"visualization_warning.json",{"warning":str(exc)})

    return {
        "requested_modes": sorted(modes), "results": results,
        "volume_core_version": volume_core.__version__,
        "volume_z_step_mm": float(volume_z_step_mm),
        "physical_spill_z_mm_estimate": physical_spill,
        "physical_common_bottom_z_mm_estimate": physical_bottom,
        "radial_profile_count": len(ray_profiles),
    }

def make_oblique_plot(mesh, groups, path, scale_to_mm, title, max_mesh_points=50000, max_group_points=9000):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    verts=np.asarray(mesh.vertices,float)
    if len(verts)>max_mesh_points:
        verts=verts[np.linspace(0,len(verts)-1,max_mesh_points,dtype=int)]
    verts=verts*scale_to_mm
    fig=plt.figure(figsize=(8.2,8.2)); ax=fig.add_subplot(111,projection="3d")
    ax.scatter(verts[:,0],verts[:,1],verts[:,2],s=.12,alpha=.07)
    for pts,color in groups:
        p=np.asarray(pts,float)
        if len(p)==0: continue
        if len(p)>max_group_points:
            p=p[np.linspace(0,len(p)-1,max_group_points,dtype=int)]
        p=p*scale_to_mm
        ax.scatter(p[:,0],p[:,1],p[:,2],s=1.0,alpha=.8)
    b=mesh.bounds*scale_to_mm; mid=b.mean(axis=0); span=(b[1]-b[0]).max(); half=span*.58
    ax.set_xlim(mid[0]-half,mid[0]+half); ax.set_ylim(mid[1]-half,mid[1]+half); ax.set_zlim(mid[2]-half,mid[2]+half)
    ax.set_xlabel("X [mm]"); ax.set_ylabel("Y [mm]"); ax.set_zlabel("Z [mm]"); ax.set_title(title); ax.view_init(elev=24,azim=-55)
    try: ax.set_box_aspect((1,1,1))
    except Exception: pass
    fig.tight_layout(); path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(path,dpi=220); plt.close(fig)


def make_axis_validation_plots(mesh, center_mm, records, full0_segments, full90_segments, vis_dir, scale):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    source_records=[r for r in records if r.get("selected_center_x_mm") is not None]

    def panel(segments, coord_idx, center_value, coord_label, filename, title):
        fig,ax=plt.subplots(figsize=(7.2,8.2))
        if len(segments):
            for seg in segments:
                ax.plot(seg[:,coord_idx]*scale, seg[:,2]*scale, linewidth=.45, alpha=.55)
        used=[r for r in source_records if r.get("used_for_axis")]
        out=[r for r in source_records if r.get("center_outlier")]
        if used:
            ax.scatter([r[f"selected_center_{coord_label.lower()}_mm"] for r in used], [r["z_mm"] for r in used], marker="o", s=24, label="used ellipse centers")
        if out:
            ax.scatter([r[f"selected_center_{coord_label.lower()}_mm"] for r in out], [r["z_mm"] for r in out], marker="x", s=40, label="center outliers")
        ax.axvline(center_value, linestyle="--", linewidth=1.2, label="final rotation axis")
        ax.set_xlabel(f"{coord_label} [mm]"); ax.set_ylabel("Z [mm]"); ax.set_title(title); ax.set_aspect("equal",adjustable="box"); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(vis_dir/filename,dpi=220); plt.close(fig)

    panel(full0_segments,0,center_mm[0],"X","axis_validation_xz.png","0° longitudinal section and horizontal-section ellipse centers")
    panel(full90_segments,1,center_mm[1],"Y","axis_validation_yz.png","90° longitudinal section and horizontal-section ellipse centers")

    fig,ax=plt.subplots(figsize=(7.2,7.2))
    used=[r for r in source_records if r.get("used_for_axis")]
    out=[r for r in source_records if r.get("center_outlier")]
    if used:
        ax.scatter([r["selected_center_x_mm"] for r in used],[r["selected_center_y_mm"] for r in used],marker="o",s=26,label="used")
    if out:
        ax.scatter([r["selected_center_x_mm"] for r in out],[r["selected_center_y_mm"] for r in out],marker="x",s=45,label="outlier")
    ax.scatter([center_mm[0]],[center_mm[1]],marker="+",s=120,label="final mean center")
    ax.set_xlabel("X [mm]"); ax.set_ylabel("Y [mm]"); ax.set_title("Horizontal-section ellipse centers in XY"); ax.set_aspect("equal",adjustable="box"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(vis_dir/"axis_centers_xy.png",dpi=220); plt.close(fig)


def make_table_images(records, summary, vis_dir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols=["id","z %","X mm","Y mm","ΔX","ΔY","Δr","RMSE","use"]
    cells=[]
    for r in records:
        cells.append([
            str(r["section_id"]), f"{100*r['z_normalized']:.1f}",
            "" if r.get("selected_center_x_mm") is None else f"{r['selected_center_x_mm']:.3f}",
            "" if r.get("selected_center_y_mm") is None else f"{r['selected_center_y_mm']:.3f}",
            "" if r.get("delta_x_mm") is None else f"{r['delta_x_mm']:+.3f}",
            "" if r.get("delta_y_mm") is None else f"{r['delta_y_mm']:+.3f}",
            "" if r.get("radial_offset_mm") is None else f"{r['radial_offset_mm']:.3f}",
            "" if r.get("selected_fit_rmse_mm") is None else f"{r['selected_fit_rmse_mm']:.3f}",
            "yes" if r.get("used_for_axis") else "no",
        ])
    fig_h=max(5.5,0.34*len(cells)+1.4)
    fig,ax=plt.subplots(figsize=(10.5,fig_h)); ax.axis("off")
    tbl=ax.table(cellText=cells,colLabels=cols,loc="center",cellLoc="center"); tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1,1.15)
    ax.set_title("Horizontal-section axis estimation table",pad=12); fig.tight_layout(); fig.savefig(vis_dir/"horizontal_sections_table.png",dpi=200,bbox_inches="tight"); plt.close(fig)

    keys=[
        ("axis_surface_mode","mode"),("axis_center_source","center source"),("attempted_sections","attempted"),("valid_source_sections","valid"),("used_for_axis","used"),
        ("axis_center_x_mm","axis X [mm]"),("axis_center_y_mm","axis Y [mm]"),("sd_x_mm","SD X [mm]"),("sd_y_mm","SD Y [mm]"),
        ("rms_radial_offset_mm","RMS Δr [mm]"),("max_radial_offset_mm","max Δr [mm]"),("axis_drift_x_mm_per_100mm","drift X /100mm"),
        ("axis_drift_y_mm_per_100mm","drift Y /100mm"),("estimated_centerline_tilt_deg","diagnostic tilt [deg]"),
    ]
    vals=[]
    for k,label in keys:
        v=summary.get(k,"")
        vals.append([label, f"{v:.6g}" if isinstance(v,float) else str(v)])
    fig,ax=plt.subplots(figsize=(7.2,5.8)); ax.axis("off"); tbl=ax.table(cellText=vals,colLabels=["statistic","value"],loc="center",cellLoc="left")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1,1.25); ax.set_title("Rotation-axis summary",pad=12); fig.tight_layout(); fig.savefig(vis_dir/"axis_summary_table.png",dpi=200,bbox_inches="tight"); plt.close(fig)


def flatten_fit(prefix: str, fit: dict | None, scale_to_mm: float) -> dict:
    out={}
    fields=["valid","reason","center_x_mm","center_y_mm","semi_major_mm","semi_minor_mm","ellipse_angle_deg","eccentricity","fit_rmse_mm","fit_rmse_relative","fit_points","fit_inlier_fraction"]
    for k in fields:
        out[f"{prefix}_{k}"] = "" if not fit or k not in fit else fit[k]
    return out


def write_horizontal_csv(path: Path, records: list[dict]) -> None:
    fields=[
        "section_id","z_input","z_mm","z_from_bottom_mm","z_normalized","contour_count","outer_area_mm2","inner_area_mm2","inner_outer_area_ratio",
        "axis_surface_mode","selected_surface","selected_center_x_mm","selected_center_y_mm","delta_x_mm","delta_y_mm","radial_offset_mm",
        "selected_fit_rmse_mm","selected_fit_rmse_relative","center_outlier","used_for_axis","exclude_reason",
        "outer_valid","outer_reason","outer_center_x_mm","outer_center_y_mm","outer_semi_major_mm","outer_semi_minor_mm","outer_ellipse_angle_deg","outer_eccentricity","outer_fit_rmse_mm","outer_fit_rmse_relative","outer_fit_points","outer_fit_inlier_fraction",
        "inner_valid","inner_reason","inner_center_x_mm","inner_center_y_mm","inner_semi_major_mm","inner_semi_minor_mm","inner_ellipse_angle_deg","inner_eccentricity","inner_fit_rmse_mm","inner_fit_rmse_relative","inner_fit_points","inner_fit_inlier_fraction",
        "outer_inner_center_offset_mm",
    ]
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in records:
            row={k:r.get(k,"") for k in fields}; w.writerow(row)


def write_axis_summary_csv(path: Path, summary: dict) -> None:
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.writer(f); w.writerow(["statistic","value","unit_or_note"])
        units={
            "axis_center_x_mm":"mm","axis_center_y_mm":"mm","sd_x_mm":"mm","sd_y_mm":"mm","mad_x_mm":"mm","mad_y_mm":"mm",
            "mean_radial_offset_mm":"mm","median_radial_offset_mm":"mm","rms_radial_offset_mm":"mm","p95_radial_offset_mm":"mm","max_radial_offset_mm":"mm",
            "axis_drift_x_mm_per_100mm":"mm/100mm","axis_drift_y_mm_per_100mm":"mm/100mm","estimated_centerline_tilt_deg":"degree",
            "center_outlier_threshold_mm":"mm",
            "mean_outer_inner_center_offset_mm":"mm",
            "median_outer_inner_center_offset_mm":"mm",
            "mean_selected_fit_rmse_mm":"mm",
            "median_selected_fit_rmse_mm":"mm",
        }
        for k,v in summary.items():
            if isinstance(v, (dict, list, tuple)):
                v = json.dumps(v, ensure_ascii=False)
            w.writerow([k,v,units.get(k,"")])


def process(input_path: Path, unit: str, step_deg: float, start_deg: float, axis_surface: str,
            z_sections: int | None, z_step_mm: float | None, center_outlier_mad: float,
            contour_spacing_mm: float, sample_spacing_mm: float, output_dir: Path | None,
            make_visualization: bool, volume_modes: list[str], single_angle_deg: float,
            volume_z_step_mm: float, min_angular_valid_fraction: float,
            section_overlay: bool = True, rim_options: dict | None = None) -> dict:
    if unit not in UNIT_SCALE_TO_MM: raise ValueError(f"Unsupported unit: {unit}")
    if contour_spacing_mm<=0 or sample_spacing_mm<=0: raise ValueError("Sampling spacings must be positive.")
    if center_outlier_mad<=0: raise ValueError("--center-outlier-mad must be positive.")
    scale=UNIT_SCALE_TO_MM[unit]
    mesh=load_geometry(input_path)
    zmin,zmax=map(float,mesh.bounds[:,2]); height_native=zmax-zmin

    full_angles=angle_series(step_deg,start_deg,180.0); ray_angles=angle_series(step_deg,start_deg,360.0)
    if output_dir is None:
        step_label=f"{step_deg:g}".replace(".","p")
        output_dir=input_path.parent/f"{input_path.stem}_RadialSections_{step_label}deg"
    full_dir=output_dir/"full_sections"; ray_dir=output_dir/"radial_half_sections"; vis_dir=output_dir/"visualization"
    axis_dir=output_dir/"axis_estimation"; horiz_dir=axis_dir/"horizontal_sections"; volume_dir=output_dir/"volume"
    for p in (full_dir,ray_dir,vis_dir,axis_dir,horiz_dir,volume_dir): p.mkdir(parents=True,exist_ok=True)

    z_values,z_sampling=horizontal_z_positions(zmin,zmax,scale,z_sections,z_step_mm)
    closure_tol_native=0.05/scale
    contour_spacing_native=contour_spacing_mm/scale
    horizontal_records=[]; horizontal_groups=[]; ellipse_groups=[]; center_pts=[]; center_colors=[]

    print(f"=== PotteryRadialSections v{__version__} ===")
    print(f"input          : {input_path}")
    print(f"unit           : {unit}")
    print(f"vertices       : {len(mesh.vertices):,}")
    print(f"faces          : {len(mesh.faces):,}")
    print(f"extents        : X={mesh.extents[0]*scale:.3f}, Y={mesh.extents[1]*scale:.3f}, Z={mesh.extents[2]*scale:.3f} mm")
    print(f"axis surface   : {axis_surface}")
    print(f"horizontal Z   : {z_sampling}")
    if axis_surface=="outer":
        print("NOTE: wall thickness is not used for axis estimation. It can be supplied later during inner-profile/volume reconstruction.")

    for i,z in enumerate(z_values, start=1):
        seg=horizontal_section_segments(mesh,z)
        contours=closed_horizontal_contours(seg,closure_tol_native)
        outer,inner=select_outer_inner(contours)
        outer_fit=inner_fit=None
        if outer is not None:
            pxy=resample_closed_polyline(np.asarray(outer["points_xy"]),contour_spacing_native)
            outer_fit=robust_ellipse_fit(pxy,scale)
        if inner is not None:
            pxy=resample_closed_polyline(np.asarray(inner["points_xy"]),contour_spacing_native)
            inner_fit=robust_ellipse_fit(pxy,scale)

        outer_area=0.0 if outer is None else float(outer["area_native2"]*scale*scale)
        inner_area=0.0 if inner is None else float(inner["area_native2"]*scale*scale)
        ratio=(inner_area/outer_area) if outer_area>0 and inner_area>0 else 0.0
        rec={
            "section_id":i,"z_input":float(z),"z_mm":float(z*scale),"z_from_bottom_mm":float((z-zmin)*scale),
            "z_normalized":float((z-zmin)/height_native),"contour_count":len(contours),"outer_area_mm2":outer_area,
            "inner_area_mm2":inner_area,"inner_outer_area_ratio":ratio,"axis_surface_mode":axis_surface,
            "outer_fit":outer_fit,"inner_fit":inner_fit,
        }
        rec.update(flatten_fit("outer",outer_fit,scale)); rec.update(flatten_fit("inner",inner_fit,scale))
        if outer_fit and outer_fit.get("valid") and inner_fit and inner_fit.get("valid"):
            rec["outer_inner_center_offset_mm"]=float(math.hypot(outer_fit["center_x_mm"]-inner_fit["center_x_mm"],outer_fit["center_y_mm"]-inner_fit["center_y_mm"]))
        else: rec["outer_inner_center_offset_mm"]=""
        horizontal_records.append(rec)

        # export raw horizontal intersection points and fitted ellipse points in native XYZ
        raw_pts=sample_segments(seg, sample_spacing_mm/scale)
        color=hsv_color(i-1,len(z_values)); horizontal_groups.append((raw_pts,color))
        stem=f"horizontal_{i:03d}_z{(z-zmin)*scale:08.3f}mm".replace(".","p")
        write_points_ply(horiz_dir/f"{stem}_intersection_points.ply",raw_pts,color)
        fit_for_display = inner_fit if axis_surface in {"inner","both"} else outer_fit
        ep=ellipse_points(fit_for_display,z)
        if len(ep):
            ellipse_groups.append((ep,color)); write_points_ply(horiz_dir/f"{stem}_fitted_ellipse_points.ply",ep,color)
        print(f"horizontal {i:02d} z={(z-zmin)*scale:8.3f} mm : contours={len(contours)} outer={'ok' if outer_fit and outer_fit.get('valid') else '-'} inner={'ok' if inner_fit and inner_fit.get('valid') else '-'}")

    center_mm,axis_summary=choose_axis_center(horizontal_records,axis_surface,center_outlier_mad)
    center_native=(center_mm[0]/scale,center_mm[1]/scale)
    source="inner" if axis_surface in {"inner","both"} else "outer"
    for r in horizontal_records:
        fit=r.get(f"{source}_fit")
        r["selected_surface"]=source
        r["selected_fit_rmse_mm"]=(fit.get("fit_rmse_mm") if fit and fit.get("valid") else None)
        r["selected_fit_rmse_relative"]=(fit.get("fit_rmse_relative") if fit and fit.get("valid") else None)
        if fit and fit.get("valid"):
            center_pts.append([fit["center_x_native"],fit["center_y_native"],r["z_input"]])
            center_colors.append(np.array([220,40,40,255],dtype=np.uint8) if r.get("center_outlier") else np.array([20,120,220,255],dtype=np.uint8))

    # additional QC summary
    both_offsets=[r["outer_inner_center_offset_mm"] for r in horizontal_records if isinstance(r.get("outer_inner_center_offset_mm"),(float,int))]
    axis_summary["mean_outer_inner_center_offset_mm"] = float(np.mean(both_offsets)) if both_offsets else None
    axis_summary["median_outer_inner_center_offset_mm"] = float(np.median(both_offsets)) if both_offsets else None
    selected_rmses=[r["selected_fit_rmse_mm"] for r in horizontal_records if r.get("used_for_axis") and r.get("selected_fit_rmse_mm") is not None]
    axis_summary["mean_selected_fit_rmse_mm"] = float(np.mean(selected_rmses)) if selected_rmses else None
    axis_summary["median_selected_fit_rmse_mm"] = float(np.median(selected_rmses)) if selected_rmses else None
    axis_summary["wall_thickness_used_for_axis"] = False
    axis_summary["wall_thickness_supported_later_for_volume_reconstruction"] = bool(axis_surface=="outer")
    axis_summary["z_sampling"] = z_sampling

    write_horizontal_csv(axis_dir/"horizontal_sections.csv",horizontal_records)
    write_axis_summary_csv(axis_dir/"axis_summary.csv",axis_summary)
    write_combined_points(axis_dir/"all_horizontal_section_points.ply",horizontal_groups)
    write_combined_points(axis_dir/"all_fitted_ellipse_points.ply",ellipse_groups)
    if center_pts:
        pc=trimesh.points.PointCloud(np.asarray(center_pts,float),colors=np.asarray(center_colors,np.uint8)); pc.export(str(axis_dir/"section_centers_points.ply"))
    write_axis_edge_ply(axis_dir/"rotation_axis_edges.ply",center_native,zmin,zmax)

    # extract longitudinal sections using estimated axis
    spacing_native=sample_spacing_mm/scale; full_cache={}; full_groups=[]; ray_groups=[]; summary_rows=[]
    for i,angle in enumerate(full_angles):
        seg=section_segments(mesh,center_native,angle); full_cache[i]=seg; opposite=(angle+180)%360
        stem=f"section_{int(round(angle))%360:03d}_{int(round(opposite))%360:03d}"
        write_segment_csv(full_dir/f"{stem}.csv",seg,center_native,angle,scale,zmin); write_edge_ply(full_dir/f"{stem}_edges.ply",seg)
        sampled=sample_segments(seg,spacing_native); color=hsv_color(i,len(full_angles)); write_points_ply(full_dir/f"{stem}_points.ply",sampled,color); full_groups.append((sampled,color))
        summary_rows.append({"type":"full_section","angle_deg":angle,"opposite_angle_deg":opposite,"segments":len(seg),"sampled_points":len(sampled)})
        print(f"full {angle:6.1f}/{opposite:6.1f} deg : {len(seg):,} segments")
    for j,ray_angle in enumerate(ray_angles):
        base_index=j%len(full_angles); seg_half=clip_segments_to_positive_half(full_cache[base_index],center_native,ray_angle)
        stem=f"ray_{int(round(ray_angle))%360:03d}"; write_segment_csv(ray_dir/f"{stem}.csv",seg_half,center_native,ray_angle,scale,zmin); write_edge_ply(ray_dir/f"{stem}_edges.ply",seg_half)
        sampled=sample_segments(seg_half,spacing_native); color=hsv_color(j,len(ray_angles)); write_points_ply(ray_dir/f"{stem}_points.ply",sampled,color); ray_groups.append((sampled,color))
        summary_rows.append({"type":"radial_half_section","angle_deg":ray_angle,"opposite_angle_deg":"","segments":len(seg_half),"sampled_points":len(sampled)})
    write_combined_points(full_dir/"all_full_section_points.ply",full_groups); write_combined_points(ray_dir/"all_radial_half_section_points.ply",ray_groups)
    overlay_summary = export_section_overlay(
        output_dir, [(a, full_cache[i]) for i, a in enumerate(full_angles)],
        center_native, scale, sample_spacing_mm, unit, enabled=section_overlay,
    )
    print(f"XY overlay / median: {overlay_summary['status']} "
          f"({overlay_summary.get('valid_sections', 0)} valid full sections)")
    if overlay_summary.get("reason"):
        print(f"median note: {overlay_summary['reason']}")
    rim_summary = export_rim_standardization(
        output_dir, [(a, full_cache[i]) for i, a in enumerate(full_angles)],
        center_native, scale, unit, spacing_mm=sample_spacing_mm, **(rim_options or {}),
    )
    print(f"right-oriented rim standardization: {rim_summary['status']} "
          f"({rim_summary.get('valid_profiles', 0)} valid half-section profiles)")
    with (output_dir/"sections_summary.csv").open("w",newline="",encoding="utf-8-sig") as f:
        fields=["type","angle_deg","opposite_angle_deg","segments","sampled_points"]; w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(summary_rows)

    print("\n=== Capacity calculation ===")
    print(f"volume mode(s) : {', '.join(volume_modes)}")
    if "single" in volume_modes or "all" in volume_modes:
        print(f"single angle   : {single_angle_deg:g} deg")
    print(f"volume Z step  : {volume_z_step_mm:g} mm")
    volume_summary=calculate_capacity_modes(
        mesh=mesh, center_native=center_native, scale_to_mm=scale,
        ray_angles=ray_angles, full_cache=full_cache, full_angles=full_angles,
        horizontal_records=horizontal_records, requested_modes=volume_modes,
        single_angle_deg=single_angle_deg, volume_z_step_mm=volume_z_step_mm,
        min_angular_valid_fraction=min_angular_valid_fraction,
        volume_dir=volume_dir, make_visualization=make_visualization,
    )
    for vr in volume_summary.get("results",[]):
        if vr.get("status")=="ok":
            print(f"{vr['method']:10s}: {vr['volume_l']:.6f} L")
        else:
            print(f"{vr['method']:10s}: FAILED - {vr.get('note','')}")

    if make_visualization:
        make_oblique_plot(mesh,horizontal_groups,vis_dir/"horizontal_sections_oblique.png",scale,"Horizontal sections used for rotation-axis estimation")
        make_oblique_plot(mesh,full_groups,vis_dir/"full_sections_oblique.png",scale,"Full longitudinal sections through estimated rotation axis")
        make_oblique_plot(mesh,ray_groups,vis_dir/"radial_half_sections_oblique.png",scale,"Radial half-sections through estimated rotation axis")
        # validation sections through final axis at exactly 0 and 90 degrees, independent of start angle
        full0=section_segments(mesh,center_native,0.0); full90=section_segments(mesh,center_native,90.0)
        make_axis_validation_plots(mesh,center_mm,horizontal_records,full0,full90,vis_dir,scale)
        make_table_images(horizontal_records,axis_summary,vis_dir)

    metadata={
        "program":"PotteryRadialSections","version":__version__,"volume_core_version":volume_core.__version__,"input":str(input_path),"input_unit":unit,
        "mesh_vertices":int(len(mesh.vertices)),"mesh_faces":int(len(mesh.faces)),
        "axis_direction":[0.0,0.0,1.0],"axis_center_xy_native":[center_native[0],center_native[1]],"axis_center_xy_mm":[center_mm[0],center_mm[1]],
        "axis_method":"horizontal-section robust ellipse centers + MAD center-outlier rejection + inlier arithmetic mean",
        "axis_surface_mode":axis_surface,"axis_center_source":source,"axis_summary":axis_summary,
        "z_sampling":z_sampling,"angle_step_deg":step_deg,"start_angle_deg":start_deg,
        "full_section_plane_angles_deg":full_angles,"radial_half_section_angles_deg":ray_angles,
        "sample_spacing_mm":sample_spacing_mm,"contour_resample_spacing_mm":contour_spacing_mm,
        "section_ply_coordinate_system":"same XYZ coordinate system and unit as input mesh",
        "wall_thickness_used_for_axis":False,
        "outer_mode_note":"If axis_surface=outer, wall thickness is intentionally not used here; single or Z-positioned thickness measurements can be used later for inner-profile reconstruction and volume calculation.",
        "tilt_note":"Center drift versus Z is diagnostic only; the model is not auto-rotated or tilt-corrected in v0.9.0.",
        "volume_reconstruction_implemented":True,
        "volume_calculation":volume_summary,
        "section_overlay":overlay_summary,
        "rim_standardization":rim_summary,
    }
    write_json(output_dir/"metadata.json",metadata)
    print("\n=== Rotation axis ===")
    print(f"X = {center_mm[0]:.6f} mm")
    print(f"Y = {center_mm[1]:.6f} mm")
    print(f"used sections = {axis_summary['used_for_axis']} / {axis_summary['attempted_sections']}")
    print(f"RMS radial center offset = {axis_summary['rms_radial_offset_mm']:.4f} mm")
    print(f"diagnostic centerline tilt = {axis_summary['estimated_centerline_tilt_deg']:.4f} deg (not corrected)")
    print(f"output dir = {output_dir}")
    return metadata


def main():
    parser=argparse.ArgumentParser(description="Estimate a pottery rotation axis, extract radial sections, and calculate vessel capacity with selectable methods.")
    parser.add_argument("input",type=Path,help="Input PLY or OBJ mesh")
    parser.add_argument("--unit",choices=["mm","cm","m"],required=True,help="Coordinate unit of input mesh")
    parser.add_argument("--axis-surface",choices=["inner","outer","both"],default="inner",help="Contour used for final axis. 'both' fits both for QC but uses inner centers (default: inner)")
    zgroup=parser.add_mutually_exclusive_group()
    zgroup.add_argument("--z-sections",type=int,default=None,help="Number of equal-height intervals; section at each interval center (default: 20)")
    zgroup.add_argument("--z-step-mm",type=float,help="Absolute horizontal-section interval in mm instead of relative equal divisions")
    parser.add_argument("--center-outlier-mad",type=float,default=3.0,help="MAD multiplier for rejecting horizontal-section center outliers (default: 3.0)")
    parser.add_argument("--contour-spacing-mm",type=float,default=0.5,help="Arc-length resampling interval used for ellipse fitting [mm] (default: 0.5)")
    parser.add_argument("--angle-step",type=float,default=30.0,help="Longitudinal radial angular interval [deg] (default: 30)")
    parser.add_argument("--start-angle",type=float,default=0.0,help="Start azimuth [deg] (default: 0)")
    parser.add_argument("--sample-spacing-mm",type=float,default=0.5,help="PLY point-cloud sampling interval [mm] (default: 0.5)")
    parser.add_argument("--volume-mode",nargs="+",choices=["all","single","optimized","angular","ellipse","none"],default=["all"],
                        help="Capacity method(s). Use one or more: all, single, optimized, angular, ellipse, none (default: all)")
    parser.add_argument("--single-angle",type=float,default=0.0,
                        help="Radial half-profile angle used by single volume mode [deg] (default: 0)")
    parser.add_argument("--volume-z-step-mm",type=float,default=0.5,
                        help="Z sampling step for longitudinal-profile volume integration [mm] (default: 0.5)")
    parser.add_argument("--min-angular-valid-fraction",type=float,default=0.75,
                        help="Minimum fraction of radial directions required at a Z cell for optimized/angular modes (default: 0.75)")
    parser.add_argument("--output-dir",type=Path,help="Output directory")
    parser.add_argument("--no-visualization",action="store_true",help="Do not generate PNG reference/QC images")
    parser.add_argument("--no-section-overlay",action="store_true",help="Disable XY section overlays and lip-registered median polyline (enabled by default)")
    parser.add_argument("--no-rim-standardization",action="store_true",help="Disable cropped right-oriented rim midlines and within-vessel registration")
    parser.add_argument("--rim-buffer-mm",type=float,help="Arc-length buffer past the horizontal stall/reversal [mm]; default: 2x local paired wall separation")
    parser.add_argument("--rim-buffer-thickness-ratio",type=float,default=2.0,help="Automatic buffer as a multiple of local paired wall separation (default: 2)")
    parser.add_argument("--rim-end-mm",type=float,help="Manually select an arc-length endpoint from the lip [mm], replacing transition+buffer detection")
    parser.add_argument("--rim-smooth-mm",type=float,default=2.0,help="Tangent smoothing and sustained direction length [mm] (default: 2)")
    parser.add_argument("--rim-x-progress-tol",type=float,default=0.1,help="Horizontal stall tolerance: dx/ds >= -tol after sustained leftward progress (default: 0.1); ignores vertical direction")
    parser.add_argument("--rim-turn-angle-deg",type=float,default=None,help="Legacy alias: uses sin(angle) as X progress tolerance; no vertical-direction test")
    parser.add_argument("--rim-points",type=int,default=129,help="Number of correspondingly sampled rim points (default: 129)")
    parser.add_argument("--rim-endpoint-weight",type=float,default=5.0,help="Endpoint weight relative to each interior point (default: 5)")
    parser.add_argument("--rim-affine-anisotropy",type=float,default=0.1,help="Bound on opposite directional stretch factors (default: 0.1)")
    parser.add_argument("--rim-affine-shear",type=float,default=0.1,help="Absolute affine shear bound (default: 0.1)")
    parser.add_argument("--rim-affine-penalty",type=float,default=1.0,help="Regularization toward similarity for stretch/shear (default: 1)")
    parser.add_argument("--rim-outlier-mad",type=float,default=3.5,help="Profile-wise robust rejection threshold for the additional inliers section model (default: 3.5)")
    parser.add_argument("--no-rim-bimodal",action="store_true",help="Disable exploratory A/B section-model screening")
    parser.add_argument("--rim-bimodal-min-profiles",type=int,default=20,help="Minimum accepted radial sides for A/B screening (default: 20; minimum 10)")
    parser.add_argument("--rim-bimodal-bic-delta",type=float,default=10.0,help="Minimum BIC improvement of the two-component first-PC mixture (default: 10)")
    parser.add_argument("--version",action="version",version=f"%(prog)s {__version__}")
    args=parser.parse_args()
    try:
        process(args.input,args.unit,args.angle_step,args.start_angle,args.axis_surface,args.z_sections,args.z_step_mm,
                args.center_outlier_mad,args.contour_spacing_mm,args.sample_spacing_mm,args.output_dir,not args.no_visualization,
                args.volume_mode,args.single_angle,args.volume_z_step_mm,args.min_angular_valid_fraction,
                not args.no_section_overlay,
                {"enabled": not args.no_rim_standardization,
                 "buffer_mm": args.rim_buffer_mm, "buffer_thickness_ratio": args.rim_buffer_thickness_ratio,
                 "end_mm": args.rim_end_mm, "smooth_mm": args.rim_smooth_mm,
                 "turn_angle_deg": args.rim_turn_angle_deg, "points": args.rim_points,
                 "x_progress_tol": args.rim_x_progress_tol,
                 "endpoint_weight": args.rim_endpoint_weight,
                 "affine_anisotropy": args.rim_affine_anisotropy,
                 "affine_shear": args.rim_affine_shear, "affine_penalty": args.rim_affine_penalty,
                 "outlier_mad": args.rim_outlier_mad, "bimodal": not args.no_rim_bimodal,
                 "bimodal_min_profiles": args.rim_bimodal_min_profiles,
                 "bimodal_bic_delta": args.rim_bimodal_bic_delta})
    except Exception as exc:
        parser.exit(1,f"ERROR: {exc}\n")

if __name__=="__main__":
    main()

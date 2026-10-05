"""Right-oriented, cropped rim midlines and within-vessel registration.

Midlines use the radial midpoint of outer/inner walls at the same input height.
They are not a medial-axis skeleton or a normal-thickness reconstruction.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import trimesh
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

from .section_overlay import project_section, registered_branches, write_polyline_ply


def arc_positions(points):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]


def sample_curve(points, positions):
    s = arc_positions(points)
    keep = np.r_[True, np.diff(s) > 1e-10]
    if keep.sum() < 2:
        raise ValueError("zero-length curve")
    return np.column_stack([np.interp(positions, s[keep], points[keep, j]) for j in range(2)])


def side_to_axis(points):
    """Keep the lip-to-axis part; require descent instead of bridging folds."""
    for i, (a, b) in enumerate(zip(points[:-1], points[1:])):
        if b[0] <= 0:
            t = a[0] / (a[0] - b[0])
            return np.vstack([points[:i + 1], a + t * (b - a)])
    raise ValueError("wall branch does not reach the axis")


def radius_at_heights(branch, heights, outer):
    if np.any(np.diff(branch[:, 1]) > 1e-5):
        raise ValueError("non-monotone wall height; automatic radial midline is ambiguous")
    # At a flat lip segment, use its wall-side extreme; the midpoint pairs caps.
    y, groups = np.unique(branch[:, 1], return_inverse=True)
    if len(y) < 2:
        raise ValueError("insufficient wall-height samples")
    x = np.full(len(y), -np.inf if outer else np.inf)
    (np.maximum.at if outer else np.minimum.at)(x, groups, branch[:, 0])
    return np.interp(heights, y, x)


def make_midline(outer_full, inner_return, side, spacing_mm):
    if side == "right":
        outer, inner = outer_full.copy(), inner_return[::-1].copy()
    elif side == "left":
        outer, inner = outer_full[::-1].copy(), inner_return.copy()
        outer[:, 0] *= -1
        inner[:, 0] *= -1
    else:
        raise ValueError("side must be right or left")
    outer, inner = side_to_axis(outer[:, :2]), side_to_axis(inner[:, :2])
    # Stop at each wall's first bottom ordinate, before the horizontal floor.
    outer = outer[:int(np.argmin(outer[:, 1])) + 1]
    inner = inner[:int(np.argmin(inner[:, 1])) + 1]
    high = min(outer[0, 1], inner[0, 1])
    low = max(outer[:, 1].min(), inner[:, 1].min())
    if high - low < 3 * spacing_mm:
        raise ValueError("wall too short for rim extraction")
    # Avoid the exact flat inner-floor ordinate, which does not describe wall radius.
    y = np.linspace(high, low, max(8, int(np.ceil((high - low) / spacing_mm)) + 1), endpoint=False)
    ro, ri = radius_at_heights(outer, y, True), radius_at_heights(inner, y, False)
    thickness = ro - ri
    if np.any(thickness <= 0):
        raise ValueError("non-positive radial wall thickness")
    midpoint = np.column_stack([(ro + ri) / 2, y])
    # Explicit lip anchor, shared by the outer and inner branches.
    midpoint[0] = outer[0]
    return midpoint, thickness


def crop_rim(midline, thickness, *, spacing_mm=.5, smooth_mm=2.,
             turn_angle_deg=10., buffer_mm=None, buffer_thickness_ratio=2., end_mm=None):
    s = arc_positions(midline)
    grid = np.linspace(0, s[-1], max(8, int(np.ceil(s[-1] / spacing_mm)) + 1))
    line = sample_curve(midline, grid)
    wall = np.interp(grid, s, thickness)
    meta = {"full_arc_length_mm": float(s[-1])}
    if end_mm is not None:
        if not 0 < end_mm <= s[-1]:
            raise ValueError("manual rim end is outside the available midline")
        end = float(end_mm)
        meta.update(selection="manual_arc_length", change_arc_length_mm=None, buffer_mm=None)
    else:
        ds = float(grid[1] - grid[0])
        smooth = gaussian_filter1d(line, max(.5, smooth_mm / (2 * ds)), axis=0, mode="nearest")
        tangent = np.gradient(smooth, axis=0)
        tangent /= np.maximum(np.linalg.norm(tangent, axis=1, keepdims=True), 1e-12)
        threshold = np.sin(np.deg2rad(turn_angle_deg))
        descending = tangent[:, 1] < -.5
        left_down = descending & (tangent[:, 0] < -threshold)
        down_or_right = descending & (tangent[:, 0] >= -threshold)
        run = max(2, int(np.ceil(smooth_mm / ds)))
        seen_left = False
        change = None
        for i in range(1, len(line) - run):
            if i >= run and np.all(left_down[i - run:i]):
                seen_left = True
            if seen_left and np.all(down_or_right[i:i + run]):
                change = i
                break
        if change is None:
            raise ValueError("no sustained left-down to down/right-down transition; specify --rim-end-mm to select the region manually")
        buffer = float(buffer_mm if buffer_mm is not None else wall[change] * buffer_thickness_ratio)
        end = float(grid[change] + buffer)
        if end > s[-1]:
            raise ValueError("insufficient midline beyond the transition for the requested buffer")
        meta.update(selection="direction_transition", change_arc_length_mm=float(grid[change]),
                    change_point_mm=line[change].tolist(), local_radial_thickness_mm=float(wall[change]), buffer_mm=buffer)
    positions = np.r_[grid[grid < end], end]
    cropped = sample_curve(midline, positions)
    if len(cropped) < 4:
        raise ValueError("cropped rim has too few points")
    meta.update(end_arc_length_mm=end, start_point_mm=cropped[0].tolist(), end_point_mm=cropped[-1].tolist())
    return cropped, meta


def fit_similarity(source, target, weights=None, scale_enabled=True):
    weights = np.ones(len(source)) if weights is None else np.asarray(weights)
    weights = weights / weights.sum()
    ca, cb = weights @ source, weights @ target
    a, b = source - ca, target - cb
    u, _, vt = np.linalg.svd((a * weights[:, None]).T @ b)
    sign = np.eye(2)
    sign[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation_row = u @ sign @ vt
    denominator = float(np.sum(weights[:, None] * a * a))
    if denominator <= 1e-12:
        raise ValueError("degenerate source curve")
    scale = float(np.sum(weights[:, None] * (a @ rotation_row) * b) / denominator) if scale_enabled else 1.
    if scale <= 0:
        raise ValueError("non-positive similarity scale")
    matrix = scale * rotation_row.T
    translation = cb - matrix @ ca
    return matrix, translation


def apply_transform(points, matrix, translation):
    return points @ matrix.T + translation


def fit_constrained_affine(source, target, weights, *, anisotropy_limit=.1, shear_limit=.1, penalty=1.):
    initial, trans = fit_similarity(source, target, weights)
    scale = float(np.sqrt(np.linalg.det(initial)))
    theta = float(np.arctan2(initial[1, 0], initial[0, 0]))
    w = weights / weights.sum()
    center = w @ source
    target_center = w @ target
    centered = source - center
    characteristic = max(float(np.sqrt(np.sum(w[:, None] * centered ** 2))), 1e-8)
    # A line alone cannot constrain deformation in its perpendicular direction.
    eigen = np.linalg.eigvalsh((centered * w[:, None]).T @ centered)
    if eigen[0] / max(eigen[-1], 1e-12) < 1e-4:
        return initial, trans, {"status": "similarity_fallback", "reason": "near-collinear midline"}

    def build(v):
        t, logscale, k, h = v[:4]
        c, s = np.cos(t), np.sin(t)
        rotation = np.array([[c, -s], [s, c]])
        stretch = np.diag([np.exp(logscale + k), np.exp(logscale - k)])
        shear = np.array([[1., h], [0., 1.]])
        return rotation @ stretch @ shear

    def residual(v):
        mapped = centered @ build(v).T + target_center + characteristic * v[4:]
        data = ((mapped - target) / characteristic * np.sqrt(w[:, None])).ravel()
        # Prior toward the similarity solution; parameters are dimensionless.
        return np.r_[data, np.sqrt(penalty) * v[2:4]]

    kmax = np.log1p(anisotropy_limit)
    v0 = np.array([theta, np.log(scale), 0., 0., 0., 0.])
    lower = [theta - np.pi / 2, np.log(scale) - np.log(2), -kmax, -shear_limit, -2, -2]
    upper = [theta + np.pi / 2, np.log(scale) + np.log(2), kmax, shear_limit, 2, 2]
    result = least_squares(residual, v0, bounds=(lower, upper), loss="soft_l1", f_scale=.1)
    matrix = build(result.x)
    translation = target_center + characteristic * result.x[4:] - matrix @ center
    return matrix, translation, {"status": "ok" if result.success else "optimizer_failed",
                                 "log_directional_stretch": float(result.x[2]), "shear": float(result.x[3]),
                                 "bound_hit": bool(np.any(result.active_mask))}


def unit_shape(points):
    center = points.mean(axis=0)
    size = float(np.linalg.norm(points - center))
    if size <= 1e-12:
        raise ValueError("zero centroid size")
    return (points - center) / size, center, size


def _write_curves(path, curves, unit_to_mm):
    with Path(path).open("w", encoding="ascii", newline="\n") as f:
        n = sum(len(p) for p in curves)
        e = sum(len(p) - 1 for p in curves)
        f.write(f"ply\nformat ascii 1.0\nelement vertex {n}\nproperty double x\nproperty double y\nproperty double z\nelement edge {e}\nproperty int vertex1\nproperty int vertex2\nend_header\n")
        for p in curves:
            for x, y in p / unit_to_mm:
                f.write(f"{x:.15g} {y:.15g} 0\n")
        offset = 0
        for p in curves:
            for i in range(len(p) - 1):
                f.write(f"{offset+i} {offset+i+1}\n")
            offset += len(p)


def export_rim_standardization(output_dir, sections, center_xy, unit_to_mm, input_unit,
                               *, enabled=True, spacing_mm=.5, smooth_mm=2., turn_angle_deg=10.,
                               buffer_mm=None, buffer_thickness_ratio=2., end_mm=None, points=129,
                               endpoint_weight=5., affine_anisotropy=.1, affine_shear=.1, affine_penalty=1.):
    out = Path(output_dir) / "rim_standardization"
    out.mkdir(parents=True, exist_ok=True)
    # Only delete files owned by this exporter. There are no user-provided profiles here.
    known = ["raw_midlines_xy.ply", "raw_midlines.csv", "transforms.json"]
    for mode in ("similarity", "affine"):
        known += [f"{mode}_midlines_xy.ply", f"standard_{mode}_midline_xy.ply", f"standard_{mode}_midline.csv", f"standard_{mode}_unit_shape.csv"]
    for name in known:
        (out / name).unlink(missing_ok=True)
    config = dict(spacing_mm=spacing_mm, smooth_mm=smooth_mm, turn_angle_deg=turn_angle_deg,
                  buffer_mm=buffer_mm, buffer_thickness_ratio=buffer_thickness_ratio,
                  end_mm=end_mm, points=points, endpoint_weight=endpoint_weight,
                  affine_anisotropy=affine_anisotropy, affine_shear=affine_shear, affine_penalty=affine_penalty)
    qa = {"enabled": enabled, "status": "disabled", "config": config, "input_unit": input_unit,
          "unit_to_mm": unit_to_mm, "profiles": [], "side_reference": "right; left reflected before fitting",
          "midline_definition": "radial midpoint of inner/outer wall at equal input height",
          "xy_frame": "X=right-oriented radius, Y=original input Z, Z=0",
          "interpretation": "within-vessel standardization; not recovery of the unfired form or a completed between-vessel Procrustes analysis"}
    if enabled:
        positive = [spacing_mm, smooth_mm, buffer_thickness_ratio, endpoint_weight,
                    affine_anisotropy, affine_shear, affine_penalty]
        if any(not np.isfinite(v) or v <= 0 for v in positive) or points < 8:
            raise ValueError("rim spacings, ratios, weights and affine constraints must be positive; points >=8")
        if not 0 < turn_angle_deg < 45:
            raise ValueError("rim turn angle must be in (0,45) degrees")
        if buffer_mm is not None and (not np.isfinite(buffer_mm) or buffer_mm < 0):
            raise ValueError("rim buffer must be finite and non-negative")
        if end_mm is not None and (not np.isfinite(end_mm) or end_mm <= 0):
            raise ValueError("manual rim end must be finite and positive")
        accepted, records = [], []
        for angle, segments in sections:
            try:
                outer, inner = registered_branches(project_section(segments, center_xy, angle) * unit_to_mm)
            except (ValueError, IndexError) as exc:
                qa["profiles"].append({"angle_deg": float(angle), "status": "excluded", "reason": str(exc)})
                continue
            for side in ("right", "left"):
                record = {"angle_deg": float(angle), "side": side}
                try:
                    midline, thickness = make_midline(outer, inner, side, spacing_mm)
                    crop, meta = crop_rim(midline, thickness, spacing_mm=spacing_mm, smooth_mm=smooth_mm,
                                          turn_angle_deg=turn_angle_deg, buffer_mm=buffer_mm,
                                          buffer_thickness_ratio=buffer_thickness_ratio, end_mm=end_mm)
                    sample = sample_curve(crop, np.linspace(0, arc_positions(crop)[-1], points))
                    accepted.append(sample)
                    record.update(status="ok", **meta)
                    records.append(record)
                except (ValueError, IndexError) as exc:
                    record.update(status="excluded", reason=str(exc))
                qa["profiles"].append(record)
        qa.update(valid_profiles=len(accepted), attempted_profiles=2 * len(sections))
        if accepted:
            _write_curves(out / "raw_midlines_xy.ply", accepted, unit_to_mm)
            with (out / "raw_midlines.csv").open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["profile_id", "angle_deg", "side", "point_id", "u", "x_input", "y_input", "x_mm", "y_mm"])
                for i, (curve, rec) in enumerate(zip(accepted, records)):
                    for j, p in enumerate(curve):
                        writer.writerow([i, rec["angle_deg"], rec["side"], j, j / (points - 1), *(p / unit_to_mm), *p])
        right = [p for p, rec in zip(accepted, records) if rec["side"] == "right"]
        if len(accepted) >= 2 and right:
            # A fixed right-side median preserves gauge (size, position, orientation).
            reference = np.median(np.stack(right), axis=0)
            weights = np.ones(points)
            weights[[0, -1]] = endpoint_weight
            transforms = []
            for mode in ("similarity", "affine"):
                aligned = []
                for index, source in enumerate(accepted):
                    if mode == "similarity":
                        matrix, trans = fit_similarity(source, reference, weights)
                        detail = {"status": "ok"}
                    else:
                        matrix, trans, detail = fit_constrained_affine(source, reference, weights,
                            anisotropy_limit=affine_anisotropy, shear_limit=affine_shear, penalty=affine_penalty)
                    if detail["status"] == "optimizer_failed":
                        matrix, trans = fit_similarity(source, reference, weights)
                        detail.update(status="similarity_fallback", reason="affine optimizer failed")
                    mapped = apply_transform(source, matrix, trans)
                    aligned.append(mapped)
                    transforms.append({"profile_id": index, "mode": mode, "angle_deg": records[index]["angle_deg"],
                                       "side": records[index]["side"], "matrix_2x2": matrix.tolist(),
                                       "translation_mm": trans.tolist(), "determinant": float(np.linalg.det(matrix)),
                                       "principal_scales": np.linalg.svd(matrix, compute_uv=False).tolist(),
                                       "rms_before_mm": float(np.sqrt(np.mean(np.sum((source-reference)**2, axis=1)))),
                                       "rms_after_mm": float(np.sqrt(np.mean(np.sum((mapped-reference)**2, axis=1)))), **detail})
                standard = np.median(np.stack(aligned), axis=0)
                deviation = np.sqrt(np.mean(np.sum((np.stack(aligned) - standard) ** 2, axis=2), axis=0))
                shape, centroid, size = unit_shape(standard)
                _write_curves(out / f"{mode}_midlines_xy.ply", aligned, unit_to_mm)
                write_polyline_ply(out / f"standard_{mode}_midline_xy.ply", np.column_stack([standard / unit_to_mm, np.zeros(points)]), closed=False)
                with (out / f"standard_{mode}_midline.csv").open("w", encoding="utf-8-sig", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["point_id", "u", "x_input", "y_input", "x_mm", "y_mm", "rms_deviation_mm"])
                    for j, p in enumerate(standard):
                        writer.writerow([j, j / (points - 1), *(p / unit_to_mm), *p, deviation[j]])
                with (out / f"standard_{mode}_unit_shape.csv").open("w", encoding="utf-8-sig", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["point_id", "u", "x_unitless", "y_unitless"])
                    for j, p in enumerate(shape):
                        writer.writerow([j, j / (points - 1), *p])
                qa[mode] = {"centroid_mm": centroid.tolist(), "centroid_size_mm": size,
                            "reference": "fixed coordinate-wise median of accepted right-side cropped midlines",
                            "rms_deviation_mm": float(np.sqrt(np.mean(deviation ** 2)))}
            (out / "transforms.json").write_text(json.dumps(transforms, ensure_ascii=False, indent=2), encoding="utf-8")
            qa["status"] = "ok"
        else:
            qa.update(status="insufficient_valid_profiles", reason="requires at least two accepted profiles and one right-side reference")
    (out / "rim_qa.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    return qa

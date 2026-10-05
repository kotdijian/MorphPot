"""Right-oriented, cropped rim midlines and within-vessel registration.

Midlines use locally paired outer/inner wall points in lip-to-end contour order.
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
    """Keep the ordered lip-to-axis part without a height-monotonicity test."""
    for i, (a, b) in enumerate(zip(points[:-1], points[1:])):
        if b[0] <= 0:
            t = a[0] / (a[0] - b[0])
            return np.vstack([points[:i + 1], a + t * (b - a)])
    raise ValueError("wall branch does not reach the axis")


def horizontal_stop(points, spacing_mm, smooth_mm, x_progress_tol):
    """First sustained stall/reversal after sustained leftward motion; ignore Y sign."""
    s = arc_positions(points)
    if s[-1] <= 2 * smooth_mm:
        raise ValueError("curve too short for a sustained horizontal transition")
    grid = np.linspace(0, s[-1], max(8, int(np.ceil(s[-1] / spacing_mm)) + 1))
    line = sample_curve(points, grid)
    ds = float(grid[1] - grid[0])
    x = gaussian_filter1d(line[:, 0], max(.5, smooth_mm / (2 * ds)), mode="nearest")
    dx_ds = np.gradient(x, grid)
    run = max(2, int(np.ceil(smooth_mm / ds)))
    seen_left = False
    for i in range(run, len(line) - run + 1):
        if np.all(dx_ds[i-run:i] < -x_progress_tol):
            seen_left = True
        if seen_left and np.all(dx_ds[i:i+run] >= -x_progress_tol):
            return float(grid[i]), line[i], float(dx_ds[i])
    raise ValueError("no sustained leftward stall/reversal; specify --rim-end-mm manually")


def _prefix(points, end):
    s = arc_positions(points)
    end = min(float(end), float(s[-1]))
    return np.vstack([points[s < end], sample_curve(points, [end])])


def ordered_wall_pairs(outer, inner, spacing_mm):
    """Local distance-minimizing DTW correspondence, preserving lip-to-end order.

    A 25% normalized arc-position band limits correspondence drift; each local
    wall is sampled at <=512 points to bound memory and runtime.
    """
    lengths = [arc_positions(p)[-1] for p in (outer, inner)]
    counts = [min(512, max(8, int(np.ceil(v / spacing_mm)) + 1)) for v in lengths]
    a, b = [sample_curve(p, np.linspace(0, length, count))
            for p, length, count in zip((outer, inner), lengths, counts)]
    n, m = len(a), len(b)
    cost = np.sum((a[:, None] - b[None, :]) ** 2, axis=2)
    cost += (.02 * max(lengths)) ** 2 * (np.linspace(0, 1, n)[:, None] - np.linspace(0, 1, m)[None, :]) ** 2
    score = np.full((n + 1, m + 1), np.inf)
    score[0, 0] = 0
    parent = np.zeros((n, m), dtype=np.uint8)
    for i in range(n):
        lo = max(0, int(np.floor((i / (n-1) - .25) * (m-1))))
        hi = min(m, int(np.ceil((i / (n-1) + .25) * (m-1))) + 1)
        for j in range(lo, hi):
            candidates = [score[i, j], score[i, j+1], score[i+1, j]]
            step = int(np.argmin(candidates))
            score[i+1, j+1] = cost[i, j] + candidates[step]
            parent[i, j] = step
    if not np.isfinite(score[n, m]):
        raise ValueError("local wall correspondence outside order-preserving band")
    pairs = []
    i, j = n-1, m-1
    while True:
        pairs.append((i, j))
        if i == 0 and j == 0:
            break
        step = parent[i, j]
        if step == 0:
            i -= 1; j -= 1
        elif step == 1:
            i -= 1
        else:
            j -= 1
        if i < 0 or j < 0:
            raise ValueError("invalid wall correspondence path")
    pairs = np.asarray(pairs[::-1])
    return a[pairs[:, 0]], b[pairs[:, 1]], counts


def make_midline(outer_full, inner_return, side, spacing_mm, *, smooth_mm=2.,
                 x_progress_tol=.1, buffer_mm=None, buffer_thickness_ratio=2.,
                 end_mm=None, diagnostics=None, source_geometry=None):
    if side == "right":
        outer, inner = outer_full.copy(), inner_return[::-1].copy()
    elif side == "left":
        outer, inner = outer_full[::-1].copy(), inner_return.copy()
        outer[:, 0] *= -1
        inner[:, 0] *= -1
    else:
        raise ValueError("side must be right or left")
    outer, inner = side_to_axis(outer[:, :2]), side_to_axis(inner[:, :2])
    meta = {"wall_pairing": "local order-preserving DTW; 25% band; at most 512 samples per wall"}
    if end_mm is None:
        so, po, _ = horizontal_stop(outer, spacing_mm, smooth_mm, x_progress_tol)
        si, pi, _ = horizontal_stop(inner, spacing_mm, smooth_mm, x_progress_tol)
        separation = float(np.linalg.norm(po - pi))
        if separation <= 1e-8 and buffer_mm is None:
            raise ValueError("zero wall separation at transition")
        provisional_buffer = float(buffer_mm if buffer_mm is not None else separation * buffer_thickness_ratio)
        # Pair a local prefix with a detection margin, then crop the central line.
        # Lower-body geometry never participates in this correspondence.
        ends = [so + provisional_buffer + 2*smooth_mm, si + provisional_buffer + 2*smooth_mm]
        if any(end - 2*smooth_mm > arc_positions(p)[-1] for end, p in zip(ends, (outer, inner))):
            raise ValueError("insufficient wall beyond the transition for the requested buffer")
        meta.update(outer_stop_arc_mm=so, inner_stop_arc_mm=si,
                    outer_stop_point_mm=po.tolist(), inner_stop_point_mm=pi.tolist())
    else:
        if any(end_mm > arc_positions(p)[-1] for p in (outer, inner)):
            raise ValueError("manual rim end exceeds an available wall branch")
        ends = [end_mm + 2*smooth_mm + 4*spacing_mm] * 2
    outer, inner = [_prefix(p, end) for p, end in zip((outer, inner), ends)]
    a, b, counts = ordered_wall_pairs(outer, inner, spacing_mm)
    midpoint = (a + b) / 2
    thickness = np.linalg.norm(a - b, axis=1)
    midpoint[0] = outer[0]
    meta.update(wall_pairing_samples=counts,
                local_wall_arc_lengths_mm=[float(arc_positions(p)[-1]) for p in (outer, inner)])
    if source_geometry is not None:
        source_geometry.update(outer=outer.copy(), inner=inner.copy(), paired_outer=a.copy(), paired_inner=b.copy())
    if diagnostics is not None:
        diagnostics.update(meta)
    return midpoint, thickness


def crop_rim(midline, thickness, *, spacing_mm=.5, smooth_mm=2.,
             turn_angle_deg=None, x_progress_tol=.1, buffer_mm=None, buffer_thickness_ratio=2., end_mm=None):
    s = arc_positions(midline)
    grid = np.linspace(0, s[-1], max(8, int(np.ceil(s[-1] / spacing_mm)) + 1))
    line = sample_curve(midline, grid)
    meta = {"full_arc_length_mm": float(s[-1])}
    if end_mm is not None:
        if not 0 < end_mm <= s[-1]:
            raise ValueError("manual rim end is outside the available midline")
        end = float(end_mm)
        meta.update(selection="manual_arc_length", change_arc_length_mm=None, buffer_mm=None)
    else:
        if turn_angle_deg is not None:
            x_progress_tol = float(np.sin(np.deg2rad(turn_angle_deg)))
        change_s, change_point, progress = horizontal_stop(midline, spacing_mm, smooth_mm, x_progress_tol)
        local_thickness = float(np.interp(change_s, s, thickness))
        buffer = float(buffer_mm if buffer_mm is not None else local_thickness * buffer_thickness_ratio)
        end = float(change_s + buffer)
        if end > s[-1]:
            raise ValueError("insufficient midline beyond the transition for the requested buffer")
        meta.update(selection="horizontal_stall_or_reversal", change_arc_length_mm=change_s,
                    change_point_mm=change_point.tolist(), x_progress_at_change=progress,
                    local_paired_wall_separation_mm=local_thickness, buffer_mm=buffer)
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


def _write_curves(path, curves, unit_to_mm, color=None):
    with Path(path).open("w", encoding="ascii", newline="\n") as f:
        n = sum(len(p) for p in curves)
        e = sum(len(p) - 1 for p in curves)
        colors = "property uchar red\nproperty uchar green\nproperty uchar blue\n" if color is not None else ""
        f.write(f"ply\nformat ascii 1.0\nelement vertex {n}\nproperty double x\nproperty double y\nproperty double z\n{colors}element edge {e}\nproperty int vertex1\nproperty int vertex2\nend_header\n")
        for p in curves:
            for x, y in p / unit_to_mm:
                rgb = " " + " ".join(map(str, color)) if color is not None else ""
                f.write(f"{x:.15g} {y:.15g} 0{rgb}\n")
        offset = 0
        for p in curves:
            for i in range(len(p) - 1):
                f.write(f"{offset+i} {offset+i+1}\n")
            offset += len(p)


def _write_verification(out, prefix, geometry, pairs, unit_to_mm):
    # Source polylines retain the original section vertices (plus prefix endpoint).
    # They include the local detection margin; connectors show the exact compared pairs.
    for key, color in (("outer", (230, 70, 50)), ("inner", (50, 100, 230))):
        _write_curves(out / f"{prefix}_source_{key}_xy.ply",
                      [g[key] for g in geometry], unit_to_mm, color)
    connectors = [np.vstack((a, b)) for pair in pairs for a, b in zip(*pair)]
    _write_curves(out / f"{prefix}_pair_connectors_xy.ply", connectors,
                  unit_to_mm, (170, 170, 170))


def export_rim_standardization(output_dir, sections, center_xy, unit_to_mm, input_unit,
                               *, enabled=True, spacing_mm=.5, smooth_mm=2., turn_angle_deg=None, x_progress_tol=.1,
                               buffer_mm=None, buffer_thickness_ratio=2., end_mm=None, points=129,
                               endpoint_weight=5., affine_anisotropy=.1, affine_shear=.1, affine_penalty=1.):
    out = Path(output_dir) / "rim_standardization"
    out.mkdir(parents=True, exist_ok=True)
    # Only delete files owned by this exporter. There are no user-provided profiles here.
    known = ["raw_midlines_xy.ply", "raw_midlines.csv", "transforms.json"]
    for mode in ("similarity", "affine"):
        known += [f"{mode}_midlines_xy.ply", f"standard_{mode}_midline_xy.ply", f"standard_{mode}_midline.csv", f"standard_{mode}_unit_shape.csv"]
    for prefix in ("raw", "similarity", "affine"):
        known += [f"{prefix}_source_outer_xy.ply", f"{prefix}_source_inner_xy.ply", f"{prefix}_pair_connectors_xy.ply"]
    known += ["raw_paired_points.csv"]
    for name in known:
        (out / name).unlink(missing_ok=True)
    config = dict(spacing_mm=spacing_mm, smooth_mm=smooth_mm, turn_angle_deg=turn_angle_deg,
                  x_progress_tol=x_progress_tol,
                  buffer_mm=buffer_mm, buffer_thickness_ratio=buffer_thickness_ratio,
                  end_mm=end_mm, points=points, endpoint_weight=endpoint_weight,
                  affine_anisotropy=affine_anisotropy, affine_shear=affine_shear, affine_penalty=affine_penalty)
    qa = {"enabled": enabled, "status": "disabled", "config": config, "input_unit": input_unit,
          "unit_to_mm": unit_to_mm, "profiles": [], "side_reference": "right; left reflected before fitting",
          "midline_definition": "midpoints of local order-preserving outer/inner wall pairs; no height monotonicity requirement",
          "xy_frame": "X=right-oriented radius, Y=original input Z, Z=0",
          "interpretation": "within-vessel standardization; not recovery of the unfired form or a completed between-vessel Procrustes analysis"}
    if enabled:
        positive = [spacing_mm, smooth_mm, buffer_thickness_ratio, endpoint_weight,
                    affine_anisotropy, affine_shear, affine_penalty]
        if any(not np.isfinite(v) or v <= 0 for v in positive) or points < 8:
            raise ValueError("rim spacings, ratios, weights and affine constraints must be positive; points >=8")
        if turn_angle_deg is not None:
            if not 0 < turn_angle_deg < 45:
                raise ValueError("legacy rim turn angle must be in (0,45) degrees")
            x_progress_tol = float(np.sin(np.deg2rad(turn_angle_deg)))
        if not np.isfinite(x_progress_tol) or not 0 < x_progress_tol < 1:
            raise ValueError("rim X progress tolerance must be in (0,1)")
        qa["config"]["effective_x_progress_tol"] = x_progress_tol
        if buffer_mm is not None and (not np.isfinite(buffer_mm) or buffer_mm < 0):
            raise ValueError("rim buffer must be finite and non-negative")
        if end_mm is not None and (not np.isfinite(end_mm) or end_mm <= 0):
            raise ValueError("manual rim end must be finite and positive")
        accepted, records, source_geometries, source_pairs = [], [], [], []
        for angle, segments in sections:
            try:
                outer, inner = registered_branches(project_section(segments, center_xy, angle) * unit_to_mm)
            except (ValueError, IndexError) as exc:
                qa["profiles"].append({"angle_deg": float(angle), "status": "excluded", "reason": str(exc)})
                continue
            for side in ("right", "left"):
                record = {"angle_deg": float(angle), "side": side}
                try:
                    geometry = {}
                    midline, thickness = make_midline(outer, inner, side, spacing_mm,
                        smooth_mm=smooth_mm, x_progress_tol=x_progress_tol, buffer_mm=buffer_mm,
                        buffer_thickness_ratio=buffer_thickness_ratio, end_mm=end_mm, diagnostics=record, source_geometry=geometry)
                    crop, meta = crop_rim(midline, thickness, spacing_mm=spacing_mm, smooth_mm=smooth_mm,
                                          x_progress_tol=x_progress_tol, buffer_mm=buffer_mm,
                                          buffer_thickness_ratio=buffer_thickness_ratio, end_mm=end_mm)
                    # Use one interpolation parameter for both walls and their midpoint.
                    positions = np.linspace(0, meta["end_arc_length_mm"], points)
                    arc = arc_positions(midline)
                    keep = np.r_[True, np.diff(arc) > 1e-10]
                    pair = [np.column_stack([np.interp(positions, arc[keep], wall[keep, k])
                                             for k in range(2)])
                            for wall in (geometry["paired_outer"], geometry["paired_inner"])]
                    sample = (pair[0] + pair[1]) / 2
                    source_geometries.append(geometry)
                    source_pairs.append(pair)
                    accepted.append(sample)
                    record.update(status="ok", **meta)
                    records.append(record)
                except (ValueError, IndexError) as exc:
                    record.update(status="excluded", reason=str(exc))
                qa["profiles"].append(record)
        qa.update(valid_profiles=len(accepted), attempted_profiles=2 * len(sections))
        if accepted:
            _write_curves(out / "raw_midlines_xy.ply", accepted, unit_to_mm, (40, 190, 70))
            _write_verification(out, "raw", source_geometries, source_pairs, unit_to_mm)
            qa["verification"] = {"source_scope": "original local wall prefixes including detection margin",
                "colors": {"outer": "red", "inner": "blue", "midline": "green", "pairs": "gray"},
                "pair_points_per_profile": points, "coordinates": "same XY frame and input unit as midlines"}
            with (out / "raw_paired_points.csv").open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["profile_id", "angle_deg", "side", "point_id", "u",
                                 "outer_x_mm", "outer_y_mm", "inner_x_mm", "inner_y_mm", "mid_x_mm", "mid_y_mm"])
                for i, (pair, curve, rec) in enumerate(zip(source_pairs, accepted, records)):
                    for j, (a, b, mid) in enumerate(zip(*pair, curve)):
                        writer.writerow([i, rec["angle_deg"], rec["side"], j, j/(points-1), *a, *b, *mid])
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
                aligned, mapped_geometry, mapped_pairs = [], [], []
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
                    mapped_geometry.append({key: apply_transform(source_geometries[index][key], matrix, trans)
                                            for key in ("outer", "inner")})
                    mapped_pairs.append([apply_transform(p, matrix, trans) for p in source_pairs[index]])
                    transforms.append({"profile_id": index, "mode": mode, "angle_deg": records[index]["angle_deg"],
                                       "side": records[index]["side"], "matrix_2x2": matrix.tolist(),
                                       "translation_mm": trans.tolist(), "determinant": float(np.linalg.det(matrix)),
                                       "principal_scales": np.linalg.svd(matrix, compute_uv=False).tolist(),
                                       "rms_before_mm": float(np.sqrt(np.mean(np.sum((source-reference)**2, axis=1)))),
                                       "rms_after_mm": float(np.sqrt(np.mean(np.sum((mapped-reference)**2, axis=1)))), **detail})
                standard = np.median(np.stack(aligned), axis=0)
                deviation = np.sqrt(np.mean(np.sum((np.stack(aligned) - standard) ** 2, axis=2), axis=0))
                shape, centroid, size = unit_shape(standard)
                _write_curves(out / f"{mode}_midlines_xy.ply", aligned, unit_to_mm, (40, 190, 70))
                _write_verification(out, mode, mapped_geometry, mapped_pairs, unit_to_mm)
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

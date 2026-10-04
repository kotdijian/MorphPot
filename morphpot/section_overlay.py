"""Axis-relative section overlays and lip-registered full-contour medians."""
from __future__ import annotations

import csv
import colorsys
import json
from pathlib import Path

import numpy as np
import trimesh


def project_section(segments, center_xy, angle_deg):
    """Native XYZ -> (signed radius, original height, 0); no scaling."""
    segments = np.asarray(segments, dtype=float).reshape(-1, 2, 3)
    angle = np.deg2rad(angle_deg)
    radial = np.array([np.cos(angle), np.sin(angle)])
    result = np.zeros_like(segments)
    result[..., 0] = (segments[..., :2] - np.asarray(center_xy)) @ radial
    result[..., 1] = segments[..., 2]
    return result


def write_polyline_ply(path, points, closed=True):
    """Write shared, ordered vertices with PLY edges (not independent pairs)."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    edges = [(i, i + 1) for i in range(len(points) - 1)]
    if closed and len(points) > 2:
        edges.append((len(points) - 1, 0))
    with Path(path).open("w", encoding="ascii", newline="\n") as f:
        f.write("ply\nformat ascii 1.0\n")
        f.write(f"element vertex {len(points)}\n")
        f.write("property double x\nproperty double y\nproperty double z\n")
        f.write(f"element edge {len(edges)}\n")
        f.write("property int vertex1\nproperty int vertex2\nend_header\n")
        for p in points:
            f.write(" ".join(f"{v:.15g}" for v in p) + "\n")
        for a, b in edges:
            f.write(f"{a} {b}\n")


def _lip_vertex(points, positive):
    """Find highest lip on one side, centering an exactly flat top edge."""
    indices = np.flatnonzero(points[:, 0] > 0 if positive else points[:, 0] < 0)
    if not len(indices):
        raise ValueError("contour does not contain both sides of the axis")
    top = float(points[indices, 1].max())
    tol = max(float(np.ptp(points[:, 1])) * 1e-9, 1e-8)
    candidates = indices[np.abs(points[indices, 1] - top) <= tol]
    target_x = float((points[candidates, 0].min() + points[candidates, 0].max()) / 2)
    # Insert an anchor on an existing flat lip segment, never bridge a gap.
    for i in candidates:
        j = (int(i) + 1) % len(points)
        if j not in candidates:
            continue
        a, b = points[i], points[j]
        if min(a[0], b[0]) < target_x < max(a[0], b[0]):
            t = (target_x - a[0]) / (b[0] - a[0])
            anchor = a + t * (b - a)
            points = np.insert(points, int(i) + 1, anchor, axis=0)
            return points, int(i) + 1
    idx = candidates[np.argmin(np.abs(points[candidates, 0] - target_x))]
    return points, int(idx)


def registered_branches(projected_mm):
    """Return right lip -> outer bottom -> left lip, then inner return.

    A single closed material contour spanning the axis is required. Disconnected
    wall fragments cannot be assigned a continuous full median automatically.
    """
    path = trimesh.load_path(np.asarray(projected_mm, dtype=float))
    # Path.discrete may omit dangling entities; check the whole graph first.
    if any(degree != 2 for _, degree in path.vertex_graph.degree):
        raise ValueError("open or branching section edges found")
    contours = []
    open_count = 0
    for discrete in path.discrete:
        p = np.asarray(discrete, dtype=float)
        if len(p) < 4 or np.linalg.norm(p[0] - p[-1]) > 1e-6:
            open_count += 1
            continue
        p = p[:-1]
        if p[:, 0].min() < 0 < p[:, 0].max():
            contours.append(p)
    if open_count or len(contours) != 1 or len(path.discrete) != 1:
        raise ValueError("requires one closed full material contour; open or multiple components found")
    p, _ = _lip_vertex(contours[0], True)
    p, _ = _lip_vertex(p, False)
    p, right = _lip_vertex(p, True)
    p = np.roll(p, -right, axis=0)
    p, left = _lip_vertex(p, False)
    first = p[:left + 1]
    second = np.vstack([p[left:], p[0]])[::-1]
    # Both candidate paths start at the positive lip and end at the negative lip.
    # The outer-bottom path is lower than the inner-floor path.
    if abs(float(first[:, 1].min() - second[:, 1].min())) <= 1e-8:
        raise ValueError("outer bottom and inner floor cannot be distinguished")
    outer, inner = (first, second) if first[:, 1].min() < second[:, 1].min() else (second, first)
    return outer, inner[::-1]


def _length(points):
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def _resample(points, count):
    distances = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    keep = np.r_[True, np.diff(distances) > 1e-12]
    if np.count_nonzero(keep) < 2:
        raise ValueError("zero-length branch")
    d, p = distances[keep], points[keep]
    t = np.linspace(0, d[-1], count)
    return np.column_stack([np.interp(t, d, p[:, j]) for j in range(3)])


def export_section_overlay(output_dir, sections, center_xy, unit_to_mm,
                           spacing_mm, input_unit, enabled=True):
    """Project all full sections; median uses only valid lip-registered contours."""
    out = Path(output_dir) / "section_overlay"
    out.mkdir(parents=True, exist_ok=True)
    # Remove only our old products so a failed/disabled rerun cannot show stale medians.
    products = ["all_full_sections_xy_points.ply", "all_full_sections_xy_edges.ply",
                "median_full_section_xy.ply", "median_full_section_xy.csv"]
    for name in products:
        (out / name).unlink(missing_ok=True)
    qa = {
        "enabled": bool(enabled), "input_unit": input_unit,
        "projection": {"x": "signed radius from axis", "y": "original input Z height", "z": 0},
        "axis_center_xy_input": list(map(float, center_xy)),
        "unit_to_mm": float(unit_to_mm), "sampling_spacing_mm": float(spacing_mm),
        "median_method": "coordinate-wise median after separate normalized arc-length registration of outer and inner branches",
        "start": "positive-side lip (highest point; midpoint for a flat lip segment)",
        "direction": "positive lip -> outer bottom -> negative lip -> inner floor -> positive lip",
        "lip_detection": "geometric highest point on each side; archaeological lip identification is not guaranteed",
        "sections": [], "status": "disabled" if not enabled else "pending",
    }
    if enabled:
        groups, projected_all, valid = [], [], []
        for index, (angle, segments) in enumerate(sections):
            projected = project_section(segments, center_xy, angle)
            projected_all.append(projected)
            pts = []
            for a, b in projected:
                n = max(2, int(np.ceil(np.linalg.norm(b - a) * unit_to_mm / spacing_mm)) + 1)
                pts.append(np.linspace(a, b, n))
            if pts:
                rgb = colorsys.hsv_to_rgb(index / max(1, len(sections)), .75, .95)
                color = np.array([round(c * 255) for c in rgb] + [255], dtype=np.uint8)
                groups.append((np.vstack(pts), color))
            record = {"angle_deg": float(angle), "segment_count": len(projected)}
            try:
                outer, inner = registered_branches(projected * unit_to_mm)
                valid.append((outer, inner))
                record.update(status="ok", positive_lip_mm=outer[0].tolist(), negative_lip_mm=outer[-1].tolist())
            except (ValueError, IndexError) as exc:
                record.update(status="excluded", reason=str(exc))
            qa["sections"].append(record)
        all_segments = np.concatenate(projected_all) if projected_all else np.empty((0, 2, 3))
        # Keep all projected segments even when they are unusable for the median.
        with (out / products[1]).open("w", encoding="ascii", newline="\n") as f:
            points = all_segments.reshape(-1, 3)
            f.write("ply\nformat ascii 1.0\n")
            f.write(f"element vertex {len(points)}\nproperty double x\nproperty double y\nproperty double z\n")
            f.write(f"element edge {len(all_segments)}\nproperty int vertex1\nproperty int vertex2\nend_header\n")
            for p in points:
                f.write(" ".join(f"{v:.15g}" for v in p) + "\n")
            for i in range(len(all_segments)):
                f.write(f"{2*i} {2*i+1}\n")
        if groups:
            points = np.vstack([p for p, _ in groups])
            colors = np.vstack([np.tile(c, (len(p), 1)) for p, c in groups])
            trimesh.points.PointCloud(points, colors=colors).export(str(out / products[0]))
        qa["valid_sections"] = len(valid)
        qa["attempted_sections"] = len(sections)
        if len(valid) >= 2:
            counts = [max(3, int(np.ceil(max(_length(v[j]) for v in valid) / spacing_mm)) + 1) for j in (0, 1)]
            medians = [np.median(np.stack([_resample(v[j], counts[j]) for v in valid]), axis=0) for j in (0, 1)]
            # Opposite lip is shared once; closure is an edge to vertex 0.
            median_mm = np.vstack([medians[0], medians[1][1:-1]])
            median_native = median_mm / unit_to_mm
            write_polyline_ply(out / products[2], median_native)
            with (out / products[3]).open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["vertex_id", "branch", "x_input", "y_input", "z_input", "x_mm", "y_mm", "z_mm"])
                for i, (native, mm) in enumerate(zip(median_native, median_mm)):
                    writer.writerow([i, "outer" if i < counts[0] else "inner", *native, *mm])
            qa.update(status="ok", median_vertices=len(median_native), closed=True,
                      positive_lip_vertex=0, negative_lip_vertex=counts[0] - 1,
                      valid_fraction=len(valid) / len(sections))
        else:
            qa.update(status="insufficient_valid_sections", reason="at least two closed full contours required; no median written")
    (out / "overlay_qa.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    return qa

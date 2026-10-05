"""Replace a highest-point lip anchor by a centreline extension intersection."""
from __future__ import annotations

import numpy as np


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def extension_tip(outer, inner, paired_outer, paired_inner, spacing_mm, smooth_mm):
    """Fit the first stable wall centreline and intersect its outward ray.

    The old common highest point is only a contour partition, not the new tip.
    The short cap-to-support bridge is interpolated, and is explicitly reported.
    No intersection means exclusion rather than silently retaining the old tip.
    """
    mid = (paired_outer + paired_inner) / 2
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(mid, axis=0), axis=1))]
    keep = np.r_[True, np.diff(arc) > 1e-10]
    arc, mid = arc[keep], mid[keep]
    a, b = paired_outer[keep], paired_inner[keep]
    separation = np.linalg.norm(a-b, axis=1)
    head = (arc <= min(.6*arc[-1], max(6*smooth_mm, 20*spacing_mm))) & (separation > 1e-8)
    if head.sum() < 3:
        raise ValueError("insufficient wall pairs for centreline-extension tip")
    width = float(np.percentile(separation[head], 75))
    start = max(1.5*width, smooth_mm, 3*spacing_mm)
    window = max(width, smooth_mm, 4*spacing_mm)
    if start + window > arc[-1]:
        raise ValueError("rim too short for centreline-extension tip support")
    grid = np.linspace(start, start+window, 16)
    support = np.column_stack([np.interp(grid, arc, mid[:, k]) for k in range(2)])
    slope = (grid-grid.mean()) @ (support-support.mean(axis=0))
    norm = float(np.linalg.norm(slope))
    if norm < 1e-10:
        raise ValueError("degenerate centreline-extension direction")
    direction = -slope/norm
    origin = support[0]
    candidates = []
    for name, contour in (("outer", outer), ("inner", inner)):
        wall_arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(contour, axis=0), axis=1))]
        for index, (v, w) in enumerate(zip(contour[:-1], contour[1:])):
            edge = w-v
            det = _cross(direction, edge)
            if abs(det) < 1e-12:
                continue
            t = _cross(v-origin, edge)/det
            u = _cross(v-origin, direction)/det
            # Restrict to the local lip, rejecting a distant body intersection.
            source_arc = wall_arc[index] + np.clip(u, 0, 1)*np.linalg.norm(edge)
            if t > 1e-8 and -1e-9 <= u <= 1+1e-9 and source_arc <= start+window+2*width:
                candidates.append((t, origin+t*direction, name, index, float(u)))
    if not candidates:
        raise ValueError("centreline extension has no local lip-contour intersection")
    distance, tip, branch, segment, fraction = min(candidates, key=lambda r: r[0])
    # Keep paired-wall correspondence beyond the stable support, replacing only
    # the curled cap portion. Both walls coincide at the true contour endpoint.
    start_a = np.array([np.interp(start, arc, a[:, k]) for k in range(2)])
    start_b = np.array([np.interp(start, arc, b[:, k]) for k in range(2)])
    rest = arc > start + 1e-10
    new_a = np.vstack([tip, start_a, a[rest]])
    new_b = np.vstack([tip, start_b, b[rest]])
    metadata = {"tip_method": "stable centreline outward extension intersected with original lip contour",
                "tip_point_mm": tip.tolist(), "old_highest_lip_mm": paired_outer[0].tolist(),
                "tip_direction_outward": direction.tolist(), "tip_support_point_mm": origin.tolist(),
                "tip_support_old_arc_mm": start, "tip_fit_window_mm": window,
                "tip_intersection_branch": branch, "tip_intersection_segment": segment,
                "tip_intersection_fraction": fraction, "tip_extension_length_mm": float(distance),
                "tip_bridge": "linear paired-wall interpolation from contour endpoint to stable support; not measured normal thickness"}
    return new_a, new_b, metadata

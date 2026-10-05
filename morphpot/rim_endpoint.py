"""Neutral geometric endpoint candidates, without archaeological vessel labels."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d


LABELS = {"radial_turn": "半径方向の停滞・反転", "horizontal_transition": "高さ方向変化の減少・水平化",
          "near_vertical": "直立傾向", "no_clear_transition": "明瞭な移行なし", "compound": "複合",
          "undetermined": "判定保留"}
END_MODES = ("auto", "radial_turn", "horizontal_transition", "manual")


def detect_candidates(points, spacing_mm=.5, smooth_mm=2., x_progress_tol=.1, horizontal_angle_deg=10.):
    p = np.asarray(points, dtype=float)[:, :2]
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    keep = np.r_[True, np.diff(arc) > 1e-10]
    result = {"candidates": [], "categories": [], "status": "ok"}
    if keep.sum() < 3 or arc[-1] <= 4*smooth_mm:
        result.update(status="too_short", categories=["undetermined"])
        return result
    grid = np.linspace(0, arc[-1], max(8, int(np.ceil(arc[-1]/spacing_mm))+1))
    line = np.column_stack([np.interp(grid, arc[keep], p[keep, j]) for j in range(2)])
    ds = float(grid[1]-grid[0])
    filtered = gaussian_filter1d(line, max(.5, smooth_mm/(2*ds)), axis=0, mode="nearest")
    tangent = np.gradient(filtered, grid, axis=0)
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1)[:, None], 1e-10)
    tr, th = tangent.T
    run = max(2, int(np.ceil(smooth_mm/ds)))
    guard = max(smooth_mm, 3*spacing_mm)
    first = max(run, int(np.ceil(guard/ds)))
    low = np.sin(np.deg2rad(horizontal_angle_deg))
    high = np.sin(np.deg2rad(horizontal_angle_deg+10))
    seen_left = seen_oblique = False
    found = set()
    for i in range(first, len(grid)-run+1):
        previous = slice(max(first, i-run), i)
        if i >= first+run:
            if np.all(tr[previous] < -x_progress_tol):
                seen_left = True
            if np.all((tr[previous] < -x_progress_tol) & (np.abs(th[previous]) >= high)):
                seen_oblique = True
        after = slice(i, i+run)
        kinds = []
        if "radial_turn" not in found and seen_left and np.all(tr[after] >= -x_progress_tol):
            kinds.append("radial_turn")
        if ("horizontal_transition" not in found and seen_oblique and
                np.all((tr[after] < -x_progress_tol) & (np.abs(th[after]) <= low))):
            kinds.append("horizontal_transition")
        for kind in kinds:
            found.add(kind)
            result["candidates"].append({"kind": kind, "label": LABELS[kind],
                "arc_length_mm": float(grid[i]), "point_mm": line[i].tolist(),
                "tangent_radial": float(tr[i]), "tangent_height": float(th[i]),
                "inclination_deg": float(np.rad2deg(np.arctan2(abs(th[i]), abs(tr[i])))),
                "sustain_length_mm": float(run*ds)})
    result["candidates"].sort(key=lambda x: x["arc_length_mm"])
    result["categories"] = [kind for kind in ("radial_turn", "horizontal_transition") if kind in found]
    if len(found) > 1:
        result["categories"].append("compound")
    probe = (grid >= guard) & (grid <= min(arc[-1], guard+max(6*smooth_mm, 10*spacing_mm)))
    vertical_fraction = float(np.mean(np.abs(tr[probe]) <= .15)) if np.any(probe) else 0.
    result["near_vertical_fraction"] = vertical_fraction
    if vertical_fraction >= .75:
        result["categories"].append("near_vertical")
    if not result["categories"]:
        result["categories"] = ["no_clear_transition"]
    return result


def select_candidate(detection, mode):
    if mode not in END_MODES:
        raise ValueError("invalid rim end mode")
    candidates = [c for c in detection["candidates"] if mode == "auto" or c["kind"] == mode]
    if not candidates:
        raise ValueError(f"no sustained endpoint candidate for {mode}; specify --rim-end-mm manually")
    return min(candidates, key=lambda c: c["arc_length_mm"])


def shape_summary(profiles):
    records = [r for r in profiles if "shape_detection" in r]
    n = len(records)
    categories = []
    for kind in LABELS:
        supporters = [r for r in records if kind in r["shape_detection"]["categories"]]
        if not supporters:
            continue
        arcs = [c["arc_length_mm"] for r in supporters for c in r["shape_detection"]["candidates"] if c["kind"] == kind]
        heights = [c["point_mm"][1] for r in supporters for c in r["shape_detection"]["candidates"] if c["kind"] == kind]
        categories.append({"category": kind, "label": LABELS[kind], "support_count": len(supporters),
            "support_fraction": len(supporters)/n,
            "supporters": [{"angle_deg": r["angle_deg"], "side": r["side"]} for r in supporters],
            "candidate_arc_median_mm": float(np.median(arcs)) if arcs else None,
            "candidate_arc_iqr_mm": float(np.subtract(*np.percentile(arcs, [75,25]))) if arcs else None,
            "candidate_height_iqr_mm": float(np.subtract(*np.percentile(heights, [75,25]))) if heights else None})
    candidates = [r for r in categories if r["category"] in ("radial_turn", "horizontal_transition")]
    ranked = sorted(candidates, key=lambda r: (-r["support_fraction"], r["candidate_arc_median_mm"], r["category"]))
    recommended = ranked[0]["category"] if ranked else None
    ambiguous = len(ranked) > 1 and ranked[0]["support_fraction"]-ranked[1]["support_fraction"] <= .15
    return {"evaluated_half_profiles": n, "categories": categories, "recommended_end_mode": recommended,
            "multiple_supported_modes": len(ranked) > 1, "ambiguous_recommendation": ambiguous,
            "selection_rule": "largest support fraction; tie chooses earlier median candidate arc",
            "status": "ok" if n else "undetermined", "interpretation": "geometric candidates; support fractions are not calibrated probabilities; no vessel-type labels",
            "arc_origin": "highest contour partition for wall screening; final selected midpoint endpoints use extension tip",
            "aggregation": "one vote per half-profile; categories may overlap; no forced single class"}


def write_endpoint_products(out, records, scale):
    rows = []
    for index, record in enumerate(records):
        for c in record.get("shape_detection", {}).get("candidates", []):
            rows.append(dict(profile_index=index, angle_deg=record["angle_deg"], side=record["side"],
                             stage="wall_screening", kind=c["kind"], point=c["point_mm"], arc=c["arc_length_mm"]))
        if record.get("status") == "ok":
            for kind, field in (("tip", "start_point_mm"), ("comparison_end", "end_point_mm")):
                rows.append(dict(profile_index=index, angle_deg=record["angle_deg"], side=record["side"], stage="final_midline",
                                 kind=kind, point=record[field], arc=0. if kind == "tip" else record["end_arc_length_mm"]))
            for c in record.get("midline_endpoint_candidates", []):
                rows.append(dict(profile_index=index, angle_deg=record["angle_deg"], side=record["side"], stage="final_midline",
                                 kind=c["kind"], point=c["point_mm"], arc=c["arc_length_mm"]))
    with (Path(out)/"endpoint_candidates.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["profile_index", "angle_deg", "side", "stage", "kind", "arc_length_mm", "x_mm", "height_mm", "x_input", "height_input"])
        for r in rows:
            writer.writerow([r["profile_index"],r["angle_deg"],r["side"],r["stage"],r["kind"],r["arc"],*r["point"],*(np.asarray(r["point"])/scale)])
    colors = {"radial_turn": (230,140,30), "horizontal_transition": (160,70,210), "tip": (30,30,30), "comparison_end": (30,190,80)}
    with (Path(out)/"endpoint_candidates_xy.ply").open("w",encoding="ascii") as f:
        f.write(f"ply\nformat ascii 1.0\nelement vertex {len(rows)}\nproperty double x\nproperty double y\nproperty double z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n")
        for r in rows:
            x,y = np.asarray(r["point"])/scale
            f.write(f"{x:.15g} {y:.15g} 0 {' '.join(map(str,colors[r['kind']]))}\n")

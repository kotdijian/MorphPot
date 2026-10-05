"""Empirical within-vessel section models; no restoration labels are inferred.

All valid profiles have equal weight. Filtering is profile-wise and auditable.
Two-mode screening is exploratory (correlated radial sections are not independent).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from .section_overlay import write_polyline_ply


def normal_frame(curve, outer_offsets):
    tangent = np.gradient(curve, axis=0)
    lengths = np.linalg.norm(tangent, axis=1)
    if np.any(lengths < 1e-10):
        raise ValueError("standard curve has a degenerate tangent")
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]]) / lengths[:, None]
    # One global sign maintains continuity rather than flipping per sample.
    if np.median(np.sum(outer_offsets[:, 1:] * normal[None, 1:], axis=2)) < 0:
        normal *= -1
    return normal


def robust_profile_scores(features, distance_floor_mm):
    center = np.median(features, axis=0)
    scale = np.maximum(1.4826*np.median(np.abs(features-center), axis=0), distance_floor_mm)
    residual = np.abs(features-center)/scale
    # A profile is the unit of rejection, so no holes are punched into its curve.
    return np.percentile(residual, 95, axis=1)


def _ray_distances(origins, normals, geometry):
    starts = np.vstack([geometry[k][:-1] for k in ("outer", "inner")])
    ends = np.vstack([geometry[k][1:] for k in ("outer", "inner")])
    edges = ends-starts
    delta = starts[None]-origins[:, None]
    denominator = normals[:, None, 0]*edges[None, :, 1]-normals[:, None, 1]*edges[None, :, 0]
    safe = np.where(np.abs(denominator) > 1e-12, denominator, np.nan)
    t = (delta[..., 0]*edges[None, :, 1]-delta[..., 1]*edges[None, :, 0])/safe
    u = (delta[..., 0]*normals[:, None, 1]-delta[..., 1]*normals[:, None, 0])/safe
    valid = (t > 1e-9) & (u >= -1e-9) & (u <= 1+1e-9)
    return np.min(np.where(valid, t, np.inf), axis=1)


def measured_distances(mids, pairs, normal, geometries=None):
    normals = np.broadcast_to(normal, mids.shape) if np.asarray(normal).ndim == 2 else np.asarray(normal)
    if normals.shape != mids.shape:
        raise ValueError("normal frames must match the individual midlines")
    outer_offsets, inner_offsets = pairs[:, 0]-mids, pairs[:, 1]-mids
    signed_outer = np.sum(outer_offsets*normals, axis=2)
    signed_inner = -np.sum(inner_offsets*normals, axis=2)
    distances = np.stack([np.abs(signed_outer), np.abs(signed_inner)], axis=1)
    measured = np.zeros_like(distances, dtype=bool)
    if geometries is not None:
        for i, geometry in enumerate(geometries):
            for side, direction in enumerate((normals[i], -normals[i])):
                values = _ray_distances(mids[i], direction, geometry)
                valid = np.isfinite(values)
                distances[i, side, valid] = values[valid]
                measured[i, side, valid] = True
    distances[:, :, 0] = 0
    measured[:, :, 0] = False  # true tip, no separate wall-distance observation
    return distances, measured, signed_outer, signed_inner


def individual_normal_frames(mids, pairs):
    """Measure thickness in each profile's local frame before model synthesis.

    A common model normal would measure an oblique chord through a differently
    oriented individual wall, inflating a parallel-wall distance by 1/cos(angle).
    """
    return np.stack([normal_frame(curve, (pair[0]-curve)[None]) for curve, pair in zip(mids, pairs)])


def screen_two_modes(features, *, min_profiles=20, bic_delta=10., min_fraction=.2):
    """Deterministic 1D Gaussian mixtures on a winsorized first PC.

    BIC advantage alone is insufficient: require two density peaks, separated
    components and substantial support in each group. This does not test an
    archaeological cause, nor all possible bimodal variation in higher PCs.
    """
    n = len(features)
    report = {"status": "insufficient_profiles", "min_profiles": min_profiles,
              "bic_delta_required": bic_delta, "minimum_group_fraction": min_fraction,
              "minimum_group_count": 5, "minimum_separation": 2.,
              "method": "winsorized first-PC, 1D 1-vs-2 Gaussian mixtures; BIC plus two density peaks",
              "interpretation": "exploratory; radial profiles of one vessel are correlated; A/B are not restoration labels"}
    if n < min_profiles:
        return None, report
    lo, hi = np.percentile(features, [5, 95], axis=0)
    centered = np.clip(features, lo, hi)
    centered -= np.median(centered, axis=0)
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    if singular[0] < 1e-8:
        report["status"] = "no_variation"
        return None, report
    axis = vt[0]
    if axis[np.argmax(np.abs(axis))] < 0:
        axis *= -1
    x = centered @ axis
    variance = max(float(np.var(x)), 1e-10)
    floor = max(variance*1e-4, 1e-10)
    single_ll = float(np.sum(-.5*(np.log(2*np.pi*variance)+(x-x.mean())**2/variance)))
    best = None
    for percentiles in ([25, 75], [15, 85], [35, 65]):
        means = np.percentile(x, percentiles)
        variances = np.full(2, variance/2)
        weights = np.full(2, .5)
        previous = -np.inf
        converged = False
        for _ in range(300):
            logs = np.log(weights)[None]-.5*(np.log(2*np.pi*variances)[None]+(x[:, None]-means)**2/variances)
            norm = logsumexp(logs, axis=1)
            ll = float(norm.sum())
            responsibilities = np.exp(logs-norm[:, None])
            masses = np.maximum(responsibilities.sum(axis=0), 1e-10)
            means = responsibilities.T @ x / masses
            variances = np.maximum(np.sum(responsibilities*(x[:, None]-means)**2, axis=0)/masses, floor)
            weights = masses/n
            if abs(ll-previous) < 1e-7*(1+abs(ll)):
                converged = True
                break
            previous = ll
        logs = np.log(weights)[None]-.5*(np.log(2*np.pi*variances)[None]+(x[:, None]-means)**2/variances)
        ll = float(logsumexp(logs, axis=1).sum())
        if best is None or ll > best[0]:
            best = (ll, means.copy(), variances.copy(), weights.copy(), logs.copy(), converged)
    ll, means, variances, weights, logs, converged = best
    order = np.argsort(means)
    means, variances, weights, logs = means[order], variances[order], weights[order], logs[:, order]
    labels = logs.argmax(axis=1)
    counts = np.bincount(labels, minlength=2)
    bic_one, bic_two = -2*single_ll+2*np.log(n), -2*ll+5*np.log(n)
    separation = float(abs(np.diff(means)[0])/np.sqrt(variances.mean()))
    grid = np.linspace(float(x.min()-np.sqrt(variance)), float(x.max()+np.sqrt(variance)), 1024)
    density = np.sum(weights[None]*np.exp(-.5*(grid[:, None]-means)**2/variances)/np.sqrt(2*np.pi*variances), axis=1)
    peaks = np.flatnonzero((density[1:-1] > density[:-2]) & (density[1:-1] > density[2:]))+1
    peaks = peaks[density[peaks] > .05*density.max()]
    passed = (converged and bic_one-bic_two >= bic_delta and separation >= 2 and
              np.all(counts >= max(5, int(np.ceil(n*min_fraction)))) and
              np.all(weights >= min_fraction) and len(peaks) == 2)
    report.update(status="two_modes" if passed else "no_supported_split", bic_one=float(bic_one),
                  bic_two=float(bic_two), bic_improvement=float(bic_one-bic_two), separation=separation,
                  counts=counts.tolist(), weights=weights.tolist(), means=means.tolist(),
                  variances=variances.tolist(), density_peaks=len(peaks), converged=converged,
                  pc1_variance_fraction=float(singular[0]**2/np.sum(singular**2)), scores=x.tolist())
    return (labels if passed else None), report


def _write_model(out, stem, mids, pairs, indices, unit_to_mm, write_curves, geometries=None):
    chosen = np.asarray(indices, dtype=int)
    curves = mids[chosen]
    standard = np.median(curves, axis=0)
    outer_offsets = pairs[chosen, 0]-curves
    normal = normal_frame(standard, outer_offsets)
    measurement_normals = individual_normal_frames(curves, pairs[chosen])
    angular_difference = np.rad2deg(np.arccos(np.clip(np.sum(measurement_normals*normal[None], axis=2), -1, 1)))
    selected_geometry = None if geometries is None else [geometries[i] for i in chosen]
    distances, measured, signed_outer, signed_inner = measured_distances(curves, pairs[chosen], measurement_normals, selected_geometry)
    quantiles = np.percentile(distances, [2.5, 25, 50, 75, 97.5], axis=0)
    outer = standard + quantiles[2, 0, :, None]*normal
    inner = standard - quantiles[2, 1, :, None]*normal
    # Close only for the visual section. The far-end connection is an artificial
    # comparison cut, not a measured surface or a complete vessel contour.
    contour = np.vstack([outer, inner[::-1][:-1]])
    write_polyline_ply(out / f"{stem}_xy.ply", np.column_stack([contour/unit_to_mm, np.zeros(len(contour))]), closed=True)
    for name, curve, color in (("outer", outer, (230, 70, 50)), ("inner", inner, (50, 100, 230)),
                               ("midline", standard, (40, 190, 70))):
        write_curves(out / f"{stem}_{name}_xy.ply", [curve], unit_to_mm, color)
    with (out / f"{stem}.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["point_id", "u", "n_profiles", "mid_x_input", "mid_y_input", "outer_x_input", "outer_y_input",
                         "inner_x_input", "inner_y_input", "mid_x_mm", "mid_y_mm", "outer_x_mm", "outer_y_mm",
                         "inner_x_mm", "inner_y_mm", "normal_x", "normal_y", "outer_ray_count", "inner_ray_count",
                         "measurement_vs_model_angle_p50_deg", "measurement_vs_model_angle_max_deg"] +
                        [f"{side}_distance_p{q}_mm" for side in ("outer", "inner") for q in ("2_5", "25", "50", "75", "97_5")])
        for j in range(len(standard)):
            writer.writerow([j, j/(len(standard)-1), len(chosen), *(standard[j]/unit_to_mm),
                             *(outer[j]/unit_to_mm), *(inner[j]/unit_to_mm), *standard[j], *outer[j], *inner[j],
                             *normal[j], int(measured[:, 0, j].sum()), int(measured[:, 1, j].sum()),
                             float(np.median(angular_difference[:,j])), float(angular_difference[:,j].max()),
                             *quantiles[:, 0, j], *quantiles[:, 1, j]])
    return {"profile_ids": chosen.tolist(), "count": len(chosen), "section_ply": f"{stem}_xy.ply",
            "distance_definition": "individual midline normal-ray intersections with its transformed original wall polylines; projection onto the same individual normal if intersection unavailable",
            "measurement_frame": "individual transformed midline normal",
            "construction_frame": "standard midline normal; distance quantiles measured before this placement",
            "measurement_vs_model_angle_p95_deg": float(np.percentile(angular_difference,95)),
            "measurement_vs_model_angle_max_deg": float(angular_difference.max()),
            "normal_ray_count": int(measured[:, :, 1:].sum()),
            "paired_projection_fallback_count": int((~measured[:, :, 1:]).sum()),
            "opposite_sign_fraction": float(np.mean(np.r_[signed_outer[:, 1:].ravel(), signed_inner[:, 1:].ravel()] < 0)),
            "quantiles": "empirical distribution quantiles; not confidence intervals",
            "cut_end": "artificial closure at comparison endpoint", "tip": "shared contour intersection; zero distance"}


def export_section_models(out, mode, aligned, mapped_pairs, records, unit_to_mm, write_curves,
                          *, geometries=None, outlier_mad=3.5, bimodal=True, bimodal_min_profiles=20, bimodal_bic_delta=10.):
    mids, pairs = np.asarray(aligned), np.asarray(mapped_pairs)
    measurement_normals = individual_normal_frames(mids, pairs)
    distances, measured, _, _ = measured_distances(mids, pairs, measurement_normals, geometries)
    widths = distances.sum(axis=1)
    floor = max(.05, .05*float(np.median(widths[:, 1:])))
    # Shape and wall separation both enter the profile-level screening.
    features = np.concatenate([mids.reshape(len(mids), -1), distances[:, :, 1:].reshape(len(mids), -1)], axis=1)
    if bimodal:
        labels, split = screen_two_modes(features, min_profiles=bimodal_min_profiles, bic_delta=bimodal_bic_delta)
    else:
        labels, split = None, {"status": "disabled"}
    scores = np.zeros(len(mids))
    groups = [np.arange(len(mids))] if labels is None else [np.flatnonzero(labels == k) for k in [0, 1]]
    filtering = []
    for group in groups:
        if len(group) < 6:
            filtering.append({"profile_ids": group.tolist(), "status": "insufficient_profiles_for_rejection"})
        else:
            scores[group] = robust_profile_scores(features[group], floor)
            filtering.append({"profile_ids": group.tolist(), "status": "evaluated"})
    inliers = np.flatnonzero(scores <= outlier_mad)
    if len(inliers) < 2:
        inliers = np.arange(len(mids))
        filter_status = "insufficient_inliers; retained all"
    else:
        filter_status = "ok"
    kept = set(inliers.tolist())
    report = {"outlier_threshold": outlier_mad, "distance_floor_mm": floor, "filter_status": filter_status,
              "measurement_frame": "individual transformed midline normal, independent of all/inliers/A/B selection",
              "filtering_groups": filtering, "bimodality": split, "models": {},
              "restoration_annotation": "unknown; all geometrically valid observed profiles eligible",
              "sampling": "equal weight per accepted radial side, no assumed independent observations",
              "profiles": [{"profile_id": i, "angle_deg": r["angle_deg"], "side": r["side"],
                            "outlier_score": float(scores[i]), "retained": i in kept,
                            "group": ("A" if labels[i] == 0 else "B") if labels is not None else None}
                           for i, r in enumerate(records)]}
    with (Path(out) / f"{mode}_section_distribution.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["profile_id", "angle_deg", "side", "point_id", "u", "outer_distance_mm",
                         "inner_distance_mm", "outer_method", "inner_method", "retained", "group",
                         "measurement_normal_x", "measurement_normal_y"])
        for i, rec in enumerate(report["profiles"]):
            for j in range(mids.shape[1]):
                methods = ["tip" if j == 0 else "normal_ray" if measured[i,k,j] else "paired_projection" for k in [0,1]]
                writer.writerow([i, rec["angle_deg"], rec["side"], j, j/(mids.shape[1]-1),
                                 *distances[i,:,j], *methods, rec["retained"], rec["group"], *measurement_normals[i,j]])
    report["distribution_csv"] = f"{mode}_section_distribution.csv"
    if "scores" in split:
        # Empirical first-PC frequency table accompanying the mixture screen.
        frequencies, edges = np.histogram(split["scores"], bins=max(8, int(np.ceil(np.sqrt(len(mids))))))
        with (Path(out) / f"{mode}_bimodal_histogram.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["pc1_lower", "pc1_upper", "count"])
            writer.writerows(zip(edges[:-1], edges[1:], frequencies))
        report["bimodality"]["histogram_csv"] = f"{mode}_bimodal_histogram.csv"
    selections = [("all", np.arange(len(mids))), ("inliers", inliers)]
    if labels is not None:
        # Group models retain all members. Outliers are reported separately;
        # within-group filtering protects a real minority mode from global rejection.
        selections += [("A", groups[0]), ("B", groups[1])]
    for label, indices in selections:
        stem = f"standard_{mode}_section_{label}"
        report["models"][label] = _write_model(Path(out), stem, mids, pairs, indices, unit_to_mm, write_curves, geometries)
    (Path(out) / f"{mode}_section_models.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def model_product_names(mode):
    names = [f"{mode}_section_models.json", f"{mode}_section_distribution.csv", f"{mode}_bimodal_histogram.csv"]
    for label in ("all", "inliers", "A", "B"):
        stem = f"standard_{mode}_section_{label}"
        names += [f"{stem}.csv", f"{stem}_xy.ply"]
        names += [f"{stem}_{name}_xy.ply" for name in ("outer", "inner", "midline")]
    return names

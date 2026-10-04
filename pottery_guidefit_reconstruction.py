#!/usr/bin/env python3
"""
PotteryGuideFitReconstruction
Version 0.1.0

Second-stage workflow for PotteryReconstruction3D.

This program explicitly separates four states:
  A. current sherd placement (immutable observation)
  B. sherd-derived vessel (from pottery_reconstruction_3d.py)
  C. guide-derived vessel (from aligned external outer/inner rings)
  D. individually fitted sherd placement against the guide-derived vessel

Key individual-fit constraints:
  * each sherd is assumed to preserve BOTH outer and inner surfaces;
  * outer and inner are fitted simultaneously with equal objective weight;
  * axial translation along the vessel axis is fixed to zero;
  * circumferential position about the vessel axis is fixed;
  * radial translation is optimized from the current position;
  * two small tilt degrees of freedom are optimized about the sherd centroid;
  * scale and sherd shape are never changed.

The current placement is always retained and exportable. Fitted geometry is a
separate diagnostic result and never overwrites the observed input placement.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable

import numpy as np

try:
    import pottery_reconstruction_3d as base
except ModuleNotFoundError as exc:
    raise RuntimeError(
        "pottery_guidefit_reconstruction.py must be in the same folder as "
        "pottery_reconstruction_3d.py"
    ) from exc

__version__ = "0.1.0"


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

@dataclass
class GuideDerivedVessel:
    profile: base.ProfileModel
    axis: base.GlobalAxisFit
    mesh_vertices_mm: np.ndarray
    mesh_faces: np.ndarray
    registration: base.GuideRegistration
    source_ring_count: int
    outer_ring_count: int
    inner_ring_count: int
    unknown_ring_count: int
    current_residual_mm: np.ndarray
    current_nearest_surface: np.ndarray
    current_fragment_quality: dict[int, dict]


@dataclass
class SurfaceSamples:
    outer_points_mm: np.ndarray
    inner_points_mm: np.ndarray
    outer_normals: np.ndarray
    inner_normals: np.ndarray
    outer_count_raw: int
    inner_count_raw: int


@dataclass
class SherdGuideFit:
    fragment_id: int
    status: str
    radial_shift_mm: float
    tilt_radial_deg: float
    tilt_tangential_deg: float
    axial_shift_mm: float
    circumferential_shift_deg: float
    outer_n: int
    inner_n: int
    outer_rmse_before_mm: float | None
    inner_rmse_before_mm: float | None
    combined_before_mm: float | None
    outer_rmse_after_mm: float | None
    inner_rmse_after_mm: float | None
    combined_after_mm: float | None
    outer_median_after_mm: float | None
    inner_median_after_mm: float | None
    matrix_mm: np.ndarray
    matrix_input_units: np.ndarray


@dataclass
class GuideFitSession:
    reconstruction: base.ReconstructionResult
    guide_rings: base.GuideRingSet
    guide_registration: base.GuideRegistration
    guide_vessel: GuideDerivedVessel
    sherd_fits: dict[int, SherdGuideFit]


# -----------------------------------------------------------------------------
# Math utilities
# -----------------------------------------------------------------------------

def _axis_angle_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    a = base.normalize(axis)
    x, y, z = a
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    C = 1.0 - c
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C],
    ], dtype=float)


def _robust_rmse(values: np.ndarray, trim_quantile: float = 0.90) -> float:
    x = np.abs(np.asarray(values, float))
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("inf")
    if len(x) >= 10 and 0.5 < trim_quantile < 1.0:
        hi = float(np.quantile(x, trim_quantile))
        x = x[x <= hi]
    if len(x) == 0:
        return float("inf")
    return float(np.sqrt(np.mean(x*x)))


def _median_abs(values: np.ndarray) -> float:
    x = np.abs(np.asarray(values, float))
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else float("nan")


def _unique_profile(z: np.ndarray, r: np.ndarray, tolerance_mm: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    """Merge rings with nearly identical axial positions by median radius."""
    z = np.asarray(z, float)
    r = np.asarray(r, float)
    good = np.isfinite(z) & np.isfinite(r)
    z, r = z[good], r[good]
    if len(z) == 0:
        return np.empty(0), np.empty(0)
    order = np.argsort(z)
    z, r = z[order], r[order]
    zz, rr = [], []
    group_z = [float(z[0])]
    group_r = [float(r[0])]
    for zi, ri in zip(z[1:], r[1:]):
        if abs(float(zi) - float(np.median(group_z))) <= tolerance_mm:
            group_z.append(float(zi)); group_r.append(float(ri))
        else:
            zz.append(float(np.median(group_z)))
            rr.append(float(np.median(group_r)))
            group_z, group_r = [float(zi)], [float(ri)]
    zz.append(float(np.median(group_z)))
    rr.append(float(np.median(group_r)))
    return np.asarray(zz), np.asarray(rr)


def _interp_no_extrap(zq: np.ndarray, z: np.ndarray, r: np.ndarray) -> np.ndarray:
    out = np.full(len(zq), np.nan, float)
    if len(z) < 2:
        return out
    sel = (zq >= z[0]) & (zq <= z[-1])
    out[sel] = np.interp(zq[sel], z, r)
    return out


# -----------------------------------------------------------------------------
# Guide-derived vessel construction
# -----------------------------------------------------------------------------

def _aligned_ring_z(ring: base.GuideRing, reg: base.GuideRegistration,
                    axis_point_mm: np.ndarray, axis_direction: np.ndarray) -> float:
    c = base.apply_homogeneous(ring.center_mm[None, :], reg.matrix_mm)[0]
    return float(np.dot(c - axis_point_mm, base.normalize(axis_direction)))


def _resolve_unknown_ring_kinds(records: list[dict], z_pair_tolerance_mm: float = 1.0) -> None:
    """Assign unknown rings conservatively when two radii occur at nearly same Z.

    Existing explicit outer/inner names are never changed. Unknown singletons remain
    unknown and are not used to build the guide-derived wall.
    """
    unknown_idx = [i for i, rec in enumerate(records) if rec["kind"] == "unknown"]
    used = set()
    for i in unknown_idx:
        if i in used:
            continue
        candidates = [j for j in unknown_idx if j != i and j not in used and
                      abs(records[j]["z_mm"] - records[i]["z_mm"]) <= z_pair_tolerance_mm]
        if not candidates:
            continue
        j = min(candidates, key=lambda k: abs(records[k]["z_mm"] - records[i]["z_mm"]))
        a, b = records[i], records[j]
        if abs(a["radius_mm"] - b["radius_mm"]) < 0.2:
            continue
        if a["radius_mm"] > b["radius_mm"]:
            a["kind"], b["kind"] = "outer", "inner"
        else:
            a["kind"], b["kind"] = "inner", "outer"
        used.add(i); used.add(j)


def build_guide_derived_vessel(result: base.ReconstructionResult,
                               guide: base.GuideRingSet,
                               reg: base.GuideRegistration,
                               z_step_mm: float = 1.0,
                               mesh_angle_step_deg: float = 10.0) -> GuideDerivedVessel:
    """Create a SECOND vessel model only from aligned external rings.

    No sherd-derived profile radii are used here. The sherd-derived global axis is
    retained because guide registration has explicitly aligned the ring-stack axis
    to it. There is no profile extrapolation outside the outer/inner guide ring
    ranges.
    """
    if z_step_mm <= 0:
        raise ValueError("z_step_mm must be > 0")
    records = []
    for ring in guide.rings:
        records.append({
            "ring": ring,
            "kind": ring.kind,
            "z_mm": _aligned_ring_z(ring, reg, result.axis.axis_point_mm, result.axis.axis_direction),
            "radius_mm": float(ring.radius_mm),
        })
    _resolve_unknown_ring_kinds(records)
    outer = [r for r in records if r["kind"] == "outer"]
    inner = [r for r in records if r["kind"] == "inner"]
    unknown = [r for r in records if r["kind"] == "unknown"]
    if len(outer) < 2 or len(inner) < 2:
        raise ValueError(
            f"Guide-derived vessel requires >=2 outer and >=2 inner rings; "
            f"got outer={len(outer)}, inner={len(inner)}, unknown={len(unknown)}"
        )
    zo, ro = _unique_profile(np.array([x["z_mm"] for x in outer]),
                             np.array([x["radius_mm"] for x in outer]))
    zi, ri = _unique_profile(np.array([x["z_mm"] for x in inner]),
                             np.array([x["radius_mm"] for x in inner]))
    if len(zo) < 2 or len(zi) < 2:
        raise ValueError("Guide rings collapse to fewer than 2 unique Z levels on outer or inner surface")

    zmin = float(min(zo[0], zi[0])); zmax = float(max(zo[-1], zi[-1]))
    zgrid = np.arange(math.floor(zmin / z_step_mm) * z_step_mm,
                      math.ceil(zmax / z_step_mm) * z_step_mm + z_step_mm * 0.5,
                      z_step_mm)
    outer_model = _interp_no_extrap(zgrid, zo, ro)
    inner_model = _interp_no_extrap(zgrid, zi, ri)
    thickness = outer_model - inner_model
    bad = np.isfinite(thickness) & (thickness <= 0)
    outer_model[bad] = np.nan
    inner_model[bad] = np.nan
    thickness = outer_model - inner_model

    # Observed fields mark only ring-supported Z nodes; model fields are linearly
    # interpolated between ring levels. No extrapolation is performed.
    outer_obs = np.full(len(zgrid), np.nan)
    inner_obs = np.full(len(zgrid), np.nan)
    support_o = np.zeros(len(zgrid), int)
    support_i = np.zeros(len(zgrid), int)
    for z0, r0 in zip(zo, ro):
        k = int(np.argmin(np.abs(zgrid - z0))); outer_obs[k] = r0; support_o[k] += 1
    for z0, r0 in zip(zi, ri):
        k = int(np.argmin(np.abs(zgrid - z0))); inner_obs[k] = r0; support_i[k] += 1

    prof = base.ProfileModel(
        z_mm=zgrid,
        inner_observed_mm=inner_obs,
        outer_observed_mm=outer_obs,
        inner_model_mm=inner_model,
        outer_model_mm=outer_model,
        thickness_model_mm=thickness,
        support_inner_points=support_i,
        support_outer_points=support_o,
        support_inner_fragments=np.zeros(len(zgrid), int),
        support_outer_fragments=np.zeros(len(zgrid), int),
        orientation_mapping="GUIDE-DERIVED: aligned external outer/inner rings; linear interpolation; no extrapolation",
    )
    mv, mf = base.build_reconstructed_mesh(prof, result.axis, angle_step_deg=mesh_angle_step_deg)
    residual, surf = base.evaluate_point_residuals(result.display_points_mm, result.axis, prof)
    q = base.fragment_fit_quality(result.display_points_mm, result.display_labels, residual)
    return GuideDerivedVessel(
        profile=prof,
        axis=result.axis,
        mesh_vertices_mm=mv,
        mesh_faces=mf,
        registration=reg,
        source_ring_count=len(records),
        outer_ring_count=len(outer),
        inner_ring_count=len(inner),
        unknown_ring_count=len(unknown),
        current_residual_mm=residual,
        current_nearest_surface=surf,
        current_fragment_quality=q,
    )


# -----------------------------------------------------------------------------
# Dual-surface sherd fitting
# -----------------------------------------------------------------------------

def _positive_normal_is_outer(mapping: str) -> bool:
    s = (mapping or "").lower()
    if s.startswith("positive radial normal = outer") or "positive radial normal = outer" in s:
        return True
    if s.startswith("negative radial normal = outer") or "negative radial normal = outer" in s:
        return False
    # Most outward-oriented triangle meshes use positive radial normal for outer.
    return True


def extract_fragment_outer_inner_samples(result: base.ReconstructionResult,
                                         fragment_id: int,
                                         compatibility_max: float = 0.35,
                                         radial_sign_min: float = 0.12,
                                         max_samples_per_surface: int = 5000,
                                         seed: int = 0) -> SurfaceSamples:
    sel = result.face_labels == int(fragment_id)
    p = result.face_centroids_mm[sel]
    n = result.face_normals[sel]
    if len(p) < 50:
        return SurfaceSamples(np.empty((0,3)), np.empty((0,3)), np.empty((0,3)), np.empty((0,3)), 0, 0)
    g = base.axis_coordinates(p, n, result.axis.axis_point_mm, result.axis.axis_direction)
    comp = g["compatibility"]
    sign = g["radial_normal_sign"]
    good = np.isfinite(comp) & np.isfinite(sign) & (comp <= compatibility_max) & (np.abs(sign) >= radial_sign_min)
    pos_outer = _positive_normal_is_outer(result.profile.orientation_mapping)
    outer_mask = good & ((sign > 0) if pos_outer else (sign < 0))
    inner_mask = good & ((sign < 0) if pos_outer else (sign > 0))
    po, no = p[outer_mask], n[outer_mask]
    pi, ni = p[inner_mask], n[inner_mask]
    no_raw, ni_raw = len(po), len(pi)
    if len(po) > max_samples_per_surface:
        idx = base.sample_indices(len(po), max_samples_per_surface, seed + fragment_id * 17 + 1)
        po, no = po[idx], no[idx]
    if len(pi) > max_samples_per_surface:
        idx = base.sample_indices(len(pi), max_samples_per_surface, seed + fragment_id * 17 + 2)
        pi, ni = pi[idx], ni[idx]
    return SurfaceSamples(po, pi, no, ni, no_raw, ni_raw)


def _fragment_frame(points_mm: np.ndarray, axis_point_mm: np.ndarray,
                    axis_direction: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centroid = np.mean(points_mm, axis=0)
    a = base.normalize(axis_direction)
    d = centroid - axis_point_mm
    d = d - np.dot(d, a) * a
    if np.linalg.norm(d) < 1e-8:
        er, et = base.orthonormal_perp_basis(a)
    else:
        er = base.normalize(d)
        et = base.normalize(np.cross(a, er))
    return centroid, er, et, a


def _fragment_transform(centroid: np.ndarray, er: np.ndarray, et: np.ndarray, a: np.ndarray,
                        radial_shift_mm: float, tilt_radial_deg: float,
                        tilt_tangential_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Rigid transform preserving centroid Z and circumferential angle.

    Tilt axes pass through the original sherd centroid. Translation is only along
    the centroid's current radial direction. Therefore centroid axial coordinate
    and circumferential direction are unchanged exactly.
    """
    Rr = _axis_angle_matrix(er, math.radians(float(tilt_radial_deg)))
    Rt = _axis_angle_matrix(et, math.radians(float(tilt_tangential_deg)))
    R = Rt @ Rr
    target_centroid = centroid + float(radial_shift_mm) * er
    t = target_centroid - R @ centroid
    M = np.eye(4, dtype=float)
    M[:3, :3] = R
    M[:3, 3] = t
    return R, M


def apply_mm_matrix(points_mm: np.ndarray, matrix_mm: np.ndarray) -> np.ndarray:
    return base.apply_homogeneous(points_mm, matrix_mm)


def _surface_residual(points_mm: np.ndarray, profile: base.ProfileModel,
                      axis: base.GlobalAxisFit, kind: str) -> np.ndarray:
    g = base.axis_coordinates(points_mm, None, axis.axis_point_mm, axis.axis_direction)
    model = profile.outer_model_mm if kind == "outer" else profile.inner_model_mm
    r = base.interp_profile(g["z"], profile.z_mm, model)
    return g["rho"] - r


def _dual_surface_metrics(po: np.ndarray, pi: np.ndarray, vessel: GuideDerivedVessel) -> dict:
    ro = _surface_residual(po, vessel.profile, vessel.axis, "outer") if len(po) else np.empty(0)
    ri = _surface_residual(pi, vessel.profile, vessel.axis, "inner") if len(pi) else np.empty(0)
    oro = _robust_rmse(ro); iri = _robust_rmse(ri)
    combined = 0.5 * oro + 0.5 * iri if np.isfinite(oro) and np.isfinite(iri) else float("inf")
    return {
        "outer_rmse": oro,
        "inner_rmse": iri,
        "combined": combined,
        "outer_median": _median_abs(ro),
        "inner_median": _median_abs(ri),
        "outer_valid": int(np.count_nonzero(np.isfinite(ro))),
        "inner_valid": int(np.count_nonzero(np.isfinite(ri))),
    }


def fit_fragment_to_guide_vessel(result: base.ReconstructionResult,
                                 vessel: GuideDerivedVessel,
                                 fragment_id: int,
                                 tilt_limit_deg: float = 5.0,
                                 radial_search_limit_mm: float = 60.0,
                                 radial_grid_step_mm: float = 1.0,
                                 min_surface_samples: int = 30,
                                 compatibility_max: float = 0.35,
                                 radial_sign_min: float = 0.12,
                                 max_samples_per_surface: int = 5000) -> SherdGuideFit:
    """Fit one sherd by simultaneous equal-weight outer+inner surface agreement."""
    samples = extract_fragment_outer_inner_samples(
        result, fragment_id,
        compatibility_max=compatibility_max,
        radial_sign_min=radial_sign_min,
        max_samples_per_surface=max_samples_per_surface,
    )
    po, pi = samples.outer_points_mm, samples.inner_points_mm
    if len(po) < min_surface_samples or len(pi) < min_surface_samples:
        I = np.eye(4)
        return SherdGuideFit(
            fragment_id, "insufficient_dual_surface", 0.0, 0.0, 0.0, 0.0, 0.0,
            len(po), len(pi), None, None, None, None, None, None, None, None,
            I, I.copy()
        )

    allp = np.vstack([po, pi])
    centroid, er, et, a = _fragment_frame(allp, result.axis.axis_point_mm, result.axis.axis_direction)
    before = _dual_surface_metrics(po, pi, vessel)

    def objective(x: np.ndarray) -> float:
        dr, tr, tt = map(float, x)
        _, M = _fragment_transform(centroid, er, et, a, dr, tr, tt)
        pox = apply_mm_matrix(po, M)
        pix = apply_mm_matrix(pi, M)
        m = _dual_surface_metrics(pox, pix, vessel)
        if not np.isfinite(m["combined"]):
            return 1e6
        # Equal surface weight is the defining condition. There is intentionally
        # no penalty anchoring radial position or tilt to the current pose.
        return float(m["combined"])

    # Current pose is an informative starting state, not a penalty term.
    grid = np.arange(-abs(radial_search_limit_mm), abs(radial_search_limit_mm) + radial_grid_step_mm * 0.5,
                     max(0.1, radial_grid_step_mm))
    vals = np.array([objective(np.array([g, 0.0, 0.0])) for g in grid], float)
    x0 = np.array([float(grid[int(np.argmin(vals))]), 0.0, 0.0])

    bounds = [(-abs(radial_search_limit_mm), abs(radial_search_limit_mm)),
              (-abs(tilt_limit_deg), abs(tilt_limit_deg)),
              (-abs(tilt_limit_deg), abs(tilt_limit_deg))]
    xbest = x0.copy()
    try:
        from scipy.optimize import minimize
        opt = minimize(objective, x0, method="Powell", bounds=bounds,
                       options={"xtol": 1e-4, "ftol": 1e-5, "maxiter": 250, "disp": False})
        if opt.success or np.isfinite(opt.fun):
            xbest = np.asarray(opt.x, float)
    except Exception:
        pass

    dr, tr, tt = map(float, xbest)
    _, M = _fragment_transform(centroid, er, et, a, dr, tr, tt)
    after = _dual_surface_metrics(apply_mm_matrix(po, M), apply_mm_matrix(pi, M), vessel)

    # Matrix expressed in the original input coordinate unit for direct use with
    # full-resolution source geometry / CloudCompare.
    s = base.unit_scale(result.input_unit)
    Min = np.eye(4)
    Min[:3, :3] = M[:3, :3]
    Min[:3, 3] = M[:3, 3] / s

    status = "ok" if np.isfinite(after["combined"]) else "no_guide_overlap"
    return SherdGuideFit(
        fragment_id=fragment_id,
        status=status,
        radial_shift_mm=dr,
        tilt_radial_deg=tr,
        tilt_tangential_deg=tt,
        axial_shift_mm=0.0,
        circumferential_shift_deg=0.0,
        outer_n=len(po), inner_n=len(pi),
        outer_rmse_before_mm=None if not np.isfinite(before["outer_rmse"]) else before["outer_rmse"],
        inner_rmse_before_mm=None if not np.isfinite(before["inner_rmse"]) else before["inner_rmse"],
        combined_before_mm=None if not np.isfinite(before["combined"]) else before["combined"],
        outer_rmse_after_mm=None if not np.isfinite(after["outer_rmse"]) else after["outer_rmse"],
        inner_rmse_after_mm=None if not np.isfinite(after["inner_rmse"]) else after["inner_rmse"],
        combined_after_mm=None if not np.isfinite(after["combined"]) else after["combined"],
        outer_median_after_mm=None if not np.isfinite(after["outer_median"]) else after["outer_median"],
        inner_median_after_mm=None if not np.isfinite(after["inner_median"]) else after["inner_median"],
        matrix_mm=M,
        matrix_input_units=Min,
    )


def fit_all_fragments(result: base.ReconstructionResult,
                      vessel: GuideDerivedVessel,
                      tilt_limit_deg: float = 5.0,
                      radial_search_limit_mm: float = 60.0,
                      radial_grid_step_mm: float = 1.0,
                      min_surface_samples: int = 30,
                      progress: Callable[[str], None] | None = None) -> dict[int, SherdGuideFit]:
    out = {}
    for cl in result.clusters:
        fit = fit_fragment_to_guide_vessel(
            result, vessel, cl.fragment_id,
            tilt_limit_deg=tilt_limit_deg,
            radial_search_limit_mm=radial_search_limit_mm,
            radial_grid_step_mm=radial_grid_step_mm,
            min_surface_samples=min_surface_samples,
        )
        out[cl.fragment_id] = fit
        if progress:
            if fit.status == "ok":
                progress(
                    f"fragment {cl.fragment_id}: combined {fit.combined_before_mm:.3f} -> "
                    f"{fit.combined_after_mm:.3f} mm / radial {fit.radial_shift_mm:+.3f} mm / "
                    f"tilt {fit.tilt_radial_deg:+.3f}, {fit.tilt_tangential_deg:+.3f} deg"
                )
            else:
                progress(f"fragment {cl.fragment_id}: {fit.status} (outer={fit.outer_n}, inner={fit.inner_n})")
    return out


def fitted_display_points(result: base.ReconstructionResult,
                          fits: dict[int, SherdGuideFit]) -> np.ndarray:
    p = result.display_points_mm.copy()
    for fid, fit in fits.items():
        if fit.status != "ok":
            continue
        sel = result.display_labels == fid
        p[sel] = apply_mm_matrix(p[sel], fit.matrix_mm)
    return p


# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

def _write_profile_csv(path: Path, vessel: GuideDerivedVessel) -> None:
    p = vessel.profile
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["z_mm","inner_ring_observed_mm","outer_ring_observed_mm",
                    "inner_model_mm","outer_model_mm","thickness_mm"])
        for row in zip(p.z_mm, p.inner_observed_mm, p.outer_observed_mm,
                       p.inner_model_mm, p.outer_model_mm, p.thickness_model_mm):
            w.writerow([base.jsonable(x) for x in row])


def export_guidefit_session(session: GuideFitSession, out: Path) -> None:
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    r = session.reconstruction
    v = session.guide_vessel
    s = base.export_scale(r.export_unit)

    # Explicitly separate the two reconstructed vessels.
    base.write_ply_mesh(out / "vessel_SHERD_DERIVED.ply",
                        r.mesh_vertices_mm * s, r.mesh_faces,
                        coordinate_unit=r.export_unit)
    base.write_ply_mesh(out / "vessel_GUIDE_DERIVED.ply",
                        v.mesh_vertices_mm * s, v.mesh_faces,
                        coordinate_unit=r.export_unit)

    _write_profile_csv(out / "profile_GUIDE_DERIVED.csv", v)

    # Current placement residuals against guide-derived vessel.
    rgb = np.zeros((len(r.display_points_mm), 3), np.uint8)
    for cl in r.clusters:
        rgb[r.display_labels == cl.fragment_id] = base.fragment_color(cl.fragment_id)
    base.write_ply_point_cloud(
        out / "current_placement_vs_GUIDE_DERIVED.ply",
        r.display_points_mm * s, rgb, coordinate_unit=r.export_unit,
        fragment_ids=r.display_labels,
        residual_mm=v.current_residual_mm,
        nearest_surface=v.current_nearest_surface,
    )

    if session.sherd_fits:
        fp = fitted_display_points(r, session.sherd_fits)
        fres, fsurf = base.evaluate_point_residuals(fp, v.axis, v.profile)
        base.write_ply_point_cloud(
            out / "fitted_placement_vs_GUIDE_DERIVED.ply",
            fp * s, rgb, coordinate_unit=r.export_unit,
            fragment_ids=r.display_labels,
            residual_mm=fres,
            nearest_surface=fsurf,
        )
        fields = [
            "fragment_id","status","radial_shift_mm","tilt_radial_deg","tilt_tangential_deg",
            "axial_shift_mm","circumferential_shift_deg","outer_n","inner_n",
            "outer_rmse_before_mm","inner_rmse_before_mm","combined_before_mm",
            "outer_rmse_after_mm","inner_rmse_after_mm","combined_after_mm",
            "outer_median_after_mm","inner_median_after_mm",
        ] + [f"mm_m{a}{b}" for a in range(4) for b in range(4)] + \
            [f"input_m{a}{b}" for a in range(4) for b in range(4)]
        with (out / "sherd_guidefit_transforms.csv").open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
            for fid in sorted(session.sherd_fits):
                q = session.sherd_fits[fid]
                row = {k: getattr(q, k) for k in fields if hasattr(q, k)}
                for a in range(4):
                    for b in range(4):
                        row[f"mm_m{a}{b}"] = float(q.matrix_mm[a,b])
                        row[f"input_m{a}{b}"] = float(q.matrix_input_units[a,b])
                w.writerow(row)

    meta = {
        "program": "PotteryGuideFitReconstruction",
        "version": __version__,
        "baseline_program": "pottery_reconstruction_3d.py",
        "baseline_version": getattr(base, "__version__", "unknown"),
        "policy": {
            "current_placement_preserved": True,
            "guide_derived_vessel_separate_from_sherd_derived_vessel": True,
            "individual_fit_outer_inner_simultaneous_equal_weight": True,
            "axial_translation_fixed_mm": 0.0,
            "circumferential_shift_fixed_deg": 0.0,
            "radial_translation_optimized": True,
            "tilt_is_small_rigid_rotation_about_sherd_centroid": True,
            "scale_change": False,
        },
        "guide_registration": base.jsonable(asdict(session.guide_registration)),
        "guide_vessel": {
            "source_ring_count": v.source_ring_count,
            "outer_ring_count": v.outer_ring_count,
            "inner_ring_count": v.inner_ring_count,
            "unknown_ring_count": v.unknown_ring_count,
        },
        "sherd_fits": {str(k): base.jsonable(asdict(x)) for k, x in session.sherd_fits.items()},
    }
    (out / "guidefit_metadata.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")


# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

def launch_gui(initial_input: Path | None = None):
    try:
        from PySide6 import QtCore, QtWidgets
        import pyvista as pv
        from pyvistaqt import QtInteractor
    except ModuleNotFoundError as exc:
        raise RuntimeError("GUI requires PySide6, pyvista and pyvistaqt") from exc

    class MainWindow(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f"PotteryGuideFitReconstruction v{__version__}")
            self.resize(1580, 960)
            self.source_path: Path | None = Path(initial_input) if initial_input else None
            self.result: base.ReconstructionResult | None = None
            self.guide: base.GuideRingSet | None = None
            self.reg: base.GuideRegistration | None = None
            self.guide_vessel: GuideDerivedVessel | None = None
            self.fits: dict[int, SherdGuideFit] = {}
            self.camera_saved = None

            central = QtWidgets.QWidget(); self.setCentralWidget(central)
            layout = QtWidgets.QHBoxLayout(central)
            controls = QtWidgets.QWidget(); controls.setMaximumWidth(460)
            form = QtWidgets.QVBoxLayout(controls); self.form = form

            self.open_btn = QtWidgets.QPushButton("1  破片群 PLY / OBJを開く")
            self.open_btn.clicked.connect(self.open_dialog); form.addWidget(self.open_btn)
            ur = QtWidgets.QHBoxLayout(); ur.addWidget(QtWidgets.QLabel("入力単位"))
            self.input_unit = QtWidgets.QComboBox(); self.input_unit.addItems(["m","mm","cm"]); ur.addWidget(self.input_unit)
            ur.addWidget(QtWidgets.QLabel("出力単位")); self.output_unit = QtWidgets.QComboBox(); self.output_unit.addItems(["m","mm","cm"]); ur.addWidget(self.output_unit)
            form.addLayout(ur)

            self.display_points = self._spin("表示点数", 20000, 800000, 220000, 20000)
            self.face_samples = self._spin("解析face sample", 20000, 1000000, 320000, 20000)
            self.cluster_voxel = self._dspin("破片分離 voxel [mm]", 0.2, 10.0, 1.5, 0.1, 2)
            self.min_frag = self._spin("最小sample点数", 30, 50000, 400, 50)

            self.analyze_btn = QtWidgets.QPushButton("2  固定配置解析（Sherd-derived）")
            self.analyze_btn.clicked.connect(self.analyze_baseline); form.addWidget(self.analyze_btn)

            form.addWidget(self._separator("External guide rings"))
            gr = QtWidgets.QHBoxLayout(); gr.addWidget(QtWidgets.QLabel("Guide PLY単位（未記録時）"))
            self.guide_unit = QtWidgets.QComboBox(); self.guide_unit.addItems(["m","mm","cm"]); gr.addWidget(self.guide_unit); form.addLayout(gr)
            self.load_guide_btn = QtWidgets.QPushButton("3  復元円 / 円弧 PLY群を読み込む")
            self.load_guide_btn.clicked.connect(self.load_guide); form.addWidget(self.load_guide_btn)
            self.align_btn = QtWidgets.QPushButton("4  Guide中心軸をアライン")
            self.align_btn.clicked.connect(self.align_guide); form.addWidget(self.align_btn)

            form.addWidget(self._separator("Guide-derived vessel"))
            self.guide_z_step = self._dspin("Guide器形 Z step [mm]", 0.2, 10.0, 1.0, 0.2, 2)
            self.guide_mesh_angle = self._dspin("Guide器形 角度step [deg]", 2.0, 30.0, 10.0, 1.0, 1)
            self.build_guide_vessel_btn = QtWidgets.QPushButton("5  Guide-derived器形を生成＋現在距離")
            self.build_guide_vessel_btn.clicked.connect(self.build_guide_vessel); form.addWidget(self.build_guide_vessel_btn)

            form.addWidget(self._separator("表示 ON / OFF"))
            self.chk_sherds = QtWidgets.QCheckBox("破片群"); self.chk_sherds.setChecked(True)
            self.chk_sherd_vessel = QtWidgets.QCheckBox("Sherd-derived vessel"); self.chk_sherd_vessel.setChecked(False)
            self.chk_guide_rings = QtWidgets.QCheckBox("Aligned guide rings"); self.chk_guide_rings.setChecked(True)
            self.chk_guide_vessel = QtWidgets.QCheckBox("Guide-derived vessel"); self.chk_guide_vessel.setChecked(True)
            self.chk_axis = QtWidgets.QCheckBox("rotation axis"); self.chk_axis.setChecked(True)
            self.chk_residual = QtWidgets.QCheckBox("Guide-derived residual color"); self.chk_residual.setChecked(True)
            for cb in [self.chk_sherds,self.chk_sherd_vessel,self.chk_guide_rings,self.chk_guide_vessel,self.chk_axis,self.chk_residual]:
                cb.stateChanged.connect(self.render_current); form.addWidget(cb)

            form.addWidget(self._separator("Individual sherd dual-surface fitting"))
            self.tilt_limit = self._dspin("tilt limit [deg]", 0.0, 20.0, 5.0, 0.5, 1)
            self.radial_limit = self._dspin("radial search limit [mm]", 5.0, 200.0, 60.0, 5.0, 1)
            self.min_surface = self._spin("outer/inner 最小sample", 10, 1000, 30, 10)
            self.fit_btn = QtWidgets.QPushButton("6  各破片をouter+inner同時fit")
            self.fit_btn.clicked.connect(self.fit_sherds); form.addWidget(self.fit_btn)
            self.show_mode = QtWidgets.QComboBox(); self.show_mode.addItems(["CURRENT placement", "FITTED placement"])
            self.show_mode.currentIndexChanged.connect(self.render_current); form.addWidget(self.show_mode)

            self.export_btn = QtWidgets.QPushButton("7  GuideFit結果を書き出す")
            self.export_btn.clicked.connect(self.export_dialog); form.addWidget(self.export_btn)

            self.table = QtWidgets.QTableWidget(0, 8)
            self.table.setHorizontalHeaderLabels(["ID","outer N","inner N","before","after","radial","tilt R","tilt T"])
            self.table.horizontalHeader().setStretchLastSection(True); self.table.setMinimumHeight(210); form.addWidget(self.table)
            self.status = QtWidgets.QPlainTextEdit(); self.status.setReadOnly(True); self.status.setMaximumBlockCount(1000); form.addWidget(self.status, 1)

            layout.addWidget(controls)
            self.plotter = QtInteractor(central); layout.addWidget(self.plotter.interactor, 1)
            self.plotter.set_background("white"); self.plotter.add_axes()
            if self.source_path:
                QtCore.QTimer.singleShot(150, self.analyze_baseline)

        def _separator(self, text):
            lab = QtWidgets.QLabel(text); f = lab.font(); f.setBold(True); lab.setFont(f); return lab
        def _spin(self,label,lo,hi,value,step=1):
            row=QtWidgets.QHBoxLayout(); row.addWidget(QtWidgets.QLabel(label)); w=QtWidgets.QSpinBox(); w.setRange(lo,hi); w.setValue(value); w.setSingleStep(step); row.addWidget(w); self.form.addLayout(row); return w
        def _dspin(self,label,lo,hi,value,step,decimals):
            row=QtWidgets.QHBoxLayout(); row.addWidget(QtWidgets.QLabel(label)); w=QtWidgets.QDoubleSpinBox(); w.setRange(lo,hi); w.setValue(value); w.setSingleStep(step); w.setDecimals(decimals); row.addWidget(w); self.form.addLayout(row); return w
        def _log(self, text):
            self.status.appendPlainText(str(text)); QtWidgets.QApplication.processEvents()

        def open_dialog(self):
            fn,_ = QtWidgets.QFileDialog.getOpenFileName(self,"土器破片3Dモデル","","Mesh (*.ply *.obj);;All files (*)")
            if fn:
                self.source_path=Path(fn); self.result=None; self.guide_vessel=None; self.fits={}; self._log(f"Input: {self.source_path}")

        def analyze_baseline(self):
            if self.source_path is None: self.open_dialog()
            if self.source_path is None: return
            try:
                self._log("--- fixed-placement baseline reconstruction ---")
                self.result = base.run_reconstruction(
                    self.source_path,
                    input_unit=self.input_unit.currentText(), export_unit_name=self.output_unit.currentText(),
                    display_points=self.display_points.value(), face_samples=self.face_samples.value(),
                    cluster_voxel_mm=self.cluster_voxel.value(), min_fragment_sample_points=self.min_frag.value(),
                    progress=self._log,
                )
                self.guide_vessel=None; self.fits={}; self.render_current(reset_camera=True)
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"解析エラー",str(exc)); self._log(f"ERROR: {exc}")

        def load_guide(self):
            start = str(self.source_path.parent) if self.source_path else str(Path.cwd())
            d = QtWidgets.QFileDialog.getExistingDirectory(self,"復元円 / 円弧 PLY群のフォルダを選択",start)
            if not d: return
            try:
                self.guide = base.load_guide_ring_folder(Path(d), fallback_unit=self.guide_unit.currentText())
                self.reg=None; self.guide_vessel=None; self.fits={}
                self._log(f"Guide rings loaded: {len(self.guide.rings)}")
                self.render_current()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"Guide読み込みエラー",str(exc)); self._log(f"GUIDE ERROR: {exc}")

        def align_guide(self):
            if self.result is None or self.guide is None:
                QtWidgets.QMessageBox.information(self,"未準備","先に固定配置解析とGuide ring読込を行ってください。")
                return
            try:
                self.reg = base.register_guide_axis_to_reconstruction(self.guide, self.result)
                self.guide_vessel=None; self.fits={}
                self._log("--- guide center-axis registration ---")
                self._log(f"radius RMSE={self.reg.radius_rmse_mm:.3f} mm / axial shift={self.reg.axial_shift_mm:.3f} mm")
                self.render_current()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"Guide alignmentエラー",str(exc)); self._log(f"ALIGN ERROR: {exc}")

        def build_guide_vessel(self):
            if self.result is None or self.guide is None or self.reg is None:
                QtWidgets.QMessageBox.information(self,"未準備","先にGuide中心軸アラインメントまで実行してください。")
                return
            try:
                self.guide_vessel = build_guide_derived_vessel(
                    self.result, self.guide, self.reg,
                    z_step_mm=self.guide_z_step.value(), mesh_angle_step_deg=self.guide_mesh_angle.value())
                self.fits={}
                v=self.guide_vessel
                self._log("--- Guide-derived vessel created ---")
                self._log(f"rings outer/inner/unknown = {v.outer_ring_count}/{v.inner_ring_count}/{v.unknown_ring_count}")
                ds=v.current_residual_mm[np.isfinite(v.current_residual_mm)]
                if len(ds): self._log(f"current residual median={np.median(ds):.3f} mm / p90={np.quantile(ds,0.9):.3f} mm")
                self._populate_table_current(); self.render_current()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"Guide-derived器形エラー",str(exc)); self._log(f"GUIDE VESSEL ERROR: {exc}")

        def fit_sherds(self):
            if self.result is None or self.guide_vessel is None:
                QtWidgets.QMessageBox.information(self,"未準備","先にGuide-derived器形を生成してください。")
                return
            try:
                self._log("--- individual sherd dual-surface fit ---")
                self._log("Constraint: outer+inner equal-weight simultaneous fit; ΔZ=0; Δtheta=0; radial free; small tilt allowed")
                self.fits = fit_all_fragments(
                    self.result, self.guide_vessel,
                    tilt_limit_deg=self.tilt_limit.value(),
                    radial_search_limit_mm=self.radial_limit.value(),
                    min_surface_samples=self.min_surface.value(), progress=self._log)
                self._populate_table_fits(); self.show_mode.setCurrentText("FITTED placement"); self.render_current()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"個別fitエラー",str(exc)); self._log(f"FIT ERROR: {exc}")

        def _populate_table_current(self):
            if self.result is None or self.guide_vessel is None: return
            self.table.setRowCount(len(self.result.clusters))
            for row, cl in enumerate(self.result.clusters):
                q=self.guide_vessel.current_fragment_quality.get(cl.fragment_id,{})
                vals=[cl.fragment_id,"","",q.get("rmse_mm",""),"","","",""]
                for c,v in enumerate(vals): self.table.setItem(row,c,QtWidgets.QTableWidgetItem(str(v)))

        def _populate_table_fits(self):
            if self.result is None: return
            self.table.setRowCount(len(self.result.clusters))
            for row,cl in enumerate(self.result.clusters):
                q=self.fits.get(cl.fragment_id)
                if q is None:
                    vals=[cl.fragment_id,"","","","","","",""]
                else:
                    vals=[q.fragment_id,q.outer_n,q.inner_n,
                          "" if q.combined_before_mm is None else f"{q.combined_before_mm:.3f}",
                          "" if q.combined_after_mm is None else f"{q.combined_after_mm:.3f}",
                          f"{q.radial_shift_mm:+.3f}",f"{q.tilt_radial_deg:+.3f}",f"{q.tilt_tangential_deg:+.3f}"]
                for c,v in enumerate(vals): self.table.setItem(row,c,QtWidgets.QTableWidgetItem(str(v)))

        def _add_polyline(self, pts, color, width=2):
            if len(pts)<2: return
            poly=pv.PolyData(np.asarray(pts,float)); poly.lines=np.hstack([[len(pts)],np.arange(len(pts),dtype=np.int64)])
            self.plotter.add_mesh(poly,color=color,line_width=width)

        def _add_rings(self, scale):
            if self.guide is None or not self.chk_guide_rings.isChecked(): return
            M=self.reg.matrix_mm if self.reg is not None else None
            for ring in self.guide.rings:
                p=base.apply_homogeneous(ring.ordered_points_mm,M) if M is not None else ring.ordered_points_mm
                if ring.arc_span_deg>=345: p=np.vstack([p,p[0]])
                color="deepskyblue" if ring.kind=="outer" else ("orange" if ring.kind=="inner" else "magenta")
                self._add_polyline(p*scale,color,2)

        def _add_mesh(self, vertices_mm, faces, scale, color, opacity):
            if len(vertices_mm)==0 or len(faces)==0: return
            f=np.hstack([np.full((len(faces),1),3,np.int64),faces.astype(np.int64)]).ravel()
            mesh=pv.PolyData(vertices_mm*scale,f)
            self.plotter.add_mesh(mesh,color=color,opacity=opacity,show_edges=True,edge_color="gray")

        def render_current(self, reset_camera=False):
            if self.result is None:
                return
            camera=self.plotter.camera_position
            self.plotter.clear(); sc=base.export_scale(self.output_unit.currentText())
            pmm=self.result.display_points_mm
            if self.show_mode.currentText()=="FITTED placement" and self.fits:
                pmm=fitted_display_points(self.result,self.fits)
            if self.chk_sherds.isChecked():
                poly=pv.PolyData(pmm*sc)
                if self.chk_residual.isChecked() and self.guide_vessel is not None:
                    res,_=base.evaluate_point_residuals(pmm,self.guide_vessel.axis,self.guide_vessel.profile)
                    poly["guide_residual_mm"]=np.nan_to_num(res,nan=-1.0)
                    finite=res[np.isfinite(res)]; hi=max(3.0,float(np.quantile(finite,0.90))) if len(finite) else 3.0
                    self.plotter.add_points(poly,scalars="guide_residual_mm",cmap="viridis",clim=[0,hi],point_size=3.0,render_points_as_spheres=False)
                else:
                    rgb=np.zeros((len(pmm),3),np.uint8)
                    for cl in self.result.clusters: rgb[self.result.display_labels==cl.fragment_id]=base.fragment_color(cl.fragment_id)
                    poly["rgb"]=rgb; self.plotter.add_points(poly,scalars="rgb",rgb=True,point_size=3.0,render_points_as_spheres=False)
            if self.chk_sherd_vessel.isChecked():
                self._add_mesh(self.result.mesh_vertices_mm,self.result.mesh_faces,sc,"lightgray",0.18)
            if self.chk_guide_vessel.isChecked() and self.guide_vessel is not None:
                self._add_mesh(self.guide_vessel.mesh_vertices_mm,self.guide_vessel.mesh_faces,sc,"wheat",0.28)
            self._add_rings(sc)
            if self.chk_axis.isChecked():
                z=self.guide_vessel.profile.z_mm if self.guide_vessel is not None else self.result.profile.z_mm
                if len(z):
                    a=self.result.axis.axis_direction;c=self.result.axis.axis_point_mm
                    zz=np.linspace(float(np.nanmin(z)-10),float(np.nanmax(z)+10),100)
                    self._add_polyline((c[None,:]+zz[:,None]*a[None,:])*sc,"black",3)
            self.plotter.add_axes()
            if reset_camera:
                self.plotter.reset_camera(); self.plotter.view_isometric()
            elif camera is not None:
                self.plotter.camera_position=camera
            self.plotter.render()

        def export_dialog(self):
            if self.result is None or self.guide is None or self.reg is None or self.guide_vessel is None:
                QtWidgets.QMessageBox.information(self,"未準備","Guide-derived器形生成まで実行してください。")
                return
            start=str(self.source_path.parent) if self.source_path else str(Path.cwd())
            d=QtWidgets.QFileDialog.getExistingDirectory(self,"GuideFit出力先",start)
            if not d:return
            out=Path(d)/f"{self.source_path.stem}_GuideFit_v010"
            try:
                session=GuideFitSession(self.result,self.guide,self.reg,self.guide_vessel,self.fits)
                export_guidefit_session(session,out)
                self._log(f"Exported: {out}")
                QtWidgets.QMessageBox.information(self,"書き出し完了",str(out))
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self,"書き出しエラー",str(exc)); self._log(f"EXPORT ERROR: {exc}")

    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    win=MainWindow(); win.show(); sys.exit(app.exec())


# -----------------------------------------------------------------------------
# Self test
# -----------------------------------------------------------------------------

def self_test():
    rng=np.random.default_rng(123)
    a=base.normalize(np.array([0.18,-0.42,0.889]))
    c=np.array([12.0,-9.0,5.0])
    e1,e2=base.orthonormal_perp_basis(a)
    zring=np.arange(-40.0,61.0,10.0)
    def ro(z): return 78.0 + 0.08*z + 3.0*np.sin(z/30.0)
    def ri(z): return ro(z)-5.0
    rings=[]
    th=np.linspace(0,2*np.pi,100,endpoint=False)
    for kind,fn in (("outer",ro),("inner",ri)):
        for i,z in enumerate(zring):
            cen=c+z*a; r=fn(z)
            pts=cen+r*np.cos(th)[:,None]*e1+r*np.sin(th)[:,None]*e2
            rings.append(base.GuideRing(Path(f"ring_{i:02d}_{kind}.ply"),"mm",kind,pts,pts,cen,a,r,0,0,360))
    guide=base.GuideRingSet(Path("."),rings,c,a,0.0)
    tax=base.GlobalAxisFit(c,a,c,a,0,0,0,1,0,0,len(rings)*100)
    # fake baseline profile only for registration target
    zgrid=np.arange(-45,66,1.0)
    prof=base.ProfileModel(zgrid,_interp_no_extrap(zgrid,zring,np.array([ri(z) for z in zring])),
                           _interp_no_extrap(zgrid,zring,np.array([ro(z) for z in zring])),
                           _interp_no_extrap(zgrid,zring,np.array([ri(z) for z in zring])),
                           _interp_no_extrap(zgrid,zring,np.array([ro(z) for z in zring])),
                           np.full(len(zgrid),5.0),np.ones(len(zgrid),int),np.ones(len(zgrid),int),
                           np.ones(len(zgrid),int),np.ones(len(zgrid),int),"positive radial normal = outer")
    fake=type("R",(),{})(); fake.axis=tax; fake.profile=prof
    reg=base.GuideRegistration("ok",1,np.eye(3),np.zeros(3),np.eye(4),0,0,0,0,0,0,len(rings),len(rings),100,-40,60,-40,60,0,0)
    # synthetic current sherd: a patch shifted radially +7mm and tilted 2deg
    zvals=np.linspace(-10,20,70); tv=np.linspace(-0.35,0.35,50)
    zz,tt=np.meshgrid(zvals,tv,indexing="ij"); zz=zz.ravel(); tt=tt.ravel()
    er=np.cos(tt)[:,None]*e1+np.sin(tt)[:,None]*e2
    po=c+zz[:,None]*a+(np.array([ro(x) for x in zz])+7.0)[:,None]*er
    pi=c+zz[:,None]*a+(np.array([ri(x) for x in zz])+7.0)[:,None]*er
    # Apply a small tilt about centroid to both surfaces to make optimization nontrivial.
    cent=np.mean(np.vstack([po,pi]),axis=0)
    Rtilt=_axis_angle_matrix(e1,math.radians(2.0)); t=cent-Rtilt@cent
    Mtilt=np.eye(4);Mtilt[:3,:3]=Rtilt;Mtilt[:3,3]=t
    po=apply_mm_matrix(po,Mtilt);pi=apply_mm_matrix(pi,Mtilt)
    # normals approximately radial with opposite inner winding
    no=er.copy();ni=-er.copy()
    P=np.vstack([po,pi]);N=np.vstack([no,ni]);L=np.ones(len(P),np.int32)
    # fake display points = face points for test
    recon=base.ReconstructionResult(Path("synthetic.ply"),"mm","mm",
        [base.FragmentCluster(1,len(P),np.mean(P,axis=0).tolist(),np.min(P,axis=0).tolist(),np.max(P,axis=0).tolist())],
        P.copy(),P.copy(),L.copy(),None,P.copy(),N.copy(),L.copy(),tax,prof,
        np.empty((0,3)),np.empty((0,3),np.int32),np.full(len(P),np.nan),np.zeros(len(P),np.uint8),{}, {})
    vessel=build_guide_derived_vessel(recon,guide,reg,z_step_mm=1.0,mesh_angle_step_deg=15)
    fit=fit_fragment_to_guide_vessel(recon,vessel,1,tilt_limit_deg=5,radial_search_limit_mm=20,min_surface_samples=20)
    if fit.status!="ok": raise SystemExit(f"SELF TEST FAILED status={fit.status}")
    if fit.combined_after_mm is None or fit.combined_before_mm is None or fit.combined_after_mm >= fit.combined_before_mm*0.35:
        raise SystemExit(f"SELF TEST FAILED improvement {fit.combined_before_mm}->{fit.combined_after_mm}")
    if abs(fit.radial_shift_mm+7.0)>1.5:
        raise SystemExit(f"SELF TEST FAILED radial shift {fit.radial_shift_mm}")
    if abs(fit.axial_shift_mm)>1e-12 or abs(fit.circumferential_shift_deg)>1e-12:
        raise SystemExit("SELF TEST FAILED fixed DOF")
    print("PotteryGuideFitReconstruction v0.1.0 SELF TEST PASSED")
    print(f"combined RMSE : {fit.combined_before_mm:.4f} -> {fit.combined_after_mm:.4f} mm")
    print(f"radial shift  : {fit.radial_shift_mm:+.4f} mm")
    print(f"tilt          : {fit.tilt_radial_deg:+.4f}, {fit.tilt_tangential_deg:+.4f} deg")
    print("axial shift   : 0.0000 mm (fixed)")
    print("theta shift   : 0.0000 deg (fixed)")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    ap=argparse.ArgumentParser(description="Guide-derived pottery vessel reconstruction and dual-surface sherd fitting")
    ap.add_argument("input",nargs="?",type=Path,help="Input sherd PLY/OBJ; omit to choose in GUI")
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--version",action="version",version=f"%(prog)s {__version__}")
    args=ap.parse_args()
    if args.self_test:
        self_test(); return
    launch_gui(args.input)


if __name__=="__main__":
    main()

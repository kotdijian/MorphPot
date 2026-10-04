#!/usr/bin/env python3
"""
PotteryVolumeCore
Version 0.6.0

Shared numerical capacity engine for profile/section based pottery volume methods.

This module does not acquire geometry. It does not know whether a profile came
from a 3-D mesh, a 2-D drawing, SVG, PDF rendering, or another measurement
source. Front-end programs are responsible for generating physical coordinates
in millimetres and for determining method-specific physical bounds such as the
spill level.

Public calculation functions return a dictionary with two keys:
    summary : scalar metadata and volume results
    profile : NumPy arrays used in the numerical integration

The voxel-fluid method is intentionally outside this module because it estimates
retained liquid from 3-D voxel-space connectivity rather than section/profile
integration.
"""
from __future__ import annotations

import argparse
import math
from typing import Iterable

import numpy as np

__version__ = "0.6.0"


def _as_1d(a: Iterable[float] | np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(a, dtype=np.float64).reshape(-1)
    if len(out) == 0:
        raise ValueError(f"{name} is empty.")
    return out


def _prepare_profile(z_mm, radius_mm, *, name: str) -> tuple[np.ndarray, np.ndarray]:
    z = _as_1d(z_mm, f"{name}.z_mm")
    r = _as_1d(radius_mm, f"{name}.radius_mm")
    if len(z) != len(r):
        raise ValueError(f"{name}: z and radius lengths differ.")
    good = np.isfinite(z) & np.isfinite(r)
    z, r = z[good], r[good]
    if len(z) < 2:
        raise ValueError(f"{name}: fewer than 2 finite profile points.")
    if np.any(r < -1e-9):
        raise ValueError(f"{name}: negative radii are not allowed.")
    r = np.maximum(r, 0.0)
    order = np.argsort(z)
    z, r = z[order], r[order]

    # Collapse duplicate/near-duplicate Z values. Median radius prevents one
    # repeated digitizing point from receiving extra numerical weight.
    out_z, out_r = [], []
    i = 0
    tol = 1e-9
    while i < len(z):
        j = i + 1
        while j < len(z) and abs(float(z[j] - z[i])) <= tol:
            j += 1
        out_z.append(float(np.mean(z[i:j])))
        out_r.append(float(np.median(r[i:j])))
        i = j
    z = np.asarray(out_z, dtype=np.float64)
    r = np.asarray(out_r, dtype=np.float64)
    if len(z) < 2 or float(np.ptp(z)) <= 0:
        raise ValueError(f"{name}: profile has no positive vertical span.")
    return z, r


def _uniform_grid(lower: float, upper: float, step_mm: float) -> np.ndarray:
    if step_mm <= 0:
        raise ValueError("z_step_mm must be positive.")
    if upper <= lower:
        raise ValueError("upper bound must be greater than lower bound.")
    internal = np.arange(lower, upper + 0.5 * step_mm, step_mm, dtype=np.float64)
    internal = internal[(internal > lower + 1e-12) & (internal < upper - 1e-12)]
    return np.r_[lower, internal, upper]


def _summary(method: str, volume_mm3: float | None, status: str = "ok", **extra) -> dict:
    row = {
        "method": method,
        "status": status,
        "core_version": __version__,
        "volume_mm3": None if volume_mm3 is None else float(volume_mm3),
        "volume_ml": None if volume_mm3 is None else float(volume_mm3 / 1000.0),
        "volume_l": None if volume_mm3 is None else float(volume_mm3 / 1_000_000.0),
    }
    row.update(extra)
    return row


def calculate_single_two_sided(
    z_left_mm,
    radius_left_mm,
    z_right_mm,
    radius_right_mm,
    *,
    z_step_mm: float = 0.5,
    bottom_z_mm: float | None = None,
    spill_z_mm: float | None = None,
    method: str = "single",
) -> dict:
    """
    Capacity from one full longitudinal section represented by opposing sides.

    The common cross-sectional area is
        A(z) = pi/2 * (r_left(z)^2 + r_right(z)^2)
    and the area-preserving equivalent radius is
        r_eq(z) = sqrt((r_left(z)^2 + r_right(z)^2)/2).

    The two profiles are linearly interpolated onto one common Z grid. If
    bottom_z_mm/spill_z_mm are explicitly supplied (e.g. from a 3-D cavity
    detector), endpoint values may be held by numpy.interp for at most the
    source front-end's intended boundary half-cell. The core records, but does
    not independently infer, the physical meaning of those bounds.
    """
    zl, rl = _prepare_profile(z_left_mm, radius_left_mm, name="left")
    zr, rr = _prepare_profile(z_right_mm, radius_right_mm, name="right")

    natural_lower = float(max(zl.min(), zr.min()))
    natural_upper = float(min(zl.max(), zr.max()))
    lower = natural_lower if bottom_z_mm is None else float(bottom_z_mm)
    upper = natural_upper if spill_z_mm is None else float(spill_z_mm)

    # Explicit bounds are allowed to extend by at most half a requested sample
    # step beyond either input profile. This supports cell-centred profiles from
    # the 3-D front end without permitting uncontrolled extrapolation.
    allowance = 0.500001 * float(z_step_mm)
    if lower < natural_lower - allowance:
        raise ValueError("bottom_z_mm lies too far below common profile coverage.")
    if upper > natural_upper + allowance:
        raise ValueError("spill_z_mm lies too far above common profile coverage.")
    if upper <= lower:
        raise ValueError("Left and right profiles have no common integration range.")

    z = _uniform_grid(lower, upper, float(z_step_mm))
    r_left = np.interp(z, zl, rl)
    r_right = np.interp(z, zr, rr)
    equivalent_radius = np.sqrt(0.5 * (r_left * r_left + r_right * r_right))
    area = 0.5 * math.pi * (r_left * r_left + r_right * r_right)
    area_left = math.pi * r_left * r_left
    area_right = math.pi * r_right * r_right

    volume_mm3 = float(np.trapezoid(area, z))
    left_volume_mm3 = float(np.trapezoid(area_left, z))
    right_volume_mm3 = float(np.trapezoid(area_right, z))
    summary = _summary(
        method,
        volume_mm3,
        bottom_z_mm=lower,
        spill_z_mm=upper,
        common_height_mm=upper - lower,
        z_step_mm=float(z_step_mm),
        integration_intervals=int(max(0, len(z) - 1)),
        integration="piecewise-linear radii resampled in Z + trapezoidal area integration",
        area_definition="A(z)=pi/2*(r_left(z)^2+r_right(z)^2)",
        equivalent_radius_definition="r_eq=sqrt((r_left^2+r_right^2)/2)",
        left_only_volume_l=float(left_volume_mm3 / 1_000_000.0),
        right_only_volume_l=float(right_volume_mm3 / 1_000_000.0),
        left_right_difference_l=float(abs(left_volume_mm3 - right_volume_mm3) / 1_000_000.0),
        left_right_difference_percent_of_mean=(
            float(abs(left_volume_mm3 - right_volume_mm3) / ((left_volume_mm3 + right_volume_mm3) * 0.5) * 100.0)
            if left_volume_mm3 + right_volume_mm3 > 0 else None
        ),
    )
    profile = {
        "z_mm": z,
        "left_radius_mm": r_left,
        "right_radius_mm": r_right,
        "equivalent_radius_mm": equivalent_radius,
        "area_mm2": area,
    }
    return {"summary": summary, "profile": profile}


def midpoint_volume(area_mm2, z_mm, z_step_mm: float,
                    lower_mm: float | None = None, upper_mm: float | None = None) -> tuple[float, int]:
    """Integrate cell-centred areas using clipped midpoint cells."""
    a = np.asarray(area_mm2, dtype=np.float64)
    z = np.asarray(z_mm, dtype=np.float64)
    if len(a) != len(z):
        raise ValueError("area_mm2 and z_mm lengths differ.")
    if z_step_mm <= 0:
        raise ValueError("z_step_mm must be positive.")
    total = 0.0
    used = 0
    for ai, zi in zip(a, z):
        if not np.isfinite(ai) or ai < 0:
            continue
        lo = float(zi - 0.5 * z_step_mm)
        hi = float(zi + 0.5 * z_step_mm)
        if lower_mm is not None:
            lo = max(lo, float(lower_mm))
        if upper_mm is not None:
            hi = min(hi, float(upper_mm))
        if hi <= lo:
            continue
        total += float(ai) * float(hi - lo)
        used += 1
    return float(total), used


def robust_radius(values, mad_k: float = 3.0) -> tuple[float, int, int]:
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan"), 0, 0
    if len(v) <= 2:
        return float(np.mean(v)), len(v), len(v)
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)))
    sigma = 1.4826 * mad
    if sigma <= 1e-12:
        use = np.ones(len(v), dtype=bool)
    else:
        use = np.abs(v - med) <= mad_k * sigma
        if np.count_nonzero(use) < max(2, len(v) // 2):
            use[:] = True
    return float(np.mean(v[use])), int(len(v)), int(np.count_nonzero(use))


def periodic_fill_radii(angles_deg, radii, min_valid: int = 3) -> np.ndarray | None:
    """Periodic linear interpolation of missing angular radii."""
    a = np.asarray(angles_deg, dtype=np.float64) % 360.0
    r = np.asarray(radii, dtype=np.float64)
    if len(a) != len(r):
        raise ValueError("angles_deg and radii lengths differ.")
    valid = np.isfinite(r)
    if int(np.count_nonzero(valid)) < min_valid:
        return None
    order = np.argsort(a)
    a_sorted, r_sorted, valid_sorted = a[order], r[order], valid[order]
    av, rv = a_sorted[valid_sorted], r_sorted[valid_sorted]
    ae = np.r_[av - 360.0, av, av + 360.0]
    re = np.r_[rv, rv, rv]
    filled_sorted = np.interp(a_sorted, ae, re)
    out = np.empty_like(filled_sorted)
    out[order] = filled_sorted
    return out


def calculate_optimized_axisymmetric(
    z_mm,
    radii_matrix_mm,
    *,
    z_step_mm: float,
    bottom_z_mm: float | None,
    spill_z_mm: float | None,
    min_valid_rays: int,
    mad_k: float = 3.0,
    method: str = "optimized",
) -> dict:
    z = _as_1d(z_mm, "z_mm")
    m = np.asarray(radii_matrix_mm, dtype=np.float64)
    if m.ndim != 2 or m.shape[1] != len(z):
        raise ValueError("radii_matrix_mm must have shape (directions, len(z_mm)).")
    r_opt = np.full(len(z), np.nan)
    n_valid = np.zeros(len(z), dtype=int)
    n_inlier = np.zeros(len(z), dtype=int)
    for k in range(len(z)):
        finite = m[:, k][np.isfinite(m[:, k])]
        if len(finite) < int(min_valid_rays):
            continue
        rr, nv, ni = robust_radius(finite, mad_k=mad_k)
        r_opt[k], n_valid[k], n_inlier[k] = rr, nv, ni
    area = math.pi * r_opt * r_opt
    if spill_z_mm is None or not np.any(np.isfinite(area)):
        return {"summary": _summary(method, None, status="failed", note="Insufficient inner radial profiles."), "profile": {}}
    vol, cells = midpoint_volume(area, z, z_step_mm, bottom_z_mm, spill_z_mm)
    summary = _summary(
        method, vol, bottom_z_mm=bottom_z_mm, spill_z_mm=spill_z_mm,
        z_step_mm=float(z_step_mm), valid_z_cells=int(cells), radial_directions=int(m.shape[0]),
        min_valid_rays=int(min_valid_rays), radius_outlier_mad_k=float(mad_k),
        note="Robust mean inner radius across directions; converted to an optimized axisymmetric body.",
    )
    return {"summary": summary, "profile": {
        "z_mm": z, "optimized_radius_mm": r_opt, "area_mm2": area,
        "valid_rays": n_valid, "radius_inliers": n_inlier,
    }}


def calculate_angular_integration(
    z_mm,
    angles_deg,
    radii_matrix_mm,
    *,
    z_step_mm: float,
    bottom_z_mm: float | None,
    spill_z_mm: float | None,
    min_valid_rays: int,
    method: str = "angular",
) -> dict:
    z = _as_1d(z_mm, "z_mm")
    angles = _as_1d(angles_deg, "angles_deg")
    m = np.asarray(radii_matrix_mm, dtype=np.float64)
    if m.ndim != 2 or m.shape != (len(angles), len(z)):
        raise ValueError("radii_matrix_mm must have shape (len(angles_deg), len(z_mm)).")
    if len(angles) < 3:
        raise ValueError("At least 3 angular directions are required.")
    area = np.full(len(z), np.nan)
    valid_counts = np.zeros(len(z), dtype=int)
    interpolated_counts = np.zeros(len(z), dtype=int)
    delta_theta = 2.0 * math.pi / len(angles)
    for k in range(len(z)):
        vals = m[:, k]
        nv = int(np.count_nonzero(np.isfinite(vals)))
        valid_counts[k] = nv
        if nv < int(min_valid_rays):
            continue
        filled = periodic_fill_radii(angles, vals, min_valid=3)
        if filled is None:
            continue
        interpolated_counts[k] = int(len(vals) - nv)
        area[k] = 0.5 * delta_theta * float(np.sum(filled * filled))
    if spill_z_mm is None or not np.any(np.isfinite(area)):
        return {"summary": _summary(method, None, status="failed", note="Insufficient angular inner profiles."), "profile": {}}
    vol, cells = midpoint_volume(area, z, z_step_mm, bottom_z_mm, spill_z_mm)
    summary = _summary(
        method, vol, bottom_z_mm=bottom_z_mm, spill_z_mm=spill_z_mm,
        z_step_mm=float(z_step_mm), valid_z_cells=int(cells), radial_directions=int(len(angles)),
        angular_step_deg=float(360.0 / len(angles)), min_valid_rays=int(min_valid_rays),
        note="Direct polar area integration A(z)=1/2 integral r(theta,z)^2 dtheta; no axial-symmetry assumption. Missing angles are periodically interpolated only when the minimum valid count is met.",
    )
    return {"summary": summary, "profile": {
        "z_mm": z, "area_mm2": area, "valid_rays": valid_counts,
        "interpolated_rays": interpolated_counts,
    }}


def calculate_ellipse_integration(
    section_z_mm,
    semi_major_mm,
    semi_minor_mm,
    *,
    z_grid_mm,
    z_step_mm: float,
    bottom_z_mm: float | None,
    spill_z_mm: float | None,
    section_ids=None,
    method: str = "ellipse",
) -> dict:
    zz = _as_1d(section_z_mm, "section_z_mm")
    a = _as_1d(semi_major_mm, "semi_major_mm")
    b = _as_1d(semi_minor_mm, "semi_minor_mm")
    if not (len(zz) == len(a) == len(b)):
        raise ValueError("ellipse sample lengths differ.")
    good = np.isfinite(zz) & np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    zz, a, b = zz[good], a[good], b[good]
    if len(zz) < 3:
        return {"summary": _summary(method, None, status="failed", note="Fewer than 3 valid inner horizontal ellipses."), "profile": {}}
    if spill_z_mm is None:
        return {"summary": _summary(method, None, status="failed", note="No spill level supplied."), "profile": {}}
    order = np.argsort(zz)
    zz, a, b = zz[order], a[order], b[order]
    area_samples = math.pi * a * b
    grid_all = _as_1d(z_grid_mm, "z_grid_mm")
    lo = float(bottom_z_mm if bottom_z_mm is not None else zz.min())
    hi = float(spill_z_mm)
    grid = grid_all[(grid_all >= lo) & (grid_all <= hi)]
    if len(grid) == 0:
        return {"summary": _summary(method, None, status="failed", note="No Z cells within ellipse integration bounds."), "profile": {}}
    area_grid = np.interp(grid, zz, area_samples)
    vol, cells = midpoint_volume(area_grid, grid, z_step_mm, lo, hi)
    summary = _summary(
        method, vol, bottom_z_mm=lo, spill_z_mm=hi, z_step_mm=float(z_step_mm),
        valid_z_cells=int(cells), horizontal_ellipse_sections=int(len(zz)),
        note="Horizontal inner-ellipse area pi*a*b interpolated along Z; intended as an independent QA method.",
    )
    ids = None
    if section_ids is not None:
        ids_arr = np.asarray(list(section_ids), dtype=object)
        if len(ids_arr) == len(good):
            ids = ids_arr[good][order]
    return {"summary": summary, "profile": {
        "section_z_mm": zz, "semi_major_mm": a, "semi_minor_mm": b,
        "ellipse_area_mm2": area_samples, "section_id": ids,
        "z_mm": grid, "area_mm2": area_grid,
    }}


def self_test() -> None:
    # exact cylinder: r=50 mm, h=100 mm
    z = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
    r = np.full_like(z, 50.0)
    one = calculate_single_two_sided(z, r, z, r, z_step_mm=0.5)
    expected = math.pi * 50.0**2 * 100.0 / 1_000_000.0
    got = one["summary"]["volume_l"]
    if abs(got - expected) > 1e-10:
        raise SystemExit(f"single self-test failed: expected {expected}, got {got}")

    # angular integration for the same cylinder at 12 directions
    zc = np.arange(0.25, 100.0, 0.5)
    m = np.full((12, len(zc)), 50.0)
    ang = calculate_angular_integration(
        zc, np.arange(0, 360, 30), m, z_step_mm=0.5,
        bottom_z_mm=0.0, spill_z_mm=100.0, min_valid_rays=9,
    )
    got_ang = ang["summary"]["volume_l"]
    if abs(got_ang - expected) > 1e-10:
        raise SystemExit(f"angular self-test failed: expected {expected}, got {got_ang}")
    print(f"PotteryVolumeCore v{__version__} SELF TEST PASSED")
    print(f"cylinder expected : {expected:.9f} L")
    print(f"single           : {got:.9f} L")
    print(f"angular          : {got_ang:.9f} L")


def main() -> None:
    parser = argparse.ArgumentParser(description="Shared numerical engine used by the pottery capacity front ends.")
    parser.add_argument("--self-test", action="store_true", help="Run numerical regression tests")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

import csv
import json

import numpy as np
import pytest

from morphpot.rim_standardization import (
    apply_transform, crop_rim, export_rim_standardization, fit_constrained_affine,
    fit_similarity, make_midline, unit_shape,
)
from morphpot.section_overlay import registered_branches


def rim_contour():
    # Everted rim -> neck -> expanding body. Include a continuous base and cavity.
    outer = np.array([[20, 40], [17, 35], [15, 30], [15, 27], [18, 22], [22, 10], [22, 0]])
    inner = outer.copy().astype(float)
    inner[:, 0] -= 2
    inner[-1, 1] = 2
    # Traverse positive outer wall, outer base, negative outer, negative inner,
    # inner floor, positive inner, then close across the positive lip.
    left_o, left_i = outer.copy(), inner.copy()
    left_o[:, 0] *= -1
    left_i[:, 0] *= -1
    p = np.vstack([outer, left_o[::-1], left_i, inner[::-1]])
    p = np.column_stack([p, np.zeros(len(p))])
    return np.stack([p, np.roll(p, -1, axis=0)], axis=1)


def source_segments(angle, scale=1):
    p = rim_contour()
    a = np.deg2rad(angle)
    s = np.zeros_like(p)
    s[..., 0] = p[..., 0] * np.cos(a) + 3
    s[..., 1] = p[..., 0] * np.sin(a) - 4
    s[..., 2] = p[..., 1]
    return s / scale


def test_both_sides_have_same_right_oriented_midline():
    outer, inner = registered_branches(rim_contour())
    a, ta = make_midline(outer, inner, "right", .25)
    b, tb = make_midline(outer, inner, "left", .25)
    np.testing.assert_allclose(a, b)
    np.testing.assert_allclose(ta, tb)
    assert (a[:, 0] > 0).all()
    np.testing.assert_allclose(a[0], [19, 40])
    assert ta[10] == pytest.approx(2)


def test_reflection_preserves_observed_left_right_asymmetry_before_registration():
    segments = rim_contour()
    segments[..., 0] = np.where(segments[..., 0] < 0, segments[..., 0] * 1.05, segments[..., 0])
    outer, inner = registered_branches(segments)
    right, _ = make_midline(outer, inner, "right", .25)
    left, _ = make_midline(outer, inner, "left", .25)
    np.testing.assert_allclose(left[:, 0], right[:, 0] * 1.05)
    np.testing.assert_allclose(left[:, 1], right[:, 1])


def test_crop_stops_after_neck_before_lower_body():
    outer, inner = registered_branches(rim_contour())
    line, thickness = make_midline(outer, inner, "right", .25)
    cropped, qa = crop_rim(line, thickness, spacing_mm=.25, smooth_mm=1., buffer_mm=3)
    assert qa["selection"] == "direction_transition"
    assert 26 < qa["change_point_mm"][1] < 33
    assert qa["end_arc_length_mm"] - qa["change_arc_length_mm"] == pytest.approx(3)
    assert cropped[-1, 1] > 20
    np.testing.assert_allclose(cropped[0], line[0])
    _, auto = crop_rim(line, thickness, spacing_mm=.25, smooth_mm=1.)
    assert auto["buffer_mm"] == pytest.approx(4)


def test_straight_or_initially_right_down_rims_require_manual_end():
    for x in [np.full(40, 10), np.linspace(10, 20, 40)]:
        p = np.column_stack([x, np.linspace(40, 0, 40)])
        with pytest.raises(ValueError, match="no sustained"):
            crop_rim(p, np.full(40, 2))
        selected, qa = crop_rim(p, np.full(40, 2), end_mm=10)
        assert qa["selection"] == "manual_arc_length"
        assert np.linalg.norm(selected[-1] - selected[0]) == pytest.approx(10)


def test_similarity_removes_only_rotation_scale_translation():
    p = np.column_stack([np.linspace(0, 2, 40), np.sin(np.linspace(0, 2, 40))])
    a = .3
    matrix = 1.3 * np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    target = apply_transform(p, matrix, np.array([4, 7]))
    fitted, t = fit_similarity(p, target)
    np.testing.assert_allclose(fitted, matrix, atol=1e-12)
    np.testing.assert_allclose(apply_transform(p, fitted, t), target, atol=1e-12)
    assert np.linalg.det(fitted) > 0


def test_affine_regularization_bounds_and_collinear_fallback():
    p = np.column_stack([np.linspace(0, 5, 80), np.sin(np.linspace(0, 5, 80))])
    target = apply_transform(p, np.array([[1.04, .03], [0, .97]]), np.array([4, 7]))
    sim, st = fit_similarity(p, target)
    a, t, qa = fit_constrained_affine(p, target, np.ones(len(p)), penalty=.0001)
    assert qa["status"] == "ok"
    assert abs(qa["log_directional_stretch"]) <= np.log1p(.1) + 1e-9
    assert abs(qa["shear"]) <= .1 + 1e-9
    assert np.linalg.det(a) > 0
    assert np.linalg.norm(apply_transform(p, a, t) - target) < np.linalg.norm(apply_transform(p, sim, st) - target)
    line = np.column_stack([np.zeros(20), np.arange(20)])
    _, _, qa = fit_constrained_affine(line, line, np.ones(len(line)))
    assert qa["status"] == "similarity_fallback"


def test_unit_shape_removes_centroid_and_size():
    p = np.array([[1, 4], [2, 1], [7, 0]], dtype=float)
    a, _, _ = unit_shape(p)
    b, _, _ = unit_shape(p * 3 + [10, -4])
    np.testing.assert_allclose(a, b)
    np.testing.assert_allclose(a.mean(axis=0), 0, atol=1e-12)
    assert np.linalg.norm(a) == pytest.approx(1)


@pytest.mark.parametrize("unit,scale", [("mm", 1), ("cm", 10), ("m", 1000)])
def test_export_records_sides_units_and_open_standard_curve(tmp_path, unit, scale):
    sections = [(a, source_segments(a, scale)) for a in [0, 30, 90]]
    qa = export_rim_standardization(tmp_path, sections, (3 / scale, -4 / scale), scale, unit,
                                   smooth_mm=1, buffer_mm=3, points=33)
    assert qa["status"] == "ok"
    assert qa["valid_profiles"] == 6
    out = tmp_path / "rim_standardization"
    lines = (out / "standard_affine_midline_xy.ply").read_text().splitlines()
    assert "element vertex 33" in lines
    assert "element edge 32" in lines  # open rim curve, no artificial closing edge
    with (out / "standard_similarity_midline.csv").open(encoding="utf-8-sig") as f:
        first = next(csv.DictReader(f))
    assert float(first["x_input"]) * scale == pytest.approx(float(first["x_mm"]))
    assert float(first["x_mm"]) == pytest.approx(19)
    transforms = json.loads((out / "transforms.json").read_text())
    assert len(transforms) == 12
    assert {r["side"] for r in transforms} == {"left", "right"}


def test_invalid_disabled_and_failed_reruns_clear_owned_outputs(tmp_path):
    args = (tmp_path, [(0, source_segments(0)), (90, source_segments(90))], (3, -4), 1, "mm")
    assert export_rim_standardization(*args, buffer_mm=3, smooth_mm=1)["status"] == "ok"
    out = tmp_path / "rim_standardization"
    assert (out / "standard_affine_midline_xy.ply").exists()
    qa = export_rim_standardization(*args, buffer_mm=1000, smooth_mm=1)
    assert qa["status"] == "insufficient_valid_profiles"
    assert not (out / "standard_affine_midline_xy.ply").exists()
    qa = export_rim_standardization(*args, enabled=False)
    assert qa["status"] == "disabled"
    assert not (out / "raw_midlines.csv").exists()
    with pytest.raises(ValueError):
        export_rim_standardization(*args, smooth_mm=-1)

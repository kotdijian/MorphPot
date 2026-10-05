import csv
import json

import numpy as np
import pytest

from morphpot.rim_standardization import (
    apply_transform, crop_rim, export_rim_standardization, fit_constrained_affine,
    fit_similarity, make_midline, unit_shape,
    horizontal_stop, sample_curve, arc_positions,
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
    # Pair separation on a sloping wall approximates normal distance, not radial 2mm.
    assert ta[10] == pytest.approx(2 / np.sqrt(1 + .6**2), abs=.15)


def test_reflection_preserves_observed_left_right_asymmetry_before_registration():
    segments = rim_contour()
    segments[..., 0] = np.where(segments[..., 0] < 0, segments[..., 0] * 1.05, segments[..., 0])
    outer, inner = registered_branches(segments)
    right, _ = make_midline(outer, inner, "right", .25)
    left, _ = make_midline(outer, inner, "left", .25)
    assert left[0, 0] == pytest.approx(right[0, 0] * 1.05)
    a = sample_curve(left, np.linspace(0, arc_positions(left)[-1], 33))
    b = sample_curve(right, np.linspace(0, arc_positions(right)[-1], 33))
    assert np.linalg.norm(a - b) > .5


def test_crop_stops_after_neck_before_lower_body():
    outer, inner = registered_branches(rim_contour())
    line, thickness = make_midline(outer, inner, "right", .25)
    cropped, qa = crop_rim(line, thickness, spacing_mm=.25, smooth_mm=1., buffer_mm=3)
    assert qa["selection"] == "horizontal_stall_or_reversal"
    assert 26 < qa["change_point_mm"][1] < 33
    assert qa["end_arc_length_mm"] - qa["change_arc_length_mm"] == pytest.approx(3)
    assert cropped[-1, 1] > 20
    np.testing.assert_allclose(cropped[0], line[0])
    _, auto = crop_rim(line, thickness, spacing_mm=.25, smooth_mm=1.)
    assert auto["buffer_mm"] == pytest.approx(2 * auto["local_paired_wall_separation_mm"])


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


def test_horizontal_stop_does_not_require_downward_motion():
    x = np.r_[np.linspace(20, 10, 41), np.full(40, 10)]
    for sign in [-1, 1]:
        y = sign * np.linspace(0, 20, len(x))
        p = np.column_stack([x, y])
        distance, point, progress = horizontal_stop(p, .25, 1, .1)
        assert point[0] == pytest.approx(10, abs=.3)
        assert distance > 10
        assert progress >= -.1


def test_short_horizontal_jitter_does_not_stop_before_sustained_stall():
    t = np.linspace(0, 40, 161)
    x = 30 - np.minimum(t, 25) + .02 * np.sin(t * 12)
    p = np.column_stack([x, t])
    _, point, _ = horizontal_stop(p, .25, 2, .1)
    assert point[0] == pytest.approx(5, abs=.3)


def test_nonmonotone_lower_body_does_not_change_local_rim():
    original = rim_contour()
    o1, i1 = registered_branches(original)
    # Insert an upward jog after the 10mm-height station on the right wall.
    # Keep outer bottom and inner floor unchanged.
    oi = np.flatnonzero((o1[:, 0] > 0) & np.isclose(o1[:, 1], 10))[0]
    ii = np.flatnonzero((i1[:, 0] > 0) & np.isclose(i1[:, 1], 10))[0]
    o2 = np.insert(o1, oi+1, o1[oi] + [0, 2, 0], axis=0)
    i2 = np.insert(i1, ii, i1[ii] + [0, 2, 0], axis=0)
    m1, _ = make_midline(o1, i1, "right", .25, buffer_mm=3)
    m2, _ = make_midline(o2, i2, "right", .25, buffer_mm=3)
    # Body modifications must not invalidate the region. Allow discretization
    # drift of stop detection due to the total branch sampling grid.
    a = sample_curve(m1, np.linspace(0, 14, 29))
    b = sample_curve(m2, np.linspace(0, 14, 29))
    np.testing.assert_allclose(a, b, atol=.2)


def _ply_vertices_edges(path):
    lines = path.read_text().splitlines()
    nv = int(next(s for s in lines if s.startswith('element vertex ')).split()[-1])
    ne = int(next(s for s in lines if s.startswith('element edge ')).split()[-1])
    start = lines.index('end_header') + 1
    return (np.array([[float(v) for v in s.split()] for s in lines[start:start+nv]]),
            np.array([[int(v) for v in s.split()] for s in lines[start+nv:start+nv+ne]]))


@pytest.mark.parametrize('scale', [1, 10, 1000])
def test_verification_pairs_are_midpoints_before_and_after_fitting(tmp_path, scale):
    qa = export_rim_standardization(tmp_path, [(a, source_segments(a, scale)) for a in [0, 90]],
        (3/scale, -4/scale), scale, 'test', smooth_mm=1, buffer_mm=3, points=33)
    assert qa['status'] == 'ok'
    out = tmp_path / 'rim_standardization'
    for prefix in ['raw', 'similarity', 'affine']:
        mid, _ = _ply_vertices_edges(out / f'{prefix}_midlines_xy.ply')
        pairs, edges = _ply_vertices_edges(out / f'{prefix}_pair_connectors_xy.ply')
        np.testing.assert_allclose((pairs[::2,:3]+pairs[1::2,:3])/2, mid[:,:3], atol=1e-12)
        np.testing.assert_array_equal(edges, np.arange(len(pairs)).reshape(-1,2))
        assert len(edges) == 4*33
        np.testing.assert_array_equal(mid[:,3:], np.tile([40,190,70], (len(mid),1)))
        for key in ['outer', 'inner']:
            source, source_edges = _ply_vertices_edges(out / f'{prefix}_source_{key}_xy.ply')
            assert len(source_edges) == len(source)-4
            assert np.all(source[:,2] == 0)
    with (out / 'raw_paired_points.csv').open(encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 4*33
    for row in rows:
        for axis in ['x','y']:
            assert float(row[f'mid_{axis}_mm']) == pytest.approx(
                (float(row[f'outer_{axis}_mm'])+float(row[f'inner_{axis}_mm']))/2)
    export_rim_standardization(tmp_path, [], (0,0), scale, 'test', enabled=False)
    assert not list(out.glob('*source*'))
    assert not list(out.glob('*pair*'))

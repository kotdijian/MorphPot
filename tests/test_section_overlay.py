import json

import numpy as np
import pytest

from morphpot.section_overlay import (
    export_section_overlay, project_section, registered_branches,
)


def contour(width=12, height=20):
    # Continuous full U-shaped material outline, with 2 mm walls and floor.
    p = np.array([[width, height, 0], [width, 0, 0], [-width, 0, 0],
                  [-width, height, 0], [-width + 2, height, 0],
                  [-width + 2, 2, 0], [width - 2, 2, 0],
                  [width - 2, height, 0]], dtype=float)
    return np.stack([p, np.roll(p, -1, axis=0)], axis=1)


def spatial_segments(profile, angle, center=(3, -4), scale=1):
    a = np.deg2rad(angle)
    s = np.empty_like(profile)
    s[..., 0] = profile[..., 0] * np.cos(a) + center[0]
    s[..., 1] = profile[..., 0] * np.sin(a) + center[1]
    s[..., 2] = profile[..., 1]
    return s / scale


def read_ply(path):
    lines = path.read_text().splitlines()
    vertices = int(next(x for x in lines if x.startswith("element vertex ")).split()[-1])
    edges = int(next(x for x in lines if x.startswith("element edge ")).split()[-1])
    start = lines.index("end_header") + 1
    p = np.array([[float(x) for x in line.split()] for line in lines[start:start + vertices]])
    e = np.array([[int(x) for x in line.split()] for line in lines[start + vertices:start + vertices + edges]])
    return p, e


def test_projection_preserves_signed_sides_height_and_native_unit():
    expected = contour()
    for angle in [0, 30, 90, 155]:
        result = project_section(spatial_segments(expected, angle), (3, -4), angle)
        np.testing.assert_allclose(result, expected, atol=1e-12)


def test_lips_and_outer_inner_paths_are_invariant_to_segment_order():
    s = contour()
    rng = np.random.default_rng(2)
    rng.shuffle(s)
    s[::2] = s[::2, ::-1]
    outer, inner = registered_branches(s)
    np.testing.assert_allclose(outer[0], [11, 20, 0])
    np.testing.assert_allclose(outer[-1], [-11, 20, 0])
    assert outer[:, 1].min() == 0
    assert inner[:, 1].min() == 2
    np.testing.assert_array_equal(outer[-1], inner[0])
    np.testing.assert_array_equal(inner[-1], outer[0])


@pytest.mark.parametrize("unit,scale", [("mm", 1), ("cm", 10), ("m", 1000)])
def test_median_resists_one_outlier_and_writes_connected_closed_polyline(tmp_path, unit, scale):
    sections = [(angle, spatial_segments(contour(width), angle, scale=scale))
                for angle, width in [(0, 12), (30, 12), (60, 30)]]
    qa = export_section_overlay(tmp_path, sections, (3 / scale, -4 / scale), scale, .5, unit)
    assert qa["status"] == "ok"
    p, edges = read_ply(tmp_path / "section_overlay/median_full_section_xy.ply")
    np.testing.assert_allclose(p[0] * scale, [11, 20, 0])
    assert (p[:, 2] == 0).all()
    assert p[:, 0].max() * scale == pytest.approx(12)
    assert p[:, 0].min() * scale == pytest.approx(-12)
    assert len(edges) == len(p)
    np.testing.assert_array_equal(edges[-1], [len(p) - 1, 0])
    assert np.array_equal(edges[:-1, 1], edges[:-1, 0] + 1)
    overlay, _ = read_ply(tmp_path / "section_overlay/all_full_sections_xy_edges.ply")
    assert overlay[:, 0].max() * scale == pytest.approx(30)
    assert qa["negative_lip_vertex"] > 0


def test_open_or_disconnected_contours_are_not_bridged():
    with pytest.raises(ValueError):
        registered_branches(contour()[:-1])
    with pytest.raises(ValueError):
        registered_branches(np.concatenate([contour(), contour() + [50, 0, 0]]))
    with pytest.raises(ValueError):
        registered_branches(np.concatenate([contour(), [[[100, 0, 0], [100, 5, 0]]]]))


def test_failed_and_disabled_reruns_remove_stale_median(tmp_path):
    sections = [(a, spatial_segments(contour(), a)) for a in [0, 90]]
    export_section_overlay(tmp_path, sections, (3, -4), 1, .5, "mm")
    median = tmp_path / "section_overlay/median_full_section_xy.ply"
    assert median.exists()
    failed = export_section_overlay(tmp_path, sections[:1], (3, -4), 1, .5, "mm")
    assert failed["status"] == "insufficient_valid_sections"
    assert not median.exists()
    disabled = export_section_overlay(tmp_path, sections, (3, -4), 1, .5, "mm", enabled=False)
    assert disabled["status"] == "disabled"
    assert not (median.parent / "all_full_sections_xy_edges.ply").exists()
    assert json.loads((median.parent / "overlay_qa.json").read_text())["enabled"] is False

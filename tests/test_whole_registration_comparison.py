import numpy as np
import unittest

from plot_whole_registration_comparison import replay, describe, measure


def test_recorded_field_preserves_guard_and_rejected_profiles():
    points = np.array([[0., 0.], [2., 10.], [5., 20.], [6., 30.]])
    record = dict(status='applied', knots_z_mm=[10., 20., 30.],
                  delta_rz_mm=[[0., 0.], [1., 2.], [3., 4.]])
    transformed = replay(points, record)
    np.testing.assert_array_equal(transformed[:2], points[:2])
    np.testing.assert_allclose(transformed[2:], points[2:] + [[1., 2.], [3., 4.]])
    np.testing.assert_array_equal(replay(points, {'status': 'rejected'}), points)
    with unittest.TestCase().assertRaises(ValueError):
        replay(points, {'status': 'missing'})


def test_paired_stations_use_same_model_and_keep_rejected_sources():
    outer = np.array([[0., 0.], [2., 10.], [5., 20.], [6., 30.]])
    inner = np.array([[6., 30.], [4., 20.], [1., 10.], [0., 1.]])
    landmarks = dict(base_outer_mm=outer[0], base_inner_mm=inner[-1],
                     lip_mm=outer[-1], outer={'arc_mm': 20.}, inner={'arc_mm': 20.})
    profile = dict(angle=0., outer=outer, inner=inner, record={'status': 'rejected'},
                   landmarks=landmarks)
    rows = measure([profile], {'outer': outer, 'inner': inner}, 1.)
    assert len(rows) > 30
    for row in rows:
        assert abs(row['original_distance_mm']) < 1e-10
        assert row['registered_distance_mm'] == row['original_distance_mm']
        assert row['registered_z_mm'] == row['raw_z_mm']
    assert np.isclose(describe([0., 2.])['std_mm'], 1.)

if __name__ == "__main__":
    test_recorded_field_preserves_guard_and_rejected_profiles()
    test_paired_stations_use_same_model_and_keep_rejected_sources()
    print("2 analytical checks passed")

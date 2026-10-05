import csv
import json

import numpy as np
import pytest

from morphpot.rim_tip import extension_tip
from morphpot.rim_models import export_section_models, screen_two_modes
from morphpot.rim_standardization import _write_curves, export_rim_standardization


def rounded_tip_geometry():
    # Highest point (20,42) is not the true end (22,40).
    outer = np.array([[20,42], [0,42]], dtype=float)
    theta = np.linspace(np.pi/2, -np.pi/2, 101)
    inner = np.vstack([np.column_stack([20+2*np.cos(theta), 40+2*np.sin(theta)]), [0,38]])
    a = np.vstack([[20,42], [20,42], np.column_stack([np.arange(19., -1, -1), np.full(20,42.)])])
    b = np.vstack([[20,42], [22,40], np.column_stack([np.arange(19., -1, -1), np.full(20,38.)])])
    return outer, inner, a, b


@pytest.mark.parametrize('angle', [0., .45, -.8])
def test_tip_follows_midline_extension_to_original_contour(angle):
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    offset = np.array([13., -7.])
    data = [p @ rotation.T + offset for p in rounded_tip_geometry()]
    a, b, qa = extension_tip(*data, .25, 1.)
    expected = np.array([22.,40.]) @ rotation.T + offset
    np.testing.assert_allclose(a[0], expected, atol=1e-10)
    np.testing.assert_allclose(b[0], expected, atol=1e-10)
    assert np.linalg.norm(a[0]-qa['old_highest_lip_mm']) > 2
    midpoint = (a+b)/2
    original_frame = (midpoint-offset) @ rotation
    np.testing.assert_allclose(original_frame[:,1], 40, atol=1e-10)
    assert np.all(np.diff(original_frame[:,0]) < 0)


def test_tip_does_not_silently_fall_back_to_highest_point():
    outer, inner, a, b = rounded_tip_geometry()
    with pytest.raises(ValueError, match='no local lip-contour intersection'):
        extension_tip(outer, np.array([[20,38],[0,38]]), a, b, .25, 1.)


def model_inputs(widths):
    t = np.linspace(0,1,33)
    mid = np.column_stack([20-10*t, 40-4*t*t])
    tangent = np.gradient(mid,axis=0)
    normal = np.column_stack([-tangent[:,1], tangent[:,0]])
    normal /= np.linalg.norm(normal,axis=1)[:,None]
    normal *= -1
    mids, pairs = [], []
    for width in widths:
        distance = np.full(len(t), width/2)
        distance[0] = 0
        mids.append(mid.copy())
        pairs.append([mid+distance[:,None]*normal, mid-distance[:,None]*normal])
    records = [dict(angle_deg=i*5, side='right' if i%2==0 else 'left') for i in range(len(widths))]
    return mids, pairs, records


@pytest.mark.parametrize('scale', [1,10,1000])
def test_standard_section_normals_distances_units_and_closed_cut(tmp_path, scale):
    mids,pairs,records = model_inputs([4]*8)
    report = export_section_models(tmp_path,'similarity',mids,pairs,records,scale,_write_curves)
    assert report['models']['all']['count'] == 8
    assert report['bimodality']['status'] == 'insufficient_profiles'
    rows = list(csv.DictReader((tmp_path/'standard_similarity_section_all.csv').open(encoding='utf-8-sig')))
    for row in rows[1:]:
        mid = np.array([float(row[f'mid_{k}_mm']) for k in ['x','y']])
        outer = np.array([float(row[f'outer_{k}_mm']) for k in ['x','y']])
        inner = np.array([float(row[f'inner_{k}_mm']) for k in ['x','y']])
        assert np.linalg.norm(outer-mid) == pytest.approx(2)
        assert np.linalg.norm(inner-mid) == pytest.approx(2)
        np.testing.assert_allclose((outer+inner)/2,mid)
        assert float(row['outer_x_input'])*scale == pytest.approx(outer[0])
    assert float(rows[0]['outer_distance_p50_mm']) == 0
    ply = (tmp_path/'standard_similarity_section_all_xy.ply').read_text().splitlines()
    assert 'element vertex 65' in ply
    assert 'element edge 65' in ply
    assert report['models']['all']['opposite_sign_fraction'] == 0


def test_profile_outlier_rejection_preserves_all_model_and_audit(tmp_path):
    mids,pairs,records = model_inputs([4]*8+[20])
    report = export_section_models(tmp_path,'affine',mids,pairs,records,1,_write_curves,bimodal=False)
    assert report['models']['all']['count'] == 9
    assert report['models']['inliers']['profile_ids'] == list(range(8))
    assert report['profiles'][-1]['retained'] is False
    assert report['profiles'][-1]['outlier_score'] > 3.5


def test_round_lip_distances_use_original_contour_not_linear_tip_bridge(tmp_path):
    outer, inner, _, _ = rounded_tip_geometry()
    xs = np.linspace(22,10,25)
    mid = np.column_stack([xs,np.full(len(xs),40.)])
    artificial = np.full(len(xs),.5)
    artificial[0] = 0
    pair = [mid+np.column_stack([np.zeros(len(xs)),artificial]),
            mid-np.column_stack([np.zeros(len(xs)),artificial])]
    records = [dict(angle_deg=0,side='right'),dict(angle_deg=0,side='left')]
    report = export_section_models(tmp_path,'similarity',[mid,mid],[pair,pair],records,1,_write_curves,
        geometries=[dict(outer=outer,inner=inner)]*2)
    rows = list(csv.DictReader((tmp_path/'standard_similarity_section_all.csv').open(encoding='utf-8-sig')))
    for row in rows[1:]:
        x = float(row['mid_x_mm'])
        expected = np.sqrt(4-(x-20)**2) if x > 20 else 2
        assert float(row['outer_distance_p50_mm']) == pytest.approx(expected, abs=.001)
        assert float(row['inner_distance_p50_mm']) == pytest.approx(expected, abs=.001)
        assert int(row['outer_ray_count']) == 2
    assert report['models']['all']['paired_projection_fallback_count'] == 0


def test_normal_intersections_follow_anisotropic_transformed_surface(tmp_path):
    outer, inner, _, _ = rounded_tip_geometry()
    xs = np.linspace(22,10,25)
    mid = np.column_stack([xs,np.full(len(xs),40.)])
    offset = np.column_stack([np.zeros(len(xs)),np.full(len(xs),2.)])
    offset[0] = 0
    matrix = np.array([[1.05,.08],[0,.95]])
    translation = np.array([-13.,7.])
    transform = lambda p: p @ matrix.T+translation
    pair = [transform(mid+offset),transform(mid-offset)]
    geometry = dict(outer=transform(outer),inner=transform(inner))
    records = [dict(angle_deg=0,side='right'),dict(angle_deg=90,side='left')]
    export_section_models(tmp_path,'affine',[transform(mid)]*2,[pair]*2,records,1,_write_curves,geometries=[geometry]*2)
    rows = list(csv.DictReader((tmp_path/'standard_affine_section_all.csv').open(encoding='utf-8-sig')))
    expected = 2*np.linalg.det(matrix)/np.linalg.norm(matrix @ np.array([1.,0.]))
    for row in rows[6:]:
        assert float(row['outer_distance_p50_mm']) == pytest.approx(expected)
        assert float(row['inner_distance_p50_mm']) == pytest.approx(expected)


def test_two_supported_modes_get_a_b_models_without_discarding_minor_group(tmp_path):
    rng = np.random.default_rng(42)
    widths = np.r_[4+rng.normal(0,.1,24), 8+rng.normal(0,.1,12)]
    mids,pairs,records = model_inputs(widths)
    report = export_section_models(tmp_path,'similarity',mids,pairs,records,1,_write_curves)
    assert report['bimodality']['status'] == 'two_modes'
    assert sorted([report['models'][k]['count'] for k in ['A','B']]) == [12,24]
    assert len(report['models']['inliers']['profile_ids']) >= 32
    assert (tmp_path/'standard_similarity_section_A_xy.ply').exists()
    assert (tmp_path/'standard_similarity_section_B_xy.ply').exists()
    assert {r['group'] for r in report['profiles']} == {'A','B'}


def test_continuous_distribution_and_single_outlier_do_not_force_two_models():
    rng = np.random.default_rng(27)
    for scores in [rng.normal(0,1,80), np.linspace(-1,1,80), np.r_[np.zeros(39),100]]:
        labels, report = screen_two_modes(scores[:,None])
        assert labels is None
        assert report['status'] != 'two_modes'


@pytest.mark.parametrize('scale', [1,10,1000])
def test_individual_normals_prevent_oblique_chord_thickening(tmp_path, scale):
    from morphpot.rim_models import measured_distances, normal_frame
    mids, pairs, geometry = [], [], []
    for angle in [-30,30]:
        theta=np.deg2rad(angle)
        rotation=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        mid=np.column_stack([np.linspace(-10,10,41),np.zeros(41)]) @ rotation.T
        normal=np.array([0,1]) @ rotation.T
        mids.append(mid); pairs.append([mid+2*normal,mid-2*normal])
        geometry.append({k:np.array([[-50,v],[50,v]]) @ rotation.T for k,v in [('outer',2),('inner',-2)]})
    mids,pairs=np.asarray(mids),np.asarray(pairs)
    common=normal_frame(np.median(mids,axis=0),pairs[:,0]-mids)
    old,*_=measured_distances(mids,pairs,common,geometry)
    assert old[:,:,1:].sum(axis=1).mean() == pytest.approx(4/np.cos(np.deg2rad(30)))
    records=[dict(angle_deg=0,side='right'),dict(angle_deg=90,side='left')]
    report=export_section_models(tmp_path,'similarity',mids,pairs,records,scale,_write_curves,geometries=geometry)
    rows=list(csv.DictReader((tmp_path/'standard_similarity_section_all.csv').open(encoding='utf-8-sig')))
    for row in rows[1:]:
        assert float(row['outer_distance_p50_mm']) == pytest.approx(2)
        assert float(row['inner_distance_p50_mm']) == pytest.approx(2)
        assert float(row['outer_y_input'])*scale == pytest.approx(2)
        assert float(row['measurement_vs_model_angle_p50_deg']) == pytest.approx(30)
    assert report['models']['all']['normal_ray_count']==160
    assert report['models']['all']['paired_projection_fallback_count']==0


def test_model_cleanup_and_invalid_config(tmp_path):
    from test_rim_standardization import source_segments
    args = (tmp_path, [(0,source_segments(0)),(90,source_segments(90))], (3,-4),1,'mm')
    qa = export_rim_standardization(*args, smooth_mm=1,buffer_mm=3)
    assert qa['status'] == 'ok'
    out = tmp_path/'rim_standardization'
    # Simulate a previous successful A/B run; removal is limited to owned names.
    (out/'standard_affine'/'standard_affine_section_A_xy.ply').write_text('stale')
    (out/'user_notes.txt').write_text('keep')
    export_rim_standardization(*args,enabled=False)
    assert not list(out.glob('*section*'))
    assert (out/'user_notes.txt').read_text() == 'keep'
    for kwargs in [dict(outlier_mad=-1),dict(bimodal_min_profiles=2),dict(bimodal_bic_delta=float('nan'))]:
        with pytest.raises(ValueError):
            export_rim_standardization(*args,**kwargs)

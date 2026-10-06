import json
import numpy as np
import pytest
from pottery_whole_dimension_validation import dimensions, thickness, run


def test_slanted_normal_thickness():
    # Parallel walls x-z = 10 and x-z = 6: normal thickness 4/sqrt(2).
    o=np.array([[10.,0.],[30.,20.]])
    i=np.array([[26.,20.],[6.,0.]])
    rows=thickness(o,i,interval=1,margin=5, tangent_window=2)
    vals=[r['thickness_mm'] for r in rows if r['status']=='ok']
    assert len(vals)>10
    assert vals==pytest.approx(np.repeat(4/np.sqrt(2),len(vals)))
    assert rows[0]['status']=='endpoint_margin'


def test_dimensions_plateau():
    o=np.array([[0.,-5.],[20.,-5.],[20.,15.],[15.,25.]])
    i=np.array([[15.,25.],[16.,15.],[0.,0.]])
    d=dimensions(o,i)
    assert d['lip_diameter_mm']==30
    assert d['maximum_diameter_mm']==40
    assert d['maximum_diameter_height_mm']==10
    assert d['vessel_height_mm']==30


def test_missing_thickness_not_filled():
    rows=thickness(np.array([[10.,0.],[10.,20.]]),np.array([[20.,20.],[20.,0.]]))
    assert all(r['thickness_mm'] is None for r in rows)
    assert any(r['status']=='no_inner_intersection' for r in rows)


def test_pipeline(tmp_path):
    root=tmp_path/'model';root.mkdir()
    o=np.array([[0.,0.],[10.,0.],[10.,20.]])
    i=np.array([[10.,20.],[8.,18.],[8.,2.],[0.,2.]])
    samples_o=np.repeat(o[None],12,axis=0);samples_i=np.repeat(i[None],12,axis=0)
    np.savez(root/'profile_distribution.npz',outer_profiles_mm=samples_o,inner_profiles_mm=samples_i,outer_median_mm=o,inner_median_mm=i,accepted_angles_deg=np.arange(12)*30)
    (root/'whole_model.json').write_text(json.dumps({'version':'0.1.0-dev'}))
    out=tmp_path/'result';report=run(root,out,plots=False)
    assert report['datasets']['training']==12
    assert (out/'thickness_comparison.csv').exists()
    assert (out/'dimension_summary.csv').exists()
    with pytest.raises(ValueError,match='empty'):run(root,out,plots=False)


def test_grazing_intersection_flag():
    # Horizontal outer normal hits an almost horizontal inner segment.
    o=np.array([[10.,0.],[10.,20.]])
    i=np.array([[9.,0.],[0.,1.]])
    rows=thickness(o,i,interval=.1,margin=0)
    bad=[r for r in rows if r['status']=='grazing_inner_intersection']
    assert bad
    assert all(r['thickness_mm'] is None and r['raw_normal_distance_mm']>0 for r in bad)

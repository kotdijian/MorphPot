import json
import numpy as np
import pytest
from pottery_whole_validation import nearest_segments,aggregate,holdout_angles,summary,validate


def test_exact_segment_distance_and_surface_branch():
    curve=np.array([[0,0],[10,0],[10,10]],float)
    d,q=nearest_segments(np.array([[5,3],[13,5],[-2,-2]],float),curve)
    assert d==pytest.approx([3,3,np.sqrt(8)])
    assert q==pytest.approx(np.array([[5,0],[10,5],[0,0]],float))


def test_disjoint_holdout_grid_rejects_training_overlap():
    angles=holdout_angles(5,0,2.5,5,0)
    assert len(angles)==36
    assert np.allclose(angles%5,2.5)
    with pytest.raises(ValueError,match='overlaps'):holdout_angles(5,0,5,5,0)
    with pytest.raises(ValueError,match='overlaps'):holdout_angles(2.5,0,2.5,5,0)


def test_aggregation_extreme_values_and_area_weighting():
    x=np.array([0]*8+[10,100],float)[:,None,None]*np.ones((1,2,2))
    assert np.all(aggregate(x,'median',.1)==0)
    assert np.all(aggregate(x,'mean',.1)==11)
    assert np.all(aggregate(x,'trimmed_mean',.1)==1.25)
    s=summary(np.array([1.,10.]),np.array([9.,1.]))
    assert s['mean_mm']==5.5
    assert s['area_weighted_mean_mm']==pytest.approx(1.9)


def test_saved_output_only_standalone_validation(tmp_path):
    root=tmp_path/'source';root.mkdir()
    o=np.array([[0,0],[10,0],[10,20]],float)
    i=np.array([[10,20],[8,20],[8,2],[0,2]],float)
    np.savez(root/'profile_distribution.npz',outer_profiles_mm=np.tile(o,(12,1,1)),inner_profiles_mm=np.tile(i,(12,1,1)),outer_median_mm=o,inner_median_mm=i,accepted_angles_deg=np.arange(12)*30)
    (root/'whole_model.json').write_text(json.dumps(dict(version='0.1.0-dev',accepted_half_sections=12,axis_xy_mm=[0,0],input_unit='mm',unit_to_mm=1)))
    out=tmp_path/'validation'
    result=validate(root,out,interval=1,plots=False)
    assert result['datasets']=={'training':12}
    import csv
    rows=list(csv.DictReader((out/'distance_summary.csv').open()))
    assert len(rows)==9
    assert max(float(r['max_mm']) for r in rows)<1e-10
    assert (out/'training/median/profiles_all_overlay_xy.ply').exists()
    with pytest.raises(ValueError,match='empty'):validate(root,out,plots=False)

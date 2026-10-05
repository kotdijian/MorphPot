import json

import numpy as np
import pytest

from morphpot.rim_endpoint import detect_candidates, shape_summary
from morphpot.rim_standardization import export_rim_standardization


def opening_contour(angle=0., scale=1.):
    outer = np.array([[30,30],[26,22],[22,14],[17,4],[14,0]],dtype=float)
    inner = np.array([[28,30],[24,22],[20,14],[15,4],[12,2]],dtype=float)
    left_o,left_i=outer.copy(),inner.copy()
    left_o[:,0]*=-1
    left_i[:,0]*=-1
    p=np.vstack([outer,left_o[::-1],left_i,inner[::-1]])
    a=np.deg2rad(angle)
    xyz=np.column_stack([p[:,0]*np.cos(a),p[:,0]*np.sin(a),p[:,1]])/scale
    return np.stack([xyz,np.roll(xyz,-1,axis=0)],axis=1)


def test_oblique_to_horizontal_needs_no_radial_reversal_or_large_angle():
    p=np.array([[30,30],[25,20],[20,10],[15,0],[0,0]],dtype=float)
    d=detect_candidates(p,.25,1.)
    assert [c['kind'] for c in d['candidates']] == ['horizontal_transition']
    c=d['candidates'][0]
    assert c['inclination_deg'] <= 10
    assert c['tangent_radial'] < 0
    assert c['point_mm'][0] > 10


def test_initial_flat_lip_does_not_end_tracking_before_oblique_wall():
    p=np.array([[30,30],[27,30],[20,15],[14,0],[0,0]],dtype=float)
    d=detect_candidates(p,.25,1.)
    assert d['candidates'][0]['kind']=='horizontal_transition'
    assert d['candidates'][0]['point_mm'][1] < 2


def test_near_vertical_and_no_clear_transition_are_neutral_candidates():
    vertical=np.column_stack([np.full(60,20),np.linspace(30,0,60)])
    straight=np.column_stack([np.linspace(30,15,60),np.linspace(30,0,60)])
    assert 'near_vertical' in detect_candidates(vertical)['categories']
    assert detect_candidates(straight)['categories']==['no_clear_transition']


@pytest.mark.parametrize('unit,scale',[('mm',1),('cm',10),('m',1000)])
def test_auto_detects_horizontal_transition_and_writes_separate_outputs(tmp_path,unit,scale):
    sections=[(a,opening_contour(a,scale)) for a in [0,45,90]]
    qa=export_rim_standardization(tmp_path,sections,(0,0),scale,unit,smooth_mm=1,buffer_mm=2,points=33)
    assert qa['status']=='ok'
    assert qa['valid_profiles']==6
    summary=qa['shape_candidates']
    assert summary['effective_end_mode']=='horizontal_transition'
    assert summary['categories'][0]['support_fraction']==1
    assert all(r['endpoint_kind']=='horizontal_transition' for r in qa['profiles'])
    out=tmp_path/'rim_standardization'
    assert (out/'shape_candidates.json').exists()
    assert (out/'endpoint_candidates_xy.ply').exists()
    for mode in ['similarity','affine']:
        sub=out/f'standard_{mode}'
        assert (sub/f'standard_{mode}_section_all_xy.ply').exists()
        assert not (out/f'standard_{mode}_midline_xy.ply').exists()
        transforms=json.loads((sub/'transforms.json').read_text())
        assert len(transforms)==6
        assert {r['mode'] for r in transforms}=={mode}


def test_explicit_wrong_policy_is_reported_and_manual_override_works(tmp_path):
    args=(tmp_path,[(0,opening_contour())],(0,0),1,'mm')
    qa=export_rim_standardization(*args,end_mode='radial_turn',smooth_mm=1,buffer_mm=2)
    assert qa['status']=='insufficient_valid_profiles'
    assert qa['shape_candidates']['recommended_end_mode']=='horizontal_transition'
    assert all(r['status']=='excluded' for r in qa['profiles'])
    qa=export_rim_standardization(*args,end_mode='radial_turn',end_mm=15,smooth_mm=1)
    assert qa['status']=='ok'
    assert qa['shape_candidates']['effective_end_mode']=='manual'
    with pytest.raises(ValueError,match='requires'):
        export_rim_standardization(*args,end_mode='manual')


def test_auto_one_policy_avoids_mixing_endpoint_meanings(tmp_path):
    from test_rim_standardization import source_segments
    opening=[(a,opening_contour(a)) for a in [0,30]]
    turning=[(90,source_segments(90))]
    # Convert the existing shifted fixture into the same axis frame.
    turning[0][1][...,0]-=3
    turning[0][1][...,1]+=4
    qa=export_rim_standardization(tmp_path,opening+turning,(0,0),1,'mm',smooth_mm=1,buffer_mm=2)
    assert qa['shape_candidates']['effective_end_mode']=='horizontal_transition'
    valid=[r for r in qa['profiles'] if r['status']=='ok']
    assert valid and {r['endpoint_kind'] for r in valid}=={'horizontal_transition'}
    assert qa['shape_candidates']['evaluated_half_profiles']==6


def test_rerun_removes_legacy_and_nested_generated_files_only(tmp_path):
    out=tmp_path/'rim_standardization'
    sub=out/'standard_affine'
    sub.mkdir(parents=True)
    for p in [out/'standard_affine_midline_xy.ply',sub/'standard_affine_midline_xy.ply',
              sub/'standard_affine_section_B_xy.ply',sub/'transforms.json']:
        p.write_text('stale')
    (sub/'user_notes.txt').write_text('keep')
    export_rim_standardization(tmp_path,[],(0,0),1,'mm',enabled=False)
    assert (sub/'user_notes.txt').read_text()=='keep'
    assert list(sub.iterdir())==[sub/'user_notes.txt']
    assert not (out/'standard_affine_midline_xy.ply').exists()


def test_summary_support_is_per_half_profile_not_probability():
    description=detect_candidates(np.array([[30,30],[15,0],[0,0]]),.25,1.)
    profiles=[dict(angle_deg=0,side='right',shape_detection=description),
              dict(angle_deg=0,side='left',shape_detection=description)]
    report=shape_summary(profiles)
    assert report['categories'][0]['support_count']==2
    assert report['categories'][0]['support_fraction']==1

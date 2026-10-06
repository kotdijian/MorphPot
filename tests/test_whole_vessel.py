import numpy as np
import pytest
from morphpot.whole_vessel import clip_positive,split_walls,aggregate_profiles,revolve,closed_sections,run


def cup_profile():
    return np.array([[0,0],[20,0],[20,40],[17,40],[17,3],[0,3],[0,0]],float)


def test_clipping_preserves_outer_inner_base_and_lip():
    full=np.array([[-20,40],[-20,0],[20,0],[20,40],[17,40],[17,3],[-17,3],[-17,40],[-20,40]],float)
    half=clip_positive(full)
    outer,inner=split_walls(half)
    assert np.allclose(outer[0],[0,0])
    assert np.allclose(inner[-1],[0,3])
    assert np.array_equal(outer[-1],inner[0])
    assert outer[:,0].max()==20


def test_known_axisymmetric_shape_preserves_cavity_floor_and_volume():
    walls=[split_walls(cup_profile())]*6
    poly,_=aggregate_profiles(walls,.1)
    mesh=revolve(poly,180,[3,-4])
    assert mesh.is_watertight and mesh.is_winding_consistent
    expected=np.pi*(20**2*40-17**2*37)
    assert mesh.volume==pytest.approx(expected,rel=.005)
    assert mesh.extents[2]==pytest.approx(40,abs=.1)
    assert len(closed_sections(mesh,[3,-4],37))==1
    outer,inner=split_walls(clip_positive(closed_sections(mesh,[3,-4],37)[0]))
    assert outer[0,1]==pytest.approx(0,abs=.1)
    assert inner[-1,1]==pytest.approx(3,abs=.1)


def test_median_resists_one_extreme_profile_and_retains_native_size():
    base=split_walls(cup_profile())
    outlier=tuple(c*np.array([1.5,1]) for c in base)
    p,samples=aggregate_profiles([base,base,outlier],.2)
    assert p[:,0].max()==pytest.approx(20)
    assert p[:,1].max()==pytest.approx(40)
    assert samples[0].shape[0]==3


def test_native_unit_equivalence_and_closed_mesh_requirement(tmp_path):
    import trimesh
    model=revolve(cup_profile(),72,[3,-4])
    stats=[]
    for unit,scale in [('mm',1),('m',1000)]:
        path=tmp_path/f'{unit}.ply';mesh=model.copy();mesh.vertices/=scale;mesh.export(path)
        stats.append(run(path,tmp_path/f'out_{unit}',unit,angle_step=30,pitch=.1,
                         axis_xy_mm=[3,-4],plots=False)['model'])
    assert stats[0]['height_mm']==pytest.approx(stats[1]['height_mm'],abs=.02)
    assert stats[0]['material_volume_mm3']==pytest.approx(stats[1]['material_volume_mm3'],rel=.002)
    broken=trimesh.Trimesh(model.vertices,model.faces[:-1],process=False)
    path=tmp_path/'broken.ply';broken.export(path)
    with pytest.raises(ValueError,match='watertight'):
        run(path,tmp_path/'bad','mm',angle_step=30,axis_xy_mm=[3,-4],plots=False,require_watertight=True)

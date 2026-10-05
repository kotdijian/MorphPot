import numpy as np
import pytest
from morphpot.rim_phase_validation import phase_indices, rebuild_subset, run_phase_experiment
from morphpot.rim_validation import load_run
from morphpot.rim_standardization import export_rim_standardization


def test_five_phases_partition_one_degree_data_once():
    records=[dict(angle_deg=a,side=side) for a in range(180) for side in ("right","left")]
    sets=[phase_indices(records,5,p) for p in range(5)]
    assert all(len(s)==72 for s in sets)
    assert sorted(i for s in sets for i in s)==list(range(360))


@pytest.mark.parametrize("step",[10,15])
def test_all_interval_phases_cover_all_dense_profiles(step):
    records=[dict(angle_deg=a,side=side) for a in range(180) for side in ("right","left")]
    sets=[phase_indices(records,step,p) for p in range(step)]
    assert all(len(s)==360//step for s in sets)
    assert sorted(i for s in sets for i in s)==list(range(360))


def test_missing_direction_not_fabricated_and_shifted_origin():
    records=[dict(angle_deg=a+.25,side="right") for a in range(180) if a!=10]
    ids=phase_indices(records,5,0,.25)
    assert len(ids)==35
    assert all(np.isclose((records[i]['angle_deg']-.25)%5,0) for i in ids)


def test_subset_refits_reference_instead_of_slicing_transformed_output(tmp_path):
    # Two complete section planes with different right-hand rim radii.
    sections=[]
    for angle,scale in [(0,1),(5,1.1)]:
        loop=np.array([[30,30],[26,22],[10,2],[-10,2],[-26,22],[-30,30],[-28,30],[-24,22],[-8,4],[8,4],[24,22],[28,30],[30,30]],float)*scale
        theta=np.deg2rad(angle)
        xyz=np.c_[loop[:,0]*np.cos(theta),loop[:,0]*np.sin(theta),loop[:,1]]
        sections.append((angle,np.stack([xyz[:-1],xyz[1:]],axis=1)))
    export_rim_standardization(tmp_path/'source',sections,(0,0),1,'mm',end_mode='horizontal_transition',bimodal=False)
    source=load_run(tmp_path/'source','similarity')
    with pytest.raises(ValueError,match="180 attempted planes"):
        run_phase_experiment(tmp_path/'source',tmp_path/'invalid',plots=False)
    ids=[i for i,r in enumerate(source['report']['profiles']) if r['angle_deg']==5]
    rebuild_subset(source,ids,tmp_path/'subset')
    subset=load_run(tmp_path/'subset','similarity')
    assert len(subset['mids'])==2
    right=next(i for i,r in enumerate(subset['transforms']) if r['side']=='right')
    assert np.allclose(subset['transforms'][right]['matrix_2x2'],np.eye(2),atol=1e-8)
    assert np.allclose(subset['transforms'][right]['translation_mm'],0,atol=1e-8)
    assert subset['qa']['source_profile_ids']==ids

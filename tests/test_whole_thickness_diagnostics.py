import numpy as np
import pytest
from pottery_whole_thickness_diagnostics import ray_hit,point_arc,measure


def test_window_rejects_remote_hit():
    curve=np.array([[5.,0.],[5.,10.]])
    hit=ray_hit(np.array([10.,8.]),np.array([-1.,0.]),curve)
    assert hit['distance']==pytest.approx(5)
    assert hit['arc']==pytest.approx(8)
    assert ray_hit(np.array([10.,8.]),np.array([-1.,0.]),curve,0,2) is None
    assert ray_hit(np.array([10.,8.]),np.array([-1.,0.]),curve,8,2)['incidence']==pytest.approx(1)


def test_window_selects_supported_second_hit():
    curve=np.array([[8.,0.],[8.,10.],[4.,10.],[4.,0.]])
    hit=ray_hit(np.array([10.,5.]),np.array([-1.,0.]),curve,19,1)
    assert hit['distance']==pytest.approx(6)
    assert hit['arc']==pytest.approx(19)


def test_point_arc():
    assert point_arc(np.array([4.,3.]),np.array([[0.,0.],[0.,10.]]))==pytest.approx(3)


def test_midline_failure_and_measurement_no_fallback():
    outer=np.array([[0.,0.],[10.,0.],[10.,20.]])
    inner=np.array([[10.,20.],[8.,18.],[8.,2.],[0.,2.]])
    rows,curves,qa=measure(outer,inner,rim_length=40)
    assert qa[0]['status']=='failed'
    assert not any(r['method']=='rim_midline' for r in rows)
    assert curves
    local=[r for r in rows if r['method']=='outer_constrained' and r['status']=='local_normal']
    assert local
    assert all(r['arc_mismatch_mm']<=8 for r in local)

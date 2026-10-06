import numpy as np
import pytest
from pottery_rim_curvature_plot import curve_metrics


def test_straight_line_has_zero_curvature_and_turning():
    curve=np.c_[np.linspace(0,30,121),np.linspace(0,12,121)]
    m=curve_metrics(curve)
    assert np.max(np.abs(m['curvature']))<1e-10
    assert np.max(np.abs(m['signed_turn']))<1e-8
    assert m['absolute_turn'][-1]<1e-8


def test_circle_curvature_and_turn_are_known():
    radius=20.;angle=np.linspace(0,np.pi/2,2001)
    curve=np.c_[radius*np.cos(angle),radius*np.sin(angle)]
    m=curve_metrics(curve)
    assert np.allclose(m['curvature'][~m['boundary']],1/radius,atol=1e-4)
    expected=np.rad2deg(m['stations'][-1]/radius)
    assert abs(m['signed_turn'][-1]-expected)<.05
    assert abs(m['absolute_turn'][-1]-expected)<.05


def test_reverse_bend_absolute_turn_not_cancelled():
    x=np.linspace(0,40,2001);curve=np.c_[x,2*np.sin(2*np.pi*x/40)]
    m=curve_metrics(curve)
    assert abs(m['signed_turn'][-1])<.1
    assert m['absolute_turn'][-1]>60


def test_rigid_transform_invariance_and_scaling():
    t=np.linspace(0,1.5,1001);curve=np.c_[20*np.cos(t),20*np.sin(t)]
    a=.8;rotation=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
    m=curve_metrics(curve);other=curve_metrics(curve@rotation.T+[37,-20])
    assert np.allclose(m['curvature'],other['curvature'],atol=1e-10)
    assert np.allclose(m['signed_turn'],other['signed_turn'],atol=1e-8)
    scaled=curve_metrics(curve*2,grid_mm=.5,smooth_mm=4)
    assert np.allclose(scaled['curvature'],m['curvature']/2,atol=1e-8)
    assert np.allclose(scaled['absolute_turn'],m['absolute_turn'],atol=1e-8)


def test_short_curve_rejected_and_boundary_flagged():
    with pytest.raises(ValueError,match='too short'):curve_metrics(np.c_[np.linspace(0,1,10),np.zeros(10)])
    m=curve_metrics(np.c_[np.linspace(0,10,41),np.zeros(41)])
    assert m['boundary'].sum()==8
    assert m['effective_smooth_mm']==2

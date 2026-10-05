import numpy as np
from morphpot.rim_validation import station_frame, ray_width, read_curves, write_overlay


def test_station_normal_width_on_original_walls():
    mid=np.c_[np.linspace(10,0,11),np.zeros(11)]
    walls={"outer":mid+[0,2],"inner":mid-[0,2]}
    points,normals=station_frame(mid,np.array([1.,3.,9.]))
    assert np.allclose(points[:,0],[9,7,1])
    assert np.allclose(ray_width(points,normals,walls),4.)


def test_missing_intersection_is_missing_not_projection():
    mid=np.c_[np.linspace(10,0,11),np.zeros(11)]
    points,normals=station_frame(mid,np.array([2.,8.]))
    walls={"outer":np.array([[10,2],[5,2]]),"inner":np.array([[10,-2],[5,-2]])}
    widths=ray_width(points,normals,walls)
    assert widths[0]==4 and np.isnan(widths[1])


def test_colored_overlay_preserves_four_disjoint_walls_and_units(tmp_path):
    walls={"outer":np.array([[10,2],[0,2]]),"inner":np.array([[10,-2],[0,-2]])}
    path=tmp_path/"overlay.ply"
    write_overlay(path,walls,walls,1000)
    curves=read_curves(path,1000)
    assert len(curves)==4
    assert np.allclose(curves[0],walls["outer"])
    assert "150 150 150" in path.read_text() and "220 50 50" in path.read_text()

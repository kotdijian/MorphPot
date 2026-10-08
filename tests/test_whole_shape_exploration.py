"""Executable analytical checks for the exploratory geometry, no pytest required."""
import numpy as np
import pottery_whole_shape_exploration as m


def test_necked_jar_detection_segmentation_and_constrained_shared_warp():
    o=np.array([[0,0],[15,0],[40,30],[65,65],[70,90],[60,120],[40,145],[41,150],[48,165]])
    i=np.array([[0,5],[15,5],[35,33],[60,65],[65,90],[55,120],[35,145],[36,150],[48,165]])[::-1]
    lm=m.scan(o,i)
    assert lm['status']=='supported',lm
    p=dict(outer=o,inner=i,landmarks=lm,angle=0.)
    t={ 'neck':(np.array(lm['outer']['point_mm'])+[.1,.1]).tolist(),'lip':(np.array(lm['lip_mm'])+[.1,.1]).tolist()}
    a,b,state=m.warp(p,t)
    assert state['status']=='applied'
    assert np.array_equal(a[:2],o[:2]) and np.array_equal(b[-2:],i[-2:])
    assert np.allclose(a[-1],t['lip']) and np.allclose(b[0],t['lip'])
    shared=dict(p,outer=np.insert(o,3,[50,100],axis=0),inner=np.insert(i,3,[50,100],axis=0))
    a,b,state=m.warp(shared,t)
    assert np.allclose(a[3],b[3]),'Inner and outer must use identical spatial field'
    large={k:(np.array(v)+[20.,0]).tolist() for k,v in t.items()}
    a,b,state=m.warp(p,large)
    assert state['status']=='rejected' and np.array_equal(a,o) and np.array_equal(b,i)
    assert m.neck(np.array([[0.,0],[20,20],[30,50],[50,100]]))['status']!='supported'
    counts,_=m.layout([p,p]);model,stacks,_=m.aggregate([p,p],'neck_segmented',counts,[100,100],t)
    assert np.allclose(model[0][-1],model[1][0])
    assert np.allclose(model[0][counts[0][0]-1],lm['outer']['point_mm'])
    m.contour_valid(np.vstack((model[0],model[1][1:],model[0][:1])))

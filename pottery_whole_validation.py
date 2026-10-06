#!/usr/bin/env python3
"""Standalone whole-vessel validation; no MorphPot package imports.

Reads WholeModel 0.1.0-dev output. Compares aggregators at fixed correspondence
and optionally samples NEW azimuths from the same original mesh. Fixed analysis
axis is intentional: holdout tests aggregation, not independent acquisition.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import trim_mean

VERSION='0.1.0-dev'


def csv_write(path,rows):
    if not rows:return
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def resample(curve, interval):
    curve=np.asarray(curve,float)
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(curve,axis=0),axis=1))]
    keep=np.r_[True,np.diff(s)>1e-10];curve=curve[keep];s=s[keep]
    if len(s)<2:raise ValueError('Degenerate profile')
    stations=np.r_[np.arange(0,s[-1],interval),s[-1]]
    return stations,np.column_stack([np.interp(stations,s,curve[:,j]) for j in [0,1]])


def nearest_segments(points,curve,chunk=128):
    """Exact point-to-segment distances; no nearest-vertex approximation."""
    a=np.asarray(curve[:-1]);d=np.diff(curve,axis=0);length2=np.sum(d*d,axis=1)
    good=length2>1e-20;a=a[good];d=d[good];length2=length2[good]
    if not len(a):raise ValueError('No model segments')
    distance=[];closest=[]
    for start in range(0,len(points),chunk):
        p=np.asarray(points[start:start+chunk]);delta=p[:,None,:]-a
        t=np.clip(np.sum(delta*d,axis=2)/length2,0,1)
        q=a+t[:,:,None]*d
        squared=np.sum((p[:,None,:]-q)**2,axis=2);index=np.argmin(squared,axis=1)
        distance.extend(np.sqrt(squared[np.arange(len(p)),index]));closest.extend(q[np.arange(len(p)),index])
    return np.asarray(distance),np.asarray(closest)


def aggregate(samples,method,trim_fraction):
    if method=='median':return np.median(samples,axis=0)
    if method=='mean':return np.mean(samples,axis=0)
    if method=='trimmed_mean':return trim_mean(samples,trim_fraction,axis=0)
    raise ValueError('Unknown aggregation')


def summary(distance,weights):
    finite=np.isfinite(distance);d=np.asarray(distance)[finite];w=np.asarray(weights)[finite]
    if not len(d):raise ValueError('No finite distances')
    return dict(n_samples=len(d),mean_mm=float(d.mean()),rms_mm=float(np.sqrt(np.mean(d*d))),
        median_mm=float(np.median(d)),p95_mm=float(np.percentile(d,95)),max_mm=float(d.max()),
        area_weighted_mean_mm=float(np.average(d,weights=w)) if w.sum()>0 else None,
        area_weighted_rms_mm=float(np.sqrt(np.average(d*d,weights=w))) if w.sum()>0 else None)


def holdout_angles(step,start,offset,training_step,training_start):
    if step<=0 or 180/step<3 or not np.isclose(180/step,round(180/step)):
        raise ValueError('Holdout step must divide 180 into at least 3 planes')
    angles=(start+offset+np.arange(int(round(180/step)))*step)%180
    remainder=(angles-training_start)%training_step
    if np.any(np.minimum(remainder,training_step-remainder)<1e-6):
        raise ValueError('Holdout azimuth overlaps the attempted training grid; choose an offset/step with no overlap')
    return angles


def clip_positive(poly):
    out=[]
    for a,b in zip(poly[:-1],poly[1:]):
        if a[0]>=0:out.append(a)
        if (a[0]>=0)!=(b[0]>=0):out.append(a+(b-a)*(-a[0]/(b[0]-a[0])))
    if len(out)<3:return np.empty((0,2))
    p=np.asarray(out);p=p[np.r_[True,np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-9]]
    if not np.allclose(p[0],p[-1],atol=1e-9,rtol=0):p=np.vstack((p,p[0]))
    return p


def split_walls(poly):
    p=poly[:-1];axis=np.flatnonzero(abs(p[:,0])<1e-7)
    if len(axis)!=2:raise ValueError('Expected two axial base endpoints')
    a,b=map(int,axis)
    if p[a,1]>p[b,1]:a,b=b,a
    f=p[np.arange(a,a+(b-a)%len(p)+1)%len(p)]
    bwd=p[np.arange(a,a-(a-b)%len(p)-1,-1)%len(p)]
    exposed=max([f,bwd],key=lambda x:float(x[:,0].max()))
    tip=int(np.argmax(exposed[:,1]))
    if tip in [0,len(exposed)-1]:raise ValueError('No supported highest lip')
    return exposed[:tip+1],exposed[tip:]


def load_holdout(mesh_path,meta,offset,step):
    import trimesh
    path=Path(mesh_path);h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(4*1024**2),b''):h.update(part)
    if h.hexdigest()!=meta['input_sha256']:raise ValueError('Input mesh SHA-256 differs from the model source')
    if path.suffix.lower()=='.ply':
        from trimesh.exchange.ply import load_ply
        with path.open('rb') as f:data=load_ply(f,fix_texture=False,skip_materials=True,prefer_color='vertex')
        mesh=trimesh.Trimesh(data['vertices'],data['faces'],process=False)
    else:mesh=trimesh.load(path,force='mesh',process=False)
    mesh.merge_vertices();mesh.remove_unreferenced_vertices();mesh.vertices*=meta['unit_to_mm']
    center=np.asarray(meta['axis_xy_mm']);setting=meta['settings']
    angles=holdout_angles(step,setting['start_angle_deg'],offset,setting['angle_step_deg'],setting['start_angle_deg'])
    profiles=[];qa=[]
    for angle in angles:
        new=[]
        try:
            radial=np.array([np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle))])
            normal=np.r_[-radial[1],radial[0],0]
            seg=trimesh.intersections.mesh_plane(mesh,normal,np.r_[center,mesh.bounds[:,2].mean()])
            if not len(seg):raise ValueError('No intersection')
            xy=np.stack(((seg[:,:,:2]-center)@radial,seg[:,:,2]),axis=-1)
            planar=trimesh.load_path(np.dstack((xy,np.zeros(xy.shape[:2]))))
            if any(d!=2 for _,d in planar.vertex_graph.degree):raise ValueError('Open or branched section')
            polys=[p[:,:2] for p in planar.discrete]
            if not polys or any(len(p)<4 or np.linalg.norm(p[0]-p[-1])>1e-6 for p in polys):raise ValueError('No complete closed contours')
            for sign,theta in [(1,angle),(-1,(angle+180)%360)]:
                halves=[clip_positive(p*np.array([sign,1])) for p in polys]
                halves=[p for p in halves if len(p)>=4]
                if len(halves)!=1:raise ValueError('Multiple/disconnected meridian contours')
                outer,inner=split_walls(halves[0]);new.append((float(theta),outer,inner))
            profiles.extend(new)
            for theta,_,_ in new:qa.append(dict(angle_deg=theta,status='accepted',reason=''))
        except ValueError as e:
            for theta in [angle,(angle+180)%360]:qa.append(dict(angle_deg=float(theta),status='excluded',reason=str(e)))
        print(f'holdout {angle:g} deg: {len(profiles)} accepted halves',flush=True)
    if len(profiles)<.8*len(angles)*2:raise ValueError('Fewer than 80% of holdout halves valid; inspect source topology')
    return profiles,qa


def write_edge_ply(path,curves,scale):
    verts=[];edges=[];colors=[]
    for p,color in curves:
        offset=len(verts);verts.extend(np.c_[p,np.zeros(len(p))]/scale);colors.extend([color]*len(p))
        edges.extend((offset+i,offset+i+1) for i in range(len(p)-1))
    with Path(path).open('w') as f:
        f.write(f'ply\nformat ascii 1.0\nelement vertex {len(verts)}\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement edge {len(edges)}\nproperty int vertex1\nproperty int vertex2\nend_header\n')
        for v,c in zip(verts,colors):f.write(' '.join(f'{x:.10g}' for x in v)+' '+' '.join(map(str,c))+'\n')
        for a,b in edges:f.write(f'{a} {b}\n')


def validate(root,output,mesh_path=None,offset=None,holdout_step=None,interval=1.,trim_fraction=.1,
             methods=('median','mean','trimmed_mean'),plots=True):
    if not np.isfinite([interval,trim_fraction]).all() or interval<=0 or not 0<trim_fraction<.5:raise ValueError('Invalid interval or trim fraction')
    root=Path(root);output=Path(output)
    if root.resolve()==output.resolve() or root.resolve() in output.resolve().parents:
        raise ValueError('Use a separate output directory outside the source model folder')
    if output.exists() and any(output.iterdir()):raise ValueError('Output directory must be empty')
    meta=json.loads((root/'whole_model.json').read_text())
    if meta.get('version')!='0.1.0-dev':raise ValueError('Unsupported WholeModel schema/version')
    with np.load(root/'profile_distribution.npz',allow_pickle=False) as data:
        outer=data['outer_profiles_mm'];inner=data['inner_profiles_mm'];angles=data['accepted_angles_deg']
        stored=[data['outer_median_mm'],data['inner_median_mm']]
    if outer.ndim!=3 or inner.ndim!=3 or outer.shape[0]!=inner.shape[0] or len(angles)!=len(outer) or outer.shape[2]!=2 or inner.shape[2]!=2:raise ValueError('Invalid profile shapes')
    if outer.shape[1]<2 or inner.shape[1]<2 or len(angles)<2:raise ValueError('Too few profiles or curve points')
    if len(angles)!=meta['accepted_half_sections'] or not np.isfinite(outer).all() or not np.isfinite(inner).all():raise ValueError('Count/nonfinite profile data mismatch')
    if not np.allclose(outer[:,-1],inner[:,0],atol=1e-7,rtol=0):raise ValueError('Outer and inner lip endpoints disagree')
    for array,median in zip([outer,inner],stored):
        if not np.allclose(np.median(array,axis=0),median,atol=1e-7,rtol=0):raise ValueError('Stored representative disagrees with sample median')
    if int(len(angles)*trim_fraction)<1 and 'trimmed_mean' in methods:raise ValueError('Too few profiles for requested trim fraction')
    models={m:[aggregate(x,m,trim_fraction) for x in [outer,inner]] for m in methods}
    train=[(float(angle),o,i) for angle,o,i in zip(angles,outer,inner)]
    groups={'training':train};holdout_qa=[]
    if mesh_path is not None:
        step=holdout_step if holdout_step is not None else meta['settings']['angle_step_deg']
        off=offset if offset is not None else meta['settings']['angle_step_deg']/2
        if not np.isfinite([step,off]).all():raise ValueError('Invalid holdout angles')
        groups['holdout'],holdout_qa=load_holdout(mesh_path,meta,off,step)
    elif offset is not None or holdout_step is not None:raise ValueError('Holdout options require --input-mesh')
    output.mkdir(parents=True,exist_ok=True)
    differences=[]
    for m,curves in models.items():
        for side,curve,base in zip(['outer','inner'],curves,stored):
            for k,(p,q) in enumerate(zip(curve,base)):
                differences.append(dict(method=m,surface=side,point_id=k,u=k/(len(curve)-1),r_mm=float(p[0]),z_mm=float(p[1]),delta_r_mm=float(p[0]-q[0]),delta_z_mm=float(p[1]-q[1]),distance_from_median_mm=float(np.linalg.norm(p-q))))
    csv_write(output/'aggregation_comparison.csv',differences)
    np.savez_compressed(output/'comparison_models.npz',**{m+'_'+surface+'_mm':p for m,curves in models.items() for surface,p in zip(['outer','inner'],curves)})
    global_rows=[]
    for group,profiles in groups.items():
        for m,curves in models.items():
            out=output/group/m;out.mkdir(parents=True,exist_ok=True)
            stations=[];per_profile=[];all_curves=[]
            for theta,o,i in profiles:
                for side,source,model in zip(['outer','inner'],[o,i],curves):
                    arc,pts=resample(source,interval)
                    distance,closest=nearest_segments(pts,model)
                    # Arc quadrature times radius approximates swept surface area.
                    ds=np.diff(arc);weights=np.r_[ds[0]/2,(ds[:-1]+ds[1:])/2,ds[-1]/2]*np.maximum(pts[:,0],0)
                    _,model_points=resample(model,interval)
                    reverse,_=nearest_segments(model_points,source)
                    per_profile.append(dict(angle_deg=theta,surface=side,**summary(distance,weights),reverse_rms_mm=float(np.sqrt(np.mean(reverse**2))),reverse_p95_mm=float(np.percentile(reverse,95)),bidirectional_max_mm=float(max(distance.max(),reverse.max()))))
                    for k,(s,p,q,d,w) in enumerate(zip(arc,pts,closest,distance,weights)):
                        stations.append(dict(angle_deg=theta,surface=side,station_id=k,arc_mm=float(s),u=float(s/arc[-1]),source_r_mm=float(p[0]),source_z_mm=float(p[1]),nearest_model_r_mm=float(q[0]),nearest_model_z_mm=float(q[1]),distance_mm=float(d),area_weight=float(w)))
                    all_curves.append((source,(150,150,150)))
            csv_write(out/'station_distances.csv',stations);csv_write(out/'profile_statistics.csv',per_profile)
            for side in ['outer','inner','combined']:
                selected=[s for s in stations if side=='combined' or s['surface']==side]
                global_rows.append(dict(dataset=group,method=m,surface=side,n_profiles=len(profiles),**summary(np.array([s['distance_mm'] for s in selected]),np.array([s['area_weight'] for s in selected]))))
            # Both source and model are open outer/inner curves: no axial closing edge.
            all_curves.extend((p,(220,30,30)) for p in curves)
            write_edge_ply(out/'profiles_all_overlay_xy.ply',all_curves,meta['unit_to_mm'])
            if plots:
                import matplotlib
                matplotlib.use('Agg')
                from matplotlib import pyplot as plt
                fig,ax=plt.subplots(figsize=(7,9))
                for p,c in all_curves:ax.plot(p[:,0],p[:,1],color=np.array(c)/255,alpha=.18 if c[0]==150 else 1,lw=.5 if c[0]==150 else 1.5)
                ax.set(xlabel='Radius [mm]',ylabel='Input Z [mm]',title=f'{group}: {m}, {len(profiles)} halves; fixed analysis axis');ax.set_aspect('equal');fig.tight_layout();fig.savefig(out/'profiles_all_overlay.png',dpi=180);plt.close(fig)
    csv_write(output/'distance_summary.csv',global_rows);csv_write(output/'holdout_section_qa.csv',holdout_qa)
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(11,5))
        for ax,side in zip(axes,['outer','inner']):
            for m,curves in models.items():ax.plot(curves[['outer','inner'].index(side)][:,0],curves[['outer','inner'].index(side)][:,1],label=m)
            ax.set(xlabel='Radius [mm]',ylabel='Input Z [mm]',title=side);ax.set_aspect('equal');ax.legend()
        fig.tight_layout();fig.savefig(output/'aggregation_models_overlay.png',dpi=180);plt.close(fig)
    report=dict(version=VERSION,source_model_sha256=hashlib.sha256((root/'whole_model.json').read_bytes()).hexdigest(),
        settings=dict(interval_mm=interval,trim_fraction=trim_fraction,methods=list(methods),holdout_offset_deg=off if mesh_path else None,holdout_angle_step_deg=step if mesh_path else None),
        datasets={g:len(p) for g,p in groups.items()},axis_xy_mm=meta['axis_xy_mm'],unit=meta['input_unit'],
        definitions=dict(reverse_distance='model to each individual source surface; separate outer/inner; sampled by same mm interval',distance='unsigned exact point-to-polyline segment distance in r,z; outer and inner evaluated separately',
            area_weights='source radius times trapezoidal arc interval; per-azimuth equal weight; not mesh-vertex count',
            trimmed_mean='coordinatewise trimmed mean at fixed correspondence; each r,z variable trimmed independently',
            holdout='new attempted-grid-disjoint azimuths from same original mesh; fixed axis estimated from original mesh'),
        limitations=['not CloudCompare signed C2M or a thickness measurement',
            'training profiles are resampled from stored arc-correspondence arrays; holdout uses raw mesh intersections',
            'holdout shares source mesh and axis, so is not independent physical validation',
            'mean/trimmed curves are comparison descriptors, not validated 3D solids; no alternative meshes exported',
            'no landmark-constrained correspondence, axis sensitivity, or phase/interval sweep implemented here',
            'lowest residual does not establish the best archaeological representative',
            'profile halves/stations are not independent statistical specimens; no confidence interval or significance test'])
    (output/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False))
    print(f'Whole-model validation written to {output}',flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('model_dir',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--input-mesh',type=Path,help='Same original mesh; SHA-256 must match whole_model.json')
    p.add_argument('--holdout-offset-deg',type=float)
    p.add_argument('--holdout-angle-step',type=float)
    p.add_argument('--interval-mm',type=float,default=1)
    p.add_argument('--trim-fraction',type=float,default=.1,help='Fraction removed at EACH end of each coordinate distribution')
    p.add_argument('--methods',nargs='+',choices=['median','mean','trimmed_mean'],default=['median','mean','trimmed_mean'])
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--version',action='version',version=VERSION)
    a=p.parse_args()
    try:validate(a.model_dir,a.output_dir,a.input_mesh,a.holdout_offset_deg,a.holdout_angle_step,a.interval_mm,a.trim_fraction,tuple(dict.fromkeys(a.methods)),not a.no_plots)
    except (ValueError,RuntimeError,FileNotFoundError,KeyError) as e:p.exit(2,f'ERROR: {e}\n')

if __name__=='__main__':main()

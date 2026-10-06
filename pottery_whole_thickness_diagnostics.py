#!/usr/bin/env python3
"""Rim thickness correspondence diagnostics, separate from whole-model synthesis.

Requires MorphPot checkout (existing rim midline implementation and validation
helpers). Coordinates are untransformed r,z in mm; edge PLY uses input units.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from pottery_whole_validation import aggregate,csv_write,resample,write_edge_ply,load_holdout,nearest_segments
from pottery_whole_dimension_validation import thickness
from morphpot.rim_standardization import make_midline,arc_positions,sample_curve
from morphpot.rim_models import normal_frame
VERSION='0.1.0-dev'


def ray_hit(origin,direction,curve,expected_arc=None,window=8.):
    starts=curve[:-1];edges=np.diff(curve,axis=0);arc=arc_positions(curve)
    cross=lambda a,b:a[...,0]*b[...,1]-a[...,1]*b[...,0]
    den=cross(direction,edges);good=abs(den)>1e-12
    t=np.full(len(den),np.nan);v=t.copy()
    t[good]=cross(starts[good]-origin,edges[good])/den[good]
    v[good]=cross(starts[good]-origin,direction)/den[good]
    hit_arc=arc[:-1]+v*np.diff(arc)
    valid=good & (t>1e-7) & (v>=-1e-9) & (v<=1+1e-9)
    if expected_arc is not None:valid &= abs(hit_arc-expected_arc)<=window
    ids=np.flatnonzero(valid)
    if not len(ids):return None
    j=ids[np.argmin(t[ids])];point=origin+t[j]*direction
    return dict(distance=float(t[j]),point=point,arc=float(hit_arc[j]),incidence=float(abs(den[j])/np.linalg.norm(edges[j])),count=len(ids))


def point_arc(point,curve):
    a=curve[:-1];d=np.diff(curve,axis=0);l2=np.sum(d*d,axis=1)
    t=np.clip(np.sum((point-a)*d,axis=1)/np.maximum(l2,1e-20),0,1)
    q=a+t[:,None]*d;j=int(np.argmin(np.sum((point-q)**2,axis=1)))
    return float(arc_positions(curve)[j]+t[j]*np.sqrt(l2[j]))


def measure(outer,inner,interval=1.,rim_length=40.,window=8.,margin=2.):
    # All correspondence coordinates run from highest lip toward the body.
    ol=outer[::-1].copy();il=inner.copy()
    # Clipped axis endpoints can retain positive floating-point roundoff.
    for wall in [ol,il]:
        if abs(wall[-1,0])<1e-7:wall[-1,0]=0.
    raw=thickness(outer,inner,interval,margin,2.)
    total=arc_positions(outer)[-1];rows=[];curves=[];qa=[]
    for r in raw:
        lip_arc=total-r['arc_mm']
        if lip_arc>rim_length:continue
        p=np.array([r['r_mm'],r['z_mm']]);n=None
        s=r['arc_mm'];before=max(0,s-1);after=min(total,s+1)
        tangent=sample_curve(outer,[after])[0]-sample_curve(outer,[before])[0]
        if np.linalg.norm(tangent)>1e-10:n=np.array([-tangent[1],tangent[0]])/np.linalg.norm(tangent)
        rawhit=np.array([r['inner_r_mm'],r['inner_z_mm']]) if r['inner_r_mm'] is not None else None
        rawarc=point_arc(rawhit,il) if rawhit is not None else None
        h=ray_hit(p,n,il,lip_arc,window) if n is not None and r['status']!='endpoint_margin' else None
        status='no_local_intersection';value=None
        if r['status']=='endpoint_margin':status='endpoint_margin'
        elif h:
            status='local_normal';value=h['distance']
            if h['incidence']<.3:status='grazing_inner_intersection';value=None
        rows.append(dict(method='outer_unrestricted',lip_arc_mm=float(lip_arc),thickness_mm=r['thickness_mm'],raw_distance_mm=r['raw_normal_distance_mm'],status=r['status'],origin_r_mm=float(p[0]),origin_z_mm=float(p[1]),outer_r_mm=float(p[0]),outer_z_mm=float(p[1]),inner_r_mm=float(rawhit[0]) if rawhit is not None else None,inner_z_mm=float(rawhit[1]) if rawhit is not None else None,outer_lip_arc_mm=float(lip_arc),inner_lip_arc_mm=rawarc,arc_mismatch_mm=abs(rawarc-lip_arc) if rawarc is not None else None,incidence_cos=r['inner_incidence_cos'],projection_fallback_mm=None))
        rows.append(dict(method='outer_constrained',lip_arc_mm=float(lip_arc),thickness_mm=value,raw_distance_mm=h['distance'] if h else None,status=status,origin_r_mm=float(p[0]),origin_z_mm=float(p[1]),outer_r_mm=float(p[0]),outer_z_mm=float(p[1]),inner_r_mm=float(h['point'][0]) if h else None,inner_z_mm=float(h['point'][1]) if h else None,outer_lip_arc_mm=float(lip_arc),inner_lip_arc_mm=h['arc'] if h else None,arc_mismatch_mm=abs(h['arc']-lip_arc) if h else None,incidence_cos=h['incidence'] if h else None,projection_fallback_mm=None))
        if rawhit is not None:curves.append((np.array([p,rawhit]),(230,140,20)))
        if h:curves.append((np.array([p,h['point']]),(30,100,230)))
    geometry={};diagnostic={}
    try:
        mid,_=make_midline(ol,il[::-1],'right',.5,end_mm=rim_length,source_geometry=geometry,diagnostics=diagnostic)
        a=geometry['paired_outer'];b=geometry['paired_inner']
        mid_arc=arc_positions(mid);positions=np.r_[np.arange(0,min(rim_length,mid_arc[-1]),interval),min(rim_length,mid_arc[-1])]
        mids=sample_curve(mid,positions)
        # Pair coordinates share the DTW midline index; interpolate by midline arc,
        # not by separate wall arc lengths.
        pairs_o=np.column_stack([np.interp(positions,mid_arc,a[:,j]) for j in [0,1]])
        pairs_i=np.column_stack([np.interp(positions,mid_arc,b[:,j]) for j in [0,1]])
        normals=normal_frame(mids,(pairs_o-mids)[None])
        for k,(p,n,po,pi) in enumerate(zip(mids,normals,pairs_o,pairs_i)):
            expected_o=point_arc(po,ol);expected_i=point_arc(pi,il)
            ho=ray_hit(p,n,ol,expected_o,window);hi=ray_hit(p,-n,il,expected_i,window)
            value=None;status='missing_wall_intersection'
            if k==0:status='tip_excluded'
            elif ho and hi:
                if min(ho['incidence'],hi['incidence'])<.3:status='grazing_wall_intersection'
                else:value=ho['distance']+hi['distance'];status='midline_normal'
            projected=float(abs(np.dot(po-p,n))+abs(np.dot(pi-p,n)))
            rows.append(dict(method='rim_midline',lip_arc_mm=float(expected_o),thickness_mm=value,raw_distance_mm=ho['distance']+hi['distance'] if ho and hi else None,status=status,origin_r_mm=float(p[0]),origin_z_mm=float(p[1]),outer_r_mm=float(ho['point'][0]) if ho else None,outer_z_mm=float(ho['point'][1]) if ho else None,inner_r_mm=float(hi['point'][0]) if hi else None,inner_z_mm=float(hi['point'][1]) if hi else None,outer_lip_arc_mm=ho['arc'] if ho else None,inner_lip_arc_mm=hi['arc'] if hi else None,arc_mismatch_mm=abs(ho['arc']-hi['arc']) if ho and hi else None,incidence_cos=min(ho['incidence'],hi['incidence']) if ho and hi else None,projection_fallback_mm=projected))
            if ho and hi:curves.append((np.array([ho['point'],hi['point']]),(20,160,70)))
        curves.append((mids,(20,160,70)))
        qa.append(dict(component='midline',status='ok',reason='',details=json.dumps(diagnostic)))
    except ValueError as e:qa.append(dict(component='midline',status='failed',reason=str(e),details=''))
    return rows,curves,qa


def plot_profile(path,outer,inner,rows,curves,title,rim_length):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    ax=axes[0]
    for p in [outer,inner]:ax.plot(p[:,0],p[:,1],color='black',lw=1)
    for p,c in curves:ax.plot(p[:,0],p[:,1],color=np.array(c)/255,alpha=.45,lw=.7)
    lip=outer[-1];ax.set_xlim(lip[0]-rim_length-10,lip[0]+15);ax.set_ylim(lip[1]-rim_length-10,lip[1]+5);ax.set_aspect('equal');ax.set(xlabel='Radius [mm]',ylabel='Input Z [mm]',title=title)
    for method,color in [('outer_unrestricted','orange'),('outer_constrained','royalblue'),('rim_midline','green')]:
        r=sorted([x for x in rows if x['method']==method],key=lambda x:x['lip_arc_mm'])
        axes[1].plot([x['lip_arc_mm'] for x in r],[x['thickness_mm'] if x['thickness_mm'] is not None else np.nan for x in r],label=method,color=color)
    axes[1].set(xlabel='Outer arc from highest lip [mm]',ylabel='Thickness [mm]',title='No missing-value fallback');axes[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)


def run(root,output,mesh=None,offset=2.5,interval=1.,rim_length=40.,window=8.,plots=True):
    root=Path(root);output=Path(output)
    if not np.isfinite([offset,interval,rim_length,window]).all() or min(interval,rim_length,window)<=0:raise ValueError('Invalid settings')
    if root.resolve()==output.resolve() or root.resolve() in output.resolve().parents:raise ValueError('Use output outside model folder')
    if output.exists() and any(output.iterdir()):raise ValueError('Output must be empty')
    meta=json.loads((root/'whole_model.json').read_text())
    if meta.get('version')!='0.1.0-dev':raise ValueError('Unsupported WholeModel version')
    with np.load(root/'profile_distribution.npz',allow_pickle=False) as d:
        outer=d['outer_profiles_mm'];inner=d['inner_profiles_mm'];angles=d['accepted_angles_deg']
    if outer.ndim!=3 or inner.ndim!=3 or outer.shape[2]!=2 or inner.shape[2]!=2 or len(outer)!=len(inner) or len(angles)!=len(outer):raise ValueError('Invalid profile shapes')
    if not np.isfinite(outer).all() or not np.isfinite(inner).all() or not np.allclose(outer[:,-1],inner[:,0],atol=1e-7,rtol=0):raise ValueError('Invalid profile coordinates')
    models={m:[aggregate(x,m,.1) for x in [outer,inner]] for m in ['median','mean','trimmed_mean'] if m!='trimmed_mean' or len(angles)>=10}
    groups={'training':list(zip(angles,outer,inner))};holdout_qa=[]
    if mesh:groups['holdout'],holdout_qa=load_holdout(mesh,meta,offset,meta['settings']['angle_step_deg'])
    output.mkdir(parents=True,exist_ok=True);allrows=[];allqa=[];summary=[];group_overlays={};group_bounds={}
    entries=[('models',m,o,i) for m,(o,i) in models.items()]
    entries += [(group,f'angle_{float(angle):07.3f}',o,i) for group,profiles in groups.items() for angle,o,i in profiles]
    for group,label,o,i in entries:
        rows,connectors,qa=measure(o,i,interval,rim_length,window)
        out=output/group/label;out.mkdir(parents=True,exist_ok=True)
        csv_write(out/'measurements.csv',rows)
        # Exact measured connectors, including invalid-but-intersecting rays.
        overlay=[(o,(60,60,60)),(i,(60,60,60))]+connectors
        write_edge_ply(out/'thickness_rays_xy.ply',overlay,meta['unit_to_mm'])
        group_overlays.setdefault(group,[]).extend(overlay)
        for wall in [o[::-1],i]:
            local=sample_curve(wall,np.linspace(0,min(rim_length,arc_positions(wall)[-1]),100))
            group_bounds.setdefault(group,[]).extend(local)
        if plots:plot_profile(out/'thickness_rays.png',o,i,rows,connectors,f'{group}: {label}',rim_length)
        allrows.extend(dict(dataset=group,profile=label,**r) for r in rows);allqa.extend(dict(dataset=group,profile=label,**q) for q in qa)
        for method in ['outer_unrestricted','outer_constrained','rim_midline']:
            r=[x for x in rows if x['method']==method];valid=[x['thickness_mm'] for x in r if x['thickness_mm'] is not None]
            summary.append(dict(dataset=group,profile=label,method=method,n_stations=len(r),n_valid=len(valid),max_thickness_mm=max(valid) if valid else None,n_arc_mismatch_gt_window=sum(x['arc_mismatch_mm'] is not None and x['arc_mismatch_mm']>window for x in r)))
        print(f'{group}/{label}: correspondence diagnostics written',flush=True)
    for group,overlay in group_overlays.items():
        write_edge_ply(output/f'{group}_all_rays_overlay_xy.ply',overlay,meta['unit_to_mm'])
        if plots:
            from matplotlib import pyplot as plt
            from matplotlib.collections import LineCollection
            fig,ax=plt.subplots(figsize=(8,7))
            for color in [(60,60,60),(230,140,20),(30,100,230),(20,160,70)]:
                lines=[seg for p,c in overlay if c==color for seg in np.stack([p[:-1],p[1:]],axis=1)]
                ax.add_collection(LineCollection(lines,colors=[np.array(color)/255],linewidths=.5,alpha=.2 if group=='models' else .06))
            bounds=np.asarray(group_bounds[group]);ax.set_xlim(bounds[:,0].min()-5,bounds[:,0].max()+5);ax.set_ylim(bounds[:,1].min()-5,bounds[:,1].max()+5);ax.set_aspect('equal');ax.set(xlabel='Radius [mm]',ylabel='Input Z [mm]',title=f'{group}: all rim profiles and measured rays')
            for color,label in [('gray','walls'),('orange','unrestricted'),('royalblue','constrained'),('green','midline')]:ax.plot([],[],color=color,label=label)
            ax.legend();fig.tight_layout();fig.savefig(output/f'{group}_all_rays_overlay.png',dpi=160);plt.close(fig)
    csv_write(output/'all_measurements.csv',allrows);csv_write(output/'measurement_summary.csv',summary);csv_write(output/'midline_qa.csv',allqa);csv_write(output/'holdout_section_qa.csv',holdout_qa)
    # Model method comparison at common recorded outer arc stations. Midline is
    # interpolated only between adjacent VALID observations; no extrapolation.
    comparison=[]
    for label in models:
        rr=[r for r in allrows if r['dataset']=='models' and r['profile']==label]
        mid=sorted([r for r in rr if r['method']=='rim_midline'],key=lambda r:r['lip_arc_mm'])
        x=np.array([r['lip_arc_mm'] for r in mid])
        for raw in [r for r in rr if r['method']=='outer_unrestricted']:
            s=raw['lip_arc_mm'];con=next(r for r in rr if r['method']=='outer_constrained' and r['lip_arc_mm']==s);mv=None
            k=int(np.searchsorted(x,s))
            if k<len(x) and abs(x[k]-s)<1e-9:mv=mid[k]['thickness_mm']
            elif 0<k<len(x) and mid[k-1]['thickness_mm'] is not None and mid[k]['thickness_mm'] is not None and x[k]>x[k-1]:mv=float(np.interp(s,x[k-1:k+1],[mid[k-1]['thickness_mm'],mid[k]['thickness_mm']]))
            comparison.append(dict(model=label,outer_lip_arc_mm=s,unrestricted_mm=raw['thickness_mm'],constrained_mm=con['thickness_mm'],midline_mm=mv,unrestricted_minus_midline_mm=raw['thickness_mm']-mv if raw['thickness_mm'] is not None and mv is not None else None,constrained_minus_midline_mm=con['thickness_mm']-mv if con['thickness_mm'] is not None and mv is not None else None))
    csv_write(output/'model_method_comparison.csv',comparison)
    report=dict(version=VERSION,source_model_sha256=hashlib.sha256((root/'whole_model.json').read_bytes()).hexdigest(),settings=dict(interval_mm=interval,rim_length_mm=rim_length,pair_window_mm=window,outer_endpoint_margin_mm=2,min_incidence_cos=.3),datasets={g:len(p) for g,p in groups.items()},definitions=dict(colors={'orange':'unrestricted outer normals','blue':'constrained outer normals','green':'rim midline and normal connectors','gray':'wall curves'},constraint='outer normals: absolute outer/inner arc-from-highest-lip mismatch <= window; midline normals: each wall hit within window of its local DTW pair arc',midline='existing make_midline manual rim length + extension tip and DTW; no similarity/affine registration; each curve measured in original r,z',fallback='paired projection saved for diagnosis only; not substituted for missing ray thickness'),limitations=['Thickness from outer and central normals describes different directions; differences are not automatically errors.','Fixed arc window is a diagnostic parameter, not established anatomical correspondence.','Manual rim length must fit both walls; midline failures reported without substitute.','Tip and highest lip are different definitions; common outer-wall arc used for comparisons.','Training uses resampled stored curves; holdout uses same mesh and fixed axis.'])
    (output/'thickness_diagnostics.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('model_dir',type=Path);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--input-mesh',type=Path);p.add_argument('--holdout-offset-deg',type=float,default=2.5);p.add_argument('--interval-mm',type=float,default=1);p.add_argument('--rim-length-mm',type=float,default=40);p.add_argument('--pair-window-mm',type=float,default=8);p.add_argument('--no-plots',action='store_true');p.add_argument('--version',action='version',version=VERSION)
    a=p.parse_args()
    try:run(a.model_dir,a.output_dir,a.input_mesh,a.holdout_offset_deg,a.interval_mm,a.rim_length_mm,a.pair_window_mm,not a.no_plots)
    except (ValueError,RuntimeError,FileNotFoundError,KeyError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()

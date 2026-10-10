#!/usr/bin/env python3
"""Exploratory necked-jar whole models. MorphPot checkout required on PYTHONPATH.
No pose transform. Provisional morphological constraints; not a classifier.
Methods: global_arc; neck_segmented; neck_segmented_warp.
Native-unit PLY; all numeric tables in mm. Original vs warped validation separate.
"""
from __future__ import annotations
import argparse,json,hashlib,csv
from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter,find_peaks
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from morphpot.asset_metadata import resolve_units
from pottery_radial_sections import load_geometry
from morphpot.whole_vessel import closed_sections,clip_positive,split_walls,revolve,resample_wall,estimate_axis
from pottery_whole_validation import nearest_segments,resample,write_edge_ply,csv_write
from pottery_whole_dimension_validation import dimensions,thickness

VERSION='0.1.0-exploration'
CONSTRAINTS=dict(detector_step_mm=.5,smoothing_window_mm=5.,neck_prominence_mm=1.,
    neck_min_upper_arc_fraction=.45,neck_max_upper_arc_fraction=.985,
    minimum_body_neck_contrast_mm=3.,minimum_lip_neck_flare_mm=1.,
    max_outer_inner_neck_height_difference_mm=8.,minimum_supported_fraction=.8,
    warp_max_displacement_mm=5.,warp_max_shear=.15,warp_max_axial_strain=.15)

def arc(p):return np.r_[0.,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]

def at(p,s):
    a=arc(p);return np.array([np.interp(s,a,p[:,j]) for j in [0,1]])

def neck(wall):
    # wall ordered from axial base to highest lip. Smooth detector only.
    a,p=resample(wall,.5);window=min(len(p)//2*2-1,11)
    if window<5:return dict(status='indeterminate',reason='too_short')
    smooth=savgol_filter(p[:,0],window,3)
    indices,properties=find_peaks(-smooth,prominence=CONSTRAINTS['neck_prominence_mm'])
    candidates=[]
    for index,prominence in zip(indices,properties['prominences']):
        u=a[index]/a[-1];body=int(np.argmax(smooth[:index+1]))
        contrast=smooth[body]-smooth[index];flare=smooth[-1]-smooth[index]
        if .45<=u<=.985 and contrast>=3 and flare>=1:
            candidates.append(dict(arc_mm=float(a[index]),u=float(u),point_mm=at(wall,a[index]).tolist(),
                                   prominence_mm=float(prominence),body_neck_contrast_mm=float(contrast),lip_flare_mm=float(flare)))
    if not candidates:return dict(status='unsupported',reason='no_upper_radial_minimum',candidates=[])
    candidates.sort(key=lambda x:x['prominence_mm'],reverse=True)
    if len(candidates)>1 and candidates[1]['prominence_mm']>.7*candidates[0]['prominence_mm']:
        return dict(status='indeterminate',reason='competing_necks',candidates=candidates)
    return dict(status='supported',**candidates[0],candidates=candidates)

def scan(o,i):
    a=neck(o);b=neck(i[::-1]);status='supported'
    if a['status']!='supported' or b['status']!='supported':status='indeterminate' if 'indeterminate' in [a['status'],b['status']] else 'unsupported'
    elif abs(a['point_mm'][1]-b['point_mm'][1])>8:status='indeterminate'
    return dict(status=status,outer=a,inner=b,lip_mm=o[-1].tolist(),base_outer_mm=o[0].tolist(),base_inner_mm=i[-1].tolist())

def extract(mesh,center,angles,group):
    profiles=[];qa=[]
    for angle in angles:
        pair=[]
        try:
            polygons=closed_sections(mesh,center,float(angle))
            for sign,theta in [(1,float(angle)),(-1,float((angle+180)%360))]:
                halves=[clip_positive(p*np.array([sign,1])) for p in polygons]
                halves=[p for p in halves if len(p)>=4]
                if len(halves)!=1:raise ValueError('multiple_or_disconnected_meridians')
                o,i=split_walls(halves[0]);pair.append(dict(angle=theta,outer=o,inner=i,landmarks=scan(o,i)))
            profiles.extend(pair)
            qa.extend(dict(dataset=group,angle_deg=p['angle'],status='accepted',reason='') for p in pair)
        except ValueError as error:
            qa.extend(dict(dataset=group,angle_deg=float(theta),status='excluded',reason=str(error)) for theta in [angle,(angle+180)%360])
        print(f'{group}: {angle:g} deg, {len(profiles)} halves',flush=True)
    return profiles,qa

def sample_parts(wall,s,count):
    lengths=arc(wall);s=float(np.clip(s,0,lengths[-1]))
    values=np.r_[np.linspace(0,s,count[0]),np.linspace(s,lengths[-1],count[1])[1:]]
    return np.column_stack([np.interp(values,lengths,wall[:,j]) for j in [0,1]])

def layout(profiles):
    counts=[];indices=[]
    for side in ['outer','inner']:
        # Common base->lip direction for detection; export inner remains lip->base.
        lengths=[arc(p[side])[-1] for p in profiles]
        s=[p['landmarks'][side]['arc_mm'] for p in profiles]
        c=[max(3,int(np.ceil(np.median(s)/.5))+1),max(3,int(np.ceil(np.median(np.array(lengths)-s)/.5))+1)]
        counts.append(c);indices.append(c[0]-1)
    return counts,indices

def warp(p,target):
    o,i=p['outer'],p['inner'];lm=p['landmarks'];height=o[-1,1]-o[0,1]
    guard=max(o[0,1],i[-1,1])+.1*height
    source=np.array([lm['outer']['point_mm'],lm['lip_mm']]);desired=np.array([target['neck'],target['lip']])
    knots=np.array([guard,source[0,1],source[1,1]])
    if np.any(np.diff(knots)<=2):return o.copy(),i.copy(),dict(status='rejected',reason='unordered_height_anchors')
    delta=np.vstack(([0.,0.],desired-source))
    # Same spatial field for BOTH walls. Identity below bottom guard.
    dr=PchipInterpolator(knots,delta[:,0],extrapolate=False);dz=PchipInterpolator(knots,delta[:,1],extrapolate=False)
    z=np.linspace(knots[0],knots[-1],500);shear=float(np.max(abs(dr.derivative()(z))));strain=float(np.max(abs(dz.derivative()(z))))
    displacement=float(np.max(np.linalg.norm(np.c_[dr(z),dz(z)],axis=1)))
    details=dict(status='applied',knots_z_mm=knots.tolist(),delta_rz_mm=delta.tolist(),maximum_displacement_mm=displacement,maximum_shear=shear,maximum_axial_strain=strain)
    if displacement>5 or shear>.15 or strain>.15:
        return o.copy(),i.copy(),dict(**{**details,'status':'rejected'},reason='deformation_limit')
    def apply(w):
        q=np.asarray(w,dtype=float).copy();active=w[:,1]>guard;zz=np.clip(w[active,1],knots[0],knots[-1]);q[active]+=np.c_[dr(zz),dz(zz)];return q
    oo,ii=apply(o),apply(i)
    if min(oo[:,0].min(),ii[:,0].min())<-1e-7:
        return o.copy(),i.copy(),dict(**{**details,'status':'rejected'},reason='crossed_axis')
    return oo,ii,details

def contour_valid(poly):
    segments=np.stack([poly[:-1],poly[1:]],axis=1);mid=segments.mean(axis=1);half=np.linalg.norm(segments[:,1]-segments[:,0],axis=1)/2
    cross=lambda x,y:x[0]*y[1]-x[1]*y[0]
    for a,b in cKDTree(mid).query_pairs(float(2*half.max()+1e-7)):
        if b==a+1 or (a==0 and b==len(segments)-1):continue
        p,q=segments[a];r,s=segments[b]
        if cross(q-p,r-p)*cross(q-p,s-p)<0 and cross(s-r,p-r)*cross(s-r,q-r)<0:raise ValueError('self_intersecting_model')

def aggregate(profiles,method,counts,global_counts,target):
    stacks=[[],[]];details=[]
    for p in profiles:
        o,i=p['outer'],p['inner'];state=dict(status='not_requested')
        if method=='neck_segmented_warp':o,i,state=warp(p,target)
        details.append(dict(angle_deg=p['angle'],**state))
        for k,w in enumerate([o,i]):
            if method=='global_arc':sample=resample_wall(w,global_counts[k])
            else:
                # Knot arc comes from RAW source correspondence; sample before
                # spatial warp to avoid reparameterizing the transformed contour.
                raw=p['outer'] if k==0 else p['inner'][::-1]
                sample=sample_parts(raw,p['landmarks'][['outer','inner'][k]]['arc_mm'],counts[k])
                if k==1:sample=sample[::-1]
                if method=='neck_segmented_warp':
                    # Evaluate same transformation field at sampled original points.
                    if state['status']=='applied':
                        knots=np.array(state['knots_z_mm']);delta=np.array(state['delta_rz_mm']);active=sample[:,1]>knots[0];zz=np.clip(sample[active,1],knots[0],knots[-1]);sample[active]+=np.c_[PchipInterpolator(knots,delta[:,0])(zz),PchipInterpolator(knots,delta[:,1])(zz)]
            stacks[k].append(sample)
    arrays=[np.stack(x) for x in stacks];median=[np.median(x,axis=0) for x in arrays]
    return median,arrays,details

def stats(v):
    v=np.asarray(v,float);return dict(n=len(v),mean_mm=float(v.mean()),rms_mm=float(np.sqrt(np.mean(v*v))),p95_mm=float(np.percentile(v,95)),max_mm=float(v.max()))

def measure_group(profiles,models,target,method,group):
    distances=[];summary=[];dimensions_rows=[];thickness_rows=[];curves=[]
    for p in profiles:
        original=[p['outer'],p['inner']];o,i,state=warp(p,target) if method=='neck_segmented_warp' else (*original,{'status':'not_requested'})
        registered=[o,i];curves.extend(registered)
        for frame,walls in [('original',original),('registered',registered)]:
            if frame=='registered' and method!='neck_segmented_warp':continue
            dimensions_rows.append(dict(dataset=group,method=method,frame=frame,angle_deg=p['angle'],warp_status=state['status'],**dimensions(*walls)))
            for side,w,model in zip(['outer','inner'],walls,models):
                s,points=resample(w,1);d,_=nearest_segments(points,model)
                distances.extend(dict(dataset=group,method=method,frame=frame,angle_deg=p['angle'],surface=side,arc_mm=float(a),source_r_mm=float(pt[0]),source_z_mm=float(pt[1]),distance_mm=float(v)) for a,pt,v in zip(s,points,d))
            for row in thickness(*walls):thickness_rows.append(dict(dataset=group,method=method,frame=frame,angle_deg=p['angle'],**row))
    for frame in ['original','registered']:
        for side in ['outer','inner','combined']:
            values=[r['distance_mm'] for r in distances if r['frame']==frame and (side=='combined' or r['surface']==side)]
            if values:summary.append(dict(dataset=group,method=method,frame=frame,surface=side,**stats(values)))
    return distances,summary,dimensions_rows,thickness_rows,curves

def run(a):
    output=a.output_dir
    if output.exists() and any(output.iterdir()):raise ValueError('Output must be empty')
    output.mkdir(parents=True,exist_ok=True);analysis=output/'analysis';validation=output/'validation';analysis.mkdir();validation.mkdir()
    unit,scale,metadata=resolve_units(a.input,a.unit);mesh=load_geometry(a.input);mesh.merge_vertices();mesh.remove_unreferenced_vertices()
    center,axisqa=estimate_axis(mesh,scale,'inner',40);mesh.vertices*=scale
    train,qa=extract(mesh,center,np.arange(0,180,a.angle_step),'training')
    np.savez_compressed(analysis/'raw_training_profiles.npz',**{f'p{k:03d}_{s}_mm':p[s] for k,p in enumerate(train) for s in ['outer','inner']},angles_deg=[p['angle'] for p in train])
    csv_write(analysis/'section_qa.csv',qa)
    (analysis/'landmarks.json').write_text(json.dumps([dict(angle_deg=p['angle'],**p['landmarks']) for p in train],indent=2))
    fig,ax=plt.subplots(figsize=(8,9))
    for p in train:
        for s in ['outer','inner']:ax.plot(*p[s].T,color='.6',lw=.4,alpha=.2)
        if p['landmarks']['status']=='supported':
            q=np.array([p['landmarks']['outer']['point_mm'],p['landmarks']['inner']['point_mm']]);ax.scatter(*q.T,s=7,c='orange')
    ax.set(aspect='equal',xlabel='Radius [mm]',ylabel='Input Z [mm]',title='Preflight: original walls / neck candidates');fig.tight_layout();fig.savefig(validation/'preflight_landmarks.png',dpi=160);plt.close(fig)
    valid=[p for p in train if p['landmarks']['status']=='supported'];requested=int(round(360/a.angle_step))
    status='proceed' if len(train)>=.8*requested and len(valid)>=.8*len(train) else ('mismatch' if len(valid)<.2*len(train) else 'indeterminate')
    report=dict(version=VERSION,category=a.category,input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),input_unit=unit,unit_to_mm=scale,axis_xy_mm=center.tolist(),axis_diagnostics=axisqa,pose_transform_applied=False,constraints=CONSTRAINTS,
                preflight=dict(status=status,attempted_halves=requested,topologically_valid=len(train),landmark_supported=len(valid)),settings=dict(angle_step_deg=a.angle_step,holdout_offset_deg=a.angle_step/2,profile_spacing_mm=.5,revolution_sections=180))
    (output/'exploration.json').write_text(json.dumps(report,indent=2))
    if status!='proceed':raise SystemExit('STOP: preflight '+status+'; inspect exploration.json and preflight_landmarks.png')
    # Same cohort across all methods. Excluded unsupported landmark profiles logged.
    holdout,hqa=extract(mesh,center,np.arange(a.angle_step/2,180,a.angle_step),'holdout');csv_write(validation/'holdout_section_qa.csv',hqa)
    np.savez_compressed(analysis/'raw_holdout_profiles.npz',**{f'p{k:03d}_{s}_mm':p[s] for k,p in enumerate(holdout) for s in ['outer','inner']},angles_deg=[p['angle'] for p in holdout])
    (analysis/'holdout_landmarks.json').write_text(json.dumps([dict(angle_deg=p['angle'],**p['landmarks']) for p in holdout],indent=2))
    held=[p for p in holdout if p['landmarks']['status']=='supported']
    counts,indices=layout(valid);global_counts=[max(3,int(np.ceil(np.median([arc(p[s])[-1] for p in valid])/.5))+1) for s in ['outer','inner']]
    target=dict(neck=np.median([p['landmarks']['outer']['point_mm'] for p in valid],axis=0).tolist(),lip=np.median([p['landmarks']['lip_mm'] for p in valid],axis=0).tolist())
    report.update(holdout=dict(topologically_valid=len(holdout),landmark_supported=len(held),attempted_halves=requested),training_angles_deg=[p['angle'] for p in valid],holdout_angles_deg=[p['angle'] for p in held],target_landmarks_mm=target,segment_counts=counts)
    summaries=[];dims=[];thick=[];models={};warp_records=[]
    for method in ['global_arc','neck_segmented','neck_segmented_warp']:
        curves,arrays,records=aggregate(valid,method,counts,global_counts,target);o,i=curves;poly=np.vstack((o,i[1:],o[:1]));contour_valid(poly)
        folder=analysis/method;folder.mkdir();vfolder=validation/method;vfolder.mkdir()
        model=revolve(poly,180,center);model.vertices/=scale;model.export(folder/'whole_representative.ply')
        models[method]=curves;dims.append(dict(dataset='model',method=method,frame='model',angle_deg='',warp_status='',**dimensions(o,i)))
        np.savez_compressed(folder/'profile_distribution.npz',outer_profiles_mm=arrays[0],inner_profiles_mm=arrays[1],outer_median_mm=o,inner_median_mm=i,angles_deg=[p['angle'] for p in valid])
        write_edge_ply(folder/'representative_profile_xy.ply',[(o,(230,30,30)),(i,(230,30,30))],scale)
        csv_write(folder/'representative_profile.csv',[dict(surface=s,point_id=k,r_mm=float(q[0]),z_mm=float(q[1])) for s,p in zip(['outer','inner'],curves) for k,q in enumerate(p)])
        for row in thickness(o,i):thick.append(dict(dataset='model',method=method,frame='model',angle_deg='',**row))
        warp_records.extend(dict(dataset='training',method=method,**r) for r in records)
        for group,profiles in [('training',valid),('holdout',held)]:
            data,summary,dimension_rows,thickness_rows,registered=measure_group(profiles,curves,target,method,group);summaries.extend(summary);dims.extend(dimension_rows);thick.extend(thickness_rows)
            csv_write(vfolder/f'{group}_station_distances.csv',data)
            for frame,walls in [('original',[p[s] for p in profiles for s in ['outer','inner']]),('registered',registered)]:
                if frame=='registered' and method!='neck_segmented_warp':continue
                write_edge_ply(vfolder/f'{group}_{frame}_all_overlay_xy.ply',[(p,(150,150,150)) for p in walls]+[(p,(230,30,30)) for p in curves],scale)
                fig,ax=plt.subplots(figsize=(7,9))
                for w in walls:ax.plot(*w.T,color='.6',lw=.4,alpha=.2)
                for w in curves:ax.plot(*w.T,color='crimson',lw=1.5)
                ax.set(aspect='equal',xlabel='Radius [mm]',ylabel='Input Z [mm]',title=f'{method}: {group} / {frame}, {len(profiles)} halves');fig.tight_layout();fig.savefig(vfolder/f'{group}_{frame}_all_overlay.png',dpi=160);plt.close(fig)
        meta=dict(version=VERSION,method=method,dimensions_mm=dimensions(o,i),watertight=model.is_watertight,winding_consistent=model.is_winding_consistent,vertices=len(model.vertices),faces=len(model.faces),unit=unit,unit_to_mm=scale)
        (folder/'model.json').write_text(json.dumps(meta,indent=2))
    # All topologically valid raw sections provide additional non-selected holdout audit.
    for method,curves in models.items():
        d,summary,_,_,_=measure_group(holdout,curves,target,'global_arc','holdout_all_raw')
        for row in summary:row['method']=method
        summaries.extend(summary)
    report['warp_training']=[r for r in warp_records if r['method']=='neck_segmented_warp']
    report['warp_holdout']=[dict(angle_deg=p['angle'],**warp(p,target)[2]) for p in held]
    csv_write(validation/'distance_summary.csv',summaries);csv_write(validation/'dimensions.csv',dims);csv_write(validation/'thickness_stations.csv',thick)
    (output/'exploration.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    fig,axes=plt.subplots(1,2,figsize=(12,7))
    colors={'global_arc':'black','neck_segmented':'royalblue','neck_segmented_warp':'crimson'}
    for ax in axes:
        for p in valid:
            for s in ['outer','inner']:ax.plot(*p[s].T,color='.7',lw=.3,alpha=.1)
        for method,curves in models.items():
            for k,p in enumerate(curves):ax.plot(*p.T,color=colors[method],lw=1.5,label=method if k==0 else None)
        ax.set(aspect='equal',xlabel='Radius [mm]',ylabel='Input Z [mm]');ax.legend(fontsize=8)
    axes[0].set_title('Whole profiles / original frame');axes[1].set_title('Neck / rim detail');axes[1].set_ylim(target['neck'][1]-10,target['lip'][1]+5);axes[1].set_xlim(target['neck'][0]-15,target['lip'][0]+15)
    fig.tight_layout();fig.savefig(validation/'method_comparison.png',dpi=170);plt.close(fig)
    # Explanatory paired figures use the same saved warp and raw source stations.
    # They do not change the model or the existing validation measurements.
    if held:
        from plot_whole_registration_comparison import run as registration_figures
        registration_figures(output)
    print('COMPLETED',output,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--category',choices=['necked_jar'],required=True);p.add_argument('--unit',choices=['auto','m','mm','cm'],default='auto');p.add_argument('--angle-step',type=float,default=5.)
    a=p.parse_args()
    if a.angle_step<=0 or 180/a.angle_step<3 or not np.isclose(180/a.angle_step,round(180/a.angle_step)):p.error('Angle step must divide 180 into >=3 planes')
    run(a)

if __name__=='__main__':main()

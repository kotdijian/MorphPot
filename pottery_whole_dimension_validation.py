#!/usr/bin/env python3
"""Validate whole-vessel dimensions and outer-normal wall thickness.

Requires sibling pottery_whole_validation.py. All measurements are in mm;
coordinates retain the input pose. No alignment or pose transform is applied.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from pottery_whole_validation import aggregate, csv_write, resample, load_holdout
VERSION='0.1.0-dev'


def dimensions(outer, inner):
    lip=outer[-1]
    maximum=float(np.max(outer[:,0]));peak=outer[np.isclose(outer[:,0],maximum,atol=1e-8,rtol=0),1]
    base=float(np.min(outer[:,1]))
    return dict(lip_diameter_mm=float(2*lip[0]),maximum_diameter_mm=2*maximum,
        maximum_diameter_height_mm=float(np.mean([peak.min(),peak.max()])-base),
        maximum_diameter_z_min_mm=float(peak.min()),maximum_diameter_z_max_mm=float(peak.max()),
        vessel_height_mm=float(max(outer[:,1].max(),inner[:,1].max())-base),base_z_mm=base)


def thickness(outer,inner,interval=1.,margin=2.,tangent_window=2.):
    arc,points=resample(outer,interval)
    original=np.r_[0,np.cumsum(np.linalg.norm(np.diff(outer,axis=0),axis=1))]
    # Local chord evaluated symmetrically in physical arc length; no dependence
    # on original vertex density. The outer wall runs axial base -> lip.
    before=np.maximum(0,arc-tangent_window/2);after=np.minimum(arc[-1],arc+tangent_window/2)
    a=np.column_stack([np.interp(before,original,outer[:,j]) for j in [0,1]])
    b=np.column_stack([np.interp(after,original,outer[:,j]) for j in [0,1]])
    tangent=b-a;norm=np.linalg.norm(tangent,axis=1)
    normals=np.column_stack([-tangent[:,1],tangent[:,0]])/np.maximum(norm[:,None],1e-20)
    segments=np.diff(inner,axis=0);starts=inner[:-1]
    rows=[]
    cross=lambda x,y:x[...,0]*y[...,1]-x[...,1]*y[...,0]
    for k,(s,p,n) in enumerate(zip(arc,points,normals)):
        status='ok';value=None;hit=None;count=0;raw=None;incidence=None
        if s<margin or arc[-1]-s<margin:status='endpoint_margin'
        elif norm[k]<1e-10:status='degenerate_tangent'
        else:
            denom=cross(n,segments);good=abs(denom)>1e-12
            t=np.full(len(denom),np.nan);v=t.copy()
            t[good]=cross(starts[good]-p,segments[good])/denom[good]
            v[good]=cross(starts[good]-p,n)/denom[good]
            hits=np.sort(t[good & (t>1e-7) & (v>=-1e-9) & (v<=1+1e-9)])
            hits=hits[np.r_[True,np.diff(hits)>1e-6]] if len(hits) else hits
            count=len(hits)
            if count:
                value=float(hits[0]);raw=value;hit=p+value*n
                candidates=np.flatnonzero(np.isfinite(t) & (abs(t-value)<1e-6) & (v>=-1e-9) & (v<=1+1e-9))
                incidence=float(max(abs(denom[j])/max(np.linalg.norm(segments[j]),1e-20) for j in candidates))
                if count>1:status='multiple_hits_first_used'
                if incidence<.3:status='grazing_inner_intersection';value=None
            else:status='no_inner_intersection'
        rows.append(dict(station_id=k,arc_mm=float(s),u=float(s/arc[-1]),r_mm=float(p[0]),z_mm=float(p[1]),
            thickness_mm=value,raw_normal_distance_mm=raw,inner_incidence_cos=incidence,status=status,intersection_count=count,
            inner_r_mm=float(hit[0]) if hit is not None else None,inner_z_mm=float(hit[1]) if hit is not None else None))
    return rows


def distribution(values):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if not len(v):return dict(n=0,mean_mm=None,median_mm=None,std_mm=None,min_mm=None,max_mm=None)
    return dict(n=len(v),mean_mm=float(v.mean()),median_mm=float(np.median(v)),std_mm=float(v.std(ddof=0)),min_mm=float(v.min()),max_mm=float(v.max()))


def run(root,output,mesh=None,offset=2.5,interval=1.,margin=2.,window=2.,plots=True):
    root=Path(root);output=Path(output)
    if not np.isfinite([interval,margin,window,offset]).all() or interval<=0 or margin<0 or window<=0:raise ValueError('Invalid measurement settings')
    if output.resolve()==root.resolve() or root.resolve() in output.resolve().parents:raise ValueError('Output must be outside source model directory')
    if output.exists() and any(output.iterdir()):raise ValueError('Output directory must be empty')
    meta=json.loads((root/'whole_model.json').read_text())
    if meta.get('version')!='0.1.0-dev':raise ValueError('Unsupported WholeModel version')
    with np.load(root/'profile_distribution.npz',allow_pickle=False) as data:
        outer=data['outer_profiles_mm'];inner=data['inner_profiles_mm'];angles=data['accepted_angles_deg']
        stored=[data['outer_median_mm'],data['inner_median_mm']]
    if outer.ndim!=3 or inner.ndim!=3 or outer.shape[2]!=2 or inner.shape[2]!=2 or len(outer)!=len(inner) or len(angles)!=len(outer) or len(angles)<2:raise ValueError('Invalid profile arrays')
    if not np.isfinite(outer).all() or not np.isfinite(inner).all() or not np.allclose(outer[:,-1],inner[:,0],atol=1e-7,rtol=0):raise ValueError('Invalid coordinates/lip endpoints')
    for samples,curve in zip([outer,inner],stored):
        if not np.allclose(np.median(samples,axis=0),curve,atol=1e-7,rtol=0):raise ValueError('Stored median mismatch')
    methods=['median','mean']+(['trimmed_mean'] if int(len(angles)*.1)>0 else [])
    models={m:[aggregate(x,m,.1) for x in [outer,inner]] for m in methods}
    groups={'training':list(zip(angles,outer,inner))};qa=[]
    if mesh:groups['holdout'],qa=load_holdout(mesh,meta,offset,meta['settings']['angle_step_deg'])
    output.mkdir(parents=True,exist_ok=True)
    model_dims=[];model_thickness={};dimension_rows=[];thickness_rows=[];comparison=[];summaries=[]
    for m,(o,i) in models.items():
        model_dims.append(dict(method=m,**dimensions(o,i)))
        model_thickness[m]=thickness(o,i,interval,margin,window)
        csv_write(output/f'{m}_model_thickness.csv',model_thickness[m])
    for group,profiles in groups.items():
        profiles_measured=[]
        for angle,o,i in profiles:
            d=dimensions(o,i);dimension_rows.append(dict(dataset=group,angle_deg=float(angle),**d))
            rows=thickness(o,i,interval,margin,window);profiles_measured.append(rows)
            thickness_rows.extend(dict(dataset=group,angle_deg=float(angle),**r) for r in rows)
        for m,(o,i) in models.items():
            md=dimensions(o,i)
            for metric in ['lip_diameter_mm','maximum_diameter_mm','maximum_diameter_height_mm','vessel_height_mm']:
                values=[r[metric] for r in dimension_rows if r['dataset']==group];stats=distribution(values)
                summaries.append(dict(dataset=group,method=m,metric=metric,model_mm=md[metric],**stats,model_minus_median_mm=md[metric]-stats['median_mm']))
            for station in model_thickness[m]:
                values=[]
                for rows in profiles_measured:
                    # Never interpolate across missing or excluded measurements.
                    u=np.array([r['u'] for r in rows]);index=int(np.searchsorted(u,station['u']))
                    if index<len(u) and abs(u[index]-station['u'])<1e-10:
                        val=rows[index]['thickness_mm']
                    elif 0<index<len(u) and rows[index-1]['thickness_mm'] is not None and rows[index]['thickness_mm'] is not None:
                        val=float(np.interp(station['u'],u[index-1:index+1],[rows[index-1]['thickness_mm'],rows[index]['thickness_mm']]))
                    else:val=None
                    if val is not None:values.append(val)
                stats=distribution(values);mv=station['thickness_mm']
                comparison.append(dict(dataset=group,method=m,station_id=station['station_id'],model_arc_mm=station['arc_mm'],u=station['u'],model_thickness_mm=mv,model_status=station['status'],**stats,model_minus_median_mm=mv-stats['median_mm'] if mv is not None and stats['n'] else None))
    csv_write(output/'model_dimensions.csv',model_dims);csv_write(output/'source_dimensions.csv',dimension_rows)
    csv_write(output/'dimension_summary.csv',summaries);csv_write(output/'source_thickness.csv',thickness_rows)
    csv_write(output/'thickness_comparison.csv',comparison);csv_write(output/'holdout_section_qa.csv',qa)
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        for group in groups:
            fig,axes=plt.subplots(2,2,figsize=(11,7))
            for ax,metric in zip(axes.flat,['lip_diameter_mm','maximum_diameter_mm','maximum_diameter_height_mm','vessel_height_mm']):
                r=[x for x in dimension_rows if x['dataset']==group];ax.scatter([x['angle_deg'] for x in r],[x[metric] for x in r],s=12,label='source halves')
                for m in methods:ax.axhline(next(x for x in model_dims if x['method']==m)[metric],label=m)
                ax.set(xlabel='Azimuth [deg]',ylabel='mm',title=metric);ax.legend(fontsize=8)
            fig.tight_layout();fig.savefig(output/f'{group}_dimensions.png',dpi=160);plt.close(fig)
            fig,axes=plt.subplots(len(methods),1,figsize=(11,3*len(methods)),squeeze=False)
            for ax,m in zip(axes.flat,methods):
                rows=[x for x in comparison if x['dataset']==group and x['method']==m]
                conv=lambda key:np.array([np.nan if r[key] is None else r[key] for r in rows])
                x=conv('model_arc_mm');ax.fill_between(x,conv('min_mm'),conv('max_mm'),alpha=.15,label='source min-max')
                ax.plot(x,conv('median_mm'),label='source median');ax.plot(x,conv('model_thickness_mm'),label='model')
                ax.set(xlabel='Model outer arc from axial base [mm]',ylabel='Normal thickness [mm]',title=f'{group}: {m}; normalized-arc correspondence');ax.legend()
            fig.tight_layout();fig.savefig(output/f'{group}_thickness.png',dpi=160);plt.close(fig)
    report=dict(version=VERSION,source_model_sha256=hashlib.sha256((root/'whole_model.json').read_bytes()).hexdigest(),datasets={k:len(v) for k,v in groups.items()},settings=dict(interval_mm=interval,endpoint_margin_mm=margin,tangent_window_mm=window,holdout_offset_deg=offset if mesh else None),
        definitions=dict(lip_diameter='2 * radius of highest lip endpoint; diameter equivalent of EACH half, not opposite-side diameter or inner aperture diameter',maximum_diameter='2 * maximum outer radius of each half',maximum_diameter_height='height above minimum outer Z; midpoint of Z range if multiple exact maxima',vessel_height='highest outer/inner Z minus minimum outer Z',thickness='first positive intersection of locally estimated inward outer normal with inner polyline; multiple intersections flagged',correspondence='equal normalized outer-wall arc from axial base to highest lip; NOT anatomical landmarks',std='population standard deviation ddof=0; equal weight per accepted half'),
        limitations=['Lip diameter is a highest-lip proxy; true opening diameter is not identified automatically.', 'Endpoint margins exclude lip and axial base thickness; missing values are not filled.', 'Inner incidence cosine below 0.3 is excluded as grazing; raw normal distance and intersection retained. Multiple-hit normals may reach a remote wall; inspect flags and geometry before interpreting extremes.', 'Training curves are stored resampled profiles, not direct mesh intersections.', 'Holdout shares mesh and fixed axis; this is not independent physical measurement.', 'Normalized arc correspondence may pair different anatomical positions.', 'Mean and trimmed profiles are descriptors, not validated solid meshes.', 'No restored/original segmentation; all accepted profiles used.'])
    (output/'dimension_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False))
    print(f'Dimension validation written to {output}')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('model_dir',type=Path);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--input-mesh',type=Path);p.add_argument('--holdout-offset-deg',type=float,default=2.5)
    p.add_argument('--interval-mm',type=float,default=1);p.add_argument('--endpoint-margin-mm',type=float,default=2);p.add_argument('--tangent-window-mm',type=float,default=2)
    p.add_argument('--no-plots',action='store_true');p.add_argument('--version',action='version',version=VERSION)
    a=p.parse_args()
    try:run(a.model_dir,a.output_dir,a.input_mesh,a.holdout_offset_deg,a.interval_mm,a.endpoint_margin_mm,a.tangent_window_mm,not a.no_plots)
    except (ValueError,RuntimeError,FileNotFoundError,KeyError) as e:p.exit(2,f'ERROR: {e}\n')
if __name__=='__main__':main()

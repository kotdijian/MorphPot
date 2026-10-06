"""Whole-vessel coordinate-median meridian; independent of rim fitting.

Complete closed mesh-plane sections are clipped at the analysis axis. Each
opposite half is one angular observation. Outer and inner contours are sampled
separately by normalized arc and aggregated without changing their source loci.
No landmark registration or local similarity/affine deformation is performed.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from morphpot.asset_metadata import resolve_units, sha256_file

VERSION = '0.1.0-dev'


def clip_positive(poly):
    """Clip a closed material polygon at r=0; preserve true axial base closure."""
    out=[]
    for a, b in zip(poly[:-1], poly[1:]):
        ina, inb = a[0] >= 0, b[0] >= 0
        if ina: out.append(a)
        if ina != inb:
            out.append(a+(b-a)*(-a[0]/(b[0]-a[0])))
    if len(out)<3: return np.empty((0,2))
    p=np.asarray(out, float)
    p=p[np.r_[True, np.linalg.norm(np.diff(p,axis=0),axis=1)>1e-9]]
    if np.linalg.norm(p[0]-p[-1])>1e-9:p=np.vstack((p,p[0]))
    return p


def split_walls(poly):
    """Exposed meridian: outer axial base -> highest lip -> inner axial base.

    Only the artificial axial clipping edge is removed. No whole-body DTW,
    pose fitting, or local deformation is used.
    """
    p=poly[:-1]
    axial=np.flatnonzero(np.abs(p[:,0])<1e-7)
    if len(axial)!=2:raise ValueError('Expected exactly two axial base endpoints')
    a,b=map(int,axial)
    if p[a,1]>p[b,1]:a,b=b,a
    forward=p[np.arange(a,a+(b-a)%len(p)+1)%len(p)]
    backward=p[np.arange(a,a-(a-b)%len(p)-1,-1)%len(p)]
    path=max([forward,backward],key=lambda c:float(np.max(c[:,0])))
    tip=int(np.argmax(path[:,1]))
    if tip==0 or tip==len(path)-1:raise ValueError('Highest lip not supported between outer and inner base')
    return path[:tip+1],path[tip:]


def resample_wall(curve,count):
    lengths=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(curve,axis=0),axis=1))]
    keep=np.r_[True,np.diff(lengths)>1e-9]
    c=curve[keep];arc=lengths[keep]
    if len(c)<2 or arc[-1]<=0:raise ValueError('Degenerate wall curve')
    u=np.linspace(0,arc[-1],count)
    return np.column_stack([np.interp(u,arc,c[:,j]) for j in range(2)])


def aggregate_profiles(walls, spacing):
    # Same normalized arc index separately for outer and inner surfaces.
    # Original individual loci remain unchanged; part lengths are not fitted.
    counts=[max(3,int(np.ceil(np.median([np.sum(np.linalg.norm(np.diff(w[k],axis=0),axis=1)) for w in walls])/spacing))+1) for k in [0,1]]
    samples=[np.stack([resample_wall(w[k],counts[k]) for w in walls]) for k in [0,1]]
    outer,inner=[np.median(x,axis=0) for x in samples]
    poly=np.vstack((outer,inner[1:],outer[:1]))
    # Reject self-crossing material contours before constructing a 3D solid.
    # This is an explicit geometric check, not automatic topology repair.
    from scipy.spatial import cKDTree
    segments=np.stack((poly[:-1],poly[1:]),axis=1)
    midpoint=segments.mean(axis=1);half=np.linalg.norm(np.diff(segments,axis=1)[:,0],axis=1)/2
    pairs=cKDTree(midpoint).query_pairs(float(2*half.max()+1e-7))
    for i,j in pairs:
        if j==i+1 or (i==0 and j==len(segments)-1):continue
        a,b=segments[i];c,d=segments[j]
        if np.any(np.maximum(a,b)<np.minimum(c,d)-1e-7) or np.any(np.maximum(c,d)<np.minimum(a,b)-1e-7):continue
        cross=lambda u,v:u[0]*v[1]-u[1]*v[0]
        v1,v2=cross(b-a,c-a),cross(b-a,d-a)
        v3,v4=cross(d-c,a-c),cross(d-c,b-c)
        if v1*v2<0 and v3*v4<0:raise ValueError('Median inner/outer meridian self-intersects; no topology repair performed')
    return poly,samples


def revolve(poly, sections, center):
    mesh=trimesh.creation.revolve(poly,sections=sections)
    mesh.merge_vertices()
    mesh.remove_unreferenced_vertices()
    mesh.fix_normals()
    mesh.vertices[:,:2]+=np.asarray(center)
    if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume<=0:
        raise ValueError('Revolved representative is not a valid closed material mesh')
    return mesh


def write_csv(path, rows):
    if not rows:return
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def write_edges(path, curves, scale):
    verts=[];edges=[];colors=[]
    for curve,color in curves:
        v=np.column_stack((curve,np.zeros(len(curve))))/scale
        offset=len(verts);verts.extend(v);colors.extend([color]*len(v))
        edges.extend((offset+i,offset+i+1) for i in range(len(v)-1))
    with Path(path).open('w') as f:
        f.write('ply\nformat ascii 1.0\nelement vertex '+str(len(verts))+'\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement edge '+str(len(edges))+'\nproperty int vertex1\nproperty int vertex2\nend_header\n')
        for v,c in zip(verts,colors):f.write(' '.join(f'{x:.10g}' for x in v)+' '+' '.join(map(str,c))+'\n')
        for a,b in edges:f.write(f'{a} {b}\n')


def estimate_axis(mesh, scale, surface, count):
    # Reuse geometric primitives only, never invoke the rim/volume pipeline.
    from pottery_radial_sections import (horizontal_section_segments, closed_horizontal_contours,
        select_outer_inner,resample_closed_polyline,robust_ellipse_fit,choose_axis_center)
    bounds=mesh.bounds[:,2]
    heights=np.linspace(bounds[0],bounds[1],count+2)[1:-1]
    records=[]
    for i,z in enumerate(heights):
        contours=closed_horizontal_contours(horizontal_section_segments(mesh,z),.05/scale)
        outer,inner=select_outer_inner(contours)
        rec=dict(section_id=i,z_mm=float(z*scale),outer_fit=None,inner_fit=None)
        for key,c in [('outer',outer),('inner',inner)]:
            if c is not None:rec[key+'_fit']=robust_ellipse_fit(resample_closed_polyline(c['points_xy'],.5/scale),scale)
        records.append(rec)
    center,summary=choose_axis_center(records,surface,3.)
    return np.asarray(center),summary


def closed_sections(mesh, center, angle):
    radial=np.array([np.cos(np.deg2rad(angle)),np.sin(np.deg2rad(angle))])
    normal=np.r_[-radial[1],radial[0],0]
    origin=np.r_[center,mesh.bounds[:,2].mean()]
    seg=trimesh.intersections.mesh_plane(mesh,normal,origin)
    if not len(seg):raise ValueError('No plane intersection')
    xy=np.stack(((seg[:,:,:2]-center)@radial,seg[:,:,2]),axis=-1)
    path=trimesh.load_path(np.dstack((xy,np.zeros(xy.shape[:2]))))
    if any(d!=2 for _,d in path.vertex_graph.degree):raise ValueError('Open or branched intersection; no repair performed')
    polys=[]
    for p in path.discrete:
        p=np.asarray(p)[:,:2]
        if len(p)<4 or np.linalg.norm(p[0]-p[-1])>1e-6:raise ValueError('Unclosed intersection')
        polys.append(p)
    if not polys:raise ValueError('No closed contours')
    return polys


def boundary_samples(poly, spacing):
    seg=np.diff(poly,axis=0);lens=np.linalg.norm(seg,axis=1)
    pieces=[]
    for a,d,length in zip(poly[:-1],seg,lens):
        # Axis clip is not an exposed surface of the 3D vessel.
        if abs(a[0])<1e-8 and abs((a+d)[0])<1e-8:continue
        n=max(1,int(np.ceil(length/spacing)))
        pieces.append(a+np.arange(n)[:,None]/n*d)
    return np.concatenate(pieces)


def run(input_path, output, unit='auto', angle_step=5., start_angle=0., pitch=.5,
        axis_xy_mm=None, axis_surface='inner', axis_sections=40, revolution_sections=180,
        max_memory_mb=1024., plots=True, require_watertight=False):
    from pottery_radial_sections import load_geometry
    if not np.isfinite([angle_step,start_angle,pitch,max_memory_mb]).all():raise ValueError('Options must be finite')
    n=180/angle_step if angle_step>0 else 0
    if n<3 or not np.isclose(n,round(n)):raise ValueError('angle-step must divide 180 and give at least 3 planes')
    if pitch<=0 or axis_sections<3 or revolution_sections<12 or max_memory_mb<=0:raise ValueError('Invalid resolution/axis/memory setting')
    output=Path(output)
    if output.exists() and any(output.iterdir()):raise ValueError('Output directory must be empty; choose a new directory to avoid stale results')
    unit,scale,metadata=resolve_units(input_path,unit)
    mesh=load_geometry(Path(input_path));original_watertight=bool(mesh.is_watertight)
    # Duplicate vertices are welded for topological QA only; positions stay fixed.
    mesh.merge_vertices();mesh.remove_unreferenced_vertices()
    if not np.isfinite(mesh.vertices).all():raise ValueError('Non-finite mesh coordinates')
    edge_counts=np.bincount(mesh.edges_unique_inverse)
    topology=dict(watertight_after_welding=bool(mesh.is_watertight),boundary_edges=int(np.count_nonzero(edge_counts==1)),nonmanifold_edges=int(np.count_nonzero(edge_counts>2)))
    if require_watertight and not mesh.is_watertight:raise ValueError('Input material mesh is not watertight after duplicate welding; no hole filling performed')
    if not mesh.is_watertight:print('Input topology warning: '+json.dumps(topology)+'; each section will be checked individually',flush=True)
    if axis_xy_mm is None:center,axisqa=estimate_axis(mesh,scale,axis_surface,axis_sections)
    else:
        center=np.asarray(axis_xy_mm,float);axisqa={'method':'explicit_axis_xy_mm'}
        if center.shape!=(2,) or not np.isfinite(center).all():raise ValueError('Invalid axis center')
    mesh.vertices*=scale
    requested=int(round(n))*2
    source_curves=[];qa=[];accepted=[];walls=[]
    for i in range(int(round(n))):
        angle=(start_angle+i*angle_step)%180
        try:
            polys=closed_sections(mesh,center,angle)
            for side,theta in [(1,angle),(-1,(angle+180)%360)]:
                halfpolys=[p*np.array([side,1]) for p in polys]
                clipped=[clip_positive(p) for p in halfpolys]
                curves=[p for p in clipped if len(p)>=4]
                if len(curves)!=1:raise ValueError('Expected one connected meridian per angular half')
                wall=split_walls(curves[0]);walls.append(wall);source_curves.extend(curves)
                accepted.append((theta,curves));qa.append(dict(angle_deg=float(theta),status='accepted',reason='',contours=len(polys),outer_length_mm=float(np.linalg.norm(np.diff(wall[0],axis=0),axis=1).sum()),inner_length_mm=float(np.linalg.norm(np.diff(wall[1],axis=0),axis=1).sum())))
        except ValueError as e:
            # Both halves of a failed plane must be excluded, including a half
            # tentatively appended before the other one failed.
            while accepted and np.isclose(accepted[-1][0]%180,angle):
                walls.pop();_,cs=accepted.pop();del source_curves[-len(cs):];qa.pop()
            for theta in [angle,(angle+180)%360]:qa.append(dict(angle_deg=float(theta),status='excluded',reason=str(e),contours=0,outer_length_mm='',inner_length_mm=''))
        print(f'plane {i+1}/{int(round(n))}: {angle:g} deg; accepted halves {len(walls)}',flush=True)
    if len(walls)<.8*requested:raise ValueError(f'Only {len(walls)}/{requested} halves valid (<80%); representative is unsupported')
    sample_counts=[max(3,int(np.ceil(np.median([np.linalg.norm(np.diff(w[k],axis=0),axis=1).sum() for w in walls])/pitch))+1) for k in [0,1]]
    estimate_mb=(len(walls)*sum(sample_counts)*8*2*3+sum(sample_counts)*revolution_sections*120)/1024**2
    if estimate_mb>max_memory_mb:raise ValueError('Profile arrays exceed memory limit')
    poly,samples=aggregate_profiles(walls,pitch)
    model=revolve(poly,revolution_sections,center)
    # Nearest-contour residuals are descriptive, not landmark correspondences.
    model_samples=boundary_samples(poly,pitch/2)
    tree=cKDTree(model_samples)
    deviations=[]
    for theta,curves in accepted:
        pts=np.concatenate([boundary_samples(p,pitch/2) for p in curves])
        distance=tree.query(pts)[0]
        deviations.append(dict(angle_deg=float(theta),n_samples=len(pts),mean_distance_mm=float(distance.mean()),rms_distance_mm=float(np.sqrt(np.mean(distance**2))),p95_distance_mm=float(np.percentile(distance,95)),max_distance_mm=float(distance.max())))
    output=Path(output)
    if output.exists() and any(output.iterdir()):raise ValueError('Output directory must be empty; choose a new directory to avoid stale results')
    output.mkdir(parents=True,exist_ok=True)
    native=model.copy();native.vertices/=scale;native.export(output/'whole_representative.ply')
    write_edges(output/'representative_profile_xy.ply',[(poly,(220,30,30))],scale)
    write_edges(output/'profiles_all_overlay_xy.ply',[(p,(150,150,150)) for p in source_curves]+[(poly,(220,30,30))],scale)
    write_csv(output/'representative_profile.csv',[dict(point_id=i,r_mm=float(p[0]),z_mm=float(p[1])) for i,p in enumerate(poly)])
    write_csv(output/'angular_deviations.csv',deviations);write_csv(output/'section_qa.csv',qa)
    np.savez_compressed(output/'profile_distribution.npz',outer_profiles_mm=samples[0],inner_profiles_mm=samples[1],
        outer_median_mm=np.median(samples[0],axis=0),inner_median_mm=np.median(samples[1],axis=0),accepted_angles_deg=[a for a,_ in accepted])
    report=dict(version=VERSION,input_sha256=sha256_file(input_path),input_unit=unit,unit_to_mm=scale,
        method='coordinate median of complete outer/inner meridians separately registered by normalized arc; axisymmetric revolution',pose_transform_applied=False,
        local_similarity_or_affine_applied=False,axis_xy_mm=center.tolist(),axis_diagnostics=axisqa,
        settings=dict(angle_step_deg=angle_step,start_angle_deg=start_angle,profile_spacing_mm=pitch,
            revolution_sections=revolution_sections,axis_sections=axis_sections,axis_surface=axis_surface),
        requested_half_sections=requested,accepted_half_sections=len(walls),input_watertight_before_welding=original_watertight,input_topology=topology,
        model=dict(watertight=bool(model.is_watertight),winding_consistent=bool(model.is_winding_consistent),
            vertices=len(model.vertices),faces=len(model.faces),height_mm=float(model.extents[2]),
            maximum_diameter_mm=float(2*poly[:,0].max()),material_volume_mm3=float(model.volume)),
        definitions=dict(lip='highest vertex of each complete exposed half-section, not the rim midline extension tip',profile_axes='X=nonnegative radius from analysis axis; Y=input Z height; PLY Z=0',
            deviation='unsigned nearest exposed meridian distance; discretized at pitch/2; not thickness or homologous point error',
            aggregation='equal weight per accepted angular half; outer axial base to highest lip and lip to inner axial base registered separately by normalized arc; coordinate median; no local transforms',
            units='PLY uses native input unit in original pose; CSV/NPZ/plots use mm'),
        limitations=['axisymmetric model suppresses angular asymmetry; deviations are retained',
            'coordinate median is not a pre-firing reconstruction or the earlier rim model',
            'no original/restored surface labels; all accepted sections participate',
            'no thickness midline extraction, part segmentation, landmark matching, or inter-individual classification',
            'normalized arc correspondence is not anatomical homology; compare with landmarks in the next stage',
            'unsigned nearest distances can match a different surface branch; not thickness errors',
            'axis positions are estimated, axis tilt and center drift are not corrected',
            'a single connected meridian meeting axis is required; handles, multiple material islands, hollow pedestals may be rejected'])
    (output/'whole_model.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    if plots:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        fig,ax=plt.subplots(figsize=(8,9))
        for p in source_curves:ax.plot(p[:,0],p[:,1],color='.6',alpha=.2,lw=.5)
        ax.plot(poly[:,0],poly[:,1],color='crimson',lw=1.5,label='Whole-vessel median profile')
        ax.set(xlabel='Radius from analysis axis [mm]',ylabel='Input Z height [mm]',title=f'{len(walls)} angular halves; profile spacing {pitch:g} mm; no local transforms')
        ax.set_aspect('equal');ax.legend();fig.tight_layout();fig.savefig(output/'profiles_all_overlay.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,4));ax.plot([d['angle_deg'] for d in sorted(deviations,key=lambda d:d['angle_deg'])],[d['rms_distance_mm'] for d in sorted(deviations,key=lambda d:d['angle_deg'])],'.-')
        ax.set(xlabel='Azimuth [deg]',ylabel='Nearest meridian RMS distance [mm]',title='Angular deviation; descriptive in-sample residual');fig.tight_layout();fig.savefig(output/'angular_deviations.png',dpi=180);plt.close(fig)
    print(json.dumps(report['model'],indent=2),flush=True)
    return report

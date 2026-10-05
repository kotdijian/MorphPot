"""Phase/interval experiment from one dense extraction, with subset refitting."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from .rim_validation import load_run, rows_csv, validate_run, compare_runs, station_frame, ray_width, write_colored_curves, read_curves
from .rim_standardization import (_write_curves, _write_verification, fit_similarity,
    fit_constrained_affine, apply_transform, arc_positions)
from .rim_models import export_section_models


def phase_indices(records, step, offset, origin=0.):
    if not isinstance(step,int) or step<=0 or 180%step:
        raise ValueError("Interval must be an integer divisor of 180")
    return [i for i,r in enumerate(records)
            if min(((float(r["angle_deg"])-origin-offset)%180)%step,
                   step-((float(r["angle_deg"])-origin-offset)%180)%step)<1e-6]


def rebuild_subset(source, indices, output, attempted_profiles=None):
    """Reuse dense extraction/cropping but refit both transforms and all models.

    This holds axis, original wall extraction and endpoint policy constant and
    isolates sampling effects in registration and model construction.
    """
    out=Path(output)/"rim_standardization";out.mkdir(parents=True,exist_ok=True)
    scale=source["scale"]; root=source["root"]
    records=[source["report"]["profiles"][i] for i in indices]
    mids=[source["raw_mids"][i] for i in indices]
    right=[m for m,r in zip(mids,records) if r["side"]=="right"]
    if len(mids)<2 or not right:raise ValueError("Subset needs >=2 profiles and a right-side reference")
    with (root/"raw_paired_points.csv").open(encoding="utf-8-sig") as f: pair_rows=list(csv.DictReader(f))
    grouped={i:[] for i in indices}
    for r in pair_rows:
        i=int(r["profile_id"])
        if i in grouped:grouped[i].append(r)
    pairs=[]
    for i in indices:
        ordered=sorted(grouped[i],key=lambda r:int(r["point_id"]))
        if len(ordered)!=len(source["raw_mids"][i]):raise ValueError("Raw pair count mismatch")
        pairs.append([np.array([[float(r[k+"_x_mm"]),float(r[k+"_y_mm"])] for r in ordered]) for k in ("outer","inner")])
    tips=read_curves(root/"raw_tip_extensions_xy.ply",scale)
    geometry=[{k:source["raw_"+k][i] for k in ("outer","inner")} for i in indices]
    for g,i in zip(geometry,indices):g["tip_extension"]=tips[i]
    _write_curves(out/"raw_midlines_xy.ply",mids,scale,(40,190,70))
    _write_verification(out,"raw",geometry,pairs,scale)
    raw_rows=[]
    for local,(i,m,p,r) in enumerate(zip(indices,mids,pairs,records)):
        for j,(a,b,c) in enumerate(zip(*p,m)):
            raw_rows.append(dict(profile_id=local,angle_deg=r["angle_deg"],side=r["side"],point_id=j,u=j/(len(m)-1),
                outer_x_mm=a[0],outer_y_mm=a[1],inner_x_mm=b[0],inner_y_mm=b[1],mid_x_mm=c[0],mid_y_mm=c[1]))
    rows_csv(out/"raw_paired_points.csv",raw_rows)
    reference=np.median(np.stack(right),axis=0); config=source["qa"]["config"]
    weights=np.ones(len(reference));weights[[0,-1]]=config["endpoint_weight"]
    qa=dict(source["qa"]);qa.update(valid_profiles=len(indices),profiles=records,
        attempted_profiles=len(indices) if attempted_profiles is None else attempted_profiles,
        source_attempted_profiles=source["qa"].get("attempted_profiles"),
        subset_source=str(root.resolve()),source_profile_ids=indices,
        extraction_policy="reuse 1-degree axis, wall extraction, tip and crop; recompute subset right reference, transforms, outlier screening and models")
    transforms=[]
    for mode in ("similarity","affine"):
        folder=out/f"standard_{mode}";folder.mkdir(exist_ok=True)
        mapped=[];mapped_pairs=[];mapped_geometry=[];mode_transforms=[]
        for local,(i,m,p,g,r) in enumerate(zip(indices,mids,pairs,geometry,records)):
            if mode=="similarity":matrix,translation=fit_similarity(m,reference,weights);detail={"status":"ok"}
            else:
                matrix,translation,detail=fit_constrained_affine(m,reference,weights,
                    anisotropy_limit=config["affine_anisotropy"],shear_limit=config["affine_shear"],penalty=config["affine_penalty"])
                if detail["status"]=="optimizer_failed":
                    matrix,translation=fit_similarity(m,reference,weights);detail.update(status="similarity_fallback")
            mapped.append(apply_transform(m,matrix,translation))
            mapped_pairs.append([apply_transform(w,matrix,translation) for w in p])
            mapped_geometry.append({k:apply_transform(w,matrix,translation) for k,w in g.items()})
            mode_transforms.append(dict(profile_id=local,source_profile_id=i,mode=mode,angle_deg=r["angle_deg"],side=r["side"],
                matrix_2x2=matrix.tolist(),translation_mm=translation.tolist(),principal_scales=np.linalg.svd(matrix,compute_uv=False).tolist(),
                determinant=float(np.linalg.det(matrix)),**detail))
        _write_curves(folder/f"{mode}_midlines_xy.ply",mapped,scale,(40,190,70))
        _write_verification(folder,mode,mapped_geometry,mapped_pairs,scale)
        report=export_section_models(folder,mode,mapped,mapped_pairs,records,scale,_write_curves,geometries=mapped_geometry,
            outlier_mad=config["outlier_mad"],bimodal=config["bimodal"],bimodal_min_profiles=config["bimodal_min_profiles"],bimodal_bic_delta=config["bimodal_bic_delta"])
        qa[mode]={"section_models":report,"output_subdir":f"standard_{mode}"}
        (folder/"transforms.json").write_text(json.dumps(mode_transforms,indent=2));transforms+=mode_transforms
    (out/"transforms.json").write_text(json.dumps(transforms,indent=2))
    (out/"rim_qa.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2))
    return Path(output)


def phase_summary(roots, output, step, mode, selection, interval_mm, plots, dense_reference=None):
    runs=[load_run(r,mode,selection) for r in roots]
    common=min(arc_positions(r["standard"])[-1] for r in runs)
    stations=np.arange(interval_mm,common+1e-9,interval_mm)
    if not len(stations):raise ValueError("No common phase stations")
    widths=[];coords=[]
    for run in runs:
        points,normals=station_frame(run["standard"],stations)
        widths.append(ray_width(points,normals,run["walls"]));coords.append(points)
    widths=np.array(widths);coords=np.array(coords);rows=[]
    for j,s in enumerate(stations):
        values=widths[:,j];values=values[np.isfinite(values)]
        # Pairwise geometry dispersion does not privilege one arbitrary phase.
        deltas=coords[:,None,j]-coords[None,:,j]
        rows.append(dict(station_mm=float(s),phase_count=len(runs),valid_thickness_phases=len(values),
            thickness_mean_mm=float(values.mean()) if len(values) else "",
            thickness_std_population_mm=float(values.std()) if len(values) else "",
            thickness_min_mm=float(values.min()) if len(values) else "",thickness_max_mm=float(values.max()) if len(values) else "",
            midline_max_pairwise_distance_mm=float(np.linalg.norm(deltas,axis=-1).max())))
    folder=Path(output)/f"summary_{step:02d}deg_{mode}_{selection}";folder.mkdir(parents=True,exist_ok=True)
    rows_csv(folder/"phase_station_statistics.csv",rows)
    walls=[r["walls"][k] for r in runs for k in ("outer","inner")]
    import colorsys
    colors=[tuple(int(c*255) for c in colorsys.hsv_to_rgb(i/len(runs),.8,.8)) for i in range(len(runs)) for _ in (0,1)]
    reference=load_run(dense_reference,mode,selection) if dense_reference is not None else None
    if reference is not None:
        walls.extend(reference["walls"][k] for k in ("outer","inner"));colors.extend([(0,0,0)]*2)
    write_colored_curves(folder/"phase_models_overlay.ply",walls,colors,runs[0]["scale"])
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,7))
        for i,r in enumerate(runs):
            for k in ("outer","inner"):ax.plot(*r["walls"][k].T,color=np.array(colors[2*i])/255,lw=1,alpha=.8,label=f"Phase {Path(roots[i]).name}" if k=="outer" else None)
        if reference is not None:
            for k in ("outer","inner"):ax.plot(*reference["walls"][k].T,color="black",lw=1.8,label="1-degree reference" if k=="outer" else None)
        ax.set(aspect="equal",xlabel="Projected radius [mm]",ylabel="Height [mm]",title=f"{step} deg / {mode}: {len(runs)} phase models")
        ax.legend(fontsize=7);fig.tight_layout();fig.savefig(folder/"phase_models_overlay.png",dpi=160);plt.close(fig)


def run_phase_experiment(root, output, *, steps=(5,10,15), phases=None, interval_mm=1., modes=("similarity","affine"), selection="all", plots=True):
    if not np.isfinite(interval_mm) or interval_mm<=0:raise ValueError("Interval must be finite and positive")
    if len(set(steps))!=len(steps) or any(not isinstance(s,int) or s<=0 or 180%s for s in steps):
        raise ValueError("Steps must be distinct positive integer divisors of 180")
    source=load_run(root,"similarity","all");out=Path(output);out.mkdir(parents=True,exist_ok=True)
    attempted=sorted({round(float(r["angle_deg"])%180,6) for r in source["qa"]["profiles"]})
    if len(attempted)!=180 or not np.allclose(np.diff(attempted),1,atol=1e-6):
        raise ValueError("Phase experiment requires all 180 attempted planes on a 1-degree grid (rim_qa.json)")
    origin=attempted[0]; records=source["report"]["profiles"];manifest=[]
    for mode in modes:
        validate_run(source["root"],out/"dense_reference_validation",interval_mm=interval_mm,mode=mode,selection=selection,plots=plots)
    for step in steps:
        offsets=list(range(step)) if phases is None else list(phases)
        if len(set(offsets))!=len(offsets) or any(o<0 or o>=step for o in offsets):raise ValueError("Phases must be distinct integers in [0,step)")
        roots=[]
        for offset in offsets:
            ids=phase_indices(records,step,offset,origin)
            run=out/f"interval_{step:02d}deg"/f"phase_{offset:02d}deg"
            rebuild_subset(source,ids,run,attempted_profiles=360//step);roots.append(run)
            manifest.append(dict(interval_deg=step,phase_deg=offset,origin_deg=origin,expected_planes=180//step,
                expected_half_profiles=360//step,accepted_half_profiles=len(ids),source_profile_ids=ids,output=str(run),
                source_angles_deg=sorted({float(records[i]["angle_deg"]) for i in ids})))
            for mode in modes:validate_run(run,run/"validation",interval_mm=interval_mm,mode=mode,selection=selection,plots=plots)
        for mode in modes:
            compare_runs([*roots,source["root"]],out/f"comparison_{step:02d}deg",interval_mm=interval_mm,mode=mode,selection=selection)
            phase_summary(roots,out,step,mode,selection,interval_mm,plots,dense_reference=source["root"])
    (out/"phase_experiment.json").write_text(json.dumps(dict(source=str(source["root"].resolve()),
        dense_reference="1-degree full model; not independent truth",intervals=list(steps),phase_policy="all integer phases" if phases is None else "specified phases",sets=manifest,
        fixed=["axis","source intersection walls","lip tip","per-profile crop","endpoint policy"],
        recomputed=["right-side median reference","similarity and affine transforms","outlier filtering","standard midline and section models"]),indent=2))
    return manifest

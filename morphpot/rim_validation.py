"""Development-only validation of saved rim models against mesh intersections.

These are geometric mesh measurements, not independent physical measurements.
No paired-point projection fallback is used in validation.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from .rim_models import _ray_distances
from .rim_standardization import arc_positions, sample_curve, _write_curves

VERSION = "0.2.2-dev"


def read_curves(path, scale):
    """Read this project's ASCII edge-PLY; enforce ordered, disjoint paths."""
    with Path(path).open(encoding="ascii") as f:
        if f.readline().strip() != "ply" or f.readline().strip() != "format ascii 1.0":
            raise ValueError(f"Expected exporter ASCII PLY: {path}")
        nv = ne = None
        for line in f:
            fields = line.split()
            if fields[:2] == ["element", "vertex"]: nv = int(fields[2])
            if fields[:2] == ["element", "edge"]: ne = int(fields[2])
            if line.strip() == "end_header": break
        if nv is None or ne is None: raise ValueError("PLY lacks vertex/edge counts")
        vertices = np.array([[float(x) for x in f.readline().split()[:2]] for _ in range(nv)])*scale
        edges = [tuple(map(int, f.readline().split())) for _ in range(ne)]
    starts = [0]
    for a, b in edges:
        if b != a+1: raise ValueError("Validation requires ordered open source paths")
    linked = {a for a, _ in edges}
    starts += [i+1 for i in range(nv-1) if i not in linked]
    return [vertices[a:b] for a, b in zip(starts, starts[1:]+[nv])]


def rows_csv(path, rows):
    if not rows: return
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def station_frame(curve, positions):
    arc = arc_positions(curve)
    if arc[-1] <= 0: raise ValueError("Degenerate curve")
    points = sample_curve(curve, positions)
    tangents = np.gradient(curve, arc, axis=0)
    tangent = np.column_stack([
        np.interp(positions, arc, tangents[:, k]) for k in range(2)])
    norm = np.linalg.norm(tangent, axis=1)
    if np.any(norm <= 1e-10): raise ValueError("Degenerate station tangent")
    tangent /= norm[:,None]
    return points, np.c_[-tangent[:,1], tangent[:,0]]


def ray_width(points, normal, geometry):
    a = _ray_distances(points, normal, geometry)
    b = _ray_distances(points, -normal, geometry)
    return np.where(np.isfinite(a)&np.isfinite(b), a+b, np.nan)


def write_overlay(path, geometry, model, scale):
    """One colored edge PLY; never create connecting edges between walls."""
    write_colored_curves(path,[geometry["outer"],geometry["inner"],model["outer"],model["inner"]],
                         [(150,150,150)]*2+[(220,50,50)]*2,scale)


def write_colored_curves(path, curves, colors, scale):
    with Path(path).open("w",encoding="ascii") as f:
        f.write(f"ply\nformat ascii 1.0\nelement vertex {sum(map(len,curves))}\nproperty double x\nproperty double y\nproperty double z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement edge {sum(len(c)-1 for c in curves)}\nproperty int vertex1\nproperty int vertex2\nend_header\n")
        for curve,color in zip(curves,colors):
            for x,y in curve/scale:f.write(f"{x:.15g} {y:.15g} 0 {' '.join(map(str,color))}\n")
        offset=0
        for curve in curves:
            for j in range(len(curve)-1):f.write(f"{offset+j} {offset+j+1}\n")
            offset+=len(curve)


def load_run(root, mode, selection="all"):
    root = Path(root)
    if root.name != "rim_standardization": root /= "rim_standardization"
    qa = json.loads((root/"rim_qa.json").read_text())
    scale = float(qa["unit_to_mm"])
    folder = root/f"standard_{mode}"
    transforms = json.loads((folder/"transforms.json").read_text())
    report = json.loads((folder/f"{mode}_section_models.json").read_text())
    model = report["models"][selection]
    stem = f"standard_{mode}_section_{selection}"
    mids = read_curves(folder/f"{mode}_midlines_xy.ply", scale)
    outer = read_curves(folder/f"{mode}_source_outer_xy.ply", scale)
    inner = read_curves(folder/f"{mode}_source_inner_xy.ply", scale)
    raw_mids = read_curves(root/"raw_midlines_xy.ply", scale)
    raw_outer = read_curves(root/"raw_source_outer_xy.ply", scale)
    raw_inner = read_curves(root/"raw_source_inner_xy.ply", scale)
    counts = [len(x) for x in (mids,outer,inner,raw_mids,raw_outer,raw_inner,transforms)]
    if len(set(counts)) != 1: raise ValueError(f"Profile ordering/count mismatch: {counts}")
    if any(t["profile_id"]!=i for i,t in enumerate(transforms)):
        raise ValueError("Transform profile IDs are not in exporter order")
    standard = read_curves(folder/f"{stem}_midline_xy.ply",scale)[0]
    walls = {k:read_curves(folder/f"{stem}_{k}_xy.ply",scale)[0] for k in ("outer","inner")}
    return dict(root=root,qa=qa,scale=scale,transforms=transforms,report=report,ids=model["profile_ids"],
                mids=mids,outer=outer,inner=inner,raw_mids=raw_mids,raw_outer=raw_outer,raw_inner=raw_inner,
                standard=standard,walls=walls,mode=mode,selection=selection)


def validate_run(root, output, *, interval_mm=1., mode="similarity", selection="all", plots=True):
    if not np.isfinite(interval_mm) or interval_mm <= 0: raise ValueError("Interval must be positive")
    run=load_run(root,mode,selection); out=Path(output)/f"standard_{mode}_{selection}"; out.mkdir(parents=True,exist_ok=True)
    standard=run["standard"]; length=arc_positions(standard)[-1]
    stations=np.arange(interval_mm,length+1e-9,interval_mm)  # lip has no two-wall thickness
    if not len(stations): raise ValueError("Interval exceeds standard midline length")
    model_points,model_normals=station_frame(standard,stations)
    model_width=ray_width(model_points,model_normals,run["walls"])
    all_walls=[run[k][i] for i in run["ids"] for k in ("outer","inner")]
    model_walls=[run["walls"][k] for k in ("outer","inner")]
    write_colored_curves(out/"profile_all_overlay.ply",all_walls+model_walls,
                         [(150,150,150)]*len(all_walls)+[(220,50,50)]*2,run["scale"])
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,7))
        for index,wall in enumerate(all_walls):
            ax.plot(*wall.T,color=".5",alpha=.35,lw=.6,label="Mesh sections" if index==0 else None)
        for index,wall in enumerate(model_walls):
            ax.plot(*wall.T,color="crimson",lw=2,label="Standard model" if index==0 else None)
        ax.scatter(*model_points.T,s=10,c="black",label=f"Stations ({interval_mm:g} mm)")
        ax.set(aspect="equal",xlabel="Projected radius [mm]",ylabel="Height [mm]",
               title=f"{mode} / {selection}: {len(run['ids'])} profiles + model")
        ax.legend();fig.tight_layout();fig.savefig(out/"profile_all_overlay.png",dpi=160);plt.close(fig)
    rows=[]; profile_stats=[]
    for i in run["ids"]:
        curve=run["mids"][i]; arc=arc_positions(curve); raw=run["raw_mids"][i]; raw_arc=arc_positions(raw)
        tr=run["transforms"][i]; geometry={k:run[k][i] for k in ("outer","inner")}
        raw_geometry={k:run["raw_"+k][i] for k in ("outer","inner")}
        ugrid=np.linspace(0,1,len(curve))
        for basis in ("absolute_mm","corresponding_u"):
            if basis=="absolute_mm": target=stations
            else:
                u=np.interp(stations,arc_positions(standard),np.linspace(0,1,len(standard)))
                target=np.interp(u,ugrid,arc)
            valid=target<=arc[-1]+1e-9
            positions=np.minimum(target,arc[-1]); points,normals=station_frame(curve,positions)
            widths=ray_width(points,normals,geometry); widths[~valid]=np.nan
            # Inverse correspondence follows the original sample index, not the
            # transformed mm scale. Raw normal is recomputed, essential for affine.
            u=np.interp(positions,arc,ugrid)
            raw_s=np.interp(u,np.linspace(0,1,len(raw)),raw_arc)
            raw_points,raw_normals=station_frame(raw,raw_s)
            raw_width=ray_width(raw_points,raw_normals,raw_geometry); raw_width[~valid]=np.nan
            for j,s in enumerate(stations):
                w=widths[j]; mw=model_width[j]
                rows.append(dict(profile_id=i,angle_deg=tr["angle_deg"],side=tr["side"],basis=basis,
                    station_mm=float(s),profile_arc_mm=float(target[j]),profile_u=float(u[j]),
                    model_thickness_mm=float(mw) if np.isfinite(mw) else "",
                    measured_transformed_mm=float(w) if np.isfinite(w) else "",
                    measured_original_mm=float(raw_width[j]) if np.isfinite(raw_width[j]) else "",
                    thickness_error_mm=float(w-mw) if np.isfinite(w+mw) else "",
                    midline_distance_mm=float(np.linalg.norm(points[j]-model_points[j])) if valid[j] else "",
                    status="outside_profile" if not valid[j] else "normal_ray" if np.isfinite(w) else "missing_intersection",
                    retained=run["report"]["profiles"][i]["retained"]))
        comparable=sample_curve(curve,np.interp(np.linspace(0,1,len(standard)),ugrid,arc))
        profile_stats.append(dict(profile_id=i,angle_deg=tr["angle_deg"],side=tr["side"],
            transformed_arc_mm=float(arc[-1]),original_arc_mm=float(raw_arc[-1]),
            midline_rms_corresponding_mm=float(np.sqrt(np.mean(np.sum((comparable-standard)**2,axis=1)))),
            tip_distance_mm=float(np.linalg.norm(curve[0]-standard[0])),
            end_distance_mm=float(np.linalg.norm(curve[-1]-standard[-1])),
            principal_scale_min=min(tr["principal_scales"]),principal_scale_max=max(tr["principal_scales"])))
        # Original intersection polylines are preserved, including detection margin.
        ply=out/f"profile_{i:03d}_overlay_xy.ply"
        write_overlay(ply,geometry,run["walls"],run["scale"])
        _write_curves(out/f"profile_{i:03d}_model_xy.ply",[run["walls"]["outer"],run["walls"]["inner"]],run["scale"],(220,50,50))
        if plots:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig,ax=plt.subplots(figsize=(7,6))
            for k in ("outer","inner"):
                ax.plot(*geometry[k].T,color=".55",lw=1,label="Mesh section" if k=="outer" else None)
                ax.plot(*run["walls"][k].T,color="crimson",lw=1.5,label="Standard model" if k=="outer" else None)
            ax.scatter(*model_points.T,s=10,c="black",label=f"Stations ({interval_mm:g} mm)")
            ax.set(aspect="equal",xlabel="Projected radius [mm]",ylabel="Height [mm]",title=f"{mode} / profile {i}: {tr['angle_deg']:g} deg {tr['side']}")
            ax.legend();fig.tight_layout();fig.savefig(out/f"profile_{i:03d}_overlay.png",dpi=150);plt.close(fig)
    rows_csv(out/"station_measurements.csv",rows); rows_csv(out/"profile_geometry.csv",profile_stats)
    stats=[]
    for basis in ("absolute_mm","corresponding_u"):
        for s in stations:
            subset=[r for r in rows if r["basis"]==basis and r["station_mm"]==s]
            values=np.array([r["measured_transformed_mm"] for r in subset if r["measured_transformed_mm"]!=""])
            errors=np.array([r["thickness_error_mm"] for r in subset if r["thickness_error_mm"]!=""])
            stats.append(dict(basis=basis,station_mm=float(s),n_requested=len(subset),n_valid=len(values),
                model_thickness_mm=subset[0]["model_thickness_mm"],
                mean_mm=float(values.mean()) if len(values) else "",min_mm=float(values.min()) if len(values) else "",
                max_mm=float(values.max()) if len(values) else "",std_population_mm=float(values.std()) if len(values) else "",
                bias_mm=float(errors.mean()) if len(errors) else "",rmse_mm=float(np.sqrt(np.mean(errors**2))) if len(errors) else ""))
    rows_csv(out/"station_statistics.csv",stats)
    meta=dict(version=VERSION,source=str(run["root"].resolve()),mode=mode,selection=selection,interval_mm=interval_mm,
        profile_count=len(run["ids"]),standard_arc_mm=float(length),
        summary_overlay="profile_all_overlay.ply",summary_png="profile_all_overlay.png" if plots else None,
        definitions={"absolute_mm":"same arc length from each transformed lip", "corresponding_u":"same normalized sample correspondence; stations spaced on standard arc", "original":"original units converted to mm at inverse sample correspondence", "std":"population ddof=0", "missing":"no projection fallback or extrapolation"},
        limitations=["mesh intersections are not independent physical ground truth", "model uses these same profiles; in-sample descriptive validation", "source walls include detection margin past comparison endpoint", "profile overlay and model PLY share coordinates and input units", "no restoration/original surface labels available"])
    (out/"validation.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2))
    return run


def compare_runs(roots, output, *, interval_mm=1.,mode="similarity",selection="all"):
    """Last run is explicit reference, never silently picked by accepted count."""
    runs=[load_run(r,mode,selection) for r in roots]; reference=runs[-1]; rows=[]
    for run in runs[:-1]:
        common=min(arc_positions(run["standard"])[-1],arc_positions(reference["standard"])[-1])
        stations=np.arange(interval_mm,common+1e-9,interval_mm)
        if not len(stations): raise ValueError("No common non-tip stations")
        a,na=station_frame(run["standard"],stations); b,nb=station_frame(reference["standard"],stations)
        wa=ray_width(a,na,run["walls"]); wb=ray_width(b,nb,reference["walls"])
        for j,s in enumerate(stations):
            rows.append(dict(run=str(run["root"]),reference=str(reference["root"]),mode=mode,station_mm=float(s),
                profile_count=len(run["ids"]),reference_profile_count=len(reference["ids"]),
                midline_distance_mm=float(np.linalg.norm(a[j]-b[j])),
                tip_aligned_midline_distance_mm=float(np.linalg.norm((a[j]-run["standard"][0])-(b[j]-reference["standard"][0]))),
                thickness_delta_mm=float(wa[j]-wb[j]) if np.isfinite(wa[j]+wb[j]) else "",
                tangent_angle_deg=float(np.rad2deg(np.arccos(np.clip(na[j]@nb[j],-1,1)))),
                endpoint_mode=run["qa"]["config"].get("end_mode"),reference_endpoint_mode=reference["qa"]["config"].get("end_mode")))
    out=Path(output);out.mkdir(parents=True,exist_ok=True); rows_csv(out/f"{mode}_{selection}_division_comparison.csv",rows)
    summary=[]
    for run in runs[:-1]:
        group=[r for r in rows if r["run"]==str(run["root"])]
        errors=[r["thickness_delta_mm"] for r in group if r["thickness_delta_mm"]!=""]
        summary.append(dict(run=str(run["root"]),reference=str(reference["root"]),
            accepted_profiles=len(run["ids"]),reference_accepted_profiles=len(reference["ids"]),
            common_stations=len(group),valid_thickness_stations=len(errors),
            midline_rms_mm=float(np.sqrt(np.mean([r["midline_distance_mm"]**2 for r in group]))),
            midline_max_mm=max(r["midline_distance_mm"] for r in group),
            thickness_rms_delta_mm=float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
            thickness_max_abs_delta_mm=max(map(abs,errors)) if errors else None,
            arc_length_delta_mm=float(arc_positions(run["standard"])[-1]-arc_positions(reference["standard"])[-1]),
            endpoint_distance_mm=float(np.linalg.norm(run["standard"][-1]-reference["standard"][-1])),
            run_config=run["qa"]["config"],reference_config=reference["qa"]["config"]))
    (out/f"{mode}_{selection}_division_summary.json").write_text(json.dumps(dict(
        comparisons=summary,reference_policy="last supplied run",limitations=[
            "A denser division is a reference, not independent truth or proof of convergence",
            "Change angular step only for sampling comparison; inspect endpoint/config differences",
            "Coordinates include sampling-dependent registration reference; tip-aligned differences also exported",
            "Inspect angular exclusion and start-angle sensitivity; no automatic unchanged threshold"]),indent=2))
    return rows

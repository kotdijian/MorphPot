"""Standalone curvature/turning validation from saved standard-model CSVs.

No mesh extraction, registration or standard model construction is repeated.
Dependencies: NumPy, SciPy and Matplotlib; no morphpot imports.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.signal import savgol_filter

VERSION = "0.1.0"


def curve_from_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if not {"mid_x_mm","mid_y_mm"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing standard midline mm columns: {path}")
        rows=list(reader)
    curve=np.array([[float(r["mid_x_mm"]),float(r["mid_y_mm"])] for r in rows])
    if len(curve)<5 or not np.all(np.isfinite(curve)):raise ValueError(f"Invalid/short midline: {path}")
    keep=np.r_[True,np.linalg.norm(np.diff(curve,axis=0),axis=1)>1e-10]
    return curve[keep]


def curve_metrics(curve, grid_mm=.25, smooth_mm=2.):
    if not np.isfinite(grid_mm) or grid_mm<=0 or not np.isfinite(smooth_mm) or smooth_mm<=0:
        raise ValueError("Grid and smoothing width must be finite and positive")
    curve=np.asarray(curve,float)
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(curve,axis=0),axis=1))]
    if len(curve)<5 or np.any(np.diff(arc)<=0):raise ValueError("Curve requires >=5 distinct ordered points")
    stations=np.arange(0,arc[-1]+1e-9,grid_mm)
    window=int(round(smooth_mm/grid_mm))+1
    if window%2==0:window+=1
    if window<5:raise ValueError("Smoothing window must span at least four grid intervals")
    if len(stations)<window+2:raise ValueError("Curve is too short for the chosen smoothing window")
    positions=np.column_stack([np.interp(stations,arc,curve[:,k]) for k in (0,1)])
    # Same physical window and polynomial degree in every model. Interpolated
    # endpoint fits are one-sided and flagged rather than silently trusted.
    first=savgol_filter(positions,window,3,deriv=1,delta=grid_mm,axis=0,mode="interp")
    second=savgol_filter(positions,window,3,deriv=2,delta=grid_mm,axis=0,mode="interp")
    speed=np.linalg.norm(first,axis=1)
    if np.any(speed<1e-9):raise ValueError("Degenerate smoothed tangent")
    theta=np.unwrap(np.arctan2(first[:,1],first[:,0]))
    curvature=(first[:,0]*second[:,1]-first[:,1]*second[:,0])/speed**3
    signed=np.rad2deg(theta-theta[0])
    absolute=np.r_[0,np.cumsum(np.abs(np.diff(np.rad2deg(theta))))]
    half=window//2;boundary=np.zeros(len(stations),bool);boundary[:half]=True;boundary[-half:]=True
    return dict(stations=stations,curvature=curvature,signed_turn=signed,absolute_turn=absolute,
                boundary=boundary,window=window,effective_smooth_mm=(window-1)*grid_mm)


def write_rows(path,rows):
    with Path(path).open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def model_csv(root,mode,selection):
    root=Path(root)
    if root.name!="rim_standardization":root/= "rim_standardization"
    return root/f"standard_{mode}"/f"standard_{mode}_section_{selection}.csv"


def create_curvature_outputs(root, output=None, *, reference_dir=None, steps=None, phases=None,
                             modes=("similarity","affine"),selection="all",grid_mm=.25,smooth_mm=2.,dpi=180):
    root=Path(root);manifest=json.loads((root/"phase_experiment.json").read_text())
    steps=list(manifest["intervals"]) if steps is None else list(steps)
    if not steps or len(set(steps))!=len(steps) or any(s<=0 for s in steps):raise ValueError("Distinct positive interval steps required")
    if phases is not None and (len(set(phases))!=len(phases) or any(p<0 for p in phases)):
        raise ValueError("Phase offsets must be distinct nonnegative integers")
    reference=Path(reference_dir) if reference_dir is not None else Path(manifest["source"])
    available={step:sorted(x["phase_deg"] for x in manifest["sets"] if x["interval_deg"]==step) for step in steps}
    inputs=[]
    for step in steps:
        selected=available[step] if phases is None else list(phases)
        if not selected or not set(selected).issubset(available[step]):raise ValueError(f"Requested phases missing at {step} degrees")
        for phase in selected:inputs.append((step,phase,root/f"interval_{step:02d}deg"/f"phase_{phase:02d}deg"))
    # Load/derive every requested mode before writing a partial output set.
    prepared={}
    for mode in modes:
        ref_path=model_csv(reference,mode,selection)
        if not ref_path.exists():raise ValueError(f"Reference model not found: {ref_path}; specify --reference-dir with the original 1-degree output")
        ref=curve_metrics(curve_from_csv(ref_path),grid_mm,smooth_mm)
        derived=[(step,phase,curve_metrics(curve_from_csv(model_csv(folder,mode,selection)),grid_mm,smooth_mm))
                 for step,phase,folder in inputs]
        n=min([len(ref["stations"]),*[len(m["stations"]) for _,_,m in derived]])
        stations=ref["stations"][:n];half=ref["window"]//2
        boundary=np.zeros(n,bool);boundary[:half]=True;boundary[-half:]=True
        interior=np.flatnonzero(~boundary)
        if not len(interior):raise ValueError("No common interior stations beyond smoothing boundary zones")
        prepared[mode]=(ref,derived,stations,boundary,interior)
    output=Path(output) if output is not None else root/"curvature_turning"
    output.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    palette=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"]
    outputs=[]
    for mode,(ref,derived,stations,boundary,interior) in prepared.items():
        n=len(stations);rows=[];summary=[]
        def add_curve(step,phase,metrics):
            k=metrics["curvature"][:n];turn=metrics["signed_turn"][:n];absolute=metrics["absolute_turn"][:n]
            peak=interior[np.argmax(np.abs(k[interior]))]
            ref_peak=float(np.max(np.abs(ref["curvature"][interior])))
            for j,s in enumerate(stations):
                rows.append(dict(interval_deg=step,phase_deg=phase,station_mm=float(s),boundary_affected=bool(boundary[j]),
                    signed_curvature_per_mm=float(k[j]),signed_turn_deg=float(turn[j]),absolute_turn_deg=float(absolute[j]),
                    curvature_delta_per_mm=float(k[j]-ref["curvature"][j]),
                    signed_turn_delta_deg=float(turn[j]-ref["signed_turn"][j]),absolute_turn_delta_deg=float(absolute[j]-ref["absolute_turn"][j])))
            summary.append(dict(interval_deg=step,phase_deg=phase,comparison_end_mm=float(stations[-1]),
                total_signed_turn_deg=float(turn[-1]),total_absolute_turn_deg=float(absolute[-1]),
                peak_abs_curvature_per_mm=float(abs(k[peak])),peak_arc_mm=float(stations[peak]),
                peak_ratio_to_reference=float(abs(k[peak])/ref_peak) if ref_peak>1e-10 else "",
                curvature_rms_delta_per_mm=float(np.sqrt(np.mean((k[interior]-ref["curvature"][interior])**2))),
                signed_turn_end_delta_deg=float(turn[-1]-ref["signed_turn"][n-1]),
                absolute_turn_end_delta_deg=float(absolute[-1]-ref["absolute_turn"][n-1])))
        add_curve(1,"reference",ref)
        for step,phase,metrics in derived:add_curve(step,phase,metrics)
        fig,axes=plt.subplots(2,3,figsize=(15,8),sharex=True)
        fields=("signed_turn","curvature","absolute_turn")
        titles=("Tangent rotation from lip","Signed curvature","Cumulative absolute turning")
        units=("Angle [deg]","Curvature [1/mm]","Angle [deg]")
        for col,(field,title,unit) in enumerate(zip(fields,titles,units)):
            axes[0,col].plot(stations,ref[field][:n],color="black",lw=1.6,label="1-degree reference")
            for index,step in enumerate(sorted(steps)):
                values=np.stack([m[field][:n] for s,_,m in derived if s==step])
                color=palette[index%len(palette)];label=f"{step} deg ({len(values)} phases)"
                for row in (0,1):
                    data=values if row==0 else values-ref[field][:n]
                    axes[row,col].plot(stations,np.median(data,axis=0),color=color,lw=1.5,label=label)
                    axes[row,col].fill_between(stations,data.min(axis=0),data.max(axis=0),color=color,alpha=.15)
            for row in (0,1):
                ax=axes[row,col];ax.set(title=title if row==0 else title+"\nDifference from reference",ylabel=unit)
                ax.axhline(0,color=".4",ls="--",lw=.6);ax.grid(alpha=.2)
                ax.axvspan(stations[0],stations[ref["window"]//2],color=".5",alpha=.07)
                ax.axvspan(stations[-(ref["window"]//2)-1],stations[-1],color=".5",alpha=.07)
                ax.ticklabel_format(axis="y",style="sci",scilimits=(-3,4),useOffset=False)
            axes[1,col].set_xlabel("Arc length from lip [mm]")
        axes[0,0].legend(fontsize=8)
        fig.suptitle(f"Curvature and turning / {mode} / {selection}",fontsize=15)
        fig.text(.5,.025,f"Grid {grid_mm:g} mm; cubic Savitzky-Golay window {ref['effective_smooth_mm']:g} mm. Line: phase median; band: min-max, not CI.\n"
                 "Gray zones: smoothing boundary effects; excluded from curvature peak/RMS. Turning totals include endpoint estimates. Reference is not ground truth.",ha="center",fontsize=9)
        fig.tight_layout(rect=(0,.09,1,.94))
        stem=f"curvature_turning_{mode}_{selection}"
        png=output/f"{stem}.png";fig.savefig(png,dpi=dpi);plt.close(fig)
        write_rows(output/f"{stem}_stations.csv",rows);write_rows(output/f"{stem}_summary.csv",summary)
        (output/f"{stem}.json").write_text(json.dumps(dict(version=VERSION,reference=str(reference),grid_mm=grid_mm,
            requested_smooth_mm=smooth_mm,effective_smooth_mm=ref["effective_smooth_mm"],polynomial_degree=3,
            comparison_end_mm=float(stations[-1]),intervals=steps,phases=available if phases is None else phases,
            definitions={"curvature":"signed (x'y''-y'x'')/(x'^2+y'^2)^(3/2); positive is counterclockwise along lip-to-body traversal",
                "signed_turn":"unwrapped tangent angle minus lip tangent angle; rotation invariant",
                "absolute_turn":"sum of absolute changes of unwrapped smoothed tangent angle; captures reverse bends",
                "peak":"maximum absolute curvature in common interior; location and ratio refer to whole-window maxima, not matched anatomical peaks"},
            limitations=["common original midline arc grid; no extrapolation", "fixed smoothing changes small-scale curvature; test window sensitivity",
                "turning totals include one-sided endpoint fits; boundary flag provided", "phase extrema depend on phase count", "no inference of bend loss from angle differences alone"]),indent=2))
        outputs.append(png)
    return outputs


def main(argv=None):
    parser=argparse.ArgumentParser(description="Generate only additional curvature/turning PNG, CSV and JSON from saved phase models")
    parser.add_argument("input_dir",type=Path,help="Existing phase_validation directory")
    parser.add_argument("--reference-dir",type=Path,help="Original 1-degree output; override manifest source path")
    parser.add_argument("--output-dir",type=Path,help="Default INPUT/curvature_turning")
    parser.add_argument("--steps",type=int,nargs="+")
    parser.add_argument("--phases",type=int,nargs="+",help="Restrict to common offsets, e.g. 0 1 2 3 4")
    parser.add_argument("--mode",choices=("similarity","affine","both"),default="both")
    parser.add_argument("--selection",choices=("all","inliers","A","B"),default="all")
    parser.add_argument("--grid-mm",type=float,default=.25)
    parser.add_argument("--smooth-mm",type=float,default=2.,help="Physical full smoothing window width; rounded to an odd sample count")
    parser.add_argument("--dpi",type=int,default=180)
    parser.add_argument("--version",action="version",version=VERSION)
    args=parser.parse_args(argv)
    if args.dpi<=0:parser.error("--dpi must be positive")
    try:
        outputs=create_curvature_outputs(args.input_dir,args.output_dir,reference_dir=args.reference_dir,steps=args.steps,
            phases=args.phases,modes=("similarity","affine") if args.mode=="both" else (args.mode,),selection=args.selection,
            grid_mm=args.grid_mm,smooth_mm=args.smooth_mm,dpi=args.dpi)
    except (OSError,ValueError,KeyError) as exc:parser.exit(2,f"ERROR: {exc}\n")
    for output in outputs:print(output)


if __name__=="__main__":main()

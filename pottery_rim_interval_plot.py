"""Plot interval sensitivity from existing comparison CSVs; no model refitting.

Standalone: this file needs only Python, NumPy and Matplotlib. Copying it outside
MorphPot is supported. The reference paths stored in CSVs are labels, not opened.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import re
import numpy as np

VERSION = "0.1.0"
METRICS = (
    ("thickness_delta_mm", "Thickness difference: interval model - reference", "Thickness difference [mm]"),
    ("midline_distance_mm", "Midline position difference", "Distance [mm]"),
    ("tip_aligned_midline_distance_mm", "Midline difference after aligning lip tips", "Distance [mm]"),
    ("tangent_angle_deg", "Midline tangent direction difference", "Angle [deg]"),
)


def load_comparisons(root, steps, mode, selection):
    root = Path(root)
    if steps is None:
        steps = sorted({int(m.group(1)) for p in root.glob("comparison_*deg")
                        if (m := re.fullmatch(r"comparison_(\d+)deg", p.name))})
    if len(steps) < 2 or len(set(steps)) != len(steps) or any(s <= 0 for s in steps):
        raise ValueError("At least two distinct positive interval steps are required")
    tables = {}; references = set(); reference_counts = set()
    for step in steps:
        path = root / f"comparison_{step:02d}deg" / f"{mode}_{selection}_division_comparison.csv"
        with path.open(encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            required = {"run", "reference", "reference_profile_count", "mode", "station_mm", *[m[0] for m in METRICS]}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"Missing comparison columns: {path}")
            rows = list(reader)
        if not rows: raise ValueError(f"Empty comparison table: {path}")
        keys = set()
        for row in rows:
            if row["mode"] != mode: raise ValueError(f"Mode mismatch in {path}")
            key = (row["run"], round(float(row["station_mm"]), 9))
            if key in keys: raise ValueError(f"Duplicate phase/station record in {path}: {key}")
            keys.add(key); references.add(row["reference"]); reference_counts.add(row["reference_profile_count"])
        tables[step] = rows
    if len(references) != 1 or len(reference_counts) != 1:
        raise ValueError("All intervals must use the same reference model and reference profile count")
    # Compare only recorded common station positions; never interpolate or extend
    # a short phase model. Available phase count may differ near comparison ends.
    common = set.intersection(*[{round(float(r["station_mm"]),9) for r in rows} for rows in tables.values()])
    if not common: raise ValueError("No common recorded stations across intervals")
    return tables, np.array(sorted(common)), next(iter(references))


def summarize_metric(rows, stations, field):
    groups = {}
    for row in rows:
        value = row[field]
        if value == "": continue
        number = float(value)
        if np.isfinite(number): groups.setdefault(round(float(row["station_mm"]),9), []).append(number)
    summaries = []
    for station in stations:
        values = np.array(groups.get(round(float(station),9),[]))
        summaries.append((float(np.median(values)),float(values.min()),float(values.max()),len(values))
                         if len(values) else (np.nan,np.nan,np.nan,0))
    return np.array(summaries)


def create_interval_plots(root, output=None, *, steps=None, modes=("similarity","affine"), selection="all", dpi=180):
    """Create only PNG files, one four-panel comparison per requested mode."""
    if dpi <= 0: raise ValueError("DPI must be positive")
    root = Path(root); output = Path(output) if output is not None else root / "interval_dependence"
    # Resolve all sources before writing, so a missing affine table does not
    # silently produce only part of an explicitly requested both-mode result.
    sources = {mode:load_comparisons(root,steps,mode,selection) for mode in modes}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True,exist_ok=True); paths=[]
    palette=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"]
    for mode,(tables,stations,reference) in sources.items():
        fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
        for ax,(field,title,ylabel) in zip(axes.flat,METRICS):
            for index,(step,rows) in enumerate(sorted(tables.items())):
                values=summarize_metric(rows,stations,field)
                color=palette[index%len(palette)]
                phase_count=len({r["run"] for r in rows})
                ax.plot(stations,values[:,0],color=color,lw=1.6,marker="o",ms=3,
                        label=f"{step} deg ({phase_count} phases)")
                ax.fill_between(stations,values[:,1],values[:,2],color=color,alpha=.15)
            ax.axhline(0,color=".35",lw=.7,ls="--")
            ax.set(title=title,ylabel=ylabel);ax.grid(alpha=.2)
            ax.ticklabel_format(axis="y",style="sci",scilimits=(-3,4),useOffset=False)
        for ax in axes[1]: ax.set_xlabel("Arc length from reference lip [mm]")
        axes[0,0].legend(fontsize=9)
        fig.suptitle(f"Interval sensitivity / {mode} / {selection}",fontsize=15)
        fig.text(.5,.025,"Line: phase median. Band: phase min-max (not a confidence interval).\n"
                 "Only common recorded stations; missing intersections omitted; no extrapolation. Reference is not ground truth.",
                 ha="center",fontsize=9)
        fig.tight_layout(rect=(0,.085,1,.94))
        path=output/f"interval_dependence_{mode}_{selection}.png"
        fig.savefig(path,dpi=dpi);plt.close(fig);paths.append(path)
    return paths


def main(argv=None):
    parser=argparse.ArgumentParser(description="Generate only interval-dependence PNGs from saved phase-validation CSVs")
    parser.add_argument("input_dir",type=Path,help="phase_validation directory containing comparison_05deg etc.")
    parser.add_argument("--output-dir",type=Path,help="PNG destination (default INPUT/interval_dependence)")
    parser.add_argument("--steps",type=int,nargs="+",help="Intervals to compare (default discover existing comparison directories)")
    parser.add_argument("--mode",choices=("similarity","affine","both"),default="both")
    parser.add_argument("--selection",choices=("all","inliers","A","B"),default="all")
    parser.add_argument("--dpi",type=int,default=180)
    parser.add_argument("--version",action="version",version=VERSION)
    args=parser.parse_args(argv)
    modes=("similarity","affine") if args.mode=="both" else (args.mode,)
    try:
        paths=create_interval_plots(args.input_dir,args.output_dir,steps=args.steps,modes=modes,selection=args.selection,dpi=args.dpi)
    except (OSError,ValueError) as exc:parser.exit(2,f"ERROR: {exc}\n")
    for path in paths:print(path)


if __name__=="__main__":main()

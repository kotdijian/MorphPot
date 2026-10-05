"""Development validation CLI, separate from release-oriented section extraction."""
import argparse
from pathlib import Path
import subprocess
import sys
from morphpot.rim_validation import validate_run, compare_runs, VERSION
from morphpot.rim_phase_validation import run_phase_experiment


def main(argv=None):
    parser=argparse.ArgumentParser(description="Development validation of rim thickness, mesh overlays and angular sampling")
    parser.add_argument("--version",action="version",version=VERSION)
    sub=parser.add_subparsers(dest="command",required=True)
    for command in ("validate","compare","sweep","phases"):
        p=sub.add_parser(command)
        p.add_argument("--output-dir",type=Path,required=True)
        p.add_argument("--interval-mm",type=float,default=1.)
        p.add_argument("--mode",choices=["similarity","affine","both"],default="both")
        p.add_argument("--selection",choices=["all","inliers","A","B"],default="all")
        p.add_argument("--no-plots",action="store_true")
        if command=="sweep":
            p.add_argument("mesh",type=Path);p.add_argument("--steps",type=float,nargs="+",default=[30,15,10,5])
        elif command=="phases":
            p.add_argument("run",type=Path,help="1-degree section output directory")
            p.add_argument("--steps",type=int,nargs="+",default=[5,10,15])
            p.add_argument("--phases",type=int,nargs="+",help="Explicit offsets; default all phases per interval")
        else: p.add_argument("runs",type=Path,nargs="+",help="Section output directories; compare uses last as reference")
    args,extra=parser.parse_known_args(argv)
    if args.interval_mm<=0: parser.error("--interval-mm must be positive")
    if args.command!="sweep" and extra:parser.error(f"Unknown arguments: {extra}")
    modes=["similarity","affine"] if args.mode=="both" else [args.mode]
    if args.command=="phases":
        run_phase_experiment(args.run,args.output_dir,steps=args.steps,phases=args.phases,interval_mm=args.interval_mm,
                             modes=modes,selection=args.selection,plots=not args.no_plots)
        print(f"Phase/interval experiment written to {args.output_dir}")
        return
    roots=getattr(args,"runs",[])
    if args.command=="sweep":
        if len(args.steps)<2 or any(s<=0 or s>180 for s in args.steps):parser.error("At least two angular steps in (0,180] required")
        if any(x in extra for x in ["--angle-step","--output-dir","--no-rim-standardization"]):parser.error("Sweep owns angle-step/output-dir and requires rim standardization")
        if extra[:1]==["--"]:extra=extra[1:]
        roots=[]
        for index,step in enumerate(args.steps):
            root=args.output_dir/f"run_{index:02d}_{step:g}deg"
            subprocess.run([sys.executable,str(Path(__file__).with_name("pottery_radial_sections.py")),str(args.mesh),
                *extra,"--angle-step",str(step),"--output-dir",str(root),"--volume-mode","none","--no-visualization"],check=True)
            roots.append(root)
    if args.command=="compare" and len(roots)<2:parser.error("Compare requires at least two runs")
    for mode in modes:
        for index,root in enumerate(roots):
            validate_run(root,args.output_dir/f"validation_{index:02d}",interval_mm=args.interval_mm,
                mode=mode,selection=args.selection,plots=not args.no_plots)
        if args.command in ("compare","sweep"):
            compare_runs(roots,args.output_dir,interval_mm=args.interval_mm,mode=mode,selection=args.selection)
    print(f"Development validation written to {args.output_dir}")


if __name__=="__main__":main()

#!/usr/bin/env python3
"""Generate an individual whole-vessel representative, separately from rim models."""
import argparse
from pathlib import Path
from morphpot.whole_vessel import VERSION,run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--unit',choices=['auto','mm','cm','m'],default='auto')
    p.add_argument('--angle-step',type=float,default=5)
    p.add_argument('--start-angle',type=float,default=0)
    p.add_argument('--profile-spacing-mm',type=float,default=.5,help='Approximate median-wall arc spacing; not voxel pitch')
    p.add_argument('--axis-xy-mm',type=float,nargs=2,help='Explicit analysis axis X Y in mm; no pose transform')
    p.add_argument('--axis-surface',choices=['inner','outer'],default='inner')
    p.add_argument('--axis-sections',type=int,default=40)
    p.add_argument('--revolution-sections',type=int,default=180)
    p.add_argument('--max-memory-mb',type=float,default=1024)
    p.add_argument('--require-watertight',action='store_true',help='Reject the input mesh if not watertight; default checks each closed plane intersection')
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--version',action='version',version=VERSION)
    a=p.parse_args()
    try:run(a.input,a.output_dir,a.unit,a.angle_step,a.start_angle,a.profile_spacing_mm,a.axis_xy_mm,
        a.axis_surface,a.axis_sections,a.revolution_sections,a.max_memory_mb,not a.no_plots,a.require_watertight)
    except (ValueError,RuntimeError) as e:p.exit(2,f'ERROR: {e}\n')

if __name__=='__main__':main()

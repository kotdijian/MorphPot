#!/usr/bin/env python3
"""Create regional comparisons from completed exploration directories (no mesh rerun).
Regions use each RAW section's neck arc and lowest 10% height; no warped reclassification.
One-way exact profile distances, evenly spaced 1 mm source stations; not mesh C2C.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

METHODS=['global_arc','neck_segmented','neck_segmented_warp']
COLORS=['black','royalblue','crimson']
def write(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def run(roots,output):
 output.mkdir(parents=True,exist_ok=True);regional=[];dimension=[]
 fig,axes=plt.subplots(len(roots),2,figsize=(13,4.5*len(roots)),squeeze=False)
 for k,root in enumerate(roots):
  source=root.name;lm={float(x['angle_deg']):x for x in json.loads((root/'analysis/holdout_landmarks.json').read_text())}
  with np.load(root/'analysis/raw_holdout_profiles.npz',allow_pickle=False) as raw:
   lengths={(float(a),side):np.linalg.norm(np.diff(raw[f'p{idx:03d}_{side}_mm'],axis=0),axis=1).sum() for idx,a in enumerate(raw['angles_deg']) for side in ['outer','inner']}
  for method in METHODS:
   groups={}
   for row in csv.DictReader((root/f'validation/{method}/holdout_station_distances.csv').open()):
    if row['frame']!='original':continue
    a=lm[float(row['angle_deg'])];side=row['surface'];s=float(row['arc_mm']);z=float(row['source_z_mm'])
    # Inner stations run from lip to base: use original length stored in NPZ.
    length=lengths[(float(row['angle_deg']),side)]
    base_s=s if side=='outer' else length-s
    guard=max(a['base_outer_mm'][1],a['base_inner_mm'][1])+.1*(a['lip_mm'][1]-a['base_outer_mm'][1])
    region='base' if z<=guard else ('neck_to_lip' if base_s>=a[side]['arc_mm'] else 'body')
    groups.setdefault(region,[]).append(float(row['distance_mm']))
   for region,v in groups.items():
    v=np.asarray(v);regional.append(dict(sample=source,method=method,region=region,n=len(v),mean_mm=v.mean(),rms_mm=np.sqrt(np.mean(v*v)),p95_mm=np.percentile(v,95),max_mm=v.max()))
  regions=['base','body','neck_to_lip'];x=np.arange(3)
  for j,(method,color) in enumerate(zip(METHODS,COLORS)):
   rows=[r for r in regional if r['sample']==source and r['method']==method];axes[k,0].bar(x+(j-1)*.23,[next(r['rms_mm'] for r in rows if r['region']==q) for q in regions],width=.23,color=color,label=method)
  axes[k,0].set(xticks=x,xticklabels=regions,ylabel='Original-frame holdout RMS [mm]',title=source);axes[k,0].legend(fontsize=8)
  rows=list(csv.DictReader((root/'validation/dimensions.csv').open()));rawrows=[r for r in rows if r['dataset']=='holdout' and r['frame']=='original' and r['method']==METHODS[0]]
  metrics=['lip_diameter_mm','maximum_diameter_mm','maximum_diameter_height_mm','vessel_height_mm']
  for j,(method,color) in enumerate(zip(METHODS,COLORS)):
   model=next(r for r in rows if r['dataset']=='model' and r['method']==method);deltas=[]
   for metric in metrics:
    v=np.array([float(r[metric]) for r in rawrows]);value=float(model[metric]);delta=value-np.median(v);deltas.append(delta)
    dimension.append(dict(sample=source,method=method,metric=metric,model_mm=value,raw_holdout_median_mm=np.median(v),raw_holdout_std_mm=v.std(),raw_holdout_min_mm=v.min(),raw_holdout_max_mm=v.max(),model_minus_raw_median_mm=delta))
   axes[k,1].bar(np.arange(4)+(j-1)*.23,deltas,width=.23,color=color,label=method)
  axes[k,1].set(xticks=np.arange(4),xticklabels=['Lip diameter','Max diameter','Max-D height','Vessel height'],ylabel='Model minus raw holdout median [mm]',title=source);axes[k,1].axhline(0,color='.5',lw=.7);axes[k,1].legend(fontsize=8)
 fig.tight_layout();fig.savefig(output/'regional_and_dimensions.png',dpi=170);plt.close(fig)
 write(output/'regional_distances.csv',regional);write(output/'dimension_comparison.csv',dimension)
 print(json.dumps(regional,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('roots',nargs='+',type=Path);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();run(a.roots,a.output_dir)

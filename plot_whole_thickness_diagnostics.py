#!/usr/bin/env python3
"""Plot model normal-ray diagnostics from exploration outputs; no mesh computation."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
METHODS=['global_arc','neck_segmented','neck_segmented_warp'];COLORS=['black','royalblue','crimson']
def run(root):
 rows=list(csv.DictReader((root/'validation/thickness_stations.csv').open()))
 fig,axes=plt.subplots(1,2,figsize=(13,6))
 for method,color in zip(METHODS,COLORS):
  a=np.load(root/f'analysis/{method}/profile_distribution.npz',allow_pickle=False)
  for side in ['outer','inner']:axes[0].plot(*a[side+'_median_mm'].T,color=color,lw=1)
  r=[q for q in rows if q['dataset']=='model' and q['method']==method];s=np.array([float(q['arc_mm']) for q in r]);v=np.array([float(q['thickness_mm']) if q['thickness_mm'] else np.nan for q in r]);axes[1].plot(s,v,color=color,label=method)
  for q in r:
   if float(q['u'])>.9 and q['inner_r_mm']:
    axes[0].plot([float(q['r_mm']),float(q['inner_r_mm'])],[float(q['z_mm']),float(q['inner_z_mm'])],color=color,alpha=.3,lw=.7)
  lip=a['outer_median_mm'][-1]
 axes[0].set(aspect='equal',xlim=(lip[0]-40,lip[0]+15),ylim=(lip[1]-35,lip[1]+5),xlabel='Radius [mm]',ylabel='Input Z [mm]',title='Rim / outer-normal ray intersections (last 10% arc)')
 axes[1].set(xlabel='Outer arc from axial base [mm]',ylabel='Normal-ray length [mm]',title='Diagnostic only: remote hits may be marked OK');axes[1].legend(fontsize=8)
 fig.tight_layout();fig.savefig(root/'validation/model_thickness_diagnostic.png',dpi=160);plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('roots',nargs='+',type=Path);a=p.parse_args()
 for root in a.roots:run(root)

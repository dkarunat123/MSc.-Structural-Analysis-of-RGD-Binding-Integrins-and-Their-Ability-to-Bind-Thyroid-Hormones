import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
base=Path(__file__).resolve().parents[1] / 'docs'
with (base/'data/species_medians.csv').open(encoding='utf8',newline='') as f: dat=list(csv.DictReader(f))
sp=['H. sapiens','S. purpuratus','N. vectensis']
lig=['Cilengitide','T3','T4','TRIAC','TETRAC','sT3','sT4','rT3','T4AM']
mat=np.array([[float(next(r['median_neg_log10_kd'] for r in dat if r['species']==s and r['ligand']==l)) for l in lig] for s in sp])
co=leaves_list(linkage(pdist(mat.T),method='complete'))
ro=leaves_list(linkage(pdist(mat),method='complete'))
mat=mat[ro,:][:,co];sp=[sp[i] for i in ro];lig=[lig[i] for i in co]
cmap=LinearSegmentedColormap.from_list('thesis',['#313695','#74ADD1','#FFFFBF','#FDAE61','#A50026'])
fig=plt.figure(figsize=(14.5,5.2),facecolor='#ffffff',dpi=165)
SPECIES_COLORS={'H. sapiens':'#396ea5','S. purpuratus':'#a641a6','N. vectensis':'#6f990c'}
ax=fig.add_axes([.105,.215,.79,.56],facecolor='#ffffff')
ax.imshow(mat,cmap=cmap,aspect='auto',vmin=mat.min(),vmax=mat.max(),interpolation='nearest')
ax.set_xticks(range(len(lig)),labels=lig,color='#3e4b61',fontsize=12)
ax.set_yticks(range(len(sp)),labels=[{'H. sapiens':'Homo sapiens','S. purpuratus':'S. purpuratus','N. vectensis':'N. vectensis'}[s] for s in sp],color='#344258',fontsize=12)
ax.tick_params(length=0,pad=13)
for i,lab in enumerate(ax.get_yticklabels()): lab.set_color(SPECIES_COLORS[sp[i]]);lab.set_fontweight('bold')
for i,sp_name in enumerate(sp): ax.add_patch(plt.Rectangle((-0.66,i-0.43),0.06,0.86,facecolor=SPECIES_COLORS[sp_name],edgecolor='none',clip_on=False))
for i in range(3):
 for j in range(9):
  v=mat[i,j];c='#172538' if 5.7<v<7.4 else '#ffffff'
  ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=10,fontweight='bold',color=c)
for s in ax.spines.values():s.set_visible(False)
ax.set_xticks(np.arange(-.5,9,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True)
ax.grid(which='minor',color='#f2f4f9',linewidth=2);ax.tick_params(which='minor',bottom=False,left=False)
fig.text(.105,.94,'PREDICTED BINDING AFFINITY · FIGURE 16',fontsize=11,color='#396ea5',fontweight='bold')
fig.text(.105,.858,'Cross-species median ligand profiles',fontsize=23,color='#26364a',fontweight='bold')
fig.text(.105,.058,'Species label colour = structure identity  ·  heatmap colour = predicted -log10(KD)  ·  higher = stronger predicted binding',fontsize=10.5,color='#526178')
cb_ax=fig.add_axes([.915,.22,.013,.55]);n=Normalize(vmin=mat.min(),vmax=mat.max());cb=plt.colorbar(plt.cm.ScalarMappable(norm=n,cmap=cmap),cax=cb_ax);cb.ax.tick_params(colors='#4f5c70',labelsize=10,length=0);cb.outline.set_visible(False)
fig.savefig(base/'assets/median-affinity-preview.png',facecolor='#ffffff',bbox_inches='tight',pad_inches=.15)
print('saved',base/'assets/median-affinity-preview.png')

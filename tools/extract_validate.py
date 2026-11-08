#!/usr/bin/env python3
"""Build audited static thesis dashboard data from the supplied workbooks.

Used to produce the November 2026 GitHub Pages dataset. This utility requires
artifact_tool, numpy and scipy in the creation environment. See docs/README.md
for the public, R-native regeneration workflow.
"""
import csv, json, math, re
from pathlib import Path
from collections import Counter, defaultdict
from artifact_tool import Blob, SpreadsheetFile
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'Data'
OUT = ROOT / 'docs' / 'data'
OUT.mkdir(parents=True, exist_ok=True)
LIGANDS=['Cilengitide','T3','T4','TRIAC','TETRAC','sT3','sT4','rT3','T4AM']
SPECIES=['H. sapiens','S. purpuratus','N. vectensis']
SHEET_SPECS=[('H. sapiens','H.S_AFReceptors',8),('S. purpuratus','S.P_AFReceptors',24),('N. vectensis','N.V_AFReceptors',16)]
HUMAN={'avb1':('αV','β1'),'avb3':('αV','β3'),'avb5':('αV','β5'),'avb6':('αV','β6'),'avb8':('αV','β8'),'a5b1':('α5','β1'),'a8b1':('α8','β1'),'aiibb3':('αIIb','β3')}
SP_ALPHA={'av3':'αV3','av2':'αV2','av1':'αV1','a9':'α9','a8':'α8','a6':'α6','aP2':'αP2','aP1':'αP1'}
SP_BETA={'b1a1':'β1A1','bg':'βG','bc':'βC'}
NV_ALPHA={'aH6':'α3','aA9':'α9','aA5':'α5','a71':'α8'}
NV_BETA={'bW0':'β2','bU5':'βn2','bD5':'β1','b23':'βb'}

def sheet_matrix(sh,cell_range):
 return sh.get_range(cell_range).values

def id_label(sp, r):
 if sp==SPECIES[0]:
  a,b=HUMAN[r]
 elif sp==SPECIES[1]:
  a,b=r.split('_');a,b=SP_ALPHA[a],SP_BETA[b.lower()]
 else:
  a,b=r.split('_');a,b=NV_ALPHA[a],NV_BETA[b]
 return a+b,a,b

def export(name,fields,rows):
 with (OUT/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

wb=SpreadsheetFile.import_xlsx(Blob.load(str(SOURCE/'BindingScores_StatisticalAnalysis.xlsx')))
sh=wb.worksheets.get_item('Sheet1')
vals=sheet_matrix(sh,'A1:J49')
assert len(vals)==49 and [str(x).strip() for x in vals[0]]==['Receptor']+LIGANDS,(len(vals),vals[0])
rows=vals[1:];receptors=[r[0] for r in rows];assert len(receptors)==len(set(receptors))==48
for r in rows:assert len(r)==10 and all(isinstance(v,(float,int)) and math.isfinite(v) for v in r[1:]),r
species_by_id={rid:sp for sp,n in zip(SPECIES,[8,24,16]) for rid in receptors[sum([8,24,16][:SPECIES.index(sp)]):sum([8,24,16][:SPECIES.index(sp)+1])]}
# Verify original Kd values against transformed score workbook. Source sheets are ligand rows × receptor columns.
raw_book=SpreadsheetFile.import_xlsx(Blob.load(str(SOURCE/'All_Species_BindingScores.xlsx')))
raw={}
for sp,sheet,count in SHEET_SPECS:
 x=sheet_matrix(raw_book.worksheets.get_item(sheet),f'A4:{chr(ord("A")+count)}13')
 assert len(x)==10 and len(x[0])==count+1
 assert sorted(v.strip() for v in [v[0] for v in x[1:]])==sorted(LIGANDS)
 for j,rid in enumerate(x[0][1:],1):
  assert species_by_id.get(rid)==sp,(sp,rid)
  for rr in x[1:]:raw[(rid,rr[0].strip())]=float(rr[j])
assert len(raw)==48*9
max_err=0.0
records=[]
for r in rows:
 rid=r[0];sp=species_by_id[rid];label,a,b=id_label(sp,rid)
 for l,v in zip(LIGANDS,r[1:]):
  kd=raw[(rid,l)]
  diff=abs(-math.log10(kd)-float(v));max_err=max(max_err,diff)
  records.append(dict(species=sp,receptor=rid,receptor_label=label,alpha_subunit=a,beta_subunit=b,ligand=l,neg_log10_kd=float(v),kd_m=kd))
assert max_err<1e-8,max_err
export('receptor_affinities.csv',list(records[0]),records)

corr_input=[]
with (ROOT/'Results/05_Species_BindingAffinity_Analysis_RESULTS/Species_LigandPair_SpearmanCorrelations.csv').open(newline='',encoding='utf-8-sig') as f:
 corr_input=list(csv.DictReader(f))
assert len(corr_input)==3*36,len(corr_input)
corr={}
max_corr_error=0.0
medians=[]
clusters={}
for sp in SPECIES:
 rec=[x for x in receptors if species_by_id[x]==sp]
 sub=np.array([[float(next(rr[1+i] for rr in rows if rr[0]==r)) for i in range(9)] for r in rec],dtype=float)
 for j,l in enumerate(LIGANDS):
  medians.append(dict(species=sp,ligand=l,n_receptors=len(rec),median_neg_log10_kd=float(np.median(sub[:,j]))))
 sp_rows=[r for r in corr_input if r['Species']==sp]
 assert len(sp_rows)==36
 for ent in sp_rows:
  i=LIGANDS.index(ent['Ligand_1']);j=LIGANDS.index(ent['Ligand_2'])
  val=float(spearmanr(sub[:,i],sub[:,j]).statistic)
  err=abs(val-float(ent['Spearman_rho']));max_corr_error=max(max_corr_error,err)
  assert err<1e-10,(sp,ent,err)
  corr[(sp,ent['Ligand_1'],ent['Ligand_2'])]=val
 # Cluster row and ligand profiles with the same method and distance as the R scripts,
 # but SCiPy leaf orientations can differ from the previously saved R pvclust object.
 row_z=linkage(pdist(sub,metric='euclidean'),method='complete')
 ligand_z=linkage(pdist(sub.T,metric='euclidean'),method='complete')
 from scipy.stats import rankdata
 ranks=rankdata(sub,axis=0)
 rho=np.corrcoef(ranks,rowvar=False)
 corr_z=linkage(squareform(1-rho,checks=False),method='complete')
 clusters[sp]={
  'figure14':{'row_ids':rec,'row_order':leaves_list(row_z).tolist(),'row_linkage':row_z.tolist(),'ligand_order':leaves_list(ligand_z).tolist(),'ligand_linkage':ligand_z.tolist()},
  'figure15':{'ligand_order':leaves_list(corr_z).tolist(),'ligand_linkage':corr_z.tolist()},
 }
 # cross-check all 9x9 values against supplied pair table
 for i,a in enumerate(LIGANDS):
  for j,b in enumerate(LIGANDS):
   if i==j: continue
   key=(sp,a,b) if (sp,a,b) in corr else (sp,b,a)
   assert abs(rho[i,j]-corr[key])<1e-10
 export(f'figure15_{sp.replace(". ","_").replace(". ","_").replace(" ","_").replace(".","")}_rho.csv', ['ligand']+LIGANDS, [dict(ligand=a,**dict(zip(LIGANDS,rho[i].tolist()))) for i,a in enumerate(LIGANDS)])
export('species_medians.csv',list(medians[0]),medians)
median_mat=np.array([[next(x['median_neg_log10_kd'] for x in medians if x['species']==sp and x['ligand']==l) for l in LIGANDS] for sp in SPECIES])
from scipy.stats import rankdata
median_rho=np.corrcoef(rankdata(median_mat,axis=1),rowvar=True)
assert abs(median_rho[0,1]-.917)<.001 and abs(median_rho[0,2]-.900)<.001 and abs(median_rho[1,2]-.983)<.001,median_rho
z_cols=linkage(pdist(median_mat.T,metric='euclidean'),method='complete')
z_rows=linkage(pdist(median_mat,metric='euclidean'),method='complete')
clusters['figure16']={'species_order':leaves_list(z_rows).tolist(),'species_linkage':z_rows.tolist(),'ligand_order':leaves_list(z_cols).tolist(),'ligand_linkage':z_cols.tolist()}
# Preserve the published species-pair statistics, for the dashboard and download.
export('species_rank_correlations.csv',['species_a','species_b','spearman_rho'],[dict(species_a=SPECIES[i],species_b=SPECIES[j],spearman_rho=float(median_rho[i,j])) for i,j in [(0,1),(0,2),(1,2)]])
with (OUT/'cluster_orders.json').open('w',encoding='utf8') as f:json.dump(clusters,f,ensure_ascii=False,indent=2)
report={
 'source_files':['All_Species_BindingScores.xlsx','BindingScores_StatisticalAnalysis.xlsx','Species_LigandPair_SpearmanCorrelations.csv'],
 'n_receptors':48,'n_ligands':9,'n_predictions':432,'n_species':3,
 'receptor_counts':dict(zip(SPECIES,[8,24,16])),
 'max_abs_neglog10_kd_difference':max_err,
 'n_spearman_ligand_pairs_verified':108,
 'max_abs_spearman_difference':max_corr_error,
 'species_spearman':[{ 'species_a':SPECIES[i], 'species_b':SPECIES[j], 'rho':float(median_rho[i,j])} for i,j in [(0,1),(0,2),(1,2)]],
 'bootstrap_annotation_status':'not included: original pvclust RDS and edge tables were not uploaded; linkage recreated from the complete numeric data with Euclidean distance and complete linkage. Node AU supports not invented.',
 'clustering_method':'SciPy hierarchical complete linkage; branch order/orientation may differ from the original R pvclust or hclust ordering',
 'calculation_notes':['Species medians are calculated on unrounded already-transformed -log10(KD) values.','No second log transform is used.','Species pair Spearman rho is calculated over nine median ligand profiles.','Raw KD values are retained independently; median(KD) differs from transform of median(-log10 KD).']
}
with (OUT/'validation_report.json').open('w',encoding='utf8') as f:json.dump(report,f,indent=2,ensure_ascii=False)
print(json.dumps(report,indent=2,ensure_ascii=False))

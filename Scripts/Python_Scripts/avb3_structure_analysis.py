#!/usr/bin/env python3
"""Independent audit of the user-provided HADDOCK receptor-ligand PDBs.

Not a docking, scoring, or hydrogen-bond energetics program. Uses only the
explicit atomic coordinates in saved complexes. Hydrogen bonds are
*geometric candidates* (N/O donor with attached H and O acceptor in opposite
partner; D-A <= 3.5 A, H-A <= 2.6 A, D-H-A angle >= 120 deg).

IMPORTANT: The prepared receptor uses *one chain A* for both subunits. In
these specific uploaded files the sequence discontinuity is between
A993 (alpha V) and A994 (beta 3); chain B is ligand UNL. Verify this mapping
for any other dataset before reuse.

Run: python avb3_structure_analysis.py
Requires pandas, numpy, scipy. Outputs three CSVs in the same folder.
"""
import re
from pathlib import Path
from math import acos, degrees
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

OUTDIR=Path(__file__).resolve().parent
AA3=dict(ALA='Ala', ARG='Arg', ASN='Asn', ASP='Asp', CYS='Cys', GLN='Gln', GLU='Glu', GLY='Gly', HIS='His', ILE='Ile', LEU='Leu', LYS='Lys', MET='Met', PHE='Phe', PRO='Pro', SER='Ser', THR='Thr', TRP='Trp', TYR='Tyr', VAL='Val')

def parse_pdb(path):
    atoms=[]; remarks={}
    for l in path.open():
        if l.startswith(('ATOM  ', 'HETATM')):
            try:
                atom=l[12:16].strip(); res=l[17:20].strip(); ch=l[21]; ri=int(l[22:26]);
                element=l[76:78].strip().upper().replace('+2','') or re.sub('[0-9]','',atom)[0]
                xyz=np.array([float(l[i:i+8]) for i in (30,38,46)],float)
            except (IndexError,ValueError):
                continue
            atoms.append({'name':atom,'res':res,'chain':ch,'num':ri,'element':element,'xyz':xyz})
        elif l.startswith('REMARK'):
            if 'buried surface area:' in l:
                remarks['HADDOCK_reported_BSA_A2']=float(l.split(':')[-1])
            if 'energies:' in l:
                vals=[v.strip() for v in l.split(':',1)[1].split(',')]
                if len(vals)>=8:
                    remarks['HADDOCK_vdw_energy']=float(vals[5]); remarks['HADDOCK_electrostatic_energy']=float(vals[6]);remarks['HADDOCK_AIR_energy']=float(vals[7]); remarks['HADDOCK_total_energy_field']=float(vals[0])
            if 'Desolvation energy:' in l: remarks['HADDOCK_desolvation_energy']=float(l.split(':')[-1])
            if 'violations.:' in l: remarks['HADDOCK_violation_count_gt_threshold']=int(l.split(':')[1].split(',')[0].strip())
    return atoms,remarks

def label(a):
    return f"{AA3.get(a['res'],a['res'])} {a['chain']}{a['num']}"

def angle(D,H,A):
    u=D-H;v=A-H
    return degrees(acos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))

def hbonds(partner1,partner2,ligand_name):
    # Explicit donor-H-bonds verified within a residue, not inferred from element proximity alone.
    candidates=[]
    for donors,acceptors,direction in ((partner1,partner2,'protein_to_ligand'),(partner2,partner1,'ligand_to_protein')):
        byres={}
        for a in donors:
            byres.setdefault((a['chain'],a['num'],a['res']),[]).append(a)
        acceptO=[a for a in acceptors if a['element']=='O']
        if not acceptO:continue
        coords=np.array([a['xyz'] for a in acceptO])
        tree=cKDTree(coords)
        for _,a_list in byres.items():
            hs=[a for a in a_list if a['element']=='H']
            for D in a_list:
                if D['element'] not in ('N','O'): continue
                bonded=[h for h in hs if 0.65 <= np.linalg.norm(D['xyz']-h['xyz'])<=1.25]
                if not bonded: continue
                for j in tree.query_ball_point(D['xyz'],3.5):
                    A=acceptO[j];d=np.linalg.norm(D['xyz']-A['xyz'])
                    for H in bonded:
                        ha=np.linalg.norm(H['xyz']-A['xyz'])
                        ang=angle(D['xyz'],H['xyz'],A['xyz'])
                        if ha<=2.6 and ang>=120:
                            p=D if direction=='protein_to_ligand' else A
                            l=A if direction=='protein_to_ligand' else D
                            candidates.append({'Ligand':ligand_name,'Direction':direction,'Receptor_residue':label(p), 'Protein_atom':p['name'], 'Ligand_atom':l['name'], 'Donor_atom':D['name'], 'Hydrogen_atom':H['name'],'Acceptor_atom':A['name'],'Donor_acceptor_A':round(d,3),'Hydrogen_acceptor_A':round(ha,3),'D_H_A_angle_degrees':round(ang,1)})
    # Avoid duplicate donor/acceptor pair when multiple H are geometrically acceptable
    candidates.sort(key=lambda x:(x['Ligand'],x['Receptor_residue'],x['Protein_atom'],x['Ligand_atom'],-x['D_H_A_angle_degrees']))
    return candidates

def analyze(path):
    ligand=re.sub(r'\(\d+\)$','',path.stem.replace('avb3_','',1))
    atoms,remarks=parse_pdb(path)
    protein=[a for a in atoms if a['chain']=='A' and a['res'] not in ('MN2','MN','MG','CA')]
    lig=[a for a in atoms if a['chain']=='B' and a['res']=='UNL']
    metals=[a for a in atoms if a['res'] in ('MN2','MN')]
    pheavy=[a for a in protein if a['element'] not in ('H','D')]
    lheavy=[a for a in lig if a['element'] not in ('H','D')]
    ptree=cKDTree([a['xyz'] for a in pheavy]);pl4={}
    close_heavy=[]
    for l in lheavy:
        for j in ptree.query_ball_point(l['xyz'],4.0):
            p=pheavy[j]
            d=np.linalg.norm(l['xyz']-p['xyz'])
            pl4[(p['num'],p['res'])]=p
            close_heavy.append((d,p,l))
    residues=sorted(pl4)
    alpha=[r for r in residues if r[0]<=993]
    beta=[r for r in residues if r[0]>=994]
    # Very short non-metal protein/ligand heavy-atom distances require manual inspection.
    veryshort=[(d,p,l) for d,p,l in close_heavy if d<2.5]
    hbs=hbonds(protein,lig,ligand)
    metal_dists=sorted([(np.linalg.norm(m['xyz']-l['xyz']),m,l) for m in metals for l in lheavy],key=lambda x:x[0])
    best=metal_dists[0]
    d,m,l=best
    result={'Ligand':ligand, **remarks,'Ligand_heavy_atoms':len(lheavy),
            'Receptor_contact_residues_4A':len(residues),'Alpha_contact_residues_4A':len(alpha), 'Beta_contact_residues_4A':len(beta),
            'Receptor_contact_residue_list':'; '.join(f'{AA3.get(rn,rn)} A{nr}' for nr,rn in residues),
            'Explicit_H_geometric_Hbond_candidates':len(hbs),
            'Distinct_receptor_residues_with_geometric_Hbond_candidates':len(set(h['Receptor_residue'] for h in hbs)),
            'Very_close_nonmetal_heavy_pairs_lt2p5A':len(veryshort),
            'Closest_Mn_residue':f'{m["chain"]}{m["num"]}', 'Closest_ligand_atom_to_Mn':l['name'], 'Closest_ligand_Mn_A':round(d,3),
            'Closest_protein_ligand_heavy_atom_A':round(min(d for d,_,_ in close_heavy),3)}
    short=[{'Ligand':ligand, 'Protein_residue':label(p),'Protein_atom':p['name'],'Ligand_atom':l['name'],'Distance_A':round(d,3)} for d,p,l in veryshort]
    return result,hbs,short

def main():
    rows=[];hb=[];cl=[]
    for f in sorted(OUTDIR.glob('avb3_*.pdb')):
        result,c,s=analyze(f)
        rows.append(result);hb+=c;cl+=s
    df=pd.DataFrame(rows);df.sort_values('Ligand').to_csv(OUTDIR/'avb3_Section3p10_structural_audit.csv',index=False)
    pd.DataFrame(hb).to_csv(OUTDIR/'avb3_Section3p10_Hbond_geometry_candidates.csv',index=False)
    pd.DataFrame(cl).to_csv(OUTDIR/'avb3_Section3p10_very_close_pairs.csv',index=False)
    print(df[['Ligand','HADDOCK_reported_BSA_A2','Receptor_contact_residues_4A','Alpha_contact_residues_4A','Beta_contact_residues_4A','Explicit_H_geometric_Hbond_candidates','Very_close_nonmetal_heavy_pairs_lt2p5A','Closest_ligand_Mn_A']].to_string(index=False))
    print('\nHBOND geometry candidates:',len(hb))
    print(pd.DataFrame(hb).to_string(index=False))
    print('\nClose pairs (<2.5 A):',len(cl))
    print(pd.DataFrame(cl).to_string(index=False))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Compare original SDFs with final HADDOCK PDBs.

Requirements: conda install -c conda-forge rdkit networkx numpy

Outputs: Ligand_Validation_RESULTS/Ligand_Validation_Results.csv and
         Ligand_Validation_RESULTS/Ligand_Validation_Stereocenters.csv

Limitations: PDB lacks reliable bond orders and partial charges; this script
compares heavy-atom *connectivity* (ignores bond orders), and compares signed
3D tetrahedral volumes at *mapped* stereocentres against reference SDF.
It does not establish protonation, ligand force-field charges, or whether
HADDOCK interpreted its ligand topology identically.
"""
from collections import Counter
from pathlib import Path
import csv
import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit import rdBase

# ---------------------------------------------------------------------------
# 1. FOLDER LOCATIONS
# ---------------------------------------------------------------------------
# Keep this script directly inside FinalThesisChecks/Ligands/.
# Its two input folders should sit beside it. Paths are based on the script's
# location, so you do not have to set a working directory manually.
ROOT = Path(__file__).resolve().parent
SDF_DIR = ROOT / 'Original_PubChemLigands'
PDB_DIR = ROOT / 'Converted_DockingLigands'
OUTPUT_DIR = ROOT / 'Ligand_Validation_RESULTS'
OUTPUT_DIR.mkdir(exist_ok=True)
# ---------------------------------------------------------------------------
# 2. MATCH EACH REFERENCE SDF TO THE ACTUAL PDB SUBMITTED TO HADDOCK
# ---------------------------------------------------------------------------
# Format: (short ligand name, original PubChem SDF, final docking PDB).
# Change a filename here only if you have renamed the corresponding file.
# The Dityrosine SDF is intentionally excluded: it is not one of the nine
# ligands in your thesis docking dataset.
PAIRS = [
    ('T3','T3_PubChem.sdf','ProtH_T3_PubChem_New.pdb'),
    ('T4','T4_PubChem.sdf','ProtH_T4_PubChem.pdb'),
    ('rT3','rT3_PubChem.sdf','ProtH_rT3_PubChem_New.pdb'),
    ('sT3','sT3_PubChem.sdf','ProtH_sT3_PubChem.pdb'),
    ('sT4','sT4_PubChem.sdf','ProtH_sT4_PubChem.pdb'),
    ('TRIAC','TRIAC_PubChem.sdf','ProtH_TRIAC_PubChem.pdb'),
    ('TETRAC','TETRAC_PubChem.sdf','ProtH_TETRAC_PubChem.pdb'),
    ('T4AM','T4AM_PubChem.sdf','ProtH_T4AM_PubChem.pdb'),
    ('Cilengitide','Cilengitide_PubChem.sdf','Lig_ProtH_Cilengitide_PubChem.pdb'),
]

# ---------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ---------------------------------------------------------------------------
def heavy_graph(mol):
    """Represent the molecule as a graph of non-hydrogen atoms and bonds.

    Node labels contain chemical element symbols. Bond orders are deliberately
    ignored, because a PDB file may not encode them faithfully. For PDB files,
    RDKit may infer bonds based on atomic coordinates (proximityBonding=True).
    """
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            graph.add_node(atom.GetIdx(), element=atom.GetSymbol())
    for bond in mol.GetBonds():
        u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if u in graph and v in graph:
            graph.add_edge(u, v)
    return graph

def signed_volume(mol, center, neighbors):
    """Calculate a signed 3D volume around a putative stereocentre.

    The *same ordered three mapped neighbours* must be used for both files.
    The volume sign indicates handedness relative to that neighbour ordering;
    it is NOT itself an R/S (CIP) assignment.
    """
    conf = mol.GetConformer()
    origin = np.array(conf.GetAtomPosition(center))
    vectors = [np.array(conf.GetAtomPosition(i)) - origin for i in neighbors]
    return float(np.linalg.det(np.stack(vectors)))

def count_elements(graph):
    """Count elements among heavy atoms (C, N, O, I, S, etc.)."""
    return dict(sorted(Counter(data['element'] for _,data in graph.nodes(data=True)).items()))

# ---------------------------------------------------------------------------
# 4. RUN THE VALIDATION
# ---------------------------------------------------------------------------
def main():
    """Compare each ligand and write machine-readable results plus a summary."""
    rows, centers = [], []
    print('RDKit version:',rdBase.rdkitVersion)
    # Examine the original PubChem structure and the submitted HADDOCK PDB
    # for each of the nine ligands.
    for name, sdf_filename, pdb_filename in PAIRS:
        sdf_path, pdb_path = SDF_DIR / sdf_filename, PDB_DIR / pdb_filename
        if not (sdf_path.is_file() and pdb_path.is_file()):
            rows.append(dict(ligand=name, sdf=sdf_filename, pdb=pdb_filename, status='MISSING_FILE', note=f'File missing: {sdf_path if not sdf_path.exists() else pdb_path}'))
            continue
        # Read the reference with full SDF bond-order and stereochemistry data.
        # Read the docking PDB without sanitization: PDB bond orders are often
        # incomplete/ambiguous, and RDKit may infer connectivity from geometry.
        refs = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        ref = refs[0] if len(refs) else None
        pdb = Chem.MolFromPDBFile(str(pdb_path), removeHs=False,
                                  sanitize=False, proximityBonding=True)
        if ref is None or pdb is None:
            rows.append(dict(ligand=name,status='PARSE_FAILURE',note='RDKit could not read one input'))
            continue
        # First compare elemental makeup and the heavy-atom bond graph.
        # Graph isomorphism allows atom numbers and names to be different.
        gr, gp = heavy_graph(ref),heavy_graph(pdb)
        same_elements = count_elements(gr) == count_elements(gp)
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            gr,gp,node_match=lambda x,y:x['element']==y['element'])
        topology_matches = matcher.is_isomorphic()
        row = dict(ligand=name,sdf=sdf_filename,pdb=pdb_filename,
                   reference_heavy_atoms=gr.number_of_nodes(),pdb_heavy_atoms=gp.number_of_nodes(),
                   reference_heavy_bonds=gr.number_of_edges(),pdb_heavy_bonds=gp.number_of_edges(),
                   elemental_composition_matches=same_elements,
                   heavy_atom_connectivity_matches=topology_matches,
                   stereocenters_checked=0,stereocenters_matching=0,
                   status='',note='')
        # If connectivity does not match, stereochemistry cannot be compared
        # reliably using the mapping. Report a failure rather than guessing.
        if not topology_matches:
            row['status']='FAIL_CONNECTIVITY'
            row['note']='PDB-inferred heavy-atom graph differs from SDF; stereo comparison cannot be trusted'
        else:
            # Match atom indices between the two files and identify assigned
            # stereocentres in the reference SDF (CIP labels: R or S).
            # NOTE: graph automorphisms can give multiple valid mappings.
            # This script uses the first mapping; ambiguous symmetric cases
            # may require additional manual review.
            mapping = matcher.mapping # reference heavy atom index -> PDB heavy atom index
            Chem.AssignStereochemistry(ref,cleanIt=True,force=True)
            stereo_atoms = [a for a in ref.GetAtoms() if a.HasProp('_CIPCode')]
            row['stereocenters_checked'] = len(stereo_atoms)
            problems=[]
            for atom in stereo_atoms:
                index = atom.GetIdx()
                neighbors = sorted(n.GetIdx() for n in atom.GetNeighbors()
                                   if n.GetAtomicNum() > 1)
                if len(neighbors) != 3:
                    problems.append(f'Cannot verify center {index} with {len(neighbors)} heavy neighbors')
                    continue
                # Three heavy neighbours + one implicit/explicit H permit a
                # handedness comparison at the mapped centre. A near-zero
                # volume is treated as indeterminate, not a stereochemical pass.
                v_ref = signed_volume(ref,index,neighbors)
                v_pdb = signed_volume(pdb,mapping[index],[mapping[x] for x in neighbors])
                same = abs(v_ref)>1e-3 and abs(v_pdb)>1e-3 and v_ref*v_pdb>0
                row['stereocenters_matching'] += int(same)
                centers.append(dict(ligand=name,reference_atom_index=index,
                                    reference_CIP=atom.GetProp('_CIPCode'),
                                    sdf_signed_volume=round(v_ref,5),
                                    pdb_signed_volume=round(v_pdb,5),
                                    configuration_preserved=same))
            if problems:
                row['status']='INCOMPLETE_STEREO';row['note']='; '.join(problems)
            elif row['stereocenters_matching']==row['stereocenters_checked']:
                row['status']='PASS'
                row['note']=('Reference has no CIP-assigned stereocenters' if not stereo_atoms
                             else '3D handedness matches for all SDF-assigned stereocenters')
            else:
                row['status']='FAIL_STEREO';row['note']='At least one reference stereocenter was inverted'
        rows.append(row)
        print(f'{name:12s} {row["status"]:20s} connectivity={topology_matches}, stereo={row["stereocenters_matching"]}/{row["stereocenters_checked"]}')
    # -----------------------------------------------------------------------
    # 5. WRITE THE TWO CSV OUTPUT FILES
    # -----------------------------------------------------------------------
    # Summary CSV: one row per ligand with status and counts.
    # Stereocentre CSV: one row per checked stereocentre with the original
    # reference CIP assignment and whether its 3D handedness was retained.
    fields=['ligand','sdf','pdb','reference_heavy_atoms','pdb_heavy_atoms',
            'reference_heavy_bonds','pdb_heavy_bonds','elemental_composition_matches',
            'heavy_atom_connectivity_matches','stereocenters_checked',
            'stereocenters_matching','status','note']
    with open(OUTPUT_DIR/'Ligand_Validation_Results.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    with open(OUTPUT_DIR/'Ligand_Validation_Stereocenters.csv','w',newline='') as f:
        fields2=['ligand','reference_atom_index','reference_CIP','sdf_signed_volume',
                 'pdb_signed_volume','configuration_preserved']
        w=csv.DictWriter(f,fieldnames=fields2);w.writeheader();w.writerows(centers)
    # -----------------------------------------------------------------------
    # 6. PRINT A READABLE END-OF-RUN SUMMARY
    # -----------------------------------------------------------------------
    counts = Counter(row.get('status', '') for row in rows)
    total = len(PAIRS)
    print('\n' + '=' * 72)
    print('LIGAND VALIDATION — OUTPUT SUMMARY')
    print('=' * 72)
    print(f'Ligands in comparison list:     {total}')
    print(f'Passed:                         {counts.get("PASS", 0)}')
    print(f'Connectivity mismatches:        {counts.get("FAIL_CONNECTIVITY", 0)}')
    print(f'Stereochemistry mismatches:     {counts.get("FAIL_STEREO", 0)}')
    print(f'Incomplete stereo assessments: {counts.get("INCOMPLETE_STEREO", 0)}')
    print(f'Missing input files:            {counts.get("MISSING_FILE", 0)}')
    print(f'Parsing failures:              {counts.get("PARSE_FAILURE", 0)}')
    print(f'Reference stereocentres checked:{sum(r.get("stereocenters_checked", 0) for r in rows):>5}')
    print(f'Stereocentres matching:         {sum(r.get("stereocenters_matching", 0) for r in rows)}')
    print('\nFILES CREATED')
    print('1. Ligand_Validation_Results.csv')
    print('   - One row per ligand: element composition, heavy-atom connectivity,')
    print('     number of stereocentres checked/matching, PASS/FAIL status, notes.')
    print('2. Ligand_Validation_Stereocenters.csv')
    print('   - One row per examined centre: reference CIP label, 3D signed volumes,')
    print('     and whether the mapped 3D configuration was preserved.')
    print(f'\nOutput directory: {OUTPUT_DIR}')
    flagged = [r for r in rows if r.get('status') != 'PASS']
    if flagged:
        print('\nLIGANDS REQUIRING REVIEW')
        for row in flagged:
            print(f' - {row["ligand"]}: {row["status"]} — {row.get("note", "")}')
    else:
        print('\nAll listed ligands passed the implemented checks.')
    print('\nIMPORTANT LIMITATIONS')
    print('- PDB connectivity can be inferred incorrectly; a mismatch is a flag for review.')
    print('- This does not verify protonation, partial charges, or HADDOCK parameters.')
    print('- PASS means checks passed against the provided SDF, not that every')
    print('  stereochemical or chemical property was independently validated.')
    print('=' * 72)

if __name__ == '__main__':
    main()
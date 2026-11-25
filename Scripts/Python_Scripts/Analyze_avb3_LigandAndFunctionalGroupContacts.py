#!/usr/bin/env python3
"""Analyze docked ligands and thyroid-hormone functional groups in human alphaVbeta3.

PROJECT
    Comparative Structural Analysis of RGD-Binding Integrins and Their Ability
    to Bind Thyroid Hormones (MSc Bioinformatics thesis).

PURPOSE
    Consolidate the former Analyze_avb3_LigandContacts.py and
    Analyze_avb3_StructuralContacts.py without changing their scientific
    distance-based calculations. One pass through each complex produces:

    A. Whole-ligand contact results for all nine ligands.
    B. Functional-group-specific results for eight thyroid-hormone compounds.
       Cilengitide is analyzed in part A only.

INPUT
    All avb3_*.pdb files in the SAME directory as this script.
    Expected HADDOCK PDB convention:
      - Integrin receptor: chain A
      - Manganese ions: chain A; residue MN2 or element Mn
      - Ligand: chain B; residue UNL
    Input filenames identify the ligands, e.g. avb3_T4.pdb.

OUTPUT
    Analyze_avb3_LigandAndFunctionalGroupContacts_RESULTS/ beside this script.
    Eight CSV files, with the same filenames and column definitions used by
    the former StructuralContacts script. Existing input PDBs are never edited.

METHOD
    Euclidean heavy-atom distances in Angstroms. Contact cutoffs are 4.0 A
    for receptor residues and ligand/Mn proximities and 3.5 A for potential
    polar atom pairs (N, O, S). Ligand functional groups are assigned using
    geometric bond inference from the PDB coordinates, not a chemical
    bond-order model. All computations use only Python's standard library.

INTERPRETATION
    A polar atom pair within 3.5 A is NOT a confirmed hydrogen bond.
    A short ligand/Mn distance is NOT proof of valid metal coordination.
    Functional-group detection depends on the deposited coordinates and the
    presence of explicit hydrogens (especially for phenolic hydroxyl groups).
"""

from pathlib import Path
import csv
import math
import re


# =============================================================================
# 0. USER-CONFIGURABLE INPUTS, CUT-OFFS, AND OUTPUT FOLDER
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILES = sorted(SCRIPT_DIR.glob("avb3_*.pdb"))
OUTPUT_DIR = SCRIPT_DIR / "Analyze_avb3_LigandAndFunctionalGroupContacts_RESULTS"

RESIDUE_CONTACT_CUTOFF = 4.0  # Receptor residue if any heavy atom is within 4 A
POLAR_CONTACT_CUTOFF = 3.5    # O/N/S atom-pair proximity (NOT a confirmed H bond)
MN_CONTACT_CUTOFF = 4.0       # Ligand/Mn proximity threshold
BOND_TOLERANCE = 0.45         # Added to the sum of covalent radii, in A

EXCLUDED_RECEPTOR_RESNAMES = {
    "UNL", "MN2", "HOH", "WAT", "H2O", "TIP", "TIP3"
}
POLAR_ELEMENTS = {"N", "O", "S"}

# Approximate covalent radii (Angstrom). Used ONLY to infer which ligand atoms
# may be connected for functional-group identification.
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
    "S": 1.05, "I": 1.39, "P": 1.07,
}

# Expected groups in thyroid-hormone derivatives. These are used for warnings,
# not to force an assignment if a functional group cannot be inferred.
EXPECTED_GROUPS = {
    "T3": {"carboxyl_group", "amine_group", "phenolic_hydroxyl",
           "ether_oxygen", "iodine_substituents"},
    "T4": {"carboxyl_group", "amine_group", "phenolic_hydroxyl",
           "ether_oxygen", "iodine_substituents"},
    "rT3": {"carboxyl_group", "amine_group", "phenolic_hydroxyl",
            "ether_oxygen", "iodine_substituents"},
    "sT3": {"sulfate_group", "carboxyl_group", "amine_group",
            "ether_oxygen", "iodine_substituents"},
    "sT4": {"sulfate_group", "carboxyl_group", "amine_group",
            "ether_oxygen", "iodine_substituents"},
    "TRIAC": {"carboxyl_group", "phenolic_hydroxyl", "ether_oxygen",
              "iodine_substituents"},
    "TETRAC": {"carboxyl_group", "phenolic_hydroxyl", "ether_oxygen",
               "iodine_substituents"},
    "T4AM": {"amine_group", "phenolic_hydroxyl", "ether_oxygen",
              "iodine_substituents"},
}


# =============================================================================
# 1. HELPERS: PDB PARSING, DISTANCES, AND CSV EXPORT
# =============================================================================

def infer_element(atom_name):
    """Infer an element when the PDB element columns (77-78) are empty."""
    letters = re.sub(r"[^A-Za-z]", "", atom_name).upper()
    for symbol in ("MN", "CL", "BR"):
        if letters.startswith(symbol):
            return symbol
    if letters.startswith("I"):
        return "I"
    return letters[:1]


def parse_pdb(path):
    """Read PDB ATOM/HETATM records; retain coordinates and atom serials."""
    atoms = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26].strip()) if line[22:26].strip() else 0
                x, y, z = (float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
                element = line[76:78].strip().upper()
            except (ValueError, IndexError):
                continue
            if not element:
                element = infer_element(name)
            elif element.startswith("MN"):
                element = "MN"
            atoms.append({
                "serial": serial, "atom": name, "resname": resname,
                "chain": chain, "resseq": resseq,
                "x": x, "y": y, "z": z, "element": element,
            })
    return atoms


def distance(a, b):
    """Return geometric distance between two atoms in Angstroms."""
    return math.sqrt(sum((a[axis] - b[axis]) ** 2 for axis in ("x", "y", "z")))


def is_hydrogen(atom):
    """Apply the hydrogen filter used in the two original scripts."""
    return atom["element"] == "H" or atom["atom"].upper().startswith("H")


def residue_label(atom):
    """Format residues as e.g. 'SER A1140', preserving docking numbering."""
    return f'{atom["resname"]} {atom["chain"]}{atom["resseq"]}'


def write_csv(filename, rows, columns):
    """Export reproducible CSV tables with explicit column order and no index."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_DIR / filename).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def readable_distance(record, index=0):
    """Round a closest-contact distance or return 'NA' if no atoms exist."""
    return round(record[index], 3) if record is not None else "NA"


# =============================================================================
# 2. IDENTIFY THYROID-HORMONE FUNCTIONAL GROUPS BY GEOMETRIC CONNECTIVITY
# =============================================================================

def infer_ligand_bonds(ligand_atoms):
    """Infer adjacency from distances and approximate radii, as in old script."""
    adjacent = {atom["serial"]: [] for atom in ligand_atoms}
    for i, a in enumerate(ligand_atoms):
        for b in ligand_atoms[i + 1:]:
            e1, e2 = a["element"], b["element"]
            if e1 == e2 == "H":
                continue
            cutoff = (COVALENT_RADII.get(e1, 0.77)
                      + COVALENT_RADII.get(e2, 0.77) + BOND_TOLERANCE)
            if distance(a, b) <= cutoff:
                adjacent[a["serial"]].append(b)
                adjacent[b["serial"]].append(a)
    return adjacent


def detect_functional_groups(ligand_atoms):
    """Return detected functional groups in the original script's group order."""
    adjacent = infer_ligand_bonds(ligand_atoms)
    by_serial = {atom["serial"]: atom for atom in ligand_atoms}
    groups = []
    assigned = set()
    carboxyl_carbons = set()

    def register(name, serials):
        groups.append({
            "name": name,
            "atoms": [by_serial[serial] for serial in sorted(serials)],
        })

    # Sulfate: sulfur connected to >=3 oxygen atoms.
    for atom in ligand_atoms:
        if atom["element"] != "S":
            continue
        oxygens = [a for a in adjacent[atom["serial"]] if a["element"] == "O"]
        if len(oxygens) >= 3:
            serials = {atom["serial"]} | {a["serial"] for a in oxygens}
            assigned.update(serials)
            register("sulfate_group", serials)

    # Carboxyl: carbon connected to >=2 oxygen atoms.
    for atom in ligand_atoms:
        if atom["element"] != "C":
            continue
        oxygens = [a for a in adjacent[atom["serial"]] if a["element"] == "O"]
        if len(oxygens) >= 2:
            serials = {atom["serial"]} | {a["serial"] for a in oxygens}
            assigned.update(serials)
            carboxyl_carbons.add(atom["serial"])
            register("carboxyl_group", serials)

    # Amine: nitrogen not directly bonded to a carboxyl carbon.
    for atom in ligand_atoms:
        if atom["element"] != "N":
            continue
        neighbors = adjacent[atom["serial"]]
        if any(a["element"] == "C" and a["serial"] in carboxyl_carbons
               for a in neighbors if a["element"] != "H"):
            continue
        serials = {atom["serial"]} | {
            a["serial"] for a in neighbors if a["element"] == "H"
        }
        assigned.update(serials)
        register("amine_group", serials)

    # Phenolic hydroxyl: unassigned oxygen attached to carbon and hydrogen.
    # Without explicit O-H coordinates this group may not be detected.
    for atom in ligand_atoms:
        if atom["element"] != "O" or atom["serial"] in assigned:
            continue
        neighbors = adjacent[atom["serial"]]
        if (any(a["element"] == "H" for a in neighbors)
                and any(a["element"] == "C" for a in neighbors)):
            serials = {atom["serial"]} | {
                a["serial"] for a in neighbors if a["element"] == "H"
            }
            assigned.update(serials)
            register("phenolic_hydroxyl", serials)

    # Ether oxygen: unassigned oxygen connected to exactly two carbons.
    for atom in ligand_atoms:
        if atom["element"] != "O" or atom["serial"] in assigned:
            continue
        carbons = [a for a in adjacent[atom["serial"]] if a["element"] == "C"]
        if len(carbons) == 2:
            assigned.add(atom["serial"])
            register("ether_oxygen", {atom["serial"]})

    # All iodine atoms are grouped into a single iodine category.
    iodine_atoms = [a for a in ligand_atoms if a["element"] == "I"]
    if iodine_atoms:
        groups.append({"name": "iodine_substituents", "atoms": iodine_atoms})

    return groups


# =============================================================================
# 3. WHOLE-LIGAND CONTACT ANALYSIS (ALL NINE LIGANDS)
# =============================================================================

def analyze_whole_ligand(ligand, lig_atoms, receptor_atoms, mn_atoms):
    """Calculate original whole-ligand summary and four contact tables."""
    mn_rows, residue_rows, polar_rows = [], [], []
    closest_mn = None

    # One nearest ligand atom is exported per Mn site within 4 A.
    # Regardless of the cutoff, the closest overall Mn pair is saved in the
    # summary (including when the nearest Mn is farther than 4 A).
    for mn in mn_atoms:
        closest_at_site = None
        for lig in lig_atoms:
            d = distance(lig, mn)
            if closest_at_site is None or d < closest_at_site[0]:
                closest_at_site = (d, lig)
            if closest_mn is None or d < closest_mn[0]:
                closest_mn = (d, mn, lig)
        if closest_at_site is not None and closest_at_site[0] <= MN_CONTACT_CUTOFF:
            d, lig = closest_at_site
            mn_rows.append({
                "Ligand": ligand, "Mn_residue": f'MN2 A{mn["resseq"]}',
                "Ligand_atom": lig["atom"], "Ligand_element": lig["element"],
                "Distance_A": round(d, 3),
            })

    # Keep the shortest heavy-atom distance for each receptor residue.
    nearest_per_residue = {}
    closest_protein = None
    for rec in receptor_atoms:
        key = (rec["chain"], rec["resseq"], rec["resname"])
        for lig in lig_atoms:
            d = distance(rec, lig)
            if closest_protein is None or d < closest_protein[0]:
                closest_protein = (d, rec, lig)
            if d <= RESIDUE_CONTACT_CUTOFF:
                old = nearest_per_residue.get(key)
                if old is None or d < old[0]:
                    nearest_per_residue[key] = (d, rec, lig)

    for d, rec, lig in sorted(nearest_per_residue.values(), key=lambda x: x[0]):
        residue_rows.append({
            "Ligand": ligand, "Receptor_residue": residue_label(rec),
            "Receptor_atom": rec["atom"], "Ligand_atom": lig["atom"],
            "Minimum_distance_A": round(d, 3),
        })

    # Record EVERY heteroatom pair <=3.5 A; these are only potential polar
    # contacts and must not automatically be interpreted as hydrogen bonds.
    ligand_polar = [a for a in lig_atoms if a["element"] in POLAR_ELEMENTS]
    receptor_polar = [a for a in receptor_atoms if a["element"] in POLAR_ELEMENTS]
    for rec in receptor_polar:
        for lig in ligand_polar:
            d = distance(rec, lig)
            if d <= POLAR_CONTACT_CUTOFF:
                polar_rows.append({
                    "Ligand": ligand, "Receptor_residue": residue_label(rec),
                    "Receptor_atom": rec["atom"], "Receptor_element": rec["element"],
                    "Ligand_atom": lig["atom"], "Ligand_element": lig["element"],
                    "Distance_A": round(d, 3),
                })

    summary = {
        "Ligand": ligand,
        "Ligand_heavy_atoms": len(lig_atoms),
        "Contact_residues_within_4A": len(nearest_per_residue),
        "Potential_polar_atom_pairs_within_3.5A": len(polar_rows),
        "Mn_sites_within_4A": len(mn_rows),
        "Closest_Mn": f'MN2 A{closest_mn[1]["resseq"]}' if closest_mn else "NA",
        "Closest_Mn_ligand_atom": closest_mn[2]["atom"] if closest_mn else "NA",
        "Closest_Mn_distance_A": readable_distance(closest_mn),
        "Closest_receptor_residue": residue_label(closest_protein[1]) if closest_protein else "NA",
        "Closest_receptor_atom": closest_protein[1]["atom"] if closest_protein else "NA",
        "Closest_receptor_ligand_atom": closest_protein[2]["atom"] if closest_protein else "NA",
        "Closest_receptor_distance_A": readable_distance(closest_protein),
    }
    return summary, mn_rows, residue_rows, polar_rows


# =============================================================================
# 4. THYROID-HORMONE FUNCTIONAL-GROUP CONTACT ANALYSIS (EIGHT LIGANDS)
# =============================================================================

def analyze_functional_groups(ligand, ligand_all, receptor_atoms, mn_atoms):
    """Calculate the four original functional-group CSV datasets."""
    summary_rows, mn_rows, residue_rows, polar_rows = [], [], [], []
    groups = detect_functional_groups(ligand_all)
    detected = {group["name"] for group in groups}
    expected = EXPECTED_GROUPS[ligand]
    if expected - detected:
        print("  WARNING: expected groups not detected: "
              + ", ".join(sorted(expected - detected)))
    if detected - expected:
        print("  NOTE: additional groups detected: "
              + ", ".join(sorted(detected - expected)))

    receptor_polar = [a for a in receptor_atoms if a["element"] in POLAR_ELEMENTS]

    for group in groups:
        name, all_atoms = group["name"], group["atoms"]
        heavy = [a for a in all_atoms if not is_hydrogen(a)]
        atom_names = ";".join(a["atom"] for a in all_atoms)

        # Every functional-group/Mn atom pair within 4 A, plus the nearest
        # group/Mn pair overall even if it is outside 4 A.
        closest_mn = None
        for ga in heavy:
            for mn in mn_atoms:
                d = distance(ga, mn)
                if closest_mn is None or d < closest_mn[0]:
                    closest_mn = (d, ga, mn)
                if d <= MN_CONTACT_CUTOFF:
                    mn_rows.append({
                        "Ligand": ligand, "Functional_group": name,
                        "Group_atom": ga["atom"],
                        "Mn_residue": f'MN2 A{mn["resseq"]}',
                        "Distance_A": round(d, 3),
                    })

        # Nearest group atom to each receptor residue within 4 A. Retain
        # group-to-receptor closest distance even if beyond that cutoff.
        per_residue = {}
        closest_receptor = None
        for ga in heavy:
            for rec in receptor_atoms:
                d = distance(ga, rec)
                if closest_receptor is None or d < closest_receptor[0]:
                    closest_receptor = (d, ga, rec)
                if d <= RESIDUE_CONTACT_CUTOFF:
                    key = (rec["chain"], rec["resseq"], rec["resname"])
                    old = per_residue.get(key)
                    if old is None or d < old[0]:
                        per_residue[key] = (d, ga, rec)

        contacts = sorted(per_residue.values(), key=lambda x: x[0])
        for d, ga, rec in contacts:
            residue_rows.append({
                "Ligand": ligand, "Functional_group": name,
                "Group_atom": ga["atom"],
                "Receptor_residue": residue_label(rec),
                "Receptor_atom": rec["atom"],
                "Minimum_distance_A": round(d, 3),
            })
        residue_list = "; ".join(
            f"{residue_label(rec)} ({d:.2f} A)" for d, ga, rec in contacts
        )

        # Potential polar pairs within 3.5 A (distance-only screen).
        group_polar = [a for a in heavy if a["element"] in POLAR_ELEMENTS]
        polar_count = 0
        for ga in group_polar:
            for rec in receptor_polar:
                d = distance(ga, rec)
                if d <= POLAR_CONTACT_CUTOFF:
                    polar_count += 1
                    polar_rows.append({
                        "Ligand": ligand, "Functional_group": name,
                        "Group_atom": ga["atom"], "Group_element": ga["element"],
                        "Receptor_residue": residue_label(rec),
                        "Receptor_atom": rec["atom"],
                        "Receptor_element": rec["element"],
                        "Distance_A": round(d, 3),
                    })

        summary_rows.append({
            "Ligand": ligand, "Functional_group": name,
            "Group_atoms": atom_names,
            "Closest_Mn": f'MN2 A{closest_mn[2]["resseq"]}' if closest_mn else "NA",
            "Closest_group_atom_to_Mn": closest_mn[1]["atom"] if closest_mn else "NA",
            "Closest_Mn_distance_A": readable_distance(closest_mn),
            "Mn_within_4A": ("Yes" if closest_mn and closest_mn[0] <= MN_CONTACT_CUTOFF
                             else "No"),
            "Closest_receptor_residue": (residue_label(closest_receptor[2])
                                         if closest_receptor else "NA"),
            "Closest_receptor_atom": (closest_receptor[2]["atom"]
                                      if closest_receptor else "NA"),
            "Closest_group_atom_to_receptor": (closest_receptor[1]["atom"]
                                               if closest_receptor else "NA"),
            "Closest_receptor_distance_A": readable_distance(closest_receptor),
            "Receptor_residues_within_4A": residue_list,
            "Potential_polar_atom_pairs_within_3.5A": polar_count,
        })

    return summary_rows, mn_rows, residue_rows, polar_rows


# =============================================================================
# 5. ANALYZE EACH HADDOCK COMPLEX
# =============================================================================

def main():
    if not INPUT_FILES:
        raise SystemExit(
            "No avb3_*.pdb files found beside this script. "
            "Place the nine prepared docking complexes in the same folder."
        )

    whole_summary, whole_mn, whole_residues, whole_polar = [], [], [], []
    group_summary, group_mn, group_residues, group_polar = [], [], [], []

    for pdb_path in INPUT_FILES:
        ligand = pdb_path.stem.replace("avb3_", "", 1)
        print(f"Analyzing {pdb_path.name} ...")
        atoms = parse_pdb(pdb_path)
        ligand_all = [a for a in atoms if a["chain"] == "B" and a["resname"] == "UNL"]
        ligand_heavy = [a for a in ligand_all if not is_hydrogen(a)]
        mn_atoms = [a for a in atoms if a["chain"] == "A" and
                    (a["resname"] == "MN2" or a["element"] == "MN")]
        receptor = [a for a in atoms if a["chain"] == "A" and
                    a["resname"] not in EXCLUDED_RECEPTOR_RESNAMES and
                    a["element"] != "MN" and not is_hydrogen(a)]

        if not ligand_heavy:
            print(f"  WARNING: no chain B/UNL ligand heavy atoms; skipping {pdb_path.name}.")
            continue
        if not receptor:
            print(f"  WARNING: no chain A receptor heavy atoms; skipping {pdb_path.name}.")
            continue
        if not mn_atoms:
            print("  WARNING: no Mn ions found; nearest-metal fields will be NA.")

        summary, mn, residues, polar = analyze_whole_ligand(
            ligand, ligand_heavy, receptor, mn_atoms
        )
        whole_summary.append(summary)
        whole_mn.extend(mn)
        whole_residues.extend(residues)
        whole_polar.extend(polar)

        # Cilengitide and other non-TH names receive whole-ligand analysis
        # only. No thyroid-hormone functional groups are assigned to them.
        if ligand not in EXPECTED_GROUPS:
            print("  Whole-ligand analysis complete; TH functional groups skipped.")
            continue

        summaries, mn, residues, polar = analyze_functional_groups(
            ligand, ligand_all, receptor, mn_atoms
        )
        group_summary.extend(summaries)
        group_mn.extend(mn)
        group_residues.extend(residues)
        group_polar.extend(polar)

    if not whole_summary:
        raise SystemExit("No valid docking complexes were analyzed; no CSV files written.")

    # ========================================================================
    # 6. EXPORT ORIGINAL EIGHT OUTPUT TABLES (SAME SCHEMA AS OLD SCRIPT)
    # ========================================================================
    write_csv("avb3_contact_summary.csv", whole_summary, [
        "Ligand", "Ligand_heavy_atoms", "Contact_residues_within_4A",
        "Potential_polar_atom_pairs_within_3.5A", "Mn_sites_within_4A",
        "Closest_Mn", "Closest_Mn_ligand_atom", "Closest_Mn_distance_A",
        "Closest_receptor_residue", "Closest_receptor_atom",
        "Closest_receptor_ligand_atom", "Closest_receptor_distance_A",
    ])
    write_csv("avb3_mn_contacts.csv", sorted(
        whole_mn, key=lambda r: (r["Ligand"], r["Distance_A"])), [
        "Ligand", "Mn_residue", "Ligand_atom", "Ligand_element", "Distance_A",
    ])
    write_csv("avb3_residue_contacts.csv", sorted(
        whole_residues, key=lambda r: (r["Ligand"], r["Minimum_distance_A"])), [
        "Ligand", "Receptor_residue", "Receptor_atom", "Ligand_atom",
        "Minimum_distance_A",
    ])
    write_csv("avb3_potential_polar_contacts.csv", sorted(
        whole_polar, key=lambda r: (r["Ligand"], r["Distance_A"])), [
        "Ligand", "Receptor_residue", "Receptor_atom", "Receptor_element",
        "Ligand_atom", "Ligand_element", "Distance_A",
    ])
    write_csv("avb3_functional_group_summary.csv", group_summary, [
        "Ligand", "Functional_group", "Group_atoms", "Closest_Mn",
        "Closest_group_atom_to_Mn", "Closest_Mn_distance_A", "Mn_within_4A",
        "Closest_receptor_residue", "Closest_receptor_atom",
        "Closest_group_atom_to_receptor", "Closest_receptor_distance_A",
        "Receptor_residues_within_4A", "Potential_polar_atom_pairs_within_3.5A",
    ])
    write_csv("avb3_functional_group_mn_contacts.csv", sorted(
        group_mn, key=lambda r: (
            r["Ligand"], r["Functional_group"], r["Distance_A"])), [
        "Ligand", "Functional_group", "Group_atom", "Mn_residue", "Distance_A",
    ])
    write_csv("avb3_functional_group_residue_contacts.csv", sorted(
        group_residues, key=lambda r: (
            r["Ligand"], r["Functional_group"], r["Minimum_distance_A"])), [
        "Ligand", "Functional_group", "Group_atom", "Receptor_residue",
        "Receptor_atom", "Minimum_distance_A",
    ])
    write_csv("avb3_functional_group_polar_contacts.csv", sorted(
        group_polar, key=lambda r: (
            r["Ligand"], r["Functional_group"], r["Distance_A"])), [
        "Ligand", "Functional_group", "Group_atom", "Group_element",
        "Receptor_residue", "Receptor_atom", "Receptor_element", "Distance_A",
    ])

    # ========================================================================
    # 7. COMPLETION REPORT AND OUTPUT SUMMARY
    # ========================================================================
    print(f"\nAnalyzed {len(whole_summary)} valid avb3_*.pdb complexes.")
    print(f"Results directory: {OUTPUT_DIR}")
    print("Created eight CSV files:")
    for filename in (
        "avb3_contact_summary.csv",
        "avb3_mn_contacts.csv",
        "avb3_residue_contacts.csv",
        "avb3_potential_polar_contacts.csv",
        "avb3_functional_group_summary.csv",
        "avb3_functional_group_mn_contacts.csv",
        "avb3_functional_group_residue_contacts.csv",
        "avb3_functional_group_polar_contacts.csv",
    ):
        print("  " + filename)
    print("\nCaution: polar contacts are distance-based, not confirmed H bonds;")
    print("         ligand-Mn proximities do not confirm coordination geometry.")


if __name__ == "__main__":
    main()

# =============================================================================
# END OF SCRIPT — OUTPUT SUMMARY
# =============================================================================
# All files in Analyze_avb3_LigandAndFunctionalGroupContacts_RESULTS/:
#   WHOLE-LIGAND ANALYSIS (all nine ligands):
#     avb3_contact_summary.csv
#     avb3_mn_contacts.csv
#     avb3_residue_contacts.csv
#     avb3_potential_polar_contacts.csv
#   FUNCTIONAL-GROUP ANALYSIS (eight TH derivatives, no Cilengitide):
#     avb3_functional_group_summary.csv
#     avb3_functional_group_mn_contacts.csv
#     avb3_functional_group_residue_contacts.csv
#     avb3_functional_group_polar_contacts.csv
# =============================================================================

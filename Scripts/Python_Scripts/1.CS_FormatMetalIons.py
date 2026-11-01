#!/usr/bin/env python3
import os
import sys
import glob

MG_RESIDUE_START = 9001   # starting residue number for Mg ions
CA_RESIDUE_START = 9101   # starting residue number for Ca ions

def process_file(input_file):
    output_file = "IonFixed_" + os.path.basename(input_file)
    mg_counter = 0
    ca_counter = 0

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if not line.startswith("HETATM"):
                outfile.write(line)
                continue

            atom_name = line[12:16].strip().upper()
            res_name  = line[17:20].strip().upper()
            element   = line[76:78].strip().upper()

            is_mg = (
                element == "MG" or
                atom_name.startswith("MG") or
                res_name == "MG"
            )

            is_ca = (
                element == "CA" or
                atom_name.startswith("CA") or
                res_name == "CA"
            )

            if is_mg or is_ca:
                atom_serial = int(line[6:11])
                chain_id = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = float(line[54:60])
                bfac = float(line[60:66])

                if is_mg:
                    new_resnum = MG_RESIDUE_START + mg_counter
                    mg_counter += 1

                    new_line = (
                        f"HETATM{atom_serial:5d} "
                        f"{'MG+2':>4} "
                        f"{'MG2':>3} {chain_id}"
                        f"{new_resnum:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          MG\n"
                    )
                    outfile.write(new_line)
                    outfile.write("TER\n")

                elif is_ca:
                    new_resnum = CA_RESIDUE_START + ca_counter
                    ca_counter += 1

                    new_line = (
                        f"HETATM{atom_serial:5d} "
                        f"{'CA+2':>4} "
                        f"{'CA2':>3} {chain_id}"
                        f"{new_resnum:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          CA\n"
                    )
                    outfile.write(new_line)
                    outfile.write("TER\n")
            else:
                outfile.write(line)

    print(f"Saved: {output_file}")

def main():
    files = sys.argv[1:]
    if not files:
        files = sorted(glob.glob("*.pdb"))

    if not files:
        print("No matching PDB files found.")
        return

    print("Files to process:", files)

    for f in files:
        process_file(f)

if __name__ == "__main__":
    main()
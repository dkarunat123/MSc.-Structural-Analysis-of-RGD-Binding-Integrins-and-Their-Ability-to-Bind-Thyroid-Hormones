#!/usr/bin/env python3
import os
import glob

# Direct RGD-contacting residues
key_residues_a = {150, 218}                       # αV (Chain A)
key_residues_b = {121, 122, 123, 214, 215, 220}  # β3 (Chain B)


def process_file(input_file):
    base_name = os.path.basename(input_file)
    folder = os.path.dirname(input_file)

    output_file = os.path.join(folder, f"renumbered_{base_name}")
    model_name = os.path.splitext(base_name)[0]
    keyres_output_file = os.path.join(folder, f"KeyRes_{model_name}.txt")

    new_lines = []
    residue_map = {}
    keyres_renumbered_A = {}
    keyres_renumbered_B = {}

    with open(input_file) as f:
        lines = f.readlines()

    # Step 1: Find last residue number in Chain A
    last_res_a = max(
        int(line[22:26])
        for line in lines
        if line.startswith("ATOM") and line[21] == "A"
    )
    res_id_counter = last_res_a

    # Step 2: Renumber Chain B → Chain A
    for line in lines:
        if not line.startswith("ATOM"):
            new_lines.append(line)
            continue

        chain = line[21]
        resnum = int(line[22:26])
        res_uid = (chain, resnum)

        if chain == "B":
            if res_uid not in residue_map:
                res_id_counter += 1
                residue_map[res_uid] = res_id_counter

                if resnum in key_residues_b:
                    keyres_renumbered_B[resnum] = res_id_counter

            new_resnum = residue_map[res_uid]

            newline = (
                line[:21] + "A" +
                f"{new_resnum:4d}" + line[26:]
            )
        else:
            if resnum in key_residues_a:
                keyres_renumbered_A[resnum] = resnum
            newline = line

        new_lines.append(newline)

    # Step 3: Save renumbered PDB
    with open(output_file, "w") as f:
        f.writelines(new_lines)

    # Step 4: Save key residue mapping
    print(f"\nKey RGD-binding residues for {input_file}:")

    with open(keyres_output_file, "w") as keyfile:
        keyfile.write("Chain A (original → new):\n")
        for orig in sorted(keyres_renumbered_A):
            print(f"αV: {orig} → {keyres_renumbered_A[orig]}")
            keyfile.write(f"αV: {orig:<5} → {keyres_renumbered_A[orig]}\n")

        keyfile.write("\nChain B (original → new in merged A):\n")
        for orig, new in sorted(keyres_renumbered_B.items()):
            print(f"β3: {orig} → {new}")
            keyfile.write(f"β3: {orig:<5} → {new}\n")

    print(f"\nRenumbered PDB: {output_file}")
    print(f"Key residues file: {keyres_output_file}")


def choose_from_list(items, prompt):
    for i, item in enumerate(items, 1):
        print(f"{i}. {item}")

    choice = input(prompt).strip()

    if choice.isdigit():
        idx = int(choice) - 1
        if 0 <= idx < len(items):
            return items[idx]

    if choice in items:
        return choice

    print("Invalid selection.")
    return None


def main():
    # Step 1: Find _PDBFiles folders
    pdb_folders = sorted([
        folder for folder in glob.glob("*_PDBFiles")
        if os.path.isdir(folder)
    ])

    if not pdb_folders:
        print("No folders ending in _PDBFiles found.")
        return

    print("Available _PDBFiles folders:")
    selected_folder = choose_from_list(
        pdb_folders,
        "\nEnter the number OR name of the folder: "
    )

    if selected_folder is None:
        return

    # Step 2: Find MNFixed file inside folder
    mnfixed_files = glob.glob(os.path.join(selected_folder, "MNFixed_*.pdb"))

    if not mnfixed_files:
        print(f"No MNFixed_*.pdb file found in {selected_folder}.")
        return

    if len(mnfixed_files) > 1:
        print("Multiple MNFixed files found, using the first one:")
    
    selected_pdb = mnfixed_files[0]
    print(f"\nUsing file: {os.path.basename(selected_pdb)}")

    # Step 3: Process file
    process_file(selected_pdb)


if __name__ == "__main__":
    main()

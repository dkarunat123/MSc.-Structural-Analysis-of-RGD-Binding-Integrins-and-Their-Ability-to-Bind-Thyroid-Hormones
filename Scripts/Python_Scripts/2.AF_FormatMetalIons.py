#!/usr/bin/env python3
import os
import glob
import re

MN_RESIDUE_START = 9001


def process_file(input_file):
    folder = os.path.dirname(input_file)
    filename = os.path.basename(input_file)
    output_file = os.path.join(folder, "MNFixed_" + filename)

    mn_counter = 0

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            is_mn = (
                line.startswith("HETATM") and
                (
                    line[76:78].strip().upper() == "MN" or
                    line[12:16].strip().upper().startswith("MN")
                )
            )

            if is_mn:
                atom_serial = int(line[6:11])
                chain_id = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occ = float(line[54:60])
                bfac = float(line[60:66])

                new_resnum = MN_RESIDUE_START + mn_counter
                mn_counter += 1

                new_line = (
                    f"HETATM{atom_serial:5d} "
                    f"{'MN+2':>4} "
                    f"{'MN2':>3} {chain_id}"
                    f"{new_resnum:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{occ:6.2f}{bfac:6.2f}          MN\n"
                )

                outfile.write(new_line)
                outfile.write("TER\n")
            else:
                outfile.write(line)

    print(f"\nSaved: {output_file}")


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


def get_models_from_ranking(selected_folder):
    ranking_files = glob.glob(os.path.join(selected_folder, "*_Model_Ranking.txt"))

    if not ranking_files:
        print(f"No Model_Ranking.txt file found in {selected_folder}.")
        return []

    ranking_file = ranking_files[0]
    print(f"\nUsing ranking file: {os.path.basename(ranking_file)}")

    ranked_models = []

    with open(ranking_file, "r") as f:
        for line in f:
            line = line.strip()

            if "_model_" not in line:
                continue

            parts = line.split()
            if not parts:
                continue

            model_name = parts[0]
            model_name = re.sub(r"\.(pdb|cif)$", "", model_name, flags=re.IGNORECASE)

            if model_name not in ranked_models:
                ranked_models.append(model_name)

    if not ranked_models:
        print("Could not find ranked models in the ranking file.")

    return ranked_models


def model_number_from_name(model_name):
    match = re.search(r"_model_(\d+)$", model_name)
    if match:
        return match.group(1)
    return None


def choose_model_from_ranking(selected_folder):
    ranked_models = get_models_from_ranking(selected_folder)

    if not ranked_models:
        return None

    best_model = ranked_models[0]
    best_model_number = model_number_from_name(best_model)

    print(f"\nBest-ranked model according to the ranking file: {best_model}.pdb")
    if best_model_number is not None:
        print(f"Best-ranked model number: {best_model_number}")

    use_best = input("\nDo you want to use this best-ranked model? Enter Y for yes or N to choose another model: ").strip().lower()

    if use_best in ["", "y", "yes"]:
        return best_model

    if use_best not in ["n", "no"]:
        print("Invalid selection.")
        return None

    available_model_numbers = []
    model_number_to_name = {}

    for model_name in ranked_models:
        model_number = model_number_from_name(model_name)
        if model_number is not None:
            available_model_numbers.append(model_number)
            model_number_to_name[model_number] = model_name

    if available_model_numbers:
        print("\nAvailable model numbers from the ranking file:")
        print(", ".join(available_model_numbers))

    model_choice = input("Enter the model number you want to use, usually between 0 and 4: ").strip()

    if model_choice in model_number_to_name:
        return model_number_to_name[model_choice]

    # Fallback: let the user enter the full model name, without .pdb
    if model_choice in ranked_models:
        return model_choice

    print("Invalid model selection.")
    return None


def main():
    pdb_folders = sorted([
        folder for folder in glob.glob("*_PDBFiles")
        if os.path.isdir(folder)
    ])

    if not pdb_folders:
        print("No folders ending in _PDBFiles were found in this directory.")
        return

    print("Available _PDBFiles folders:")
    selected_folder = choose_from_list(
        pdb_folders,
        "\nEnter the number OR name of the folder: "
    )

    if selected_folder is None:
        return

    selected_model = choose_model_from_ranking(selected_folder)

    if selected_model is None:
        return

    selected_pdb_path = os.path.join(selected_folder, selected_model + ".pdb")

    if not os.path.exists(selected_pdb_path):
        print(f"Selected PDB file not found: {selected_pdb_path}")
        return

    print(f"\nSelected model: {selected_model}.pdb")

    process_file(selected_pdb_path)


if __name__ == "__main__":
    main()

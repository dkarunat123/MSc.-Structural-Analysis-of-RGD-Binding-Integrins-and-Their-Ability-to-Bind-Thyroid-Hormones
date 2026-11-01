#!/usr/bin/env python3

import subprocess
import shutil
import os
import glob
import math


R = 0.0019872041  # kcal/mol/K
T = 298.15        # Kelvin

LIGAND_ORDER = [
    "Cilengitide",
    "T3",
    "T4",
    "TRIAC",
    "TETRAC",
    "sT3",
    "sT4",
    "rT3",
    "T4AM",
]

LIGAND_ALIASES = {
    "TERAC": "TETRAC",   # catches accidental TERAC typo
}

TARGET_FOLDER_PREFIXES = ("SA_", "SU_", "HU_")


def check_prodigy_environment():
    if not shutil.which("prodigy_lig"):
        print("\nERROR: PRODIGY-LIG was not found.\n")
        print("Activate the correct environment first:")
        print("    source ~/downloads/prodigy_env/bin/activate\n")
        return False
    return True


def show_current_directory_contents():
    print("\nCurrent directory:")
    print(f"    {os.getcwd()}\n")

    for item in sorted(os.listdir(".")):
        print(f"[DIR]  {item}" if os.path.isdir(item) else f"[FILE] {item}")

    print()


def extract_electrostatic_energy(pdb_file):
    try:
        with open(pdb_file, "r") as file:
            for line in file:
                if line.startswith("REMARK energies:"):
                    values = line.replace("REMARK energies:", "").strip().split(",")
                    return values[6].strip()
    except:
        return None

    return None


def convert_dg_to_kd(dg_value):
    try:
        return math.exp(float(dg_value) / (R * T))
    except:
        return None


def format_kd(kd_molar):
    return "ERROR" if kd_molar is None else f"{kd_molar:.2e}"


def get_job_name(pdb_file):
    return os.path.splitext(os.path.basename(pdb_file))[0]


def get_ligand_name(job_name):
    raw_ligand = job_name.split("_")[-1]
    return LIGAND_ALIASES.get(raw_ligand, raw_ligand)


def sort_results(results):
    return sorted(
        results,
        key=lambda row: (
            LIGAND_ORDER.index(get_ligand_name(row[0]))
            if get_ligand_name(row[0]) in LIGAND_ORDER
            else len(LIGAND_ORDER)
        )
    )


def run_prodigy_lig(pdb_file, chain_parts, electrostatic_energy):
    cmd = ["prodigy_lig", "-i", pdb_file, "-c", chain_parts[0], chain_parts[1]]

    if electrostatic_energy:
        cmd.extend(["-e", electrostatic_energy])

    result = subprocess.run(cmd, capture_output=True, text=True)

    dg_prediction = None

    for line in result.stdout.splitlines():
        parts = line.split()

        if len(parts) >= 3 and parts[0] == get_job_name(pdb_file):
            dg_prediction = parts[1]

    return dg_prediction, result.stderr


def process_pdb_files(pdb_files, chain_parts):
    results = []

    for pdb_file in pdb_files:
        job_name = get_job_name(pdb_file)
        elec = extract_electrostatic_energy(pdb_file)

        if elec:
            print(f"{job_name}: found electrostatic energy {elec}")
        else:
            print(f"{job_name}: no electrostatic energy found, running without -e")

        dg_prediction, stderr = run_prodigy_lig(pdb_file, chain_parts, elec)

        if dg_prediction:
            kd_formatted = format_kd(convert_dg_to_kd(dg_prediction))
            results.append((job_name, dg_prediction, kd_formatted))
            print(f"{job_name}: DGprediction = {dg_prediction}, Kd (M) = {kd_formatted}")
        else:
            results.append((job_name, "ERROR", "ERROR"))
            print(f"{job_name}: ERROR reading DGprediction")

        if stderr:
            print("Warnings/Errors:")
            print(stderr)

    return sort_results(results)


def print_kd_list_to_terminal(folder_name, results):
    print("\n" + "-" * 60)
    print(f"Copy/paste Kd(M) list for: {folder_name}")
    print("-" * 60)
    print("Kd(M)")

    for job_name, dg_prediction, kd_formatted in results:
        print(kd_formatted)

    print("-" * 60 + "\n")


def write_results(output_file, results):
    with open(output_file, "w") as out:
        out.write("Job name(s)\tDGprediction (Kcal/mol)\tKd (M)\n")

        for job_name, dg_prediction, kd_formatted in results:
            out.write(f"{job_name}\t{dg_prediction}\t{kd_formatted}\n")

        out.write("\nKd(M)\n")

        for job_name, dg_prediction, kd_formatted in results:
            out.write(f"{kd_formatted}\n")


def get_target_folders(root_directory):
    target_folders = []

    for item in sorted(os.listdir(root_directory)):
        full_path = os.path.join(root_directory, item)

        if os.path.isdir(full_path) and item.startswith(TARGET_FOLDER_PREFIXES):
            pdb_files = sorted(glob.glob(os.path.join(full_path, "*.pdb")))

            if pdb_files:
                target_folders.append((full_path, pdb_files))

    return target_folders


def main():
    print("\n=== PRODIGY-LIG Batch Binding Affinity Calculator ===\n")

    if not check_prodigy_environment():
        return

    print("Choose input mode:")
    print("1. Single PDB file")
    print("2. Multiple PDB files entered manually")
    print("3. All PDB files in one directory")
    print("4. All PDB files in SA_/SU_/HU_ folders only")

    mode = input("\nEnter 1, 2, 3, or 4: ").strip()

    show_current_directory_contents()

    use_default_chains = input(
        "\nUse default chain specification A B:UNL for all files? (Y/N): "
    ).strip().upper()

    if use_default_chains == "Y" or use_default_chains == "":
        chains = "A B:UNL"
    else:
        chains = input("Enter chain specification: ").strip()

    chain_parts = chains.split()

    if len(chain_parts) != 2:
        print("\nERROR: chain specification must look like:")
        print("    A B:UNL")
        return

    if mode == "1":
        pdb_files = [input("\nEnter the PDB file: ").strip()]
        pdb_files = [pdb for pdb in pdb_files if os.path.isfile(pdb)]

        if not pdb_files:
            print("\nERROR: No valid PDB files found.")
            return

        output_file = input(
            "\nEnter output txt filename, or press Enter for Prodigy_Lig_Results.txt: "
        ).strip() or "Prodigy_Lig_Results.txt"

        results = process_pdb_files(pdb_files, chain_parts)
        results = sort_results(results)
        write_results(output_file, results)
        print_kd_list_to_terminal("single_file", results)

    elif mode == "2":
        print("\nEnter PDB filenames one at a time.")
        print("Press Enter on a blank line when finished.\n")

        pdb_files = []

        while True:
            pdb_file = input("PDB file: ").strip()

            if pdb_file == "":
                break

            pdb_files.append(pdb_file)

        pdb_files = [pdb for pdb in pdb_files if os.path.isfile(pdb)]

        if not pdb_files:
            print("\nERROR: No valid PDB files found.")
            return

        output_file = input(
            "\nEnter output txt filename, or press Enter for Prodigy_Lig_Results.txt: "
        ).strip() or "Prodigy_Lig_Results.txt"

        results = process_pdb_files(pdb_files, chain_parts)
        results = sort_results(results)
        write_results(output_file, results)
        print_kd_list_to_terminal("manual_files", results)

    elif mode == "3":
        directory = input("\nEnter directory name/path: ").strip() or "."
        pdb_files = sorted(glob.glob(os.path.join(directory, "*.pdb")))

        if not pdb_files:
            print("\nERROR: No valid PDB files found.")
            return

        folder_name = os.path.basename(os.path.abspath(directory))
        output_file = os.path.join(directory, f"{folder_name}_Prodigy_Lig_Results.txt")

        results = process_pdb_files(pdb_files, chain_parts)
        results = sort_results(results)
        write_results(output_file, results)
        print_kd_list_to_terminal(folder_name, results)

    elif mode == "4":
        root_directory = input(
            "\nEnter parent directory, or press Enter for current directory: "
        ).strip() or "."

        target_folders = get_target_folders(root_directory)

        if not target_folders:
            print("\nERROR: No SA_/SU_/HU_ folders containing PDB files were found.")
            return

        print(f"\nFound {len(target_folders)} target folder(s):\n")

        for folder, pdb_files in target_folders:
            print(f"  {os.path.basename(folder)}  ({len(pdb_files)} PDB files)")

        confirm = input(
            "\nContinue and run PRODIGY-LIG on these folders? (Y/N): "
        ).strip().upper()

        if confirm != "Y":
            print("\nCancelled.")
            return

        for folder, pdb_files in target_folders:
            folder_name = os.path.basename(os.path.abspath(folder))
            output_file = os.path.join(folder, f"{folder_name}_Prodigy_Lig_Results.txt")

            print("\n" + "=" * 80)
            print(f"Processing folder: {folder}")
            print(f"Output file: {output_file}")
            print("=" * 80 + "\n")

            results = process_pdb_files(pdb_files, chain_parts)
            results = sort_results(results)
            write_results(output_file, results)
            print_kd_list_to_terminal(folder_name, results)

    else:
        print("\nERROR: Please enter 1, 2, 3, or 4.")
        return

    print("\nDone!")

    print("\nIf you are finished using PRODIGY-LIG, type:")
    print("    deactivate")


if __name__ == "__main__":
    main()
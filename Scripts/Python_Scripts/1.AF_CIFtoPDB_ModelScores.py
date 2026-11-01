from Bio.PDB import MMCIFParser, PDBIO
import zipfile
import os
import glob
import json
import re


def convert_cif_to_pdb(cif_path, output_path):
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    try:
        structure = parser.get_structure("model", cif_path)
        io.set_structure(structure)
        io.save(output_path)
        print(f"Converted: {os.path.basename(cif_path)} → {os.path.basename(output_path)}")
    except Exception as e:
        print(f"Failed to convert {cif_path}: {e}")


def extract_plddt_from_pdb(pdb_file):
    plddt_scores = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_scores.append(b_factor)
                except ValueError:
                    continue

    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    else:
        return None


def extract_scores_from_json(json_file):
    try:
        with open(json_file) as f:
            data = json.load(f)
            iptm = data.get("iptm", None)
            ptm = data.get("ptm", None)
            return iptm, ptm
    except Exception:
        return None, None


def main():
    zip_files = glob.glob("*.zip")

    if not zip_files:
        print("No zip files found in this directory.")
        return

    print("Available zip files:")
    for i, z in enumerate(zip_files, 1):
        print(f"{i}. {z}")

    choice = input("\nEnter the number OR name of the zip file: ").strip()

    if choice.isdigit():
        idx = int(choice) - 1
        if idx < 0 or idx >= len(zip_files):
            print("Invalid selection.")
            return
        zip_name = zip_files[idx]
    else:
        zip_name = choice

    if not os.path.exists(zip_name):
        print("Zip file not found.")
        return

    base_name = zip_name.replace(".zip", "")
    extract_folder = base_name + "_extracted"
    output_folder = base_name + "_PDBFiles"

    os.makedirs(extract_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)

    print("\nExtracting zip file...")
    with zipfile.ZipFile(zip_name, "r") as zip_ref:
        zip_ref.extractall(extract_folder)

    print("Converting AlphaFold model CIF files to PDB...\n")

    converted_models = []

    for root, dirs, files in os.walk(extract_folder):
        for file in files:
            if file.endswith(".cif") and "_model_" in file:
                cif_path = os.path.join(root, file)
                pdb_name = file.replace(".cif", ".pdb")
                output_path = os.path.join(output_folder, pdb_name)

                convert_cif_to_pdb(cif_path, output_path)
                converted_models.append((file, output_path))

    print(f"\nDone. Model PDB files saved in: {output_folder}")

    results = []

    for cif_file, pdb_path in converted_models:
        match = re.match(r"(.+?)_model_(\d+)\.cif", cif_file)

        if not match:
            continue

        model_base_name = match.group(1)
        model_index = match.group(2)

        model_tag = f"{model_base_name}_model_{model_index}"
        json_name = f"{model_base_name}_summary_confidences_{model_index}.json"

        json_path = None

        for root, dirs, files in os.walk(extract_folder):
            if json_name in files:
                json_path = os.path.join(root, json_name)
                break

        avg_plddt = extract_plddt_from_pdb(pdb_path)

        if json_path is not None:
            iptm, ptm = extract_scores_from_json(json_path)
        else:
            print(f"JSON file not found for {model_tag}; ipTM and pTM will be N/A")
            iptm, ptm = None, None

        results.append({
            "Model": model_tag,
            "pLDDT": round(avg_plddt, 2) if avg_plddt is not None else "N/A",
            "ipTM": round(iptm, 3) if iptm is not None else "N/A",
            "pTM": round(ptm, 3) if ptm is not None else "N/A"
        })

    if not results:
        print("\nNo model score results found.")
        return

    results = sorted(
        results,
        key=lambda x: x["pLDDT"] if isinstance(x["pLDDT"], float) else 0,
        reverse=True
    )

    print("\nModel Comparison Summary:")
    header = f"{'Model':<45} {'pLDDT':>10} {'ipTM':>10} {'pTM':>10}"
    separator = "-" * 80

    print(header)
    print(separator)

    for r in results:
        print(f"{r['Model']:<45} {r['pLDDT']:>10} {r['ipTM']:>10} {r['pTM']:>10}")

    output_txt = os.path.join(output_folder, f"{base_name}_Model_Ranking.txt")

    with open(output_txt, "w") as f:
        f.write("Model Comparison Summary:\n")
        f.write(header + "\n")
        f.write(separator + "\n")

        for r in results:
            line = f"{r['Model']:<45} {r['pLDDT']:>10} {r['ipTM']:>10} {r['pTM']:>10}\n"
            f.write(line)

    print(f"\nRanking saved to: {output_txt}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
FindAFScores.py

Run this script from the directory that contains one subfolder per receptor.

Workflow
--------
1. Lists all subfolders in the current directory.
2. Lets you choose one receptor or all receptors.
3. Lets you choose model numbers interactively or from a text file.
4. For each chosen receptor/model, finds these matching original AF3 files:
       *_summary_confidences_<model>.json
       *_full_data_<model>.json
       *_model_<model>.cif
5. Extracts mean pLDDT, ipTM, pTM, and median PAE.
6. Saves one PAE heatmap inside each receptor folder.
7. Creates/updates one cumulative AlphaFold_Scores.csv in the current directory.

Text-file format
----------------
One folder and model number per line, separated by a comma, tab, or spaces:

    HS_avb1,0
    HS_avb3,1
    NV_a3b2,0

Blank lines and lines beginning with # are ignored.
"""

from __future__ import annotations

import csv
import json
import re
import shlex
import sys
from pathlib import Path
from typing import Any

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    sys.exit(
        "ERROR: matplotlib/numpy are required. On Ubuntu install with:\n"
        "    sudo apt update\n"
        "    sudo apt install python3-matplotlib python3-numpy"
    )

SUMMARY_RE = re.compile(
    r"^(?P<prefix>.+)_summary_confidences_(?P<model>[0-4])\.json$"
)

CSV_FIELDS = [
    "Receptor folder",
    "AlphaFold prefix",
    "Model number",
    "Mean pLDDT",
    "ipTM",
    "pTM",
    "Median PAE (Å)",
    "PAE heatmap",
]


def list_receptor_folders(base: Path) -> list[Path]:
    excluded = {"__pycache__", ".git", "AlphaFold_Score_Results"}
    return sorted(
        p for p in base.iterdir()
        if p.is_dir() and not p.name.startswith(".") and p.name not in excluded
    )


def choose_folders(folders: list[Path]) -> list[Path]:
    print("\nReceptor subfolders found:\n")
    for i, folder in enumerate(folders, 1):
        print(f"  {i:>2}. {folder.name}")

    print(
        "\nWhich receptor(s) would you like to process?\n"
        "Enter a folder number, the exact folder name, or 'all'."
    )

    while True:
        choice = input("Selection: ").strip()

        if choice.lower() in {"all", "a", "*"}:
            return folders

        if choice.isdigit():
            i = int(choice)
            if 1 <= i <= len(folders):
                return [folders[i - 1]]

        for folder in folders:
            if choice == folder.name:
                return [folder]

        print("Invalid selection. Enter a listed number/name or 'all'.")


def choose_model_method() -> str:
    print(
        "\nHow would you like to choose the AlphaFold model(s)?\n"
        "  1. Enter model numbers interactively\n"
        "  2. Read model numbers from a text file"
    )
    while True:
        choice = input("Selection [1/2]: ").strip()
        if choice == "1":
            return "interactive"
        if choice == "2":
            return "file"
        print("Please enter 1 or 2.")


def ask_model(folder_name: str) -> int | None:
    while True:
        value = input(
            f"Model for {folder_name} [0-4, or 'skip']: "
        ).strip().lower()
        if value in {"skip", "s"}:
            return None
        if value in {"0", "1", "2", "3", "4"}:
            return int(value)
        print("Enter 0, 1, 2, 3, 4, or 'skip'.")


def parse_mapping_file(path: Path) -> dict[str, int]:
    mapping: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = [p.strip() for p in line.split(",")] if "," in line else line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Line {line_number}: expected 'folder,model', got {raw.rstrip()!r}"
                )

            folder, model_text = parts
            if model_text not in {"0", "1", "2", "3", "4"}:
                raise ValueError(
                    f"Line {line_number}: model must be 0-4, got {model_text!r}"
                )
            mapping[folder] = int(model_text)

    if not mapping:
        raise ValueError("The mapping file contained no usable entries.")
    return mapping


def ask_mapping_file(base: Path) -> dict[str, int]:
    while True:
        raw = input("Model-selection text file: ").strip().strip('"').strip("'")
        path = Path(raw).expanduser()
        if not path.is_absolute():
            path = base / path
        path = path.resolve()

        if not path.is_file():
            print(f"File not found: {path}")
            continue

        try:
            return parse_mapping_file(path)
        except (OSError, ValueError) as exc:
            print(f"Could not read mapping file: {exc}")


def build_plan(
    folders: list[Path], method: str, base: Path
) -> list[tuple[Path, int]]:
    plan: list[tuple[Path, int]] = []

    if method == "interactive":
        for folder in folders:
            model = ask_model(folder.name)
            if model is not None:
                plan.append((folder, model))
        return plan

    mapping = ask_mapping_file(base)
    for folder in folders:
        model = mapping.get(folder.name)
        if model is None:
            print(f"WARNING: no model listed for {folder.name}; skipping.")
            continue
        plan.append((folder, model))
    return plan


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def split_cif_line(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = "#"
    return list(lexer)


def read_cif_loop(
    lines: list[str], start: int, columns: int
) -> tuple[list[list[str]], int]:
    rows: list[list[str]] = []
    tokens: list[str] = []
    i = start

    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped or stripped == "#":
            i += 1
            if tokens:
                break
            continue
        if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
            break

        tokens.extend(split_cif_line(lines[i]))
        while len(tokens) >= columns:
            rows.append(tokens[:columns])
            tokens = tokens[columns:]
        i += 1

    return rows, i


def extract_mean_plddt(cif_path: Path) -> float:
    lines = cif_path.read_text(encoding="utf-8", errors="replace").splitlines()
    metric_names: dict[str, str] = {}
    global_values: list[tuple[str, float]] = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        if line == "loop_":
            headers: list[str] = []
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                headers.append(lines[j].strip())
                j += 1

            if headers and headers[0].startswith("_ma_qa_metric."):
                rows, i = read_cif_loop(lines, j, len(headers))
                idx = {h: n for n, h in enumerate(headers)}
                id_col = idx.get("_ma_qa_metric.id")
                name_col = idx.get("_ma_qa_metric.name")
                if id_col is not None and name_col is not None:
                    for row in rows:
                        if len(row) == len(headers):
                            metric_names[row[id_col]] = row[name_col].strip("'\"")
                continue

            if headers and headers[0].startswith("_ma_qa_metric_global."):
                rows, i = read_cif_loop(lines, j, len(headers))
                idx = {h: n for n, h in enumerate(headers)}
                metric_col = idx.get("_ma_qa_metric_global.metric_id")
                value_col = idx.get("_ma_qa_metric_global.metric_value")
                if metric_col is not None and value_col is not None:
                    for row in rows:
                        if len(row) == len(headers):
                            try:
                                global_values.append((row[metric_col], float(row[value_col])))
                            except ValueError:
                                pass
                continue

        if line.startswith("_ma_qa_metric_global.metric_id"):
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                metric_id = parts[1].strip("'\"")
                for k in range(i + 1, min(i + 10, len(lines))):
                    candidate = lines[k].strip()
                    if candidate.startswith("_ma_qa_metric_global.metric_value"):
                        value_parts = candidate.split(maxsplit=1)
                        if len(value_parts) == 2:
                            try:
                                global_values.append(
                                    (metric_id, float(value_parts[1].strip("'\"")))
                                )
                            except ValueError:
                                pass
                        break
        i += 1

    for metric_id, value in global_values:
        if "plddt" in metric_names.get(metric_id, "").lower():
            return value

    if global_values:
        return global_values[0][1]

    raise ValueError(f"Could not find global mean pLDDT in {cif_path.name}")


def find_file_triplets(
    folder: Path, model: int
) -> list[tuple[str, Path, Path, Path]]:
    triplets: list[tuple[str, Path, Path, Path]] = []

    for summary in sorted(folder.glob(f"*_summary_confidences_{model}.json")):
        match = SUMMARY_RE.match(summary.name)
        if not match:
            continue

        prefix = match.group("prefix")
        full_data = folder / f"{prefix}_full_data_{model}.json"
        cif = folder / f"{prefix}_model_{model}.cif"

        if full_data.is_file() and cif.is_file():
            triplets.append((prefix, summary, full_data, cif))

    return triplets


def choose_triplet(
    folder: Path, model: int, triplets: list[tuple[str, Path, Path, Path]]
) -> tuple[str, Path, Path, Path]:
    if len(triplets) == 1:
        return triplets[0]

    print(f"Multiple complete model-{model} file sets found in {folder.name}:")
    for i, triplet in enumerate(triplets, 1):
        print(f"  {i}. {triplet[0]}")

    while True:
        choice = input("Choose file-set number or exact prefix: ").strip()
        if choice.isdigit():
            i = int(choice)
            if 1 <= i <= len(triplets):
                return triplets[i - 1]
        for triplet in triplets:
            if choice == triplet[0]:
                return triplet
        print("Invalid choice.")


def chain_boundaries(chain_ids: list[str]) -> tuple[list[int], list[tuple[str, int, int]]]:
    if not chain_ids:
        return [], []

    boundaries: list[int] = []
    spans: list[tuple[str, int, int]] = []
    start = 0
    current = chain_ids[0]

    for i, chain in enumerate(chain_ids[1:], 1):
        if chain != current:
            spans.append((current, start, i))
            boundaries.append(i)
            start = i
            current = chain
    spans.append((current, start, len(chain_ids)))
    return boundaries, spans


def save_pae_heatmap(
    pae: np.ndarray, chain_ids: list[str], output: Path, title: str
) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(
        pae,
        origin="upper",
        interpolation="nearest",
        aspect="equal",
        vmin=0,
        vmax=31.75,
    )

    boundaries, spans = chain_boundaries(chain_ids)
    for boundary in boundaries:
        ax.axhline(boundary - 0.5, linewidth=0.8)
        ax.axvline(boundary - 0.5, linewidth=0.8)

    if spans and len(spans) <= 20:
        centres = [(start + end - 1) / 2 for _, start, end in spans]
        labels = [chain for chain, _, _ in spans]
        ax.set_xticks(centres)
        ax.set_yticks(centres)
        ax.set_xticklabels(labels, rotation=90)
        ax.set_yticklabels(labels)

    ax.set_xlabel("Scored residue/token")
    ax.set_ylabel("Aligned residue/token")
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Predicted aligned error (Å)")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def process_one(
    base: Path,
    folder: Path,
    model: int,
    prefix: str,
    summary_path: Path,
    full_data_path: Path,
    cif_path: Path,
) -> dict[str, Any]:
    summary = read_json(summary_path)
    full_data = read_json(full_data_path)

    pae_raw = full_data.get("pae")
    if not isinstance(pae_raw, list) or not pae_raw:
        raise ValueError(f"No PAE matrix found in {full_data_path.name}")

    pae = np.asarray(pae_raw, dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE matrix is not square in {full_data_path.name}")

    chain_ids = full_data.get("token_chain_ids")
    if not isinstance(chain_ids, list) or len(chain_ids) != pae.shape[0]:
        chain_ids = [""] * pae.shape[0]

    mean_plddt = extract_mean_plddt(cif_path)
    median_pae = float(np.nanmedian(pae))

    heatmap = folder / f"{prefix}_model_{model}_PAE.png"
    save_pae_heatmap(pae, chain_ids, heatmap, f"{folder.name}, model {model}")

    return {
        "Receptor folder": folder.name,
        "AlphaFold prefix": prefix,
        "Model number": model,
        "Mean pLDDT": round(mean_plddt, 3),
        "ipTM": summary.get("iptm"),
        "pTM": summary.get("ptm"),
        "Median PAE (Å)": round(median_pae, 3),
        "PAE heatmap": str(heatmap.relative_to(base)),
    }


def load_existing(csv_path: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
    rows: dict[tuple[str, str, str], dict[str, Any]] = {}
    if not csv_path.is_file():
        return rows

    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                row.get("Receptor folder", ""),
                row.get("AlphaFold prefix", ""),
                row.get("Model number", ""),
            )
            rows[key] = row
    return rows


def save_csv(rows: dict[tuple[str, str, str], dict[str, Any]], csv_path: Path) -> None:
    ordered = sorted(
        rows.values(),
        key=lambda row: (
            str(row["Receptor folder"]),
            int(row["Model number"]),
            str(row["AlphaFold prefix"]),
        ),
    )

    temp = csv_path.with_suffix(".csv.tmp")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(ordered)
    temp.replace(csv_path)


def row_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row["Receptor folder"]),
        str(row["AlphaFold prefix"]),
        str(row["Model number"]),
    )


def main() -> int:
    base = Path.cwd()
    csv_path = base / "AlphaFold_Scores.csv"

    print("\nAlphaFold Score Extractor")
    print("=========================")
    print(f"Current directory: {base}")

    folders = list_receptor_folders(base)
    if not folders:
        print("ERROR: no receptor subfolders found in the current directory.")
        return 1

    selected = choose_folders(folders)
    method = choose_model_method()
    plan = build_plan(selected, method, base)

    if not plan:
        print("Nothing selected for processing.")
        return 1

    print("\nProcessing plan:")
    for folder, model in plan:
        print(f"  {folder.name}: model {model}")

    existing = load_existing(csv_path)
    successes = 0
    skipped = 0
    errors = 0

    total = len(plan)
    for index, (folder, model) in enumerate(plan, 1):
        print(f"\n[{index}/{total}] Processing {folder.name}, model {model}...")
        triplets = find_file_triplets(folder, model)

        if not triplets:
            print(f"  SKIPPED: no complete model-{model} file triplet found.")
            skipped += 1
            continue

        prefix, summary, full_data, cif = choose_triplet(folder, model, triplets)
        print("  Using:")
        print(f"    {summary.name}")
        print(f"    {full_data.name}")
        print(f"    {cif.name}")

        try:
            row = process_one(base, folder, model, prefix, summary, full_data, cif)
            existing[row_key(row)] = row
            save_csv(existing, csv_path)
            successes += 1
            print(
                "  SUCCESS\n"
                f"    Mean pLDDT: {row['Mean pLDDT']}\n"
                f"    ipTM: {row['ipTM']}\n"
                f"    pTM: {row['pTM']}\n"
                f"    Median PAE: {row['Median PAE (Å)']} Å\n"
                f"    CSV updated: {csv_path.name}"
            )
        except Exception as exc:
            errors += 1
            print(f"  ERROR: {exc}")

    print("\nFinished")
    print("========")
    print(f"Successful: {successes}")
    print(f"Skipped:    {skipped}")
    print(f"Errors:     {errors}")
    print(f"Cumulative CSV: {csv_path}")

    return 0 if successes else 1


if __name__ == "__main__":
    raise SystemExit(main())

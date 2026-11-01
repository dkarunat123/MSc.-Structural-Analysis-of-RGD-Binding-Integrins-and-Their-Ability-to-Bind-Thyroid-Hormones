"""
compare_af_crystal_batch_interactive.py

Interactive PyMOL tool for comparing multiple already-loaded crystal and
AlphaFold receptor structures in one session.

FEATURES
--------
- Lists all currently loaded PyMOL objects.
- Lets the user select:
    1. A crystal-structure object
    2. Its corresponding AlphaFold object
- Calculates:
    - whole-structure unrefined C-alpha RMSD
    - whole-structure refined C-alpha RMSD
    - number and percentage of C-alpha atoms removed during refinement
    - optional user-defined residue-region RMSDs
- Appends/updates results in one cumulative CSV:
      RMSD_Comparisons.csv
- After each comparison, asks whether the user wants to compare another pair.
- Does not require PyMOL or the script to be restarted.
- Assumes HADDOCK preprocessing may have merged each receptor into one
  continuous chain, so no chain-by-chain comparison is performed.

LOAD IN PYMOL
-------------
    run C:/full/path/to/compare_af_crystal_batch_interactive.py

START BATCH INTERACTIVE MODE
----------------------------
    compare_af_crystal_batch

DIRECT SINGLE-COMPARISON MODE
-----------------------------
    compare_af_crystal crystal_object, alphafold_object

OPTIONAL MATCHING REGION
------------------------
    compare_af_crystal crystal_object, alphafold_object, \
        crystal_region="crystal_object and resi 1-790", \
        af_region="alphafold_object and resi 1-790"

The two region selections must describe corresponding biological regions.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from pymol import cmd


CSV_FIELDS = [
    "Crystal object",
    "AlphaFold object",
    "Comparison",
    "Crystal selection",
    "AlphaFold selection",
    "Cycles requested",
    "Unrefined RMSD (A)",
    "Atoms before refinement",
    "Refined RMSD (A)",
    "Atoms after refinement",
    "Atoms removed",
    "Atoms removed (%)",
    "Aligned residues",
    "Alignment score",
]


def _safe_name(text: str) -> str:
    return "".join(
        char if char.isalnum() or char in {"_", "-"} else "_"
        for char in text
    )


def _loaded_objects() -> list[str]:
    return list(cmd.get_object_list("all"))


def _selection_exists(selection: str) -> bool:
    try:
        return cmd.count_atoms(selection) > 0
    except Exception:
        return False


def _align(
    mobile_selection: str,
    target_selection: str,
    cycles: int,
    alignment_object: str,
) -> dict[str, float | int]:
    result = cmd.align(
        mobile_selection,
        target_selection,
        cycles=int(cycles),
        object=alignment_object,
        transform=1,
        reset=1,
    )

    return {
        "refined_rmsd": float(result[0]),
        "refined_atom_count": int(result[1]),
        "cycles": int(result[2]),
        "unrefined_rmsd": float(result[3]),
        "unrefined_atom_count": int(result[4]),
        "alignment_score": float(result[5]),
        "aligned_residue_count": int(result[6]),
    }


def _result_row(
    crystal_object: str,
    af_object: str,
    comparison: str,
    crystal_selection: str,
    af_selection: str,
    cycles_requested: int,
    result: dict[str, float | int],
) -> dict[str, Any]:
    atoms_before = int(result["unrefined_atom_count"])
    atoms_after = int(result["refined_atom_count"])
    atoms_removed = max(atoms_before - atoms_after, 0)

    if atoms_before > 0:
        atoms_removed_percent = (atoms_removed / atoms_before) * 100
    else:
        atoms_removed_percent = 0.0

    return {
        "Crystal object": crystal_object,
        "AlphaFold object": af_object,
        "Comparison": comparison,
        "Crystal selection": crystal_selection,
        "AlphaFold selection": af_selection,
        "Cycles requested": cycles_requested,
        "Unrefined RMSD (A)": round(
            float(result["unrefined_rmsd"]), 4
        ),
        "Atoms before refinement": atoms_before,
        "Refined RMSD (A)": round(
            float(result["refined_rmsd"]), 4
        ),
        "Atoms after refinement": atoms_after,
        "Atoms removed": atoms_removed,
        "Atoms removed (%)": round(atoms_removed_percent, 2),
        "Aligned residues": result["aligned_residue_count"],
        "Alignment score": round(
            float(result["alignment_score"]), 4
        ),
    }


def _row_key(row: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(row["Crystal object"]),
        str(row["AlphaFold object"]),
        str(row["Comparison"]),
        str(row["Crystal selection"]),
        str(row["AlphaFold selection"]),
    )


def _update_cumulative_csv(
    output_path: Path,
    new_rows: list[dict[str, Any]],
) -> None:
    indexed_rows: dict[
        tuple[str, str, str, str, str],
        dict[str, Any],
    ] = {}

    if output_path.is_file():
        with output_path.open(
            "r",
            newline="",
            encoding="utf-8",
        ) as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                indexed_rows[_row_key(row)] = row

    for row in new_rows:
        indexed_rows[_row_key(row)] = row

    rows = sorted(
        indexed_rows.values(),
        key=lambda row: (
            str(row["Crystal object"]),
            str(row["AlphaFold object"]),
            str(row["Comparison"]),
        ),
    )

    temporary_path = output_path.with_suffix(".csv.tmp")

    with temporary_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CSV_FIELDS,
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary_path.replace(output_path)


def _print_loaded_objects(objects: list[str] | None = None) -> None:
    if objects is None:
        objects = _loaded_objects()

    print("\nCurrently loaded PyMOL objects:")

    if not objects:
        print("  None")
        return

    for index, object_name in enumerate(objects, start=1):
        print(f"  {index}. {object_name}")


def _qt_choose_object(
    title: str,
    label: str,
    objects: list[str],
) -> str | None:
    try:
        from pymol.Qt import QtWidgets

        parent = QtWidgets.QApplication.activeWindow()
        selected, accepted = QtWidgets.QInputDialog.getItem(
            parent,
            title,
            label,
            objects,
            0,
            False,
        )

        if accepted and selected:
            return str(selected)

        return None

    except Exception:
        return None


def _console_choose_object(
    prompt: str,
    objects: list[str],
) -> str | None:
    _print_loaded_objects(objects)

    while True:
        try:
            choice = input(
                f"\n{prompt}\n"
                "Enter an object number, exact object name, or 'cancel': "
            ).strip()
        except EOFError:
            print(
                "\nInteractive console input is unavailable. "
                "Use direct command mode instead."
            )
            return None

        if choice.lower() in {"cancel", "c", "q", "quit"}:
            return None

        if choice.isdigit():
            index = int(choice)
            if 1 <= index <= len(objects):
                return objects[index - 1]

        if choice in objects:
            return choice

        print("Invalid selection. Please try again.")


def _choose_object(
    title: str,
    label: str,
    objects: list[str],
) -> str | None:
    selected = _qt_choose_object(title, label, objects)

    if selected is not None:
        return selected

    return _console_choose_object(label, objects)


def _ask_region() -> tuple[str, str] | None:
    try:
        from pymol.Qt import QtWidgets

        parent = QtWidgets.QApplication.activeWindow()

        answer = QtWidgets.QMessageBox.question(
            parent,
            "Optional region comparison",
            "Would you like to compare a matching residue region "
            "in addition to the whole structure?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No,
        )

        if answer != QtWidgets.QMessageBox.Yes:
            return ("", "")

        crystal_region, accepted = QtWidgets.QInputDialog.getText(
            parent,
            "Crystal-region selection",
            "Enter the PyMOL selection for the crystal region:",
        )
        if not accepted:
            return None

        af_region, accepted = QtWidgets.QInputDialog.getText(
            parent,
            "AlphaFold-region selection",
            "Enter the matching PyMOL selection for the AlphaFold region:",
        )
        if not accepted:
            return None

        return str(crystal_region).strip(), str(af_region).strip()

    except Exception:
        try:
            answer = input(
                "\nCompare a matching residue region too? [y/N]: "
            ).strip().lower()
        except EOFError:
            return ("", "")

        if answer not in {"y", "yes"}:
            return ("", "")

        crystal_region = input(
            "Crystal-region PyMOL selection: "
        ).strip()
        af_region = input(
            "Matching AlphaFold-region PyMOL selection: "
        ).strip()

        return crystal_region, af_region


def _ask_continue() -> bool:
    try:
        from pymol.Qt import QtWidgets

        parent = QtWidgets.QApplication.activeWindow()
        answer = QtWidgets.QMessageBox.question(
            parent,
            "Continue RMSD comparisons",
            "Would you like to compare another crystal/AlphaFold pair?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.Yes,
        )
        return answer == QtWidgets.QMessageBox.Yes

    except Exception:
        try:
            answer = input(
                "\nCompare another crystal/AlphaFold pair? [Y/n]: "
            ).strip().lower()
        except EOFError:
            return False

        return answer not in {"n", "no", "q", "quit"}


def compare_af_crystal(
    crystal_object: str = "",
    af_object: str = "",
    crystal_region: str = "",
    af_region: str = "",
    refined_cycles: int = 5,
    output_csv: str = "RMSD_Comparisons.csv",
) -> list[dict[str, Any]] | None:
    """
    Compare one loaded crystal structure and one loaded AlphaFold model.

    Direct use:
        compare_af_crystal crystal_object, alphafold_object
    """
    if not crystal_object or not af_object:
        print(
            "\nUsage:\n"
            "  compare_af_crystal crystal_object, alphafold_object\n\n"
            "Batch interactive mode:\n"
            "  compare_af_crystal_batch\n"
        )
        _print_loaded_objects()
        return None

    loaded = _loaded_objects()

    if crystal_object not in loaded:
        raise ValueError(
            f"Crystal object {crystal_object!r} is not loaded."
        )

    if af_object not in loaded:
        raise ValueError(
            f"AlphaFold object {af_object!r} is not loaded."
        )

    if crystal_object == af_object:
        raise ValueError(
            "The crystal and AlphaFold objects must be different."
        )

    crystal_ca = (
        f"({crystal_object}) and polymer.protein and name CA"
    )
    af_ca = (
        f"({af_object}) and polymer.protein and name CA"
    )

    if not _selection_exists(crystal_ca):
        raise ValueError(
            f"No protein C-alpha atoms were found in {crystal_object!r}."
        )

    if not _selection_exists(af_ca):
        raise ValueError(
            f"No protein C-alpha atoms were found in {af_object!r}."
        )

    refined_cycles = int(refined_cycles)
    rows: list[dict[str, Any]] = []

    print(
        f"\nComparing:\n"
        f"  Crystal:   {crystal_object}\n"
        f"  AlphaFold: {af_object}"
    )

    print("\n=== Whole-structure comparison ===")

    unrefined = _align(
        mobile_selection=af_ca,
        target_selection=crystal_ca,
        cycles=0,
        alignment_object=(
            f"{_safe_name(af_object)}_to_"
            f"{_safe_name(crystal_object)}_whole_unrefined"
        ),
    )

    print(
        f"Unrefined C-alpha RMSD: "
        f"{unrefined['unrefined_rmsd']:.3f} A\n"
        f"Matched C-alpha atoms: "
        f"{unrefined['unrefined_atom_count']}\n"
        f"Aligned residues: "
        f"{unrefined['aligned_residue_count']}"
    )

    rows.append(
        _result_row(
            crystal_object=crystal_object,
            af_object=af_object,
            comparison="Whole structure - unrefined",
            crystal_selection=crystal_ca,
            af_selection=af_ca,
            cycles_requested=0,
            result=unrefined,
        )
    )

    refined = _align(
        mobile_selection=af_ca,
        target_selection=crystal_ca,
        cycles=refined_cycles,
        alignment_object=(
            f"{_safe_name(af_object)}_to_"
            f"{_safe_name(crystal_object)}_whole_refined"
        ),
    )

    print(
        f"\nRefined alignment:\n"
        f"  RMSD before refinement: "
        f"{refined['unrefined_rmsd']:.3f} A\n"
        f"  Atoms before refinement: "
        f"{refined['unrefined_atom_count']}\n"
        f"  RMSD after refinement: "
        f"{refined['refined_rmsd']:.3f} A\n"
        f"  Atoms after refinement: "
        f"{refined['refined_atom_count']}\n"
        f"  Atoms removed: "
        f"{refined['unrefined_atom_count'] - refined['refined_atom_count']}\n"
        f"  Atoms removed (%): "
        f"{((refined['unrefined_atom_count'] - refined['refined_atom_count']) / refined['unrefined_atom_count'] * 100) if refined['unrefined_atom_count'] else 0:.2f}%"
    )

    rows.append(
        _result_row(
            crystal_object=crystal_object,
            af_object=af_object,
            comparison="Whole structure - refined",
            crystal_selection=crystal_ca,
            af_selection=af_ca,
            cycles_requested=refined_cycles,
            result=refined,
        )
    )

    if bool(crystal_region) != bool(af_region):
        raise ValueError(
            "Provide both crystal_region and af_region, "
            "or leave both blank."
        )

    if crystal_region and af_region:
        crystal_region_ca = (
            f"({crystal_region}) and polymer.protein and name CA"
        )
        af_region_ca = (
            f"({af_region}) and polymer.protein and name CA"
        )

        if not _selection_exists(crystal_region_ca):
            raise ValueError(
                "The crystal-region selection contains no C-alpha atoms."
            )

        if not _selection_exists(af_region_ca):
            raise ValueError(
                "The AlphaFold-region selection contains no C-alpha atoms."
            )

        print("\n=== User-defined region comparison ===")

        region_unrefined = _align(
            mobile_selection=af_region_ca,
            target_selection=crystal_region_ca,
            cycles=0,
            alignment_object=(
                f"{_safe_name(af_object)}_to_"
                f"{_safe_name(crystal_object)}_region_unrefined"
            ),
        )

        print(
            f"Region unrefined RMSD: "
            f"{region_unrefined['unrefined_rmsd']:.3f} A\n"
            f"Matched C-alpha atoms: "
            f"{region_unrefined['unrefined_atom_count']}"
        )

        rows.append(
            _result_row(
                crystal_object=crystal_object,
                af_object=af_object,
                comparison="User-defined region - unrefined",
                crystal_selection=crystal_region_ca,
                af_selection=af_region_ca,
                cycles_requested=0,
                result=region_unrefined,
            )
        )

        region_refined = _align(
            mobile_selection=af_region_ca,
            target_selection=crystal_region_ca,
            cycles=refined_cycles,
            alignment_object=(
                f"{_safe_name(af_object)}_to_"
                f"{_safe_name(crystal_object)}_region_refined"
            ),
        )

        print(
            f"\nRegion refined alignment:\n"
            f"  RMSD before refinement: "
            f"{region_refined['unrefined_rmsd']:.3f} A\n"
            f"  RMSD after refinement: "
            f"{region_refined['refined_rmsd']:.3f} A\n"
            f"  Atoms after refinement: "
            f"{region_refined['refined_atom_count']}\n"
            f"  Atoms removed: "
            f"{region_refined['unrefined_atom_count'] - region_refined['refined_atom_count']}\n"
            f"  Atoms removed (%): "
            f"{((region_refined['unrefined_atom_count'] - region_refined['refined_atom_count']) / region_refined['unrefined_atom_count'] * 100) if region_refined['unrefined_atom_count'] else 0:.2f}%"
        )

        rows.append(
            _result_row(
                crystal_object=crystal_object,
                af_object=af_object,
                comparison="User-defined region - refined",
                crystal_selection=crystal_region_ca,
                af_selection=af_region_ca,
                cycles_requested=refined_cycles,
                result=region_refined,
            )
        )

    output_path = Path(output_csv).expanduser().resolve()
    _update_cumulative_csv(output_path, rows)

    print(f"\nCumulative CSV updated:\n  {output_path}")

    return rows


def compare_af_crystal_batch() -> None:
    """
    Repeatedly select crystal/AlphaFold pairs and append results to one CSV.
    """
    comparison_count = 0

    print(
        "\n=== Batch AlphaFold-vs-crystal RMSD mode ===\n"
        "Select each crystal/AlphaFold pair in turn.\n"
        "Results will be added to RMSD_Comparisons.csv.\n"
    )

    while True:
        objects = _loaded_objects()

        if len(objects) < 2:
            print(
                "At least two PyMOL objects must be loaded before "
                "starting a comparison."
            )
            _print_loaded_objects(objects)
            return

        crystal_object = _choose_object(
            "Select crystal structure",
            "Choose the crystal-structure object:",
            objects,
        )

        if crystal_object is None:
            print("Batch comparison finished.")
            break

        af_candidates = [
            object_name
            for object_name in objects
            if object_name != crystal_object
        ]

        af_object = _choose_object(
            "Select AlphaFold structure",
            "Choose the corresponding AlphaFold-model object:",
            af_candidates,
        )

        if af_object is None:
            print("Current comparison cancelled.")
            if not _ask_continue():
                break
            continue

        regions = _ask_region()

        if regions is None:
            print("Current comparison cancelled.")
            if not _ask_continue():
                break
            continue

        crystal_region, af_region = regions

        try:
            result = compare_af_crystal(
                crystal_object=crystal_object,
                af_object=af_object,
                crystal_region=crystal_region,
                af_region=af_region,
            )

            if result:
                comparison_count += 1

        except Exception as error:
            print(f"\nERROR: {error}")

        if not _ask_continue():
            break

    print(
        f"\nBatch session complete.\n"
        f"Successful receptor pairs compared: {comparison_count}\n"
        f"Cumulative results file:\n"
        f"  {Path('RMSD_Comparisons.csv').resolve()}"
    )


def compare_af_crystal_help() -> None:
    print(
        "\nBatch interactive mode:\n"
        "  compare_af_crystal_batch\n\n"
        "Direct single-comparison mode:\n"
        "  compare_af_crystal crystal_object, alphafold_object\n\n"
        "Optional region example:\n"
        '  compare_af_crystal crystal_object, alphafold_object, \\\n'
        '      crystal_region="crystal_object and resi 1-790", \\\n'
        '      af_region="alphafold_object and resi 1-790"\n'
    )
    _print_loaded_objects()


cmd.extend("compare_af_crystal", compare_af_crystal)
cmd.extend("compare_af_crystal_batch", compare_af_crystal_batch)
cmd.extend("compare_af_crystal_help", compare_af_crystal_help)

print(
    "\nLoaded batch AlphaFold-vs-crystal RMSD tool.\n\n"
    "Start the repeated comparison workflow with:\n"
    "  compare_af_crystal_batch\n\n"
    "Show help with:\n"
    "  compare_af_crystal_help\n"
)

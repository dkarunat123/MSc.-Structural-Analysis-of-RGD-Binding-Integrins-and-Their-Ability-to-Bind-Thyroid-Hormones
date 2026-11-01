from pymol import cmd

# Metal ion residue/element names to recognize.
# Add more names here later if needed.
METAL_SELECTION = "(resn MN or resn MN2 or elem MN or resn MG or resn MG2 or elem MG or resn CA or resn CA2 or elem CA)"


def _chain_resi_selection(chain, resi):
    """Builds a PyMOL selection for chain:resi, including blank-chain cases."""
    if chain == "":
        return f"(chain '' and resi {resi})"
    return f"(chain {chain} and resi {resi})"


def label_metal_ions():
    """
    Finds MN/MN2, MG/MG2, and CA/CA2 ions, prints them, shows them as spheres,
    and labels the actual atoms in the original object.

    This does NOT create separate ion objects like CA2_9101 or MG2_9001.
    """

    metal_ions = []

    cmd.iterate(
        METAL_SELECTION,
        "metal_ions.append((chain, resi, resn))",
        space={"metal_ions": metal_ions}
    )

    unique_metal_ions = []
    seen = set()

    for chain, resi, resn in metal_ions:
        key = (chain, resi, resn)
        if key not in seen:
            seen.add(key)
            unique_metal_ions.append((chain, resi, resn))

    if not unique_metal_ions:
        print("No MN/MN2, MG/MG2, or CA/CA2 ions found.")
        return

    print("\nMN/MN2, MG/MG2, and CA/CA2 ions found:")
    for chain, resi, resn in unique_metal_ions:
        print(f"{resn}  {chain}:{resi}")

    # Show and label the actual ion atoms directly, without creating new objects.
    cmd.show("spheres", METAL_SELECTION)
    cmd.color("red", METAL_SELECTION)
    cmd.set("sphere_scale", 0.45, METAL_SELECTION)
    cmd.label(METAL_SELECTION, 'resn + "_" + resi')

    # Make the labels much easier to see.
    # Some label settings are global/object-level in PyMOL, so setting them
    # globally avoids the pink Setting-Warning messages.
    cmd.set("label_size", 30)
    cmd.set("label_color", "white")
    cmd.set("label_outline_color", "black")
    cmd.set("label_shadow_mode", 1)

    # This one can be applied to the ion selection. It offsets labels away
    # from the ion spheres so they are not buried inside them.
    cmd.set("label_position", [3.0, 2.0, 2.5], METAL_SELECTION)

    # Make metal spheres a bit smaller so the labels do not get buried.
    

    print("\nLabeled the actual ion atoms directly with larger, offset labels. No separate ion objects were created.")


def clear_metal_ion_labels():
    """Removes labels from MN/MN2, MG/MG2, and CA/CA2 ions."""
    cmd.label(METAL_SELECTION, "")
    print("Removed metal ion labels.")


def find_metal_neighbors(*args, _self=None):
    """
    Examples:

    find_metal_neighbors B:9001
    find_metal_neighbors B:9001,B:9101,B:9102
    find_metal_neighbors B:9001 B:9101 B:9102
    find_metal_neighbors B:9001,B:9101,B:9102 5.0
    find_metal_neighbors B:9001 B:9101 B:9102 5.0

    This recognizes MN/MN2, MG/MG2, and CA/CA2 ions.
    It creates only temporary selections, not separate ion objects.
    """

    if len(args) == 0:
        print("Usage: find_metal_neighbors chain:resi,... [distance]")
        return

    try:
        distance = float(args[-1])
        metal_args = args[:-1]
    except ValueError:
        distance = 4.0
        metal_args = args

    metal_parts = []
    for arg in metal_args:
        metal_parts.extend([x.strip() for x in arg.split(",") if x.strip()])

    selection_parts = []
    input_metal_resis = []

    for metal in metal_parts:
        chain, resi = metal.split(":")
        input_metal_resis.append(resi)
        selection_parts.append(
            f"({_chain_resi_selection(chain, resi)} and {METAL_SELECTION})"
        )

    metal_selection = "(" + " or ".join(selection_parts) + ")"

    cmd.select("target_metal", metal_selection)

    cmd.select(
        "target_metal_neighbors",
        f"byres polymer.protein within {distance} of target_metal"
    )

    cmd.show("spheres", "target_metal")
    cmd.color("red", "target_metal")
    cmd.set("sphere_scale", 0.45, "target_metal")
    cmd.label("target_metal", 'resn + "_" + resi')

    # Make selected-ion labels visible instead of buried inside the spheres.
    # Global label settings avoid PyMOL object/state-level warning messages.
    cmd.set("label_size", 30)
    cmd.set("label_color", "white")
    cmd.set("label_outline_color", "black")
    cmd.set("label_shadow_mode", 1)
    cmd.set("label_position", [3.0, 2.0, 2.5], "target_metal")
    

    cmd.show("sticks", "target_metal_neighbors")

    print("\nSelected metal ions:")
    print(metal_selection)

    print(f"\nProtein residues within {distance} Å:")
    cmd.iterate(
        "target_metal_neighbors and name CA",
        "print(resn, resi, chain)"
    )

    stored_resis = []
    cmd.iterate(
        "target_metal_neighbors and name CA",
        "stored_resis.append(resi)",
        space={"stored_resis": stored_resis}
    )

    final_resis = stored_resis + input_metal_resis

    print("\nResidue numbers only, including selected metal ions at the end:")
    print(",".join(final_resis))


# New command names
cmd.extend("find_metal_neighbors", find_metal_neighbors)
cmd.extend("label_metal_ions", label_metal_ions)
cmd.extend("clear_metal_ion_labels", clear_metal_ion_labels)

# Old command names still work too.
cmd.extend("find_mn_neighbors", find_metal_neighbors)
cmd.extend("list_mn_ions", label_metal_ions)
cmd.extend("list_metal_ions", label_metal_ions)

# Automatically label ions when the script is run.
label_metal_ions()

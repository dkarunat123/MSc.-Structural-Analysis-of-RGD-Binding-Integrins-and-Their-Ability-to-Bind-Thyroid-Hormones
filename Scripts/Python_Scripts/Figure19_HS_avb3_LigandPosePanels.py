from pathlib import Path
from pymol import cmd

# ============================================================
# FIGURE 19
# Structural basis of ligand-specific affinity differences
# in the human AlphaFold alphaVbeta3 model
#
# IMPORTANT:
# This script assumes these objects are ALREADY loaded in PyMOL:
#
#   avb3_T4
#   avb3_TETRAC
#   avb3_sT4
#   avb3_T4AM
#
# It does NOT load PDB files.
# It does NOT change directories.
# It saves one PNG per panel in a dedicated results folder.
#
# Temporary copies are created so that the original structures
# are not modified.
#
# Panels:
#   A = T4
#   B = TETRAC
#   C = sT4
#   D = T4AM
#
# Three MIDAS-associated Mn2+ ions are displayed:
#   A9004
#   A9006
#   A9007
#
# A9007 is the ligand-proximal Mn2+ identified in the
# structural-contact analysis.
# ============================================================

# CONFIGURATION: PyMOL must already contain the four named receptor-ligand
# complexes. Run this file from PyMOL with: run /path/to/this_script.py
# Files are saved next to this script, regardless of the PyMOL working folder.
SAVE_PANEL_PNGS = True
PANEL_WIDTH_PX = 1800
PANEL_HEIGHT_PX = 1800
PANEL_DPI = 300

OUTPUT_DIR = (
    Path(__file__).resolve().parent
    / "Figure19_HumanAlphaVBeta3_LigandPosePanels_RESULTS"
)
if SAVE_PANEL_PNGS:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# 1. PANEL INFORMATION
# ------------------------------------------------------------

PANELS = {

    "A_T4": {

        "source":
            "avb3_T4",

        "object":
            "fig19_T4",

        "ligand_colour":
            "fig19_orange",

        # T4 carboxyl group
        "functional_group_atoms":
            "C9+O2+O3",

        # Carboxyl oxygen closest to Mn2+ A9007
        "metal_atom":
            "O3",

        # Relevant receptor residues
        "pocket_residues":
            "1140+1142+1234+1239",

        "draw_mn_distance":
            True,
    },


    "B_TETRAC": {

        "source":
            "avb3_TETRAC",

        "object":
            "fig19_TETRAC",

        "ligand_colour":
            "fig19_lightgreen",

        # TETRAC carboxyl group
        "functional_group_atoms":
            "C14+O3+O4",

        # Carboxyl oxygen closest to Mn2+ A9007
        "metal_atom":
            "O3",

        "pocket_residues":
            "1140+1142+1234+1239",

        "draw_mn_distance":
            True,
    },


    "C_sT4": {

        "source":
            "avb3_sT4",

        "object":
            "fig19_sT4",

        "ligand_colour":
            "fig19_blue",

        # sT4 sulfate group
        "functional_group_atoms":
            "S1+O3+O5+O6+O7",

        # Sulfate oxygen closest to Mn2+ A9007
        "metal_atom":
            "O7",

        "pocket_residues":
            "1140+1142+1234+1239",

        "draw_mn_distance":
            True,
    },


    "D_T4AM": {

        "source":
            "avb3_T4AM",

        "object":
            "fig19_T4AM",

        "ligand_colour":
            "fig19_lavender",

        # T4AM lacks a carboxyl group.
        # The amine is highlighted for orientation.
        "functional_group_atoms":
            "N1",

        "metal_atom":
            None,

        # T4AM-specific nearby receptor residues
        "pocket_residues":
            "1144+1145+1148",

        "draw_mn_distance":
            False,
    },
}


# ------------------------------------------------------------
# 2. CHECK THAT REQUIRED OBJECTS ARE ALREADY LOADED
# ------------------------------------------------------------

loaded_objects = set(
    cmd.get_object_list()
)

missing_objects = []

for info in PANELS.values():

    if info["source"] not in loaded_objects:

        missing_objects.append(
            info["source"]
        )


if missing_objects:

    raise RuntimeError(
        "These required PyMOL objects are not loaded: "
        + ", ".join(missing_objects)
    )


print(
    "All required structures are already loaded."
)

print(
    "Using: "
    + ", ".join(
        info["source"]
        for info in PANELS.values()
    )
)


# ------------------------------------------------------------
# 3. DEFINE CUSTOM COLOURS
# ------------------------------------------------------------

# Ligand carbon colours used in Figure 19
cmd.set_color(
    "fig19_orange",
    [0.95, 0.60, 0.10]
)

cmd.set_color(
    "fig19_lightgreen",
    [0.65, 0.90, 0.20]
)

cmd.set_color(
    "fig19_blue",
    [0.10, 0.50, 0.95]
)

cmd.set_color(
    "fig19_lavender",
    [0.75, 0.75, 0.95]
)


# Integrin colours
cmd.set_color(
    "fig19_receptor_grey",
    [0.78, 0.78, 0.78]
)

cmd.set_color(
    "fig19_pocket_grey",
    [0.48, 0.48, 0.48]
)


# Mn2+ colour
cmd.set_color(
    "fig19_metal_grey",
    [0.28, 0.28, 0.28]
)


# ------------------------------------------------------------
# 4. REMOVE OLD FIGURE-19 TEMPORARY OBJECTS
#
# Original avb3_* objects are NOT deleted.
# ------------------------------------------------------------

for name in [

    "fig19_T4",
    "fig19_TETRAC",
    "fig19_sT4",
    "fig19_T4AM",

]:

    cmd.delete(
        name
    )


for name in [

    "fig19_T4_Mn_distance",
    "fig19_TETRAC_Mn_distance",
    "fig19_sT4_Mn_distance",

]:

    cmd.delete(
        name
    )


# ------------------------------------------------------------
# 5. CREATE TEMPORARY FIGURE COPIES
# ------------------------------------------------------------

for info in PANELS.values():

    cmd.create(
        info["object"],
        info["source"]
    )


print(
    "Temporary Figure 19 copies created."
)


# ------------------------------------------------------------
# 6. DISABLE ORIGINAL LOADED STRUCTURES
#
# This prevents the original coloured receptors from appearing
# underneath the grey Figure-19 copies.
#
# They remain loaded and completely unchanged.
# ------------------------------------------------------------

for info in PANELS.values():

    cmd.disable(
        info["source"]
    )


print(
    "Original avb3_* objects hidden."
)


# ------------------------------------------------------------
# 7. ALIGN THE FIGURE COPIES
#
# Align receptor C-alpha atoms so that all four panels use
# directly comparable orientations.
# ------------------------------------------------------------

REFERENCE = "fig19_T4"


for info in PANELS.values():

    obj = info["object"]

    if obj == REFERENCE:
        continue

    cmd.align(

        f"{obj} and polymer.protein and name CA",

        f"{REFERENCE} and polymer.protein and name CA",

        cycles=0
    )


# ------------------------------------------------------------
# 8. GLOBAL DISPLAY SETTINGS
# ------------------------------------------------------------

cmd.bg_color(
    "white"
)

cmd.set(
    "ray_opaque_background",
    1
)

cmd.set(
    "antialias",
    2
)

cmd.set(
    "orthoscopic",
    1
)

cmd.set(
    "depth_cue",
    0
)

cmd.set(
    "specular",
    0.25
)

cmd.set(
    "shininess",
    20
)

cmd.set(
    "cartoon_fancy_helices",
    1
)

cmd.set(
    "cartoon_transparency",
    0.55
)

cmd.set(
    "stick_radius",
    0.18
)

cmd.set(
    "sphere_scale",
    0.38
)

cmd.set(
    "dash_width",
    2.5
)

cmd.set(
    "dash_length",
    0.25
)

cmd.set(
    "dash_gap",
    0.20
)

cmd.set(
    "dash_round_ends",
    1
)

cmd.set(
    "label_size",
    18
)

cmd.set(
    "label_color",
    "black"
)

cmd.set(
    "label_digits",
    2
)


# ------------------------------------------------------------
# 9. LIGAND COLOURING FUNCTION
#
# Carbon scaffold = ligand-specific colour
# O = red
# N = blue
# S = yellow
# I = purple
# H = white
# ------------------------------------------------------------

def colour_ligand(
    selection,
    carbon_colour
):

    cmd.color(
        carbon_colour,
        f"({selection}) and elem C"
    )

    cmd.color(
        "red",
        f"({selection}) and elem O"
    )

    cmd.color(
        "blue",
        f"({selection}) and elem N"
    )

    cmd.color(
        "yellow",
        f"({selection}) and elem S"
    )

    cmd.color(
        "purple",
        f"({selection}) and elem I"
    )

    cmd.color(
        "white",
        f"({selection}) and elem H"
    )


# ------------------------------------------------------------
# 10. STYLE EACH FIGURE OBJECT
# ------------------------------------------------------------

for panel_name, info in PANELS.items():

    obj = info["object"]


    # --------------------------------------------------------
    # Selections
    # --------------------------------------------------------

    receptor = (
        f"{obj} and polymer.protein"
    )


    ligand = (
        f"{obj} and chain B and resn UNL"
    )


    # All three MIDAS-associated metal ions
    midas_metals = (
        f"{obj} and resn MN2 "
        f"and resi 9004+9006+9007"
    )


    # A9007 is the ligand-proximal metal
    interacting_metal = (
        f"{obj} and resn MN2 "
        f"and resi 9007"
    )


    pocket = (
        f"{obj} and chain A and "
        f"resi {info['pocket_residues']}"
    )


    functional_group = (
        f"{ligand} and "
        f"name {info['functional_group_atoms']}"
    )


    # --------------------------------------------------------
    # Hide everything first
    # --------------------------------------------------------

    cmd.hide(
        "everything",
        obj
    )


    # --------------------------------------------------------
    # Integrin receptor
    #
    # Entire receptor is grey.
    # --------------------------------------------------------

    cmd.show(
        "cartoon",
        receptor
    )

    cmd.color(
        "fig19_receptor_grey",
        receptor
    )


    # --------------------------------------------------------
    # Ligand
    # --------------------------------------------------------

    cmd.show(
        "sticks",
        ligand
    )

    colour_ligand(
        ligand,
        info["ligand_colour"]
    )

    cmd.set(
        "stick_radius",
        0.22,
        ligand
    )


    # Slightly emphasize the functional group
    cmd.set(
        "stick_radius",
        0.30,
        functional_group
    )


    # --------------------------------------------------------
    # Relevant receptor residues
    #
    # These remain GREY as part of the integrin.
    # They are simply slightly darker than the cartoon.
    # --------------------------------------------------------

    cmd.show(
        "sticks",
        pocket
    )

    cmd.color(
        "fig19_pocket_grey",
        pocket
    )

    cmd.set(
        "stick_radius",
        0.20,
        pocket
    )


    # Label one CA atom per selected residue
    cmd.label(

        f"({pocket}) and name CA",

        '"%s%s" % (resn,resi)'
    )


    # --------------------------------------------------------
    # Three MIDAS-associated Mn2+ ions
    # --------------------------------------------------------

    cmd.show(
        "spheres",
        midas_metals
    )

    cmd.color(
        "fig19_metal_grey",
        midas_metals
    )

    cmd.set(
        "sphere_scale",
        0.45,
        midas_metals
    )


    # Label ONLY the ligand-proximal metal to reduce clutter
    cmd.label(
        interacting_metal,
        '"Mn2+"'
    )


    # --------------------------------------------------------
    # Ligand -> Mn2+ distance
    #
    # Only A9007 is used for the measurement line.
    #
    # T4AM receives no line because no ligand atom was
    # within 4 A of an Mn2+ ion in the contact analysis.
    # --------------------------------------------------------

    if info["draw_mn_distance"]:

        ligand_atom = (

            f"{ligand} and "
            f"name {info['metal_atom']}"
        )


        distance_name = (
            f"{obj}_Mn_distance"
        )


        cmd.distance(

            distance_name,

            ligand_atom,

            interacting_metal
        )


        cmd.set(
            "dash_color",
            "black",
            distance_name
        )


        cmd.set(
            "label_color",
            "black",
            distance_name
        )


# ------------------------------------------------------------
# 11. CREATE A COMMON CAMERA VIEW
#
# The view includes:
#   ligand
#   relevant pocket residues
#   all three MIDAS-associated Mn2+ ions
# ------------------------------------------------------------

reference_region = (

    "("

    "fig19_T4 and chain A and "
    "resi 1140+1141+1142+1144+1145+1148+1234+1239"

    ") or ("

    "fig19_T4 and chain B and resn UNL"

    ") or ("

    "fig19_T4 and resn MN2 "
    "and resi 9004+9006+9007"

    ")"
)


# Disable all Figure-19 copies temporarily
for info in PANELS.values():

    cmd.disable(
        info["object"]
    )


# Show reference panel
cmd.enable(
    "fig19_T4"
)


cmd.orient(
    reference_region
)


cmd.zoom(
    reference_region,
    5.0
)


COMMON_VIEW = (
    cmd.get_view()
)


# ------------------------------------------------------------
# 12. CREATE FOUR PYMOL SCENES
# ------------------------------------------------------------

for panel_name, info in PANELS.items():


    # --------------------------------------------------------
    # Keep original objects hidden
    # --------------------------------------------------------

    for other_info in PANELS.values():

        cmd.disable(
            other_info["source"]
        )


    # --------------------------------------------------------
    # Hide all Figure-19 copies
    # --------------------------------------------------------

    for other_info in PANELS.values():

        cmd.disable(
            other_info["object"]
        )


    # --------------------------------------------------------
    # Hide all Figure-19 distance objects
    # --------------------------------------------------------

    for distance_obj in [

        "fig19_T4_Mn_distance",
        "fig19_TETRAC_Mn_distance",
        "fig19_sT4_Mn_distance",

    ]:

        if distance_obj in cmd.get_names(
            "objects"
        ):

            cmd.disable(
                distance_obj
            )


    # --------------------------------------------------------
    # Enable the current receptor-ligand complex
    # --------------------------------------------------------

    cmd.enable(
        info["object"]
    )


    # --------------------------------------------------------
    # Enable its Mn-distance line when applicable
    # --------------------------------------------------------

    if info["draw_mn_distance"]:

        distance_name = (
            f"{info['object']}_Mn_distance"
        )

        cmd.enable(
            distance_name
        )


    # --------------------------------------------------------
    # Apply identical camera orientation
    # --------------------------------------------------------

    cmd.set_view(
        COMMON_VIEW
    )


    # --------------------------------------------------------
    # Store scene
    # --------------------------------------------------------

    cmd.scene(
        panel_name,
        "store"
    )

    # Export the exact scene view as a high-resolution, ray-traced PNG.
    # This does not modify the original avb3_* PDB objects.
    if SAVE_PANEL_PNGS:
        output_file = OUTPUT_DIR / f"Figure19_{panel_name}.png"
        cmd.png(
            str(output_file),
            width=PANEL_WIDTH_PX,
            height=PANEL_HEIGHT_PX,
            dpi=PANEL_DPI,
            ray=1,
        )
        print(f"Saved: {output_file}")


    print(
        f"Created scene: {panel_name}"
    )


# ------------------------------------------------------------
# 13. RETURN TO PANEL A
# ------------------------------------------------------------

cmd.scene(
    "A_T4",
    "recall"
)


# ------------------------------------------------------------
# 14. COMPLETION MESSAGE
# ------------------------------------------------------------

print()

print(
    "============================================================"
)

print(
    "FIGURE 19 SETUP COMPLETE"
)

print(
    "============================================================"
)

print()

print(
    "Original avb3_* objects remain loaded but are hidden."
)

print()

print(
    "Integrin:"
)

print(
    "  grey cartoon"
)

print(
    "  nearby residues shown as darker grey sticks"
)

print()

print(
    "Ligands:"
)

print(
    "  T4     = orange"
)

print(
    "  TETRAC = light green"
)

print(
    "  sT4    = blue"
)

print(
    "  T4AM   = lavender"
)

print()

print(
    "MIDAS-associated Mn2+ ions shown:"
)

print(
    "  A9004"
)

print(
    "  A9006"
)

print(
    "  A9007"
)

print()

print(
    "Distance lines:"
)

print(
    "  T4     -> A9007"
)

print(
    "  TETRAC -> A9007"
)

print(
    "  sT4    -> A9007"
)

print(
    "  T4AM   -> no distance line"
)

print()

print(
    "Scenes created:"
)

print(
    "  A_T4"
)

print(
    "  B_TETRAC"
)

print(
    "  C_sT4"
)

print(
    "  D_T4AM"
)

print()

print(
    "Use the PyMOL Scenes menu to switch between panels."
)

# OUTPUT SUMMARY
# With SAVE_PANEL_PNGS = True, four images are written to:
# Figure19_HumanAlphaVBeta3_LigandPosePanels_RESULTS/
#   Figure19_A_T4.png
#   Figure19_B_TETRAC.png
#   Figure19_C_sT4.png
#   Figure19_D_T4AM.png
# The script also stores the four original scene names in PyMOL.
# Figures are single-panel images, not an assembled multi-panel figure.
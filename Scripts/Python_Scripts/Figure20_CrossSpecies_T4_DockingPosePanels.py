"""Figure 20: cross-species visualization of docked T4 binding poses.

Run inside PyMOL, NOT with a plain Python interpreter:
    run /absolute/path/Figure20_CrossSpecies_T4_DockingPosePanels.py

Required PyMOL objects (loaded beforehand):
    HS_avb3_T4      Homo sapiens alphaVbeta3 with docked T4
    SP_aP1_bG_T4    Strongylocentrotus purpuratus alphaP1betaG with T4
    NV_a5_bn2_T4   Nematostella vectensis alpha5betan2 with T4

Expected molecular structure: receptor protein plus ligand UNL, chain B;
Mn2+ ions have element Mn or residue name MN/MN2.

Operations:
  * Create independent copies, keeping source coordinates intact.
  * Align the two non-human receptor copies to the human receptor using C-alpha.
  * Identify each protein residue having a HEAVY ATOM within 4.0 A of a
    HEAVY ATOM of ligand UNL (not of other ions or heteroatoms).
  * Show receptor cartoons and contacting residues in species-specific colours;
    show T4 in orange and Mn2+ ions as dark spheres.
  * Produce three independent, uniformly oriented PNG panels and PyMOL scenes.

Geometric 4 A contacts are proximity measures, not validated interactions.
All 8 input Mn2+ ions may exist in each model; only the visible pocket region
is interpreted in the thesis figure. Metal placement is not experimental proof
of conservation.
"""

from pathlib import Path
from pymol import cmd

# =============================================================================
# 0. USER-ADJUSTABLE SETTINGS
# =============================================================================

# Set False to use the old script's individually reoriented cameras.
# True is recommended for meaningful visual comparison across aligned models.
USE_COMMON_CAMERA_VIEW = True

# Set False if you only need PyMOL scenes and do not want PNG outputs.
SAVE_PANEL_PNGS = True
PANEL_WIDTH_PX = 1800
PANEL_HEIGHT_PX = 1800
PANEL_DPI = 300

# Label the selected contact residues. The earlier script displayed their
# sticks but did not add labels, even though the thesis figure caption mentions
# labelled residues.
LABEL_CONTACT_RESIDUES = True
CONTACT_CUTOFF_A = 4.0

OUTPUT_DIR = (
    Path(__file__).resolve().parent
    / "Figure20_CrossSpecies_T4_DockingPosePanels_RESULTS"
)
if SAVE_PANEL_PNGS:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Ordered by their panel positions in Figure 20.
PANELS = [
    {"source": "HS_avb3_T4", "copy": "fig20_HS", "scene": "A_HS",
     "output": "Figure20_A_HS_T4.png", "colour": "fig20_hs_blue"},
    {"source": "SP_aP1_bG_T4", "copy": "fig20_SP", "scene": "B_SP",
     "output": "Figure20_B_SP_T4.png", "colour": "fig20_sp_purple"},
    {"source": "NV_a5_bn2_T4", "copy": "fig20_NV", "scene": "C_NV",
     "output": "Figure20_C_NV_T4.png", "colour": "fig20_nv_green"},
]

# =============================================================================
# 1. VERIFY REQUIRED PYMOL INPUTS BEFORE CHANGING THE SESSION
# =============================================================================

loaded = set(cmd.get_object_list())
missing = [p["source"] for p in PANELS if p["source"] not in loaded]
if missing:
    raise RuntimeError(
        "Load the following PDB objects in PyMOL before running Figure 20: "
        + ", ".join(missing)
    )

for panel in PANELS:
    source = panel["source"]
    ligand = f"({source} and chain B and resn UNL and not elem H)"
    protein = f"({source} and polymer.protein and not elem H)"
    if cmd.count_atoms(ligand) == 0:
        raise RuntimeError(
            f"{source}: missing ligand heavy atoms in chain B, residue UNL."
        )
    if cmd.count_atoms(protein) == 0:
        raise RuntimeError(f"{source}: missing protein heavy atoms.")

# =============================================================================
# 2. RECREATE TEMPORARY FIGURE OBJECTS (THE ORIGINAL INPUTS ARE NOT ALTERED)
# =============================================================================

for panel in PANELS:
    cmd.delete(panel["copy"])  # Safe if object does not exist.
    cmd.create(panel["copy"], panel["source"])

# The original script used PyMOL's default align() outlier-refinement behavior.
# State cycles=5 explicitly to reproduce it, rather than introducing a new
# superposition method during code cleanup.
for panel in PANELS[1:]:
    alignment = cmd.align(
        f"{panel['copy']} and polymer.protein and name CA",
        "fig20_HS and polymer.protein and name CA",
        cycles=5,
    )
    print(
        f"Aligned {panel['source']} to HS_avb3_T4; "
        f"refined C-alpha RMSD = {alignment[0]:.3f} A"
    )

# =============================================================================
# 3. COLOUR PALETTE AND VIEW SETTINGS
# =============================================================================

# Match the original species and T4 colours.
cmd.set_color("fig20_hs_blue", [0.27, 0.49, 0.74])
cmd.set_color("fig20_sp_purple", [0.67, 0.24, 0.66])
cmd.set_color("fig20_nv_green", [0.53, 0.68, 0.10])
cmd.set_color("fig20_t4_orange", [1.00, 0.55, 0.00])
cmd.set_color("fig20_mn_dark", [0.20, 0.20, 0.20])

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("orthoscopic", 1)
cmd.set("antialias", 2)
cmd.set("depth_cue", 0)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_flat_sheets", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("sphere_scale", 0.32)
cmd.set("stick_radius", 0.18)
cmd.set("label_font_id", 7)
cmd.set("label_color", "black")
cmd.set("label_size", 15)

# =============================================================================
# 4. STYLE EACH RECEPTOR / LIGAND COMPLEX
# =============================================================================

def style_panel(obj_name, receptor_colour):
    """Show ligand UNL and protein residues within a 4 A HEAVY-ATOM cutoff.

    This explicitly excludes ions and other non-protein heteroatoms from
    the ligand selection, unlike the earlier broad 'not polymer.protein'
    expression. It matches the definition used for the CSV contact analysis.
    """
    protein = f"({obj_name} and polymer.protein)"
    protein_heavy = f"({obj_name} and polymer.protein and not elem H)"
    ligand = f"({obj_name} and chain B and resn UNL)"
    ligand_heavy = f"({obj_name} and chain B and resn UNL and not elem H)"
    metals = (
        f"({obj_name} and "
        "(elem Mn or resn MN or resn MN2 or resn MN+))"
    )
    nearby = (
        f"byres (({protein_heavy}) within {CONTACT_CUTOFF_A:.1f} "
        f"of ({ligand_heavy}))"
    )

    cmd.hide("everything", obj_name)

    cmd.show("cartoon", protein)
    cmd.color(receptor_colour, protein)
    cmd.set("cartoon_transparency", 0.35, protein)

    cmd.show("sticks", ligand)
    cmd.color("fig20_t4_orange", ligand)
    cmd.set("stick_radius", 0.20, ligand)

    cmd.show("spheres", metals)
    cmd.color("fig20_mn_dark", metals)

    cmd.show("sticks", nearby)
    cmd.color(receptor_colour, nearby)

    if LABEL_CONTACT_RESIDUES:
        # One label per residue using the C-alpha atom (not every atom).
        cmd.label(f"({nearby}) and name CA", '"%s%s" % (resn,resi)')

    residue_count = cmd.count_atoms(f"({nearby}) and name CA")
    print(f"{obj_name}: {residue_count} protein residues within 4.0 A of T4.")


for panel in PANELS:
    style_panel(panel["copy"], panel["colour"])

# =============================================================================
# 5. CREATE COMPARABLE CAMERA VIEWS, SCENES AND EXPORTED PANEL IMAGES
# =============================================================================

# First reproduce the original human-panel viewing procedure.
cmd.disable("all")
cmd.enable("fig20_HS")
cmd.orient("fig20_HS and not solvent")
cmd.zoom(
    "fig20_HS and (not polymer.protein or "
    "(polymer.protein within 8 of (not polymer.protein and not solvent)))",
    2.5,
)
cmd.turn("y", 10)
cmd.turn("x", -10)
COMMON_VIEW = cmd.get_view()

for index, panel in enumerate(PANELS):
    # Showing only one copied receptor ensures that no other open PyMOL
    # structures contaminate the rendered image.
    cmd.disable("all")
    cmd.enable(panel["copy"])

    if index == 0 or USE_COMMON_CAMERA_VIEW:
        cmd.set_view(COMMON_VIEW)
    else:
        # Original script: individual reorientation of SP/NV after alignment.
        cmd.orient(f"{panel['copy']} and not solvent")
        cmd.zoom(
            f"{panel['copy']} and (not polymer.protein or "
            "(polymer.protein within 8 of (not polymer.protein and not solvent)))",
            2.5,
        )
        cmd.turn("y", 10)
        cmd.turn("x", -10)

    cmd.scene(panel["scene"], "store")

    if SAVE_PANEL_PNGS:
        output_path = OUTPUT_DIR / panel["output"]
        cmd.png(
            str(output_path),
            width=PANEL_WIDTH_PX,
            height=PANEL_HEIGHT_PX,
            dpi=PANEL_DPI,
            ray=1,
        )
        print(f"Saved {output_path}")

# Leave the three aligned copies available for visual inspection.
# Scene switching recreates the separately formatted views.
cmd.disable("all")
for panel in PANELS:
    cmd.enable(panel["copy"])
cmd.orient("fig20_HS and polymer.protein")
cmd.zoom("all", 8)

print("\nFigure 20 preparation complete.")
print("Temporary PyMOL objects:", ", ".join(p["copy"] for p in PANELS))
print("PyMOL scenes:", ", ".join(p["scene"] for p in PANELS))
print("Output folder:", OUTPUT_DIR if SAVE_PANEL_PNGS else "PNG export disabled")
print(
    "Interpretation: contact residues are selected with a 4.0 A heavy-atom "
    "cutoff; these contacts do not establish binding mechanisms."
)

# =============================================================================
# END OF SCRIPT / OUTPUT SUMMARY
# =============================================================================
# Figure20_CrossSpecies_T4_DockingPosePanels_RESULTS/
#   Figure20_A_HS_T4.png
#   Figure20_B_SP_T4.png
#   Figure20_C_NV_T4.png
# PyMOL scenes: A_HS, B_SP, C_NV.
# =============================================================================
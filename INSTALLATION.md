# How to add the upgrade kit to your existing GitHub repository

The ZIP is **not a replacement for your full existing thesis repository**. It is a set of files to merge into the root of your current repository, retaining all of your existing data, docking results, structures and figures.

## Step 1 — Merge the files

1. Make a backup of your existing `README.md`, or use Git history to retain the previous version.
2. Extract the ZIP **into your repository root** (the folder that already contains `Data/`, `Docking/`, `Figures/` and `Scripts/`). Merge directories instead of replacing them wholesale. The kit includes numerical `Data/` files, four R scripts, the Figure 15 correlation output, and the `docs/` static dashboard.
3. Confirm `Data/OMA_Subunit_Relationships.xlsx` still exists in the root-level `Data/`. The uploaded files did not include this workbook, but it was visible in your original repository screenshots.
4. Decide whether to keep the replaced main `README.md` as-is or customize its introduction and layout to your preferred style.
5. Check `Scripts/R_Scripts/05_Species_BindingAffinity_Analysis.R` has a **single `.R` extension**; the uploaded copy had `.R..R`.

## Step 2 — Verify the local dashboard

Open `docs/index.html` in a browser. The initial page has all 432 values embedded as static JavaScript data. It does not require a development server or CDN. The dashboard allows you to click Figures 14, 15, and 16, change species or ordering, inspect scores and export SVG/CSV.

The uploaded source workbooks and Figure 15 correlations were audited numerically. Saved AU bootstrap data were **not** provided; the dashboard deliberately omits those labels. See `docs/README.md`.

## Step 3 — Commit your changes

From your repository root, check and commit the files you intend to publish:

```bash
git status
git add README.md CITATION.cff .gitignore INSTALLATION.md REPOSITORY_AUDIT.md docs/ Data/ Scripts/R_Scripts/ Results/
git commit -m "Add interactive thesis heatmaps and validated result data"
git push origin main
```

If you prefer to use GitHub's web interface, upload the files to their **matching repository-relative paths**; make sure the `docs/` directory and its subfolders remain intact. Avoid GitHub's "Upload files" action placing all nested files in the root.

## Step 4 — Enable GitHub Pages

Open **Settings → Pages** in the GitHub repository. Under **Build and deployment**, choose **Deploy from a branch**, branch `main`, folder `/docs`; save. Once GitHub reports a published site, open:

`https://dkarunat123.github.io/MSc.-Structural-Analysis-of-RGD-Binding-Integrins-and-Their-Ability-to-Bind-Thyroid-Hormones/`

The main README already links to this intended URL. If the repository was renamed, update both the README and site navigation.

## Step 5 — Final scientific/reproducibility checks

- Keep your original input workbooks in `Data/` and the complete selected structure/docking outputs in their designated folders.
- For an exact recreation of Figure 14's red AU annotations, supply or publish the three provenance-matched `pvclust` `.rds` files and companion edge/cluster CSVs listed in `docs/README.md`.
- Ensure the final deposited thesis PDF does not include template prompts or broken table-of-contents bookmarks before tagging a public release.
- Do not add a blanket LICENSE until redistribution permissions for third-party scientific structures, model outputs and reused figures have been checked.

## Rebuilding the embedded browser data after workbook changes

The site is currently complete with the values supplied here. If you **change the source workbooks**, rerun the bundled `tools/extract_validate.py`, followed by `tools/make_datajs.py` and `tools/make_preview.py` in a Python environment containing `artifact_tool`, NumPy, SciPy and Matplotlib. Alternatively, regenerate the standardized outputs from the R workflow and adapt the site data serialization to those resulting CSVs.

Unlike the static site, rebuilding from `.xlsx` requires the extraction dependencies. The main scientific analysis does not depend on the browser dashboard.

## Thesis-matched structure colours (v3)

The dashboard and README now apply the Figure 20 model identity colours: human AlphaFold blue (`#396EA5`), crystal-reference grey (`#9EA5AF`), sea urchin AlphaFold magenta (`#A641A6`), and *N. vectensis* AlphaFold olive-green (`#6F990C`). The three AlphaFold species colours are approximated from the supplied screenshot; human crystal-reference grey follows your stated convention. Change the `--human`, `--crystal`, `--urchin` and `--anemone` CSS variables in `docs/index.html` and the corresponding `MODEL_COLORS`/`MODEL_LIGHT` constants in `docs/assets/app.js` if you have the original PyMOL RGB values.

These colours identify receptor model/source only; the affinity and correlation gradients remain exactly the same as in v2.

## Light pastel site theme (v4)

This version retains **the original saturated thesis structure colours**, unchanged from v3:

- Human AlphaFold blue: `#396EA5`
- Human crystal-reference grey: `#9EA5AF`
- Sea urchin AlphaFold magenta: `#A641A6`
- *N. vectensis* AlphaFold olive-green: `#6F990C`

Only the **dashboard background, cards, controls, text contrast and static previews** have been made lighter. The page uses a very pale blue/lilac wash (`#F7F8FC`), near-white surfaces and subtle borders. Heatmap cell values **still use the original scientific continuous gradients**, not the species identity colours. The update is applied in `docs/index.html`, `docs/assets/app.js`, the structure colour SVG and updated PNG previews.

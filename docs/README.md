# Interactive results atlas: data provenance and deployment

This directory is a self-contained static site implementing interactive companion views of thesis Figures **14, 15 and 16**. Open `index.html` locally or enable GitHub Pages from **Settings → Pages → Deploy from branch → main → /docs**. The public URL is expected to follow `https://dkarunat123.github.io/<repository-name>/` when Pages is configured.

## Visual design and scientific colour meaning

The model-identity colours are used throughout the thesis to distinguish the four model types: **human AlphaFold = blue `#396EA5`**, **human crystal structures = grey `#9EA5AF`**, **sea urchin AlphaFold = magenta `#A641A6`**, and **sea anemone AlphaFold = olive-green `#6F990C`**. The digital hex values are visual approximations of the thesis-wide colour convention, rather than verified original PyMOL RGB settings; the crystal grey is a neutral approximation. Structure colours appear on labels, strips, legends and accents. **Cell fill remains a continuous quantitative map**: Figures 14/16 represent `-log10(KD)` and Figure 15 represents ρ.

Grey only identifies crystal structure controls in the structure colour key. The Figure 14-16 plots themselves use the AlphaFold model affinity data, so the dashboard does not fabricate a crystal-structure comparison.

## Directory files

- `index.html`: accessible static layout and figure controls.
- `assets/app.js`: no-third-party-runtime JavaScript rendering, interactive SVG heatmaps, inspection, CSV and SVG export.
- `assets/data.js`: audited numerical data included in the site so it works offline.
- `assets/structure-colour-key.svg`: reusable legend for the thesis structural colour identities.
- `assets/median-affinity-preview.png`: Figure 16 static preview used in the main README.
- `data/receptor_affinities.csv`: complete **48 × 9** long-format table of 432 predicted values, with species, historical receptor ID, readable heterodimer name, alpha and beta subunit, source raw KD and transformed score.
- `data/species_medians.csv`: 3 × 9 unrounded species medians from transformed affinity scores.
- `data/figure15_*_rho.csv`: 9 × 9 Spearman ligand-correlation matrix for each species.
- `data/species_rank_correlations.csv`: three species-pair Spearman coefficients across nine ligand medians.
- `data/cluster_orders.json`: recomputed hierarchical cluster merges and leaf orders.
- `data/validation_report.json`: quantitative audit and limitations.

## Validation performed

The supplied `BindingScores_StatisticalAnalysis.xlsx` contains 48 unique receptor rows and nine numeric transformed affinity columns. All **432** `-log10(KD)` values agree with the corresponding raw KD values in `All_Species_BindingScores.xlsx` (maximum absolute discrepancy **8.88 × 10⁻¹⁶**). All **108** ligand-pair correlations agree with the supplied `Species_LigandPair_SpearmanCorrelations.csv` (maximum absolute discrepancy **5.00 × 10⁻¹⁶**). The independently recalculated cross-species Spearman coefficients are **0.9166667**, **0.9000000**, **0.9833333**.

No filtering, Z-scoring, imputation or second log transformation was applied. Species medians are calculated directly in transformed-score space. Note that the original workbook's `S.P_AFReceptors` tab has an inaccurate *N. vectensis* description in its note row; the species was identified using the worksheet's correct name, receptor IDs and contents, as in Script 02.

## Figure accuracy and limits

**Figure 14:** Receptor rows are clustered using Euclidean/complete linkage. Ligand columns use Euclidean/complete linkage as in the R `pvclust` workflow. In the dashboard, SciPy regenerated the tree from the numerical data. The original saved bootstrap `.rds` objects were **not uploaded**, so node-wise approximately unbiased (AU) support values are **not displayed**. The Figure 14 R script combines the `pvclust` ligand tree with its exact AU edges, which cannot be reconstructed simply from the score matrix without rerunning the bootstrap.

**Figure 15:** Pairwise Spearman correlations were recalculated across receptors separately in each species. Rows and columns use the same correlation-based hierarchical clustering (1 - ρ, complete linkage); the value scale is fixed from -1 to +1. The **correlations match** the supplied CSV. Left/right tree orientation may differ from `pheatmap` output.

**Figure 16:** Median -log10(KD) values are calculated separately for nine ligands across 8/24/16 models. Species and ligand axes use Euclidean/complete linkage and no standardization. The values and species-pair Spearman rankings reproduce the thesis. Interactive tree orientation can differ from R due to independent tie-handling.

These predictions do **not** demonstrate experimentally validated hormone binding or receptor activity. Modelled heterodimers are not independent biological samples.

## R-native reproduction

From the **GitHub repository root**, the original scripts' relative paths `Data/` and `Results/` resolve as follows:

```bash
Rscript Scripts/R_Scripts/02_Prepare_Integrin_BindingData.R
# Optional only when the saved, provenance-matched pvclust RDS files are absent:
Rscript Scripts/R_Scripts/04_LigandCluster_BootstrapAnalysis.R
Rscript Scripts/R_Scripts/05_Species_BindingAffinity_Analysis.R
Rscript Scripts/R_Scripts/06_CrossSpecies_LigandProfile_Comparison.R
```

`Script 02` additionally requires `Data/OMA_Subunit_Relationships.xlsx`, which was shown in the user's existing repository and is not bundled here. `Script 04` uses one million replicates at each of ten multiscale bootstrap sample sizes *for every species*. `Script 05` stops if the correct saved bootstrap files are unavailable. The `05_Species_BindingAffinity_Analysis.R..R` uploaded filename was normalized to `05_Species_BindingAffinity_Analysis.R` in the upgrade kit; the script body was left unchanged.

### To restore the original AU annotations

Copy the provenance-matched output files produced by `04_LigandCluster_BootstrapAnalysis.R`:

- `Step12A_Ligand_pvclust_H_sapiens_1000000.rds`, and similarly `S_purpuratus` and `N_vectensis`
- `Step12A_Ligand_pvclust_Edges_<species>_1000000.csv`
- `Step12A_Ligand_pvclust_Clusters_<species>_1000000.csv`

These are required to place the exact published AU labels on the correct dendrogram merges, as explained in Script 05. The current interactive site intentionally makes no claims about individual cluster AU percentages.

## Data dictionary

| Column | Meaning |
|---|---|
| `species` | Abbreviated scientific species name (`H. sapiens`, `S. purpuratus`, `N. vectensis`) |
| `receptor` | Original computational heterodimer identifier used in docking and scripts |
| `receptor_label` | Readable alpha/beta designation in the final thesis nomenclature |
| `alpha_subunit`, `beta_subunit` | Readable constituent subunit names |
| `ligand` | One of the nine submitted compounds |
| `neg_log10_kd` | Predicted `-log10(KD)`, with raw KD in M; higher = stronger predicted binding |
| `kd_m` | Original predicted dissociation constant in molar units, not an experimentally measured KD |

In particular, for *N. vectensis*, the historical codes `aH6 → α3`, `aA9 → α9`, `aA5 → α5`, `a71 → α8`, and `bW0 → β2`, `bU5 → βn2`, `bD5 → β1`, `b23 → βb` follow the mapping in `02_Prepare_Integrin_BindingData.R`.

## Updating the dashboard

After updating the original R dataset, run Script 02 and regenerate the CSV files. The fully static `assets/data.js` is a serialized copy of those outputs and must be regenerated when inputs change. For this delivery, the data extraction and audit utility is included in the kit under `tools/extract_validate.py` (requires `artifact_tool`, NumPy and SciPy in that environment). The public site itself does **not** require these Python dependencies.

No backend server, cookies, authentication, external CDN or usage telemetry are required by the dashboard. The GitHub Pages deployment is simply the contents of `docs/`.

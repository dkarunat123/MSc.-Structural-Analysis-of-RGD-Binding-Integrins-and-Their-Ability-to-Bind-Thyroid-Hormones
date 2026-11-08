# Repository completion audit · 8 November 2026

**Evidence:** Thesis PDF; repository screenshots supplied in the conversation; complete uploaded binding-score workbooks, Spearman CSV, and R Scripts 02/04/05/06. This is a **partial visibility audit** rather than a recursive Git checkout. An item not visible in screenshots is **unverified**, not necessarily absent.

## Completed in this upgrade

- [x] Converted 48 × 9 receptor-level transformed scores to 432 machine-readable records with canonical species and human-readable receptor names.
- [x] Matched every one of the 432 transformed values against the source raw KD workbook; maximum error `8.88e-16`.
- [x] Matched all 108 within-species Spearman ligand-pair values against the supplied CSV; maximum error `5.00e-16`.
- [x] Regenerated Figure 16 medians and three rank correlations, matching the thesis.
- [x] Added interactive, exportable Figure 14/15/16 heatmaps with recreated clustering; omitted unprovided AU support labels.
- [x] Added README preview, GitHub Pages-ready `docs/`, source CSVs, `CITATION.cff` and `.gitignore`.
- [x] Normalized double extension of uploaded `05_Species_BindingAffinity_Analysis.R..R` **in the kit** to `.R`.
- [x] Bundled supplied source score workbooks, Scripts 02/04/05/06, and correlation CSV with original scripts unchanged.

## High-priority remaining repository checks

| Priority | File or task | Why it matters | Visibility/status |
|---|---|---|---|
| 1 | `Data/OMA_Subunit_Relationships.xlsx` | Required by Script 02. | Visible in screenshots, **not uploaded here**; must be present when running the bundled R workflow. |
| 1 | Original saved `pvclust` output per species: `.rds`, `Edges_*.csv`, `Clusters_*.csv` | Required to reproduce **exact Figure 14 AU labels** without the costly bootstrap rerun. | **Not provided** in these uploads. |
| 1 | Complete `All_Species_AlphaFoldScores.xlsx` | Source model-selection/confidence values for all candidate AlphaFold structures. | Only `All_Species_AlphaFoldScores_ScoresOnly.xlsx` visible in screenshots. |
| 1 | Verify `Results/`, `Structures/`, `Supplementary/` on `main` | The prior README described all three, but original screenshots showed only four folders at root. | **Unverified** on current branch; audit actual folder tree rather than assume missing. |
| 2 | `Ligand_Validation_Results.csv`, `Ligand_Validation_Stereocenters.csv` | Thesis Methods §2.6 validation evidence. | Not visible in screenshots. |
| 2 | Main figures 9–20 as PDF/SVG/PNG, supplementary figures 1–4 and supplementary Table 1 | Accessible viewing and traceability without PyMOL `.pse` files. | Only `.pse` sessions and a ZIP visible in `Figures/`; check `Supplementary/` first. |
| 2 | Rename `HS_CSandAF_RMSDScores.csv.xlsx` if it is actually an `.xlsx` workbook | Align references in thesis and analysis scripts. | Double extension visible in screenshot. |
| 2 | `03_Human_CS_vs_AF_BindingValidation.R`, `07_OMA_BindingAffinity_Comparisons.R`, `01_Integrin_SequenceDistance_Analysis.R`, `Check_Rejected_SeaUrchin_Sequences.R` | Other exact named stages of the reported R workflow. | May be inside `Scripts/R_Scripts/`; not independently verified. |
| 2 | Detailed `Scripts/Python_Scripts/` inventory and `Docking/` contents | Trace structure preparation, ligand validation and docking workflows. | Folder listings in screenshots are partial. |
| 2 | Confirm full 48-receptor selected structures, ligand inputs and cluster output mappings | Reproduce structural figures and docking calculations. | Selected folder lists visible but not recursively audited. |
| 3 | `LICENSE`, Python and R dependency manifests, reproducibility tests | Usability and reuse; avoid attributing third-party files to the wrong licenses. | Not visible; license must be chosen deliberately. |
| 3 | Final thesis PDF | The uploaded 86-page PDF appears to retain front-matter template text and broken TOC references. | Correct source Word document and regenerate PDF before marking as final. |

## Workflow path issue worth fixing

The uploaded R scripts use relative paths such as `Data/BindingScores_StatisticalAnalysis.xlsx` and `Results/02_Prepare_Integrin_BindingData_RESULTS/`. Their header comments say to set the working directory to `R_Scripts`, apparently reflecting an earlier local folder layout. In the GitHub repository shown, `Data/` is **at the root**. Therefore, when scripts live in `Scripts/R_Scripts/`, run them **from the repository root**, e.g.:

```bash
Rscript Scripts/R_Scripts/02_Prepare_Integrin_BindingData.R
```

The scripts themselves use relative paths and do not invoke `setwd()` in the uploaded versions. A collaborator running them from inside `Scripts/R_Scripts/` will otherwise see missing data. Keep the documentation consistent with your actual directory tree.

## Recommended further enhancements

1. **Original AU-supported dendrograms** inside the dashboard after the `.rds` bootstrap outputs are supplied. These should be aligned with each ligand cluster using its saved `pvclust` edge, not guessed from branch order or thesis image.
2. **3D molecular viewer** (3Dmol.js or Mol*) to rotate Figure 19/20 αVβ3, αP1βG and α5βn2 complexes, after selecting and validating exact PDB files for presentation.
3. **Figure lookup pages** connecting each thesis figure and table to the source file, script, exact output and method in one index.
4. **Citable release** (`v1.0-thesis`) after completeness/licensing review, and optionally a DOI through a suitable public research repository.

## Scientific cautions for presentation

- Non-human heterodimers were **modeled combinations**, not confirmed physiological receptors.
- A Spearman correlation of predicted ligand rankings does not establish conserved physiological binding or signaling.
- Different Figure 14 species panels use independent colour scales by default; use the toggle in the interactive site to compare their magnitudes on a common scale.
- Figure 16 is based on medians in transformed `−log10(KD)` space. Figure 18 requires medians from **raw KD**, which can yield different numerical values for even-sized groups.
- The interactive site is a faithful *data companion* to the thesis figures, not a claim of pixel-identical recreation of R output or of biological validation.

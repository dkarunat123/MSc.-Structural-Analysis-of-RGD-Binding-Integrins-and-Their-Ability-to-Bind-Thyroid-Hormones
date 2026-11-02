# Comparative Structural Analysis of RGD-Binding Integrins and Their Ability to Bind Thyroid Hormones

**MSc Bioinformatics Thesis**  
**University of Guelph | Heyland Lab**  
**Researcher:** Dilma Karunathilake

## Project Overview

This repository contains research data, analysis scripts, molecular structures, docking outputs, figures, and supplementary materials associated with my MSc Bioinformatics thesis at the University of Guelph.

The study investigates the potential for thyroid hormone recognition among RGD-binding integrins and compares predicted binding properties across three metazoan species:

- *Homo sapiens* (human)
- *Strongylocentrotus purpuratus* (purple sea urchin)
- *Nematostella vectensis* (starlet sea anemone)

## Abstract

Integrin αVβ3 mediates non-genomic thyroid hormone signaling, but whether this capacity extends to other integrins or predates vertebrates remains unclear. This study investigated predicted thyroid hormone binding in eight human RGD-binding integrins and 40 candidate heterodimers from the sea urchin *Strongylocentrotus purpuratus* and sea anemone *Nematostella vectensis*. Integrin ectodomains were modeled using AlphaFold 3, docked with thyroid hormone derivatives using HADDOCK, and evaluated using PRODIGY-LIG alongside sequence and orthology analyses. All eight human receptors and multiple non-human models exhibited structurally plausible interactions, with human αVβ8 showing particularly favorable predicted affinities. Across species, ligand-affinity rankings were strongly correlated (Spearman ρ = 0.900–0.983), with TRIAC and TETRAC generally favored and T4AM consistently weakest. Conserved binding-pocket architecture and similar affinities among β1-related integrins suggest an evolutionarily ancient structural compatibility with thyroid hormones. These computational findings identify candidate receptors but do not establish physiological binding or signaling.

## Research Objectives

1. Investigate predicted thyroid hormone interactions with the eight human RGD-binding integrins.
2. Identify sea urchin integrin heterodimers with structural features compatible with thyroid hormone binding.
3. Compare ligand-binding architectures and predicted affinities among evolutionarily related integrins across the three species.

## Computational Approach

The analyses cover **48 integrin heterodimer models** (8 human, 24 sea urchin, and 16 sea anemone) and **nine ligands**: T3, T4, rT3, TRIAC, TETRAC, sT3, sT4, T4AM, and the RGD-mimetic peptide cilengitide. The workflow includes integrin sequence and orthology analyses, structural modelling, docking, affinity estimation, statistical comparisons, and visualization of binding poses.

### Key Software and Versions

| Tool | Version used | Purpose |
|---|---|---|
| AlphaFold | 3 | Prediction of integrin extracellular heterodimer structures |
| HADDOCK | 2.4 | Protein–ligand molecular docking |
| PRODIGY-LIG | 1.1.4 | Prediction of binding free energies |
| PyMOL | 3.2.0 | Structural comparison and visualization |
| UCSF Chimera | 1.19 | Structure and ligand preparation |
| Open Babel | 3.1.1 | Ligand file conversion |
| RDKit | 2026.09.1 | Validation of ligand structures and stereochemistry |
| MUSCLE | 5.2 | Multiple sequence alignment |
| R | 4.6.1 | Sequence, statistical, and visualization analyses |
| ape (R package) | 5.8-1 | Sequence distances and bootstrap analyses |
| pheatmap (R package) | 1.0.13 | Heatmap visualization |
| Phobius | 1.01 | Prediction of signal peptides and membrane topology |
| DeepTMHMM | 1.0 | Identification of transmembrane regions |
| SignalP | 6.0 | Additional signal-peptide assessment |

Additional resources include the **OMA** database for orthology relationships, **PubChem** for ligand structures, and the **Protein Data Bank (PDB)** for experimental integrin reference structures. The thesis does not specify versions for all additional databases or software packages.

## Repository Contents

| Folder | Contents |
|---|---|
| `Data/` | Sequence information, receptor metadata, AlphaFold confidence scores, predicted binding-affinity datasets, and OMA relationships |
| `Docking/` | Docking inputs, selected integrin–ligand complexes, and affinity-prediction outputs |
| `Figures/` | Figures and visualizations prepared for the thesis |
| `Results/` | Processed data, statistical analysis outputs, correlations, clustering results, and related summaries |
| `Scripts/` | Python and R scripts for structure preparation, sequence and evolutionary analyses, docking post-processing, statistics, and figures |
| `Structures/` | AlphaFold-generated integrin structures and experimentally determined reference structures |
| `Supplementary/` | Supplementary figures, tables, validation outputs, and supporting materials |

Individual folder README files provide further details where available.

## Interpretation and Scope

- Binding affinities and ligand poses generated in this study are **computational predictions**, not direct experimental measurements.
- The *S. purpuratus* and *N. vectensis* heterodimers are **candidate subunit combinations**; their physiological pairings have not been established experimentally.
- Similar predicted ligand preferences and binding-pocket architectures are consistent with potential evolutionary conservation, but do not independently establish receptor-mediated thyroid hormone signaling.

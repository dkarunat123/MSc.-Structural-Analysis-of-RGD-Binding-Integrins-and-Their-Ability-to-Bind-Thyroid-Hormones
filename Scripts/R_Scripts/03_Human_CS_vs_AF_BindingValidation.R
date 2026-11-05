
# ==============================================================================
# SCRIPT 03: HUMAN CRYSTAL STRUCTURE VS ALPHAFOLD BINDING VALIDATION
# ==============================================================================
#
# PROJECT:
# Comparative Structural Analysis of RGD-Binding Integrins and Their Ability
# to Bind Thyroid Hormones
#
# PURPOSE:
# Compare predicted ligand-binding affinity rankings obtained by docking the
# same nine ligands to (1) experimentally determined human integrin structures
# and (2) the corresponding AlphaFold 3-predicted receptor models.
#
# This is a VALIDATION OF DOCKING-PROFILE CONSISTENCY BETWEEN STRUCTURAL MODELS.
# The binding affinities in both columns are computational predictions; the
# crystal-structure column does NOT represent experimental binding measurements.
#
# INPUT:
# Results/02_Prepare_Integrin_BindingData_RESULTS/Human_CS_vs_AF_BindingData.csv
#   Generated automatically by 02_Prepare_Integrin_BindingData.R.
#   Each row contains one receptor-ligand pair and its predicted -log10(KD)
#   using the crystal structure and the AlphaFold model.
#
# DESIGN:
#   Five human integrin receptor pairs, each evaluated with nine ligands:
#     alphaVbeta1   (PDB 8W30)
#     alphaVbeta8   (PDB 8TCF)
#     alphaIIbbeta3 (PDB 8T2V)
#     alphaVbeta3   (PDB 1L5G)
#     alpha5beta1   (PDB 4WK4)
#
# STATISTICAL APPROACH:
#   Spearman rank correlation (stats::cor.test)
#   method = "spearman", exact = FALSE, alternative = "two.sided"
#   n = 9 ligand pairs for each integrin
#
# Spearman rho tests whether ligands predicted to bind relatively strongly
# to one structural model also tend to rank strongly in the other model.
# A positive rho indicates agreement in ligand RANK ORDER, not necessarily
# agreement in the absolute predicted affinities or ligand docking poses.
#
# Asymptotic p-values reproduce the settings of the original
# Human_CS_vs_AF_Spearman.R analysis. No multiple-testing correction was
# applied to these five exploratory/model-benchmarking comparisons.
#
# Tied scores receive average ranks, as in R's Spearman implementation.
# The alpha5beta1 crystal profile contains one tied pair of ligand scores.
#
# THESIS REFERENCES:
# Methods: Section 2.8.1
# Results: Sections 3.4.1-3.4.2; Table 5
#
# WORKING DIRECTORY:
# Set the working directory to R_Scripts manually before running this script.
# The script does not change the working directory or the source CSV.
#
# OUTPUTS:
# Results/03_Human_CS_vs_AF_BindingValidation_RESULTS/
# ==============================================================================


# ==============================================================================
# 0. INPUT AND OUTPUT DIRECTORIES
# ==============================================================================

# Read the standardized CSV directly from the results of Script 02.
# This is an intermediate analysis input, so it need not be manually copied
# into Data/. The original Excel workbooks remain in Data/ for Script 02.

INPUT_FILE <- file.path(
  "Results",
  "02_Prepare_Integrin_BindingData_RESULTS",
  "Human_CS_vs_AF_BindingData.csv"
)

# Give this analysis its own results folder; create the folder once.

OUTPUT_DIR <- file.path(
  "Results",
  "03_Human_CS_vs_AF_BindingValidation_RESULTS"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Confirm the prepared input file exists.

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found: ",
    INPUT_FILE,
    "\nRun 02_Prepare_Integrin_BindingData.R first to create this CSV."
  )
}

# No external R packages are required. read.csv(), cor.test(), and rank()
# are all available in a standard R installation.


# ==============================================================================
# 1. DEFINE THE EXPECTED HUMAN RECEPTORS AND LIGANDS
# ==============================================================================

# Historical receptor IDs are retained so that the profiles can be matched
# directly to the identifiers used in docking and in Script 02.
#
# The display labels and PDB IDs correspond to the five comparisons in Table 5.

RECEPTORS <- data.frame(
  Receptor = c(
    "avb1", "avb8", "aiibb3", "avb3", "a5b1"
  ),
  
  Integrin = c(
    "alphaVbeta1",
    "alphaVbeta8",
    "alphaIIbbeta3",
    "alphaVbeta3",
    "alpha5beta1"
  ),
  
  Crystal_PDB = c(
    "8W30", "8TCF", "8T2V", "1L5G", "4WK4"
  ),
  
  CrystalReceptor = c(
    "avb1_8W30",
    "avb8_8TCF",
    "aIIbb3_8T2V",
    "avb3_1L5G",
    "a5b1_4WK4"
  ),
  
  stringsAsFactors = FALSE
)

# Ligand order is used for reproducible exports.
# Correlation coefficients do not depend on the order of correctly
# matched ligand pairs.

LIGANDS <- c(
  "Cilengitide",
  "T3",
  "T4",
  "TRIAC",
  "TETRAC",
  "sT3",
  "sT4",
  "rT3",
  "T4AM"
)


# ==============================================================================
# 2. IMPORT AND VALIDATE THE PREPARED BINDING DATA
# ==============================================================================

# Import the output generated by Script 02.
#
# check.names = FALSE prevents R from modifying the source column labels.

dat <- utils::read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Define the columns required for this analysis.

required_columns <- c(
  "Receptor",
  "CrystalReceptor",
  "Crystal_PDB",
  "Ligand",
  "Crystal_NegLog10KD",
  "AlphaFold_NegLog10KD"
)

# Verify that all required columns are present.

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# Restrict subsequent calculations to the required columns.

dat <- dat[, required_columns, drop = FALSE]


# ------------------------------------------------------------------------------
# 2.1 STANDARDIZE AND CHECK IDENTIFIERS
# ------------------------------------------------------------------------------

# Remove accidental leading or trailing whitespace from text identifiers.
#
# This does not alter binding-affinity scores or biological identities.

for (column in c(
  "Receptor",
  "CrystalReceptor",
  "Crystal_PDB",
  "Ligand"
)) {
  
  dat[[column]] <- trimws(dat[[column]])
  
  if (
    anyNA(dat[[column]]) ||
    any(!nzchar(dat[[column]]))
  ) {
    stop(
      "Missing or empty identifier in column: ",
      column
    )
  }
}


# ------------------------------------------------------------------------------
# 2.2 CHECK DATASET COMPLETENESS
# ------------------------------------------------------------------------------

# The analysis requires:
#
#   5 human integrin pairs
#   9 ligands per integrin
#   45 paired binding-affinity measurements
#
# The two structural models for each integrin must be compared using
# the same nine ligands.

if (nrow(dat) != 45L) {
  stop(
    "Expected 45 receptor-ligand pairs; found ",
    nrow(dat),
    "."
  )
}

# Confirm that all five expected human integrins are present.

if (!setequal(
  unique(dat$Receptor),
  RECEPTORS$Receptor
)) {
  stop(
    "Input receptor IDs do not match the five expected human integrins."
  )
}

# Check for duplicated receptor-ligand combinations.
#
# Duplicate combinations would incorrectly count an observation
# more than once when calculating the correlations.

if (anyDuplicated(
  dat[, c("Receptor", "Ligand")]
)) {
  stop(
    "Duplicate receptor-ligand pairs found in the input data."
  )
}


# ------------------------------------------------------------------------------
# 2.3 VALIDATE NUMERICAL BINDING-AFFINITY SCORES
# ------------------------------------------------------------------------------

# Binding-affinity scores are already expressed as -log10(KD).
#
# They were calculated and validated against the original raw KD values
# in Script 02.
#
# IMPORTANT:
# These values must NOT be logarithmically transformed again.
#
# Higher -log10(KD) values represent stronger predicted binding.

score_columns <- c(
  "Crystal_NegLog10KD",
  "AlphaFold_NegLog10KD"
)

for (column in score_columns) {
  
  if (
    !is.numeric(dat[[column]]) ||
    any(!is.finite(dat[[column]]))
  ) {
    stop(
      "Column must contain finite numeric affinities: ",
      column
    )
  }
}


# ------------------------------------------------------------------------------
# 2.4 VALIDATE LIGAND AND CRYSTAL-STRUCTURE ASSIGNMENTS
# ------------------------------------------------------------------------------

# For each integrin:
#
#   - Verify exactly nine ligands.
#   - Confirm that all expected ligands are represented.
#   - Confirm that the correct experimental PDB structure is assigned.
#
# This prevents data from different receptor models being compared
# accidentally.

for (i in seq_len(nrow(RECEPTORS))) {
  
  info <- RECEPTORS[i, ]
  
  block <- dat[
    dat$Receptor == info$Receptor,
    ,
    drop = FALSE
  ]
  
  if (
    nrow(block) != length(LIGANDS) ||
    !setequal(block$Ligand, LIGANDS)
  ) {
    stop(
      "Missing, duplicated, or unexpected ligands for ",
      info$Receptor
    )
  }
  
  if (
    !all(block$Crystal_PDB == info$Crystal_PDB) ||
    !all(block$CrystalReceptor == info$CrystalReceptor)
  ) {
    stop(
      "Crystal-structure identifier mismatch for ",
      info$Receptor
    )
  }
}

message(
  "Validated ",
  nrow(dat),
  " measurements across ",
  nrow(RECEPTORS),
  " human integrins and ",
  length(LIGANDS),
  " ligands."
)


# ==============================================================================
# 3. CALCULATE CRYSTAL-VERSUS-ALPHAFOLD SPEARMAN CORRELATIONS
# ==============================================================================

# PURPOSE:
#
# Evaluate whether the two structural representations of the same
# integrin produce similar predicted ligand-binding affinity rankings.
#
# A Spearman rank correlation measures the monotonic association
# between two variables.
#
# In this study, the two variables are:
#
#   X = Predicted -log10(KD) from docking to the experimental structure
#   Y = Predicted -log10(KD) from docking to the AlphaFold structure
#
# Both values are predicted binding affinities, not experimentally
# measured ligand dissociation constants.
#
# INTERPRETATION:
#
#   rho = +1: Perfect agreement in ligand ranking
#   rho =  0: No monotonic rank association
#   rho = -1: Perfect reversal of ligand ranking
#
# Positive correlations indicate that ligands predicted to bind
# relatively strongly to one structural model also tend to bind
# relatively strongly to the other model.
#
# Correlation does not establish agreement in absolute affinity
# magnitudes or docking poses.
#
# Statistical settings:
#
#   method = "spearman":
#     Calculate Spearman's rank correlation coefficient.
#
#   exact = FALSE:
#     Use asymptotic p-value estimation, matching the original
#     analysis and thesis Methods section.
#
#   alternative = "two.sided":
#     Test for either positive or negative rank association.
#
# Five separate correlations are calculated, one for each integrin.
# Each correlation uses nine paired ligand-binding predictions.
#
# The p-values are reported without multiple-testing correction,
# matching the original descriptive benchmarking analysis.


# Create containers for statistical summaries and ligand-level rank tables.

results <- vector(
  "list",
  nrow(RECEPTORS)
)

ligand_rank_tables <- vector(
  "list",
  nrow(RECEPTORS)
)


# ------------------------------------------------------------------------------
# 3.1 PERFORM EACH RECEPTOR COMPARISON
# ------------------------------------------------------------------------------

for (i in seq_len(nrow(RECEPTORS))) {
  
  info <- RECEPTORS[i, ]
  
  # Extract all nine ligand observations for the current integrin.
  
  block <- dat[
    dat$Receptor == info$Receptor,
    ,
    drop = FALSE
  ]
  
  # Arrange ligands in a consistent order for exporting.
  
  block <- block[
    match(LIGANDS, block$Ligand),
    ,
    drop = FALSE
  ]
  
  # Extract the paired binding-affinity values.
  
  crystal <- block$Crystal_NegLog10KD
  
  alphafold <- block$AlphaFold_NegLog10KD
  
  
  # A correlation cannot be calculated if either vector has
  # no variation among its ligand scores.
  
  if (
    length(unique(crystal)) < 2L ||
    length(unique(alphafold)) < 2L
  ) {
    stop(
      "Insufficient variation in scores for ",
      info$Integrin
    )
  }
  
  
  # --------------------------------------------------------------------------
  # 3.1.1 SPEARMAN RANK CORRELATION TEST
  # --------------------------------------------------------------------------
  
  # Calculate Spearman rho and its associated two-sided p-value.
  #
  # The nine ligands are paired by their identifiers, ensuring
  # that each crystal-derived score is compared with the
  # corresponding AlphaFold-derived score.
  
  test <- stats::cor.test(
    crystal,
    alphafold,
    method = "spearman",
    exact = FALSE,
    alternative = "two.sided"
  )
  
  
  # --------------------------------------------------------------------------
  # 3.1.2 SAVE THE STATISTICAL RESULTS
  # --------------------------------------------------------------------------
  
  # Store:
  #
  #   Integrin: Human integrin identity
  #   Crystal_PDB: Experimental structure identifier
  #   N_Ligands: Number of paired ligand observations
  #   Spearman_rho: Rank correlation coefficient
  #   P_value: Two-sided asymptotic p-value
  #
  # Full numerical precision is retained here.
  # Rounding is performed only for the separate thesis table.
  
  results[[i]] <- data.frame(
    Integrin = info$Integrin,
    Crystal_PDB = info$Crystal_PDB,
    N_Ligands = length(crystal),
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # 3.1.3 CALCULATE INDIVIDUAL LIGAND RANKS
  # --------------------------------------------------------------------------
  
  # Rank the nine ligands according to predicted binding affinity.
  #
  # Because higher -log10(KD) indicates stronger predicted binding,
  # ranks are calculated using NEGATIVE score values.
  #
  # Rank 1 = strongest predicted binder
  # Rank 9 = weakest predicted binder
  #
  # ties.method = "average":
  #   Ligands with identical predicted affinities receive
  #   the average of their occupied rank positions.
  #
  # This is particularly relevant for alpha5beta1, whose crystal
  # profile contains a tied affinity value.
  
  crystal_rank <- rank(
    -crystal,
    ties.method = "average"
  )
  
  af_rank <- rank(
    -alphafold,
    ties.method = "average"
  )
  
  
  # Export the original paired scores alongside their ranks.
  #
  # Rank_Difference indicates how much a ligand's relative position
  # changes between crystal-structure and AlphaFold docking.
  #
  # This provides additional interpretive detail beyond rho alone.
  
  ligand_rank_tables[[i]] <- data.frame(
    Integrin = info$Integrin,
    Crystal_PDB = info$Crystal_PDB,
    Ligand = block$Ligand,
    Crystal_NegLog10KD = crystal,
    AlphaFold_NegLog10KD = alphafold,
    Crystal_Rank = crystal_rank,
    AlphaFold_Rank = af_rank,
    Rank_Difference = crystal_rank - af_rank,
    stringsAsFactors = FALSE
  )
  
  
  # Display the result for this receptor in the R console.
  
  message(
    sprintf(
      "%-16s (%s): rho = %.3f, p = %.6g",
      info$Integrin,
      info$Crystal_PDB,
      unname(test$estimate),
      test$p.value
    )
  )
  
}


# ------------------------------------------------------------------------------
# 3.2 COMBINE ALL FIVE RECEPTOR COMPARISONS
# ------------------------------------------------------------------------------

# Combine the individual receptor results into complete tables.

spearman_results <- do.call(
  rbind,
  results
)

ligand_rank_results <- do.call(
  rbind,
  ligand_rank_tables
)

# Remove inherited row numbering from the data frames.

rownames(spearman_results) <- NULL
rownames(ligand_rank_results) <- NULL


# ==============================================================================
# 4. EXPORT NUMERICAL CORRELATION RESULTS AND LIGAND RANKS
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 SPEARMAN CORRELATION RESULTS
# ------------------------------------------------------------------------------

# Export the five Spearman correlations and associated p-values.
#
# These columns match the original Human_CS_vs_AF_Spearman.R analysis.
#
# Unlike the original script, full numerical precision is retained
# in this file so the statistical estimates are not prematurely rounded.

utils::write.csv(
  spearman_results,
  
  file.path(
    OUTPUT_DIR,
    "Human_CS_vs_AF_Spearman.csv"
  ),
  
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 4.2 INDIVIDUAL LIGAND RANKING DATA
# ------------------------------------------------------------------------------

# Export all 45 ligand observations with their:
#
#   - Predicted crystal-derived affinity scores
#   - Predicted AlphaFold-derived affinity scores
#   - Crystal-model ligand ranks
#   - AlphaFold-model ligand ranks
#   - Differences between the two rankings
#
# These data make it easier to determine which ligands contribute
# to disagreements between the models.

utils::write.csv(
  ligand_rank_results,
  
  file.path(
    OUTPUT_DIR,
    "Human_CS_vs_AF_LigandRanks.csv"
  ),
  
  row.names = FALSE
)


# ==============================================================================
# 5. FORMAT THE RESULTS FOR THESIS TABLE 5
# ==============================================================================

# The previous sections retain full-precision numerical values.
#
# This section creates a separate formatted table suitable for
# comparison against Table 5 in the thesis.
#
# Reporting conventions:
#
#   Spearman rho: three decimal places
#   p < 0.001:    "<0.001"
#   Other p-values: three or four decimal places
#
# IMPORTANT:
# Formatting does not change the statistical calculations.
#
# The preformatting of p to four significant digits reproduces
# the original script's reporting convention.

format_table_p <- function(p) {
  
  if (p < 0.001) {
    return("<0.001")
  }
  
  p <- signif(
    p,
    digits = 4
  )
  
  if (p >= 0.1) {
    sprintf("%.3f", p)
  } else {
    sprintf("%.4f", p)
  }
  
}


# Construct the reporting table in the original receptor order.
#
# Greek-letter integrin labels are included for readability.

table_5 <- data.frame(
  Integrin = c(
    "αVβ1",
    "αVβ8",
    "αIIbβ3",
    "αVβ3",
    "α5β1"
  ),
  
  Crystal_PDB = spearman_results$Crystal_PDB,
  
  Spearman_rho = sprintf(
    "%.3f",
    spearman_results$Spearman_rho
  ),
  
  P_value = vapply(
    spearman_results$P_value,
    format_table_p,
    character(1)
  ),
  
  stringsAsFactors = FALSE
)


# Export the formatted thesis table.

utils::write.csv(
  table_5,
  
  file.path(
    OUTPUT_DIR,
    "Table5_Human_CS_vs_AF_Spearman.csv"
  ),
  
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Print the formatted Table 5 in the R console.

cat(
  "\nHuman crystal-vs-AlphaFold validation (Table 5):\n"
)

print(
  table_5,
  row.names = FALSE
)


# ==============================================================================
# 6. ANALYSIS COMPLETION
# ==============================================================================

# Confirm successful completion and report the output location.

message(
  "Completed human crystal-vs-AlphaFold rank comparisons."
)

message(
  "Results saved to: ",
  OUTPUT_DIR
)


# ==============================================================================
# END OF SCRIPT
# ==============================================================================
#
# OUTPUT SUMMARY
#
# All outputs are saved to:
#
# Results/03_Human_CS_vs_AF_BindingValidation_RESULTS/
#
# Full-precision statistical comparison results:
#   Human_CS_vs_AF_Spearman.csv
#
# Per-ligand affinity scores and ranks for auditing correlations:
#   Human_CS_vs_AF_LigandRanks.csv
#
# Formatted results corresponding to thesis Table 5:
#   Table5_Human_CS_vs_AF_Spearman.csv
#
# Total: 3 output files
#
# ==============================================================================

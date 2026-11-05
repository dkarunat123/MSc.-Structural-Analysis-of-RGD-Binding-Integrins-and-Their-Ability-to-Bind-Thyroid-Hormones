
# ==============================================================================
# SCRIPT 06: CROSS-SPECIES LIGAND-AFFINITY PROFILE COMPARISON
# ==============================================================================
#
# PROJECT:
# Comparative Structural Analysis of RGD-Binding Integrins and Their Ability
# to Bind Thyroid Hormones
#
# PURPOSE:
# Compare overall predicted ligand-binding preferences across Homo sapiens,
# Strongylocentrotus purpuratus, and Nematostella vectensis. For each species,
# the median -log10(KD) across its modelled integrin heterodimers is calculated
# for each of nine ligands, producing one nine-ligand profile per species.
#
# The three species profiles are compared using (1) a clustered heatmap of
# median affinity values and (2) Spearman rank correlations of the nine ligand
# medians. The clustered median heatmap corresponds to thesis Figure 16.
#
# INPUTS (outputs of earlier R scripts; do not copy these into Data/):
#   Results/02_Prepare_Integrin_BindingData_RESULTS/
#     BindingAffinity_NegLog10KD_Wide.csv
#     Receptor_Metadata_OMA.csv
#
#   Results/05_Species_BindingAffinity_Analysis_RESULTS/
#     Species_Ligand_DescriptiveStatistics.csv
#
# STUDY DESIGN:
#   H. sapiens:      8 modelled receptors x 9 ligands
#   S. purpuratus:  24 modelled receptors x 9 ligands
#   N. vectensis:   16 modelled receptors x 9 ligands
#   Total:          48 modelled receptors, 432 predicted affinity scores
#
# STATISTICAL METHODS (following the original Script2B.R, Analysis D):
#   - Summarize each ligand within each species using the MEDIAN predicted
#     -log10(KD) value across that species' receptor models.
#   - Cluster the resulting 3 x 9 species-by-ligand median matrix with
#     Euclidean distance and complete-linkage hierarchical clustering.
#   - Calculate pairwise Spearman rank correlations between the THREE species
#     profiles, each based on NINE corresponding ligand medians.
#   - Cluster the optional Spearman similarity heatmap using distance 1-rho
#     and complete linkage, on a shared colour scale from -1 to +1.
#   - Rank the nine ligands within each species by their median affinity.
#
# SCIENTIFIC INTERPRETATION:
#   - All binding affinities are computational predictions, not experimentally
#     measured KD values.
#   - Input scores are ALREADY -log10(KD), with KD in molar units. Higher
#     values indicate stronger predicted binding; do not log-transform again.
#   - Species profiles summarize DIFFERENT sets of 8, 24, and 16 modelled
#     heterodimers. They are not matched receptor measurements or biological
#     animal replicates. Medians reduce sensitivity to unusual receptor scores.
#   - Spearman rho assesses similarity in ligand RANKINGS across species,
#     not equivalence of absolute KD values or biological receptor function.
#   - No formal species-level significance tests are performed. In particular,
#     modelled heterodimers are not treated as independent biological samples.
#   - This analysis does not establish orthology; OMA-based comparisons are
#     analyzed separately in Script 07.
#
# THESIS REFERENCES:
#   Methods: Section 2.8.4
#   Results: Sections 3.7.1 and 3.7.2; Figure 16
#
# WORKING DIRECTORY:
#   Manually set R's working directory to R_Scripts before running.
#   Do not change the original data files or their directory structure.
#   This script never calls setwd().
#
# OUTPUT DIRECTORY:
#   Results/06_CrossSpecies_LigandProfile_Comparison_RESULTS/
# ==============================================================================


# ==============================================================================
# 0. INPUT AND OUTPUT DIRECTORIES; REQUIRED PACKAGE
# ==============================================================================

# Read prepared, validated data directly from earlier scripts' results folders.
# There is no need to rename or copy inputs after downloading the GitHub repo.

PREPARED_DIR <- file.path(
  "Results", "02_Prepare_Integrin_BindingData_RESULTS"
)

SPECIES_RESULTS_DIR <- file.path(
  "Results", "05_Species_BindingAffinity_Analysis_RESULTS"
)

INPUT_SCORES <- file.path(
  PREPARED_DIR, "BindingAffinity_NegLog10KD_Wide.csv"
)

INPUT_METADATA <- file.path(
  PREPARED_DIR, "Receptor_Metadata_OMA.csv"
)

INPUT_STATISTICS <- file.path(
  SPECIES_RESULTS_DIR, "Species_Ligand_DescriptiveStatistics.csv"
)


# Create one results directory for all outputs of Script 06.

OUTPUT_DIR <- file.path(
  "Results", "06_CrossSpecies_LigandProfile_Comparison_RESULTS"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# pheatmap is used to reproduce the two original clustered heatmaps.
# All numerical calculations and CSV exports use base R / standard packages.

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  
  stop(
    "The 'pheatmap' package is required. Install it using ",
    "install.packages('pheatmap')."
  )
  
}


# Stop with a clear message when a preceding analysis has not been run.

for (path in c(INPUT_SCORES, INPUT_METADATA, INPUT_STATISTICS)) {
  
  if (!file.exists(path)) {
    
    stop(
      "Missing required input: ", path,
      "\nRun Scripts 02 and 05 first, retaining their results directories."
    )
    
  }
  
}


# Preserve a single explicit ligand definition throughout the repository.

LIGANDS <- c(
  "Cilengitide", "T3", "T4", "TRIAC", "TETRAC",
  "sT3", "sT4", "rT3", "T4AM"
)


# Species names and expected receptor counts come from the completed design.
# The same order is used for the exported tables; dendrograms may reorder them.

SPECIES <- data.frame(
  Name = c(
    "H. sapiens",
    "S. purpuratus",
    "N. vectensis"
  ),
  N_Receptors = c(8L, 24L, 16L),
  stringsAsFactors = FALSE
)


# Helper to export tables without R's default row-number index column.

export_csv <- function(x, filename) {
  
  utils::write.csv(
    x,
    file.path(OUTPUT_DIR, filename),
    row.names = FALSE,
    na = ""
  )
  
}



# ==============================================================================
# 1. IMPORT AND VALIDATE THE 48 x 9 INPUT AFFINITY MATRIX
# ==============================================================================

# The prepared score table has receptors in rows and ligands in columns.
# Its values were checked against the corresponding raw KD values in Script 02.

scores <- utils::read.csv(
  INPUT_SCORES,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Metadata supplies the species associated with each receptor identifier.

metadata <- utils::read.csv(
  INPUT_METADATA,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Script 05 provides independently generated within-species medians for
# cross-validation. These are not assumed to be correct without checking.

previous_statistics <- utils::read.csv(
  INPUT_STATISTICS,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Fail early if required columns or study dimensions have changed.

if (!identical(names(scores), c("Receptor", LIGANDS))) {
  stop("Unexpected receptor/ligand columns in Script 02 affinity data.")
}

if (!all(c("Receptor", "Species") %in% names(metadata))) {
  stop("Script 02 metadata must contain Receptor and Species columns.")
}

if (nrow(scores) != 48L || nrow(metadata) != 48L) {
  stop("Expected exactly 48 receptors in both Script 02 input files.")
}


# Confirm that receptor identifiers are complete, unique and consistent.
#
# Without these checks, a species label could accidentally be associated
# with the wrong integrin heterodimer.

if (
  anyNA(scores$Receptor) ||
  anyNA(metadata$Receptor) ||
  any(!nzchar(scores$Receptor)) ||
  any(!nzchar(metadata$Receptor)) ||
  anyDuplicated(scores$Receptor) ||
  anyDuplicated(metadata$Receptor) ||
  !setequal(scores$Receptor, metadata$Receptor)
) {
  
  stop("Missing, duplicated, or inconsistent receptor identifiers.")
  
}


# Match by identifier rather than by row number. If a collaborator reorders a
# CSV, a receptor's species assignment must still remain correct.

metadata <- metadata[
  match(scores$Receptor, metadata$Receptor),
  ,
  drop = FALSE
]


if (
  anyNA(metadata$Species) ||
  !setequal(unique(metadata$Species), SPECIES$Name)
) {
  
  stop("Unexpected or missing species in receptor metadata.")
  
}


# Count receptor models separately for each species and confirm that
# the complete set of 48 models is represented.

observed_counts <- as.integer(
  table(
    factor(
      metadata$Species,
      levels = SPECIES$Name
    )
  )
)

if (!identical(observed_counts, SPECIES$N_Receptors)) {
  stop("Expected 8 human, 24 sea-urchin, and 16 sea-anemone receptors.")
}


# Confirm that all nine affinity columns contain numerical values.

if (!all(
  vapply(
    scores[, LIGANDS, drop = FALSE],
    is.numeric,
    logical(1)
  )
)) {
  
  stop("Every ligand-affinity column must be numeric.")
  
}


# Convert to a numeric receptor x ligand matrix, preserving receptor IDs.

affinity_matrix <- as.matrix(
  scores[, LIGANDS, drop = FALSE]
)

rownames(affinity_matrix) <- scores$Receptor

storage.mode(affinity_matrix) <- "double"


# Check for missing or non-finite predicted affinities.

if (any(!is.finite(affinity_matrix))) {
  stop("Missing, infinite, or nonnumeric predicted affinity values found.")
}


message(
  "Validated ",
  nrow(affinity_matrix),
  " receptor models x ",
  ncol(affinity_matrix),
  " ligands (",
  length(affinity_matrix),
  " values)."
)



# ==============================================================================
# 2. CALCULATE ONE MEDIAN NINE-LIGAND PROFILE PER SPECIES
# ==============================================================================

# A median is calculated separately for each ligand within each species.
#
# For example, the human T4 median comes from eight predicted human receptor
# affinities for T4. The analogous medians use 24 sea-urchin or 16 sea-anemone
# model scores, respectively.
#
# This yields a 3 x 9 matrix:
#
#   Rows    = species
#   Columns = ligands
#   Cells   = median predicted -log10(KD) across that species' receptors
#
# WHY USE THE MEDIAN?
#
# Individual integrin heterodimers sometimes exhibit unusually high or low
# predicted binding affinities. The median summarizes the centre of each
# species-specific ligand distribution while reducing sensitivity to these
# extreme receptor scores.
#
# IMPORTANT:
# We take the median of the ALREADY transformed score values.
# We do not calculate the median raw KD and then transform it; these are not
# generally equivalent for even-sized receptor groups.

median_matrix <- matrix(
  NA_real_,
  nrow = nrow(SPECIES),
  ncol = length(LIGANDS),
  dimnames = list(
    SPECIES$Name,
    LIGANDS
  )
)


# Calculate the nine median values independently for each species.

for (i in seq_len(nrow(SPECIES))) {
  
  species_name <- SPECIES$Name[i]
  
  # Identify the heterodimer models belonging to this species.
  
  receptor_ids <- metadata$Receptor[
    metadata$Species == species_name
  ]
  
  # Extract the corresponding receptor x ligand affinity matrix.
  
  species_matrix <- affinity_matrix[
    receptor_ids,
    LIGANDS,
    drop = FALSE
  ]
  
  if (nrow(species_matrix) != SPECIES$N_Receptors[i]) {
    stop("Unexpected receptor count for ", species_name)
  }
  
  # Calculate the median separately for each of the nine ligand columns.
  #
  # apply(..., 2, median) applies the median function to columns.
  # The resulting nine values constitute the species median profile.
  
  median_matrix[i, ] <- apply(
    species_matrix,
    2,
    stats::median
  )
  
}


if (any(!is.finite(median_matrix))) {
  stop("Species median affinity matrix contains invalid values.")
}


# ------------------------------------------------------------------------------
# 2.1 EXPORT THE SPECIES x LIGAND MEDIAN MATRIX
# ------------------------------------------------------------------------------

# Retain all digits of precision in the CSV.
# Rounding is only performed when displaying values for interpretation.

median_wide <- data.frame(
  Species = rownames(median_matrix),
  median_matrix,
  check.names = FALSE,
  row.names = NULL
)

export_csv(
  median_wide,
  "Figure16_CrossSpecies_MedianAffinity_Matrix.csv"
)


# ------------------------------------------------------------------------------
# 2.2 EXPORT THE MEDIAN PROFILES IN LONG FORMAT
# ------------------------------------------------------------------------------

# Long-format data make it easier for other researchers to reproduce
# plots or import the results into other statistical software.
#
# Each row represents one species-ligand combination.
#
# There are 3 species x 9 ligands = 27 observations in this summary table.

median_long <- do.call(
  rbind,
  lapply(seq_len(nrow(SPECIES)), function(i) {
    
    data.frame(
      Species = SPECIES$Name[i],
      Ligand = LIGANDS,
      N_Receptors = SPECIES$N_Receptors[i],
      Median_NegLog10KD = unname(median_matrix[i, ]),
      stringsAsFactors = FALSE
    )
    
  })
)

rownames(median_long) <- NULL

export_csv(
  median_long,
  "Figure16_CrossSpecies_MedianAffinity_Long.csv"
)



# ==============================================================================
# 3. VALIDATE THE MEDIANS AGAINST SCRIPT 05 DESCRIPTIVE STATISTICS
# ==============================================================================

# The species-specific median scores were already calculated in Script 05.
#
# However, this script independently recalculates the medians from the
# receptor-level data generated by Script 02.
#
# We compare the two sets of results to ensure consistency between scripts.
#
# If the medians disagree, the script stops rather than generating
# cross-species comparisons based on inconsistent values.

expected_columns <- c(
  "Species",
  "Ligand",
  "N_Receptors",
  "Median"
)

if (!all(expected_columns %in% names(previous_statistics))) {
  stop("Script 05 descriptive statistics lack required columns.")
}


# Check for missing values, duplicate combinations or unexpected dimensions.

if (
  nrow(previous_statistics) != nrow(median_long) ||
  anyNA(previous_statistics$Species) ||
  anyNA(previous_statistics$Ligand) ||
  anyDuplicated(previous_statistics[, c("Species", "Ligand")]) ||
  !is.numeric(previous_statistics$Median) ||
  !is.numeric(previous_statistics$N_Receptors)
) {
  
  stop("Invalid/incomplete descriptive statistics from Script 05.")
  
}


# Construct unique identifiers for each species-ligand combination.

key <- function(species, ligand) {
  paste(species, ligand, sep = " || ")
}

previous_keys <- key(
  previous_statistics$Species,
  previous_statistics$Ligand
)

current_keys <- key(
  median_long$Species,
  median_long$Ligand
)


# Match records by species and ligand rather than assuming identical row order.

matched_rows <- match(
  current_keys,
  previous_keys
)

if (
  anyNA(matched_rows) ||
  !setequal(previous_keys, current_keys)
) {
  
  stop("Species-ligand identities differ between Scripts 05 and 06.")
  
}

previous_matched <- previous_statistics[
  matched_rows,
  ,
  drop = FALSE
]


# Ensure the receptor counts and numerical medians are valid.

if (
  anyNA(previous_matched$Median) ||
  any(!is.finite(previous_matched$Median)) ||
  anyNA(previous_matched$N_Receptors) ||
  any(previous_matched$N_Receptors != median_long$N_Receptors)
) {
  
  stop("Script 05 median values or receptor counts are invalid.")
  
}


# Calculate the largest absolute difference between the two sets of medians.
#
# Small differences caused by floating-point representation are acceptable.
# A tolerance of 1e-9 is used; larger discrepancies stop the script.

median_difference <- max(
  abs(
    previous_matched$Median -
      median_long$Median_NegLog10KD
  )
)

if (median_difference > 1e-9) {
  
  stop(
    "Medians do not agree with Script 05; maximum difference = ",
    signif(median_difference, 8)
  )
  
}


message(
  "All 27 median affinity values agree with Script 05 (max difference: ",
  format(median_difference, scientific = TRUE),
  ")."
)



# ==============================================================================
# 4. GENERATE FIGURE 16: CLUSTERED CROSS-SPECIES MEDIAN HEATMAP
# ==============================================================================

# THESIS OUTPUT:
#
# Figure 16 - Cross-species comparison of median predicted ligand-binding
# affinity profiles.
#
# The heatmap contains:
#
#   Rows:    H. sapiens, S. purpuratus, N. vectensis
#   Columns: Nine ligands
#   Values:  Median predicted -log10(KD) within each species
#
# Both species and ligands are hierarchically clustered using Euclidean
# distance and complete linkage.
#
# EUCLIDEAN DISTANCE:
#
# Euclidean distance measures differences in the numerical affinity profiles.
# Two species with similar median affinities across the nine ligands will
# have a relatively small distance.
#
# Similarly, two ligands with comparable median affinities across species
# will tend to cluster together.
#
# COMPLETE LINKAGE:
#
# The distance between two clusters is determined by the greatest pairwise
# distance between their members.
#
# IMPORTANT:
#
# No standardization, Z-scoring, or rescaling is performed.
# Heatmap colours represent the actual median predicted -log10(KD) values.
#
# The original blue -> yellow -> red palette from Script2A/Script2B is used
# for visual consistency with the thesis analysis.


# ------------------------------------------------------------------------------
# 4.1 DEFINE THE ORIGINAL AFFINITY COLOUR PALETTE
# ------------------------------------------------------------------------------

affinity_colours <- grDevices::colorRampPalette(
  c(
    "#313695",
    "#74ADD1",
    "#FFFFBF",
    "#FDAE61",
    "#A50026"
  )
)(100)


# ------------------------------------------------------------------------------
# 4.2 GENERATE AND SAVE FIGURE 16
# ------------------------------------------------------------------------------

pheatmap::pheatmap(
  median_matrix,
  
  # Independently cluster the three species and nine ligands.
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  
  clustering_method = "complete",
  
  # Display original median affinity values, not standardized scores.
  
  scale = "none",
  color = affinity_colours,
  
  # Remove cell borders and adjust text appearance.
  
  border_color = NA,
  fontsize = 10,
  angle_col = 45,
  
  main = "Median ligand-binding profiles across species",
  
  # Export the thesis figure as a vector PDF.
  
  filename = file.path(
    OUTPUT_DIR,
    "Figure16_CrossSpecies_MedianAffinity_Heatmap.pdf"
  ),
  
  width = 9,
  height = 5
)



# ==============================================================================
# 5. COMPARE SPECIES MEDIAN PROFILES USING SPEARMAN CORRELATION
# ==============================================================================

# PURPOSE:
#
# Determine whether the three species exhibit similar RELATIVE preferences
# among the nine ligands.
#
# This is different from comparing the absolute predicted affinities.
#
# For example, one species could have generally higher affinity scores
# than another but still rank the nine ligands in a similar order.
#
# Spearman rank correlation quantifies this rank-order similarity.
#
# INTERPRETATION:
#
#   rho = +1: Identical ligand rankings
#   rho =  0: No monotonic association
#   rho = -1: Opposite ligand rankings
#
# Positive correlations therefore indicate that ligands with relatively
# high median affinities in one species also tend to have high median
# affinities in another species.
#
# IMPORTANT:
#
# This is NOT the within-species ligand-to-ligand analysis of Figure 15.
#
# Figure 15 compares ligand columns across individual receptor models.
#
# Script 06 compares species rows across nine corresponding median values.
#
# Each comparison uses the same nine ligand identities, ensuring that
# corresponding medians are correctly paired.
#
# No p-values are calculated. The original Script2B.R Analysis D reported
# descriptive Spearman coefficients, not formal species-level tests.


# ------------------------------------------------------------------------------
# 5.1 CALCULATE THE THREE-SPECIES SPEARMAN CORRELATION MATRIX
# ------------------------------------------------------------------------------

# The median matrix has:
#
#   Rows    = species
#   Columns = ligands
#
# Since cor() normally compares columns, transpose the matrix using t().
#
# The result is a 3 x 3 correlation matrix comparing the three species.

species_spearman <- stats::cor(
  t(median_matrix),
  method = "spearman",
  use = "everything"
)


# Validate matrix dimensions and numerical values.

if (
  !identical(dim(species_spearman), c(3L, 3L)) ||
  any(!is.finite(species_spearman)) ||
  !isTRUE(
    all.equal(
      species_spearman,
      t(species_spearman),
      tolerance = 1e-12
    )
  )
) {
  
  stop("Cross-species Spearman correlation matrix is invalid.")
  
}


# ------------------------------------------------------------------------------
# 5.2 EXPORT THE COMPLETE CORRELATION MATRIX
# ------------------------------------------------------------------------------

# Preserve the full numerical precision of the correlation coefficients.
#
# The diagonal contains 1 because each species profile is perfectly
# correlated with itself.

correlation_wide <- data.frame(
  Species = rownames(species_spearman),
  species_spearman,
  check.names = FALSE,
  row.names = NULL
)

export_csv(
  correlation_wide,
  "CrossSpecies_MedianProfile_Spearman_Matrix.csv"
)


# ------------------------------------------------------------------------------
# 5.3 EXPORT THE THREE UNIQUE SPECIES COMPARISONS
# ------------------------------------------------------------------------------

# With three species, there are three distinct pairwise comparisons:
#
#   H. sapiens vs S. purpuratus
#   H. sapiens vs N. vectensis
#   S. purpuratus vs N. vectensis
#
# Reverse-order duplicates are unnecessary because Spearman correlation
# is symmetric.

species_pairs <- utils::combn(
  SPECIES$Name,
  2
)

correlation_pairs <- data.frame(
  Species_1 = species_pairs[1, ],
  Species_2 = species_pairs[2, ],
  
  N_Ligands = rep(
    length(LIGANDS),
    ncol(species_pairs)
  ),
  
  Spearman_rho = vapply(
    seq_len(ncol(species_pairs)),
    
    function(i) {
      
      species_spearman[
        species_pairs[1, i],
        species_pairs[2, i]
      ]
      
    },
    
    numeric(1)
  ),
  
  stringsAsFactors = FALSE
)

export_csv(
  correlation_pairs,
  "CrossSpecies_MedianProfile_Spearman_Pairs.csv"
)


# Display the results in the R console for immediate inspection.

cat(
  "\nCross-species Spearman correlations of ligand median profiles:\n"
)

print(
  round(species_spearman, 3)
)



# ==============================================================================
# 6. EXPORT THE CROSS-SPECIES SPEARMAN CORRELATION HEATMAP
# ==============================================================================

# The original Script2B.R Analysis D4 also generated a small heatmap showing
# correlations among the three species-level median profiles.
#
# We retain this as a SUPPORTING visualization.
#
# Figure 16 itself is the absolute-affinity median heatmap in Section 4.
#
# The correlation heatmap clusters species according to:
#
#   Distance = 1 - Spearman rho
#
# Under this transformation:
#
#   rho = +1 corresponds to distance 0
#   rho =  0 corresponds to distance 1
#   rho = -1 corresponds to distance 2
#
# Complete-linkage clustering is used, matching the original analysis.
#
# The colour scale is fixed from -1 to +1.


# ------------------------------------------------------------------------------
# 6.1 DEFINE THE CORRELATION COLOUR PALETTE
# ------------------------------------------------------------------------------

correlation_colours <- grDevices::colorRampPalette(
  c(
    "#313695",
    "#74ADD1",
    "#FFFFBF",
    "#F46D43",
    "#A50026"
  )
)(100)


# ------------------------------------------------------------------------------
# 6.2 CALCULATE CORRELATION-BASED DISTANCES
# ------------------------------------------------------------------------------

species_distance <- stats::as.dist(
  1 - species_spearman
)


# ------------------------------------------------------------------------------
# 6.3 GENERATE AND SAVE THE SUPPORTING CORRELATION HEATMAP
# ------------------------------------------------------------------------------

pheatmap::pheatmap(
  species_spearman,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  clustering_distance_rows = species_distance,
  clustering_distance_cols = species_distance,
  
  clustering_method = "complete",
  
  color = correlation_colours,
  
  breaks = seq(
    -1,
    1,
    length.out = 101
  ),
  
  border_color = NA,
  
  display_numbers = TRUE,
  number_format = "%.2f",
  
  fontsize = 11,
  
  main = "Cross-species similarity of median ligand-preference profiles",
  
  filename = file.path(
    OUTPUT_DIR,
    "CrossSpecies_MedianProfile_Spearman_Heatmap.pdf"
  ),
  
  width = 7,
  height = 6
)



# ==============================================================================
# 7. RANK LIGANDS BY MEDIAN PREDICTED AFFINITY WITHIN EACH SPECIES
# ==============================================================================

# PURPOSE:
#
# Provide a straightforward numerical comparison of predicted ligand
# preferences within each species.
#
# This helps identify ligands that consistently exhibit relatively
# strong or weak predicted affinities across the three species.
#
# Ranking is based on the median predicted -log10(KD) values:
#
#   Rank 1 = strongest median predicted binding
#   Rank 9 = weakest median predicted binding
#
# Because HIGHER -log10(KD) indicates stronger binding, the negative
# median values are supplied to rank().
#
# Exact ties, if present, receive their average rank.
#
# The rankings use FULL-PRECISION medians rather than rounded values.
# This is important because two ligands with similar displayed medians
# may still have slightly different underlying values.


# ------------------------------------------------------------------------------
# 7.1 CALCULATE THE SPECIES-SPECIFIC LIGAND RANKS
# ------------------------------------------------------------------------------

ligand_ranks <- median_long

ligand_ranks$Median_Rank <- NA_real_


for (i in seq_len(nrow(SPECIES))) {
  
  # Identify the nine ligand records for the current species.
  
  rows <- which(
    ligand_ranks$Species == SPECIES$Name[i]
  )
  
  # Rank ligands from highest to lowest median predicted affinity.
  
  ligand_ranks$Median_Rank[rows] <- rank(
    -ligand_ranks$Median_NegLog10KD[rows],
    ties.method = "average"
  )
  
}


# ------------------------------------------------------------------------------
# 7.2 ARRANGE AND EXPORT THE RANK TABLE
# ------------------------------------------------------------------------------

# Arrange records by species, then from strongest to weakest ligand.

ligand_ranks <- ligand_ranks[
  order(
    match(ligand_ranks$Species, SPECIES$Name),
    ligand_ranks$Median_Rank,
    ligand_ranks$Ligand
  ),
  ,
  drop = FALSE
]

rownames(ligand_ranks) <- NULL

export_csv(
  ligand_ranks,
  "CrossSpecies_Ligand_MedianRanks.csv"
)



# ==============================================================================
# 8. EXPORT QUALITY-CONTROL AND ANALYSIS SUMMARY
# ==============================================================================

# The following table records dataset sizes, statistical settings, and
# the agreement between the independently calculated medians and Script 05.
#
# These values describe the computational dataset.
#
# The receptor counts are NOT biological replicate sample sizes.

qc <- data.frame(
  
  Check = c(
    "Total_receptors",
    "H_sapiens_receptors",
    "S_purpuratus_receptors",
    "N_vectensis_receptors",
    "Ligands_per_species",
    "Total_receptor_ligand_scores",
    "Species_median_affinity_values",
    "Independent_species_pairs_compared",
    "Paired_ligand_medians_per_correlation",
    "Maximum_median_difference_vs_script05",
    "Figure16_clustering_distance",
    "Figure16_clustering_linkage",
    "Species_profile_correlation_method"
  ),
  
  Value = c(
    nrow(affinity_matrix),
    observed_counts,
    length(LIGANDS),
    length(affinity_matrix),
    length(median_matrix),
    nrow(correlation_pairs),
    length(LIGANDS),
    format(
      median_difference,
      digits = 15,
      scientific = TRUE
    ),
    "Euclidean",
    "complete",
    "Spearman"
  ),
  
  stringsAsFactors = FALSE
)


# Export the quality-control summary.

export_csv(
  qc,
  "CrossSpecies_Profile_ValidationSummary.csv"
)


# Confirm successful completion.

message(
  "\nCross-species ligand-profile analysis completed successfully."
)

message(
  "All results saved to: ",
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
# Results/06_CrossSpecies_LigandProfile_Comparison_RESULTS/
#
# Thesis Figure 16 - clustered species median affinity heatmap:
#   Figure16_CrossSpecies_MedianAffinity_Heatmap.pdf
#
# Numerical species x ligand median affinity matrix (3 x 9):
#   Figure16_CrossSpecies_MedianAffinity_Matrix.csv
#
# Long-form table of all 27 species-ligand median scores:
#   Figure16_CrossSpecies_MedianAffinity_Long.csv
#
# Full-precision Spearman similarity matrix among the three species:
#   CrossSpecies_MedianProfile_Spearman_Matrix.csv
#
# Three unique pairwise cross-species Spearman coefficients:
#   CrossSpecies_MedianProfile_Spearman_Pairs.csv
#
# Supporting heatmap of the three species-level Spearman correlations:
#   CrossSpecies_MedianProfile_Spearman_Heatmap.pdf
#
# Rankings of the nine ligands by median affinity within each species:
#   CrossSpecies_Ligand_MedianRanks.csv
#
# Input validation, dataset counts, and analysis settings:
#   CrossSpecies_Profile_ValidationSummary.csv
#
# Total: 8 output files (2 PDF figures + 6 CSV tables).
#
# This script uses NO bootstrap calculations and NO formal species-level
# significance tests. OMA-specific comparisons are performed in Script 07.
# ==============================================================================

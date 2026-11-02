
# ==============================================================================
# SCRIPT 02: PREPARATION OF INTEGRIN BINDING-AFFINITY AND OMA DATA
# ==============================================================================
#
# PROJECT:
# Comparative Structural Analysis of RGD-Binding Integrins and Their Ability
# to Bind Thyroid Hormones
#
# PURPOSE:
# Create validated, reusable datasets for the subsequent statistical analyses.
# This script does NOT conduct clustering, correlations, hypothesis tests,
# or generate figures. Those analyses are performed in Scripts 03-07.
#
# INPUT FILES (place all three in R_Scripts/Data/):
#   BindingScores_StatisticalAnalysis.xlsx : 48 receptors x 9 ligands;
#                                           values ALREADY equal -log10(KD).
#   All_Species_BindingScores.xlsx         : original predicted KD (molar),
#                                           plus human crystal/AF comparisons.
#   OMA_Subunit_Relationships.xlsx         : reviewed alpha- and beta-subunit
#                                           OMA evidence and interpretations.
#
# BIOLOGICAL DESIGN:
#   H. sapiens       :  8 integrin heterodimers
#   S. purpuratus    : 24 integrin heterodimers
#   N. vectensis     : 16 integrin heterodimers
#   Total            : 48 heterodimers x 9 ligands = 432 measurements
#
# CRITICAL DISTINCTIONS:
#   - Higher -log10(KD) indicates stronger predicted binding.
#   - Lower raw KD (M) indicates stronger predicted binding.
#   - RAW KD is preserved for Figure 18. Median(KD) must NOT be replaced with
#     10^(-median(-log10(KD))); these differ for some even-sized groups.
#   - OMA supports relationships between INDIVIDUAL alpha or beta subunits,
#     not orthology of entire alpha-beta receptor heterodimers.
#   - OMA analysis groups below are curated categories based on the supplied
#     OMA evidence. The groups are NOT new OMA algorithm outputs.
#
# WORKING DIRECTORY:
#   Set the working directory manually to R_Scripts before running this file.
#
# OUTPUT DIRECTORY:
#   Results/02_Prepare_Integrin_BindingData_RESULTS/
# ==============================================================================


# ==============================================================================
# 0. CONFIGURATION AND REQUIRED PACKAGE
# ==============================================================================

INPUT_DIR <- "Data"
OUTPUT_DIR <- file.path("Results", "02_Prepare_Integrin_BindingData_RESULTS")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# readxl reads Excel workbooks without changing or overwriting the originals.
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

score_file <- file.path(INPUT_DIR, "BindingScores_StatisticalAnalysis.xlsx")
raw_file   <- file.path(INPUT_DIR, "All_Species_BindingScores.xlsx")
oma_file   <- file.path(INPUT_DIR, "OMA_Subunit_Relationships.xlsx")

for (path in c(score_file, raw_file, oma_file)) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}

# The nine ligand labels and their reporting order are kept consistent across
# all later analyses. 'sT4 ' in the workbook has a trailing space; trimws()
# removes that formatting inconsistency without changing the ligand identity.
LIGANDS <- c("Cilengitide", "T3", "T4", "TRIAC", "TETRAC",
             "sT3", "sT4", "rT3", "T4AM")
SPECIES <- c("H. sapiens", "S. purpuratus", "N. vectensis")

# The source spreadsheets have been independently checked to contain 48
# receptors, 9 ligands, and positive KD measurements throughout.
N_EXPECTED <- c(8L, 24L, 16L)

# Small floating-point differences are allowed when checking log values.
# The provided score workbook and raw KD workbook agree at ~1e-15.
LOG_TOLERANCE <- 1e-8

# Export CSV files consistently without an extra unnamed row-number column.
export_csv <- function(x, filename) {
  utils::write.csv(x, file.path(OUTPUT_DIR, filename), row.names = FALSE,
                   na = "")
}

# Stop immediately if names are missing or duplicated; mismatched names would
# otherwise silently misassign affinities to receptors or ligands.
check_names <- function(actual, expected, description) {
  if (anyNA(actual) || any(!nzchar(actual)) || anyDuplicated(actual) ||
      !setequal(actual, expected) || length(actual) != length(expected)) {
    stop(description, " do not match the expected names.\nMissing: ",
         paste(setdiff(expected, actual), collapse = ", "),
         "\nUnexpected: ", paste(setdiff(actual, expected), collapse = ", "))
  }
}

# Pivot a numeric receptor-by-ligand matrix to a long table. The matrix is
# transposed for vectorization so each receptor's nine values stay together.
matrix_to_long <- function(x, species_by_receptor, value_name) {
  out <- data.frame(
    Receptor = rep(rownames(x), each = ncol(x)),
    Species  = rep(species_by_receptor, each = ncol(x)),
    Ligand   = rep(colnames(x), times = nrow(x)),
    Value    = as.vector(t(x)),
    stringsAsFactors = FALSE
  )
  names(out)[4] <- value_name
  out
}


# ==============================================================================
# 1. IMPORT THE PRECOMPUTED -log10(KD) RECEPTOR MATRIX
# ==============================================================================

# This file contains the FINAL 48 x 9 score matrix used by the earlier
# statistical-analysis scripts. Its values must NOT be logarithmically
# transformed again.
score <- as.data.frame(readxl::read_excel(score_file, sheet = "Sheet1"),
                       stringsAsFactors = FALSE)
names(score) <- trimws(names(score))

if (!"Receptor" %in% names(score)) stop("Missing 'Receptor' column.")
check_names(names(score)[names(score) != "Receptor"], LIGANDS,
            "Score-sheet ligand columns")
if (nrow(score) != sum(N_EXPECTED)) {
  stop("Expected 48 receptors in the score workbook; found ", nrow(score))
}

score$Receptor <- trimws(as.character(score$Receptor))
if (anyNA(score$Receptor) || any(!nzchar(score$Receptor)) ||
    anyDuplicated(score$Receptor)) stop("Invalid or duplicate receptor IDs.")

# The original analysis assumes that this particular worksheet is ordered as
# 8 human, 24 sea-urchin, and 16 sea-anemone receptors. We retain that known
# layout here, and validate the receptor names against the three RAW KD sheets
# below so that an accidental row reordering is detected.
receptor_species <- rep(SPECIES, times = N_EXPECTED)
receptors <- score$Receptor
score_matrix <- as.matrix(score[, LIGANDS, drop = FALSE])
if (!all(vapply(score[, LIGANDS, drop = FALSE], is.numeric, logical(1)))) {
  stop("At least one transformed-score ligand column is not numeric.")
}
rownames(score_matrix) <- receptors
storage.mode(score_matrix) <- "double"
if (any(!is.finite(score_matrix))) stop("Invalid/missing transformed scores.")


# ==============================================================================
# 2. IMPORT ORIGINAL, UNTRANSFORMED KD VALUES (MOLAR UNITS)
# ==============================================================================

# Each raw species worksheet has:
#   Excel row 4: Ligand header + integrin receptor identifiers
#   Excel rows 5-13: Nine ligands and their predicted KD (M)
#
# The headers on one worksheet include an inaccurate descriptive note about
# species, so the sheet name, expected receptor identifiers, and dimensions
# (rather than that note) are used to establish species identity.
read_raw_kd <- function(sheet_name, cell_range, species_name) {
  dat <- as.data.frame(readxl::read_excel(
    raw_file, sheet = sheet_name, range = cell_range,
    .name_repair = "minimal"
  ), stringsAsFactors = FALSE)
  names(dat) <- trimws(names(dat))
  
  if (names(dat)[1] != "Ligand" || nrow(dat) != length(LIGANDS)) {
    stop("Unexpected KD layout in worksheet: ", sheet_name)
  }
  dat$Ligand <- trimws(as.character(dat$Ligand))
  check_names(dat$Ligand, LIGANDS, paste(sheet_name, "ligands"))
  
  raw_receptors <- names(dat)[-1]
  expected_receptors <- receptors[receptor_species == species_name]
  check_names(raw_receptors, expected_receptors,
              paste(sheet_name, "receptors"))
  if (!all(vapply(dat[-1], is.numeric, logical(1)))) {
    stop("Nonnumeric KD column in worksheet: ", sheet_name)
  }
  
  # Original sheet orientation: ligand rows x receptor columns.
  # Transpose to the common receptor rows x ligand columns orientation.
  kd <- t(as.matrix(dat[, -1, drop = FALSE]))
  rownames(kd) <- raw_receptors
  colnames(kd) <- dat$Ligand
  kd <- kd[expected_receptors, LIGANDS, drop = FALSE]
  storage.mode(kd) <- "double"
  if (any(!is.finite(kd)) || any(kd <= 0)) {
    stop("KD values must all be finite and strictly positive in ", sheet_name)
  }
  kd
}

human_kd <- read_raw_kd("H.S_AFReceptors", "A4:I13", "H. sapiens")
urchin_kd <- read_raw_kd("S.P_AFReceptors", "A4:Y13", "S. purpuratus")
anemone_kd <- read_raw_kd("N.V_AFReceptors", "A4:Q13", "N. vectensis")

raw_kd_matrix <- rbind(human_kd, urchin_kd, anemone_kd)
raw_kd_matrix <- raw_kd_matrix[receptors, LIGANDS, drop = FALSE]

# A critical cross-workbook validation: every source raw KD must produce the
# corresponding supplied -log10(KD) score. If this fails, analyses should not
# proceed because the workbooks would disagree about the same docking result.
log_difference <- abs((-log10(raw_kd_matrix)) - score_matrix)
max_log_difference <- max(log_difference)
if (!is.finite(max_log_difference) || max_log_difference > LOG_TOLERANCE) {
  stop("Raw KD and transformed scores disagree. Maximum |difference| = ",
       signif(max_log_difference, 8), ".")
}

# Export BOTH units. The raw values must remain raw for the Figure 18 medians.
score_wide <- data.frame(Receptor = receptors, score_matrix,
                         check.names = FALSE, row.names = NULL)
raw_wide <- data.frame(Receptor = receptors, raw_kd_matrix,
                       check.names = FALSE, row.names = NULL)
export_csv(score_wide, "BindingAffinity_NegLog10KD_Wide.csv")
export_csv(raw_wide, "BindingAffinity_RawKD_Wide.csv")
export_csv(matrix_to_long(score_matrix, receptor_species, "NegLog10KD"),
           "BindingAffinity_NegLog10KD_Long.csv")
export_csv(matrix_to_long(raw_kd_matrix, receptor_species, "KD_M"),
           "BindingAffinity_RawKD_Long.csv")


# ==============================================================================
# 3. MAP HISTORICAL RECEPTOR IDS TO BIOLOGICAL ALPHA AND BETA SUBUNITS
# ==============================================================================

# The spreadsheet uses original computational-model IDs as receptor names.
# These are preserved unchanged for exact matching to the docking results.
# For biological interpretation, the alpha/beta subunits are also recorded
# under the FINAL subunit nomenclature used in the thesis.
#
# Especially important N. vectensis conversions:
#   aH6 -> alpha3   aA9 -> alpha9   aA5 -> alpha5   a71 -> alpha8
#   bW0 -> beta2    bU5 -> betan2  bD5 -> beta1    b23 -> betabeta

human_labels <- data.frame(
  Receptor     = c("avb1", "avb3", "avb5", "avb6", "avb8",
                   "a5b1", "a8b1", "aiibb3"),
  AlphaSubunit = c(rep("alphaV", 5), "alpha5", "alpha8", "alphaIIb"),
  BetaSubunit  = c("beta1", "beta3", "beta5", "beta6", "beta8",
                   "beta1", "beta1", "beta3"),
  stringsAsFactors = FALSE
)

sp_alpha <- c(av3 = "alphaV3", av2 = "alphaV2", av1 = "alphaV1",
              a9 = "alpha9", a8 = "alpha8", a6 = "alpha6",
              aP2 = "alphaP2", aP1 = "alphaP1")
sp_beta <- c(b1a1 = "beta1A1", bg = "betaG", bc = "betaC")
nv_alpha <- c(aH6 = "alpha3", aA9 = "alpha9", aA5 = "alpha5", a71 = "alpha8")
nv_beta <- c(bW0 = "beta2", bU5 = "betan2", bD5 = "beta1", b23 = "betabeta")

identify_subunits <- function(receptor, species) {
  if (species == "H. sapiens") {
    i <- match(receptor, human_labels$Receptor)
    if (is.na(i)) stop("Unknown human receptor: ", receptor)
    return(c(human_labels$AlphaSubunit[i], human_labels$BetaSubunit[i]))
  }
  parts <- strsplit(receptor, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2L) stop("Unexpected receptor ID: ", receptor)
  if (species == "S. purpuratus") {
    alpha <- unname(sp_alpha[parts[1]])
    beta <- unname(sp_beta[tolower(parts[2])])
  } else if (species == "N. vectensis") {
    alpha <- unname(nv_alpha[parts[1]])
    beta <- unname(nv_beta[parts[2]])
  } else stop("Unexpected species: ", species)
  if (length(alpha) != 1L || length(beta) != 1L ||
      is.na(alpha) || is.na(beta)) stop("Unknown subunit code: ", receptor)
  c(alpha, beta)
}

subunits <- t(vapply(seq_along(receptors), function(i) {
  identify_subunits(receptors[i], receptor_species[i])
}, character(2)))

metadata <- data.frame(
  Receptor = receptors,
  Species = receptor_species,
  AlphaSubunit = subunits[, 1],
  BetaSubunit = subunits[, 2],
  stringsAsFactors = FALSE
)


# ==============================================================================
# 4. READ AND ATTACH THE ORIGINAL OMA SUBUNIT EVIDENCE
# ==============================================================================

# OMA_Subunit_Relationships.xlsx is the data source for the human relationships,
# evidence types (Pair, HOG, OMA Group), and interpretation notes.
# Each row describes one NONHUMAN subunit. The human subunits are the
# comparison reference, not additional nonhuman OMA lookup entries.
#
# Unicode identifiers in the OMA workbook use Greek alpha/beta. Convert only
# the lookup keys to readable ASCII ('alpha', 'beta') for matching; retain the
# original relationship and evidence TEXT unmodified in the exported tables.

as_ascii_key <- function(x) {
  x <- gsub("\u03b1", "alpha", x, fixed = TRUE)
  gsub("\u03b2", "beta", x, fixed = TRUE)
}

read_oma_sheet <- function(sheet_name, type) {
  x <- as.data.frame(readxl::read_excel(oma_file, sheet = sheet_name),
                     stringsAsFactors = FALSE)
  if (ncol(x) != 4L) stop("Expected four columns in OMA sheet ", sheet_name)
  names(x) <- c("OriginalOMAName", "HumanRelationship", "OMAEvidence",
                "OMAInterpretation")
  x$OriginalOMAName <- trimws(as.character(x$OriginalOMAName))
  x$OMAKey <- as_ascii_key(x$OriginalOMAName)
  if (anyDuplicated(x$OMAKey) || anyNA(x$OMAKey)) {
    stop("Duplicate/missing subunit in OMA sheet: ", sheet_name)
  }
  x$SubunitType <- type
  x
}

alpha_oma <- read_oma_sheet("\u03b1-subunit relationships", "Alpha")
beta_oma <- read_oma_sheet("\u03b2-subunit relationships", "Beta")

# ------------------------------------------------------------------------------
# 4.1 CURATED ANALYSIS GROUPS BASED ON THE SUPPLIED OMA WORKBOOK
# ------------------------------------------------------------------------------
# These names retain the group definitions from the original Script2A C9
# analysis. They organize receptor comparisons; they are not formal OMA HOG
# identifiers, and subunits with only HOG-level evidence are NOT elevated to
# pairwise orthologs.
#
# 'Not_evaluable' differs from 'No_OMA_supported_group': the former had no
# corresponding OMA/UniProt entry; the latter was evaluated without a useful
# human grouping. Human alphaIIb has a separate special category.

alpha_groups <- c(
  "S. purpuratus alphaV1" = "No_OMA_supported_group",
  "S. purpuratus alphaV2" = "alphaV_5_8_related",
  "S. purpuratus alphaV3" = "No_OMA_supported_group",
  "S. purpuratus alphaP1" = "alpha4_9_related",
  "S. purpuratus alphaP2" = "alphaV_5_8_related",
  "S. purpuratus alpha6"  = "alpha6_7_related",
  "S. purpuratus alpha8"  = "alphaV_5_8_related",
  "S. purpuratus alpha9"  = "alpha4_9_related",
  "N. vectensis alpha3"    = "alpha3_6_7_related",
  "N. vectensis alpha5"    = "Not_evaluable",
  "N. vectensis alpha8"    = "No_OMA_supported_group",
  "N. vectensis alpha9"    = "Not_evaluable"
)
beta_groups <- c(
  "S. purpuratus beta1A1" = "beta1_related",
  "S. purpuratus betaG"   = "beta1_related",
  "S. purpuratus betaC"   = "beta1_related", # HOG-only: weaker evidence
  "N. vectensis betabeta"  = "beta1_related",
  "N. vectensis beta1"     = "beta1_related",
  "N. vectensis beta2"     = "beta1_related",
  "N. vectensis betan2"    = "beta6_related"
)

# Require exact agreement between subunits in the uploaded OMA source and the
# curated group lookup. This catches changes to the workbook or transcription
# errors instead of silently assigning an unreviewed group.
check_names(alpha_oma$OMAKey, names(alpha_groups), "OMA alpha subunits")
check_names(beta_oma$OMAKey, names(beta_groups), "OMA beta subunits")
alpha_oma$AlphaGroup <- unname(alpha_groups[alpha_oma$OMAKey])
beta_oma$BetaGroup   <- unname(beta_groups[beta_oma$OMAKey])

# Export the subunit-level evidence and the explicit group classifications,
# providing a traceable bridge from the OMA workbook to receptor metadata.
export_csv(alpha_oma, "OMA_Alpha_Subunit_Evidence.csv")
export_csv(beta_oma, "OMA_Beta_Subunit_Evidence.csv")

# Build keys for each receptor's two individual subunits. Human receptors
# receive a descriptive 'Human reference subunit' rather than a fabricated
# nonhuman OMA hit.
metadata$AlphaOMAKey <- paste(metadata$Species, metadata$AlphaSubunit)
metadata$BetaOMAKey  <- paste(metadata$Species, metadata$BetaSubunit)
alpha_match <- match(metadata$AlphaOMAKey, alpha_oma$OMAKey)
beta_match  <- match(metadata$BetaOMAKey, beta_oma$OMAKey)
nonhuman <- metadata$Species != "H. sapiens"
if (anyNA(alpha_match[nonhuman]) || anyNA(beta_match[nonhuman])) {
  stop("A nonhuman receptor subunit could not be matched to the OMA workbook.")
}

metadata$AlphaGroup <- unname(alpha_oma$AlphaGroup[alpha_match])
metadata$BetaGroup  <- unname(beta_oma$BetaGroup[beta_match])
metadata$AlphaEvidence <- alpha_oma$OMAEvidence[alpha_match]
metadata$BetaEvidence  <- beta_oma$OMAEvidence[beta_match]
metadata$AlphaHumanRelationship <- alpha_oma$HumanRelationship[alpha_match]
metadata$BetaHumanRelationship  <- beta_oma$HumanRelationship[beta_match]
metadata$AlphaOMAInterpretation <- alpha_oma$OMAInterpretation[alpha_match]
metadata$BetaOMAInterpretation  <- beta_oma$OMAInterpretation[beta_match]

is_human <- metadata$Species == "H. sapiens"
metadata$AlphaGroup[is_human] <- ifelse(
  metadata$AlphaSubunit[is_human] == "alphaIIb",
  "alphaIIb_no_broad_group", "alphaV_5_8_related"
)
metadata$BetaGroup[is_human] <- ifelse(
  metadata$BetaSubunit[is_human] == "beta1", "beta1_related",
  ifelse(metadata$BetaSubunit[is_human] == "beta6", "beta6_related",
         "No_cross_species_group")
)
for (field in c("AlphaEvidence", "BetaEvidence",
                "AlphaHumanRelationship", "BetaHumanRelationship",
                "AlphaOMAInterpretation", "BetaOMAInterpretation")) {
  metadata[[field]][is_human] <- "Human reference subunit"
}

# OMA keys are intermediate join identifiers, not extra biological results.
metadata$AlphaOMAKey <- NULL
metadata$BetaOMAKey <- NULL

# Keep the six original C9 metadata field names for compatibility with
# prior Figure 17/18 scripts, and retain the more detailed source notes.
metadata <- metadata[, c(
  "Receptor", "Species", "AlphaSubunit", "BetaSubunit",
  "AlphaGroup", "BetaGroup", "AlphaEvidence", "BetaEvidence",
  "AlphaHumanRelationship", "BetaHumanRelationship",
  "AlphaOMAInterpretation", "BetaOMAInterpretation"
)]
export_csv(metadata, "Receptor_Metadata_OMA.csv")


# ==============================================================================
# 5. PREPARE THE ORIGINAL C9 OMA-ANNOTATED RECEPTOR TABLE
# ==============================================================================

# Three selected ligand/group columns were used in the original cross-species
# OMA analysis. The two grouped scores are arithmetic means in -log10(KD)
# SPACE for each receptor, not mean raw KD values.
#
# Keeping these familiar column names makes the output usable by scripts
# originally written around C9_OMA_Annotated_Receptor_Binding_Table.csv.

oma_binding <- metadata
oma_binding$TRIAC_TETRAC <- rowMeans(score_matrix[, c("TRIAC", "TETRAC"),
                                                  drop = FALSE])
oma_binding$T3_T4 <- rowMeans(score_matrix[, c("T3", "T4"), drop = FALSE])
oma_binding$T4AM <- score_matrix[, "T4AM"]
export_csv(oma_binding, "C9_OMA_Annotated_Receptor_Binding_Table.csv")


# ==============================================================================
# 6. PREPARE HUMAN CRYSTAL-STRUCTURE VS ALPHAFOLD VALIDATION DATA
# ==============================================================================

# The "H.S_CS_vs_ AF" sheet contains five VERTICALLY STACKED three-column
# blocks (not two side-by-side blocks as assumed by an older analysis script).
# Each block contains one header row and the same nine ligands below it.
# This export lets Script 03 calculate Spearman correlations directly from a
# standardized table while preserving the exact supplied -log10(KD) values.
#
# The five experimental PDB controls are 8W30, 8TCF, 8T2V, 1L5G, and 4WK4.

cs_kd <- as.data.frame(readxl::read_excel(
  raw_file, sheet = "H.S_CSReceptors", range = "A4:F13",
  .name_repair = "minimal"
), stringsAsFactors = FALSE)
names(cs_kd) <- trimws(names(cs_kd))
cs_kd$Ligand <- trimws(as.character(cs_kd$Ligand))
check_names(cs_kd$Ligand, LIGANDS, "Human crystal ligands")

block_starts <- c(1L, 13L, 25L, 37L, 49L)
block_receptor <- c("avb1", "avb8", "aiibb3", "avb3", "a5b1")
block_crystal <- c("avb1_8W30", "avb8_8TCF", "aIIbb3_8T2V",
                   "avb3_1L5G", "a5b1_4WK4")
block_pdb <- c("8W30", "8TCF", "8T2V", "1L5G", "4WK4")

validation_parts <- lapply(seq_along(block_starts), function(i) {
  first <- block_starts[i]
  x <- as.data.frame(readxl::read_excel(
    raw_file, sheet = "H.S_CS_vs_ AF",
    range = paste0("A", first, ":C", first + 9L),
    .name_repair = "minimal"
  ), stringsAsFactors = FALSE)
  names(x) <- trimws(names(x))
  if (nrow(x) != 9L || !identical(names(x),
                                  c("Ligand", block_crystal[i], block_receptor[i]))) {
    stop("Unexpected human crystal/AF block at Excel row ", first)
  }
  x$Ligand <- trimws(as.character(x$Ligand))
  check_names(x$Ligand, LIGANDS, "Human validation ligands")
  if (!is.numeric(x[[2]]) || !is.numeric(x[[3]])) {
    stop("Human validation profiles must be numeric (row ", first, ").")
  }
  
  # Cross-check against independent score and raw CS KD worksheets.
  af_expected <- score_matrix[block_receptor[i], x$Ligand]
  cs_idx <- match(x$Ligand, cs_kd$Ligand)
  cs_expected <- -log10(as.numeric(cs_kd[[block_crystal[i]]][cs_idx]))
  if (any(!is.finite(cs_expected)) ||
      max(abs(as.numeric(x[[2]]) - cs_expected)) > LOG_TOLERANCE ||
      max(abs(as.numeric(x[[3]]) - af_expected)) > LOG_TOLERANCE) {
    stop("Human validation data do not match original source values for ",
         block_receptor[i])
  }
  
  data.frame(
    Receptor = block_receptor[i],
    CrystalReceptor = block_crystal[i],
    Crystal_PDB = block_pdb[i],
    Ligand = x$Ligand,
    Crystal_NegLog10KD = as.numeric(x[[2]]),
    AlphaFold_NegLog10KD = as.numeric(x[[3]]),
    stringsAsFactors = FALSE
  )
})
human_validation <- do.call(rbind, validation_parts)
export_csv(human_validation, "Human_CS_vs_AF_BindingData.csv")


# ==============================================================================
# 7. OUTPUT QUALITY-CONTROL SUMMARY
# ==============================================================================

# Counts are checks on dataset completeness, NOT biological replicate sizes.
# Each row of the affinity matrices represents a computational receptor model.
qc <- data.frame(
  Check = c("Total_receptors", "Human_receptors", "Sea_urchin_receptors",
            "Sea_anemone_receptors", "Ligands", "Affinity_measurements",
            "Human_validation_pairs", "Human_validation_measurements",
            "OMA_alpha_nonhuman_subunits", "OMA_beta_nonhuman_subunits",
            "Max_abs_log_score_discrepancy"),
  Value = c(nrow(score_matrix), N_EXPECTED[1], N_EXPECTED[2], N_EXPECTED[3],
            ncol(score_matrix), length(score_matrix), length(block_starts),
            nrow(human_validation), nrow(alpha_oma), nrow(beta_oma),
            max_log_difference),
  stringsAsFactors = FALSE
)
export_csv(qc, "Data_Preparation_QC.csv")

message("Preparation complete: ", nrow(score_matrix), " receptors x ",
        ncol(score_matrix), " ligands (", length(score_matrix), " values).")
message("Raw KD / -log10(KD) agreement; maximum discrepancy: ",
        format(max_log_difference, scientific = TRUE))
message("Validated OMA evidence for ", nrow(alpha_oma), " alpha and ",
        nrow(beta_oma), " beta nonhuman subunits.")
message("Human CS/AF validation profiles: ", length(block_starts),
        " receptor pairs (", nrow(human_validation), " ligand measurements).")
message("All prepared data saved in: ", OUTPUT_DIR)


# ==============================================================================
# END OF SCRIPT
# ==============================================================================
#
# OUTPUT SUMMARY
#
# All outputs are saved to:
#
# Results/02_Prepare_Integrin_BindingData_RESULTS/
#
# Binding-affinity matrices (-log10 KD):
#   BindingAffinity_NegLog10KD_Wide.csv
#   BindingAffinity_NegLog10KD_Long.csv
#
# Original binding-affinity matrices (raw KD in molar units):
#   BindingAffinity_RawKD_Wide.csv
#   BindingAffinity_RawKD_Long.csv
#
# OMA subunit relationship evidence:
#   OMA_Alpha_Subunit_Evidence.csv
#   OMA_Beta_Subunit_Evidence.csv
#
# Receptor metadata and OMA-annotated binding data:
#   Receptor_Metadata_OMA.csv
#   C9_OMA_Annotated_Receptor_Binding_Table.csv
#
# Human crystal-structure vs AlphaFold validation data (Table 5):
#   Human_CS_vs_AF_BindingData.csv
#
# Data preparation and quality-control summary:
#   Data_Preparation_QC.csv
#
# Total: 10 output files
#
# These prepared datasets are used by Scripts 03-07 for subsequent
# statistical analyses, binding-affinity comparisons, and figure generation.
#
# ==============================================================================


# ==============================================================
# Check_Rejected_SeaUrchin_Sequences.R
#
# PURPOSE:
# Investigate why three S. purpuratus integrin sequences
# were excluded from the final sequence dataset.
#
# Local pairwise sequence alignments are used to:
#   1. Identify the closest retained integrin sequence.
#   2. Determine which portion of that sequence aligns
#      with the excluded sequence.
#   3. Calculate amino-acid sequence identity and coverage.
#   4. Visualize the aligned regions to assess whether
#      the excluded sequences represent partial proteins.
#
# NOTE:
# The retained sequences have already been trimmed to their
# extracellular domains. Therefore, missing sequence coverage
# does not necessarily indicate a missing full-length domain.
# Domain-level conclusions require additional annotation.
#
# INPUT FILES:
#   Data/Alpha_MSA.fasta
#   Data/Beta_MSA.fasta
#   Data/Rejected_SeaUrchinSeqs.txt
#
# OUTPUT DIRECTORY:
#   Results/Check_Rejected_SeaUrchin_Sequences_RESULTS/
#
# WORKING DIRECTORY:
#   Run this script with R_Scripts as the working directory.
# ==============================================================


# ==============================================================
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# ==============================================================

# Biostrings: Read and manipulate amino-acid sequences.
# pwalign: Perform pairwise amino-acid sequence alignments.
# ggplot2: Generate the alignment coverage figure.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("Biostrings", "pwalign")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(Biostrings)
library(pwalign)
library(ggplot2)


# ==============================================================
# 2. DEFINE INPUT AND OUTPUT DIRECTORIES
# ==============================================================

# Input FASTA files are stored in R_Scripts/Data.
data_dir <- "Data"

# All outputs are saved in a dedicated results folder.
out_dir <- file.path(
  "Results",
  "Check_Rejected_SeaUrchin_Sequences_RESULTS"
)

# Create the folder automatically if it does not exist.
dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Define input file paths.
alpha_file <- file.path(data_dir, "Alpha_MSA.fasta")
beta_file <- file.path(data_dir, "Beta_MSA.fasta")

rejected_file <- file.path(
  data_dir,
  "Rejected_SeaUrchinSeqs.txt"
)

# Stop with an informative error if an input is missing.
input_files <- c(alpha_file, beta_file, rejected_file)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\n\nCheck getwd() and your Data folder."
  )
}


# ==============================================================
# 3. READ THE AMINO-ACID SEQUENCES
# ==============================================================

# The alpha and beta FASTA files contain retained sequences.
# They are already multiple-sequence alignments and therefore
# contain gap characters.

alpha <- readAAStringSet(alpha_file)
beta <- readAAStringSet(beta_file)

# The rejected file contains three additional sequences
# excluded from the original dataset.
rejected <- readAAStringSet(
  rejected_file,
  format = "fasta"
)

# Remove existing alignment gaps before new pairwise alignment.
# Keeping old MSA gaps would interfere with local alignment.

remove_gaps <- function(sequences) {
  
  cleaned <- gsub(
    "[-.]",
    "",
    as.character(sequences)
  )
  
  output <- AAStringSet(cleaned)
  
  # Preserve original FASTA identifiers.
  names(output) <- names(sequences)
  
  return(output)
}

alpha <- remove_gaps(alpha)
beta <- remove_gaps(beta)

# Check that the rejected sequences have the expected accessions.
expected_ids <- c(
  "XP_030845946",
  "XP_030845945",
  "XP_030836103"
)

for (id in expected_ids) {
  
  if (!any(grepl(id, names(rejected), fixed = TRUE))) {
    stop("Rejected sequence not found: ", id)
  }
}


# ==============================================================
# 4. IDENTIFY SEQUENCE TYPE AND SPECIES
# ==============================================================

# The FASTA identifiers supplied for the excluded sequences
# distinguish beta (SP_b) and alpha (SP_a) subunits.
#
# Their NCBI accession numbers are also used here so that
# the analysis does not depend entirely on header formatting.

get_sequence_type <- function(sequence_name) {
  
  if (grepl("XP_030836103", sequence_name, fixed = TRUE)) {
    return("alpha")
  }
  
  if (grepl("XP_030845946|XP_030845945", sequence_name)) {
    return("beta")
  }
  
  stop("Unrecognized excluded sequence: ", sequence_name)
}

# Retained reference names may use different abbreviations.
# This function recognizes commonly used sea urchin prefixes.
#
# If your FASTA headers use a different naming convention,
# adjust this function to match them.

is_sea_urchin <- function(sequence_name) {
  
  grepl(
    "^(SP|Spur|S_purpuratus|S\\.purpuratus)[_\\-]",
    sequence_name,
    ignore.case = TRUE
  )
}


# ==============================================================
# 5. PREPARE THE LOCAL ALIGNMENT PROCEDURE
# ==============================================================

# Local alignment identifies the most similar region between
# two sequences without requiring end-to-end matching.
#
# This is useful for investigating truncated proteins, since
# an excluded sequence may correspond to only part of a
# retained integrin extracellular domain.
#
# BLOSUM62 is used as the amino-acid substitution matrix.
# Gap penalties discourage excessive insertions and deletions.

align_sequences <- function(query, reference) {
  
  pwalign::pairwiseAlignment(
    pattern = query,
    subject = reference,
    type = "local",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5
  )
}


# ==============================================================
# 6. COMPARE EACH EXCLUDED SEQUENCE WITH RETAINED SEQUENCES
# ==============================================================

results <- list()
alignment_records <- list()
plot_records <- list()

for (i in seq_along(rejected)) {
  
  query_name <- names(rejected)[i]
  query <- rejected[[i]]
  
  sequence_type <- get_sequence_type(query_name)
  
  # Alpha fragments are compared with alpha references.
  # Beta fragments are compared with beta references.
  references <- if (sequence_type == "alpha") {
    alpha
  } else {
    beta
  }
  
  # Keep only S. purpuratus reference sequences.
  keep <- vapply(
    names(references),
    is_sea_urchin,
    logical(1)
  )
  
  references <- references[keep]
  
  if (length(references) == 0) {
    stop(
      "No retained S. purpuratus ",
      sequence_type,
      " sequences found.\n",
      "Inspect the FASTA headers and update is_sea_urchin()."
    )
  }
  
  # Align the excluded sequence against every retained
  # sea urchin sequence of the corresponding subunit type.
  
  candidate_alignments <- lapply(
    seq_along(references),
    function(j) {
      align_sequences(query, references[[j]])
    }
  )
  
  # Retrieve the local alignment scores.
  alignment_scores <- vapply(
    candidate_alignments,
    function(aln) as.numeric(score(aln)),
    numeric(1)
  )
  
  # Select the best-scoring local alignment.
  best_index <- which.max(alignment_scores)
  
  best_alignment <- candidate_alignments[[best_index]]
  reference <- references[[best_index]]
  
  reference_name <- names(references)[best_index]
  
  # ----------------------------------------------------------
  # Extract aligned coordinates and sequence lengths
  # ----------------------------------------------------------
  
  query_length <- length(query)
  reference_length <- length(reference)
  
  query_start <- start(pattern(best_alignment))
  query_end <- end(pattern(best_alignment))
  
  reference_start <- start(subject(best_alignment))
  reference_end <- end(subject(best_alignment))
  
  aligned_query <- as.character(
    alignedPattern(best_alignment)
  )
  
  aligned_reference <- as.character(
    alignedSubject(best_alignment)
  )
  
  # Count actual aligned residues, excluding gaps.
  aligned_query_residues <- nchar(
    gsub("-", "", aligned_query)
  )
  
  aligned_reference_residues <- nchar(
    gsub("-", "", aligned_reference)
  )
  
  # ----------------------------------------------------------
  # Calculate sequence identity and coverage
  # ----------------------------------------------------------
  
  # PID1 uses aligned residue pairs as the denominator.
  identity_percent <- as.numeric(
    pid(best_alignment, type = "PID1")
  )
  
  # Coverage of the retained extracellular reference.
  reference_coverage <- 100 *
    aligned_reference_residues / reference_length
  
  # Coverage of the excluded sequence.
  query_coverage <- 100 *
    aligned_query_residues / query_length
  
  # Approximate number of residues outside the aligned
  # interval at either end of the retained reference.
  #
  # These are unaligned flanking residues, not independently
  # confirmed missing biological domains.
  
  reference_n_flank <- reference_start - 1
  reference_c_flank <- reference_length - reference_end
  
  # ----------------------------------------------------------
  # Store the summary statistics
  # ----------------------------------------------------------
  
  results[[i]] <- data.frame(
    Excluded_Sequence = query_name,
    Subunit_Type = sequence_type,
    Best_Retained_Match = reference_name,
    Excluded_Length_aa = query_length,
    Retained_ECD_Length_aa = reference_length,
    Query_Start = query_start,
    Query_End = query_end,
    Reference_Start = reference_start,
    Reference_End = reference_end,
    Reference_N_Flank_aa = reference_n_flank,
    Reference_C_Flank_aa = reference_c_flank,
    Identity_Percent = round(identity_percent, 2),
    Reference_Coverage_Percent = round(reference_coverage, 2),
    Query_Coverage_Percent = round(query_coverage, 2),
    Local_Alignment_Score = alignment_scores[best_index]
  )
  
  # Store the original alignment for the detailed report.
  alignment_records[[i]] <- list(
    excluded = query_name,
    retained = reference_name,
    alignment = best_alignment
  )
  
  # Convert alignment coordinates into a percentage of the
  # retained extracellular reference sequence.
  plot_records[[i]] <- data.frame(
    Excluded_Sequence = query_name,
    Reference = reference_name,
    Start_Percent = 100 *
      (reference_start - 1) / reference_length,
    End_Percent = 100 *
      reference_end / reference_length
  )
  
  cat(
    "Processed:", query_name,
    "| Best match:", reference_name, "\n"
  )
}


# ==============================================================
# 7. COMBINE AND EXPORT ALIGNMENT STATISTICS
# ==============================================================

summary_df <- do.call(rbind, results)
coverage_df <- do.call(rbind, plot_records)

# Export a CSV with one row per rejected sequence.
summary_file <- file.path(
  out_dir,
  "Rejected_Sequence_Alignment_Summary.csv"
)

write.csv(
  summary_df,
  summary_file,
  row.names = FALSE
)


# ==============================================================
# 8. EXPORT THE DETAILED ALIGNMENTS
# ==============================================================

# Capture the alignment objects as plain text.
# This provides a record of the residue-level alignments.

alignment_file <- file.path(
  out_dir,
  "Rejected_Sequence_Alignments.txt"
)

capture.output({
  
  cat("REJECTED SEA URCHIN INTEGRIN SEQUENCE ALIGNMENTS\n")
  cat("Alignment method: Local pairwise alignment\n")
  cat("Substitution matrix: BLOSUM62\n")
  cat("Gap opening: 10\n")
  cat("Gap extension: 0.5\n\n")
  
  for (item in alignment_records) {
    
    cat("================================================\n")
    cat("Excluded:", item$excluded, "\n")
    cat("Retained:", item$retained, "\n")
    cat("================================================\n\n")
    
    print(item$alignment)
    
    cat("\n\n")
  }
  
}, file = alignment_file)


# ==============================================================
# 9. GENERATE THE ALIGNMENT COVERAGE FIGURE
# ==============================================================

# Grey bars represent the full retained extracellular
# reference sequences.
#
# Blue bars indicate the regions spanned by the strongest
# local alignment with each rejected sequence.
#
# The figure does not display individual alignment gaps or
# annotated protein domains.
#
# Short labels are used for readability in the thesis figure.
# Full FASTA identifiers are preserved in the CSV and TXT files.

# Match short labels to the excluded accession numbers.
# This avoids depending on the order of the FASTA entries.

short_labels <- c(
  "XP_030845946" = "XP_030845946 vs βG",
  "XP_030845945" = "XP_030845945 vs β1A1",
  "XP_030836103" = "XP_030836103 vs αV2"
)

# Extract accession numbers from the excluded sequence names.
coverage_df$Accession <- sub(
  "^.*?(XP_[0-9]+).*$",
  "\\1",
  coverage_df$Excluded_Sequence
)

# Assign a short figure label to each comparison.
coverage_df$Label <- unname(
  short_labels[coverage_df$Accession]
)

# Stop if any accession was not recognized.
if (any(is.na(coverage_df$Label))) {
  stop("One or more excluded accessions lack figure labels.")
}

# Preserve the original top-to-bottom order in the figure.
coverage_df$Label <- factor(
  coverage_df$Label,
  levels = rev(unique(coverage_df$Label))
)

coverage_plot <- ggplot(
  coverage_df,
  aes(y = Label)
) +
  
  # Complete retained extracellular reference.
  geom_segment(
    aes(
      x = 0,
      xend = 100,
      yend = Label
    ),
    linewidth = 9,
    colour = "grey80",
    lineend = "butt"
  ) +
  
  # Region covered by the local alignment.
  geom_segment(
    aes(
      x = Start_Percent,
      xend = End_Percent,
      yend = Label
    ),
    linewidth = 9,
    colour = "#377CA6",
    lineend = "butt"
  ) +
  
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    title = "Alignment Coverage of Excluded Sea Urchin Integrins",
    subtitle = "Comparison with the closest retained extracellular sequence",
    x = "Position along retained extracellular sequence (%)",
    y = NULL,
    caption = paste(
      "Grey: retained extracellular reference.",
      "Blue: local alignment span."
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.caption = element_text(hjust = 0)
  )


# ==============================================================
# 10. SAVE FIGURES
# ==============================================================

pdf_file <- file.path(
  out_dir,
  "Rejected_Sequence_Coverage.pdf"
)

png_file <- file.path(
  out_dir,
  "Rejected_Sequence_Coverage.png"
)

# PDF is suitable for scalable supplementary figures.
ggsave(
  filename = pdf_file,
  plot = coverage_plot,
  width = 10,
  height = 4.5
)

# PNG is useful for presentations and document previews.
ggsave(
  filename = png_file,
  plot = coverage_plot,
  width = 10,
  height = 4.5,
  dpi = 300
)


# ==============================================================
# 11. PRINT RESULTS TO THE R CONSOLE
# ==============================================================

cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n\n")

print(summary_df)

print(coverage_plot)

cat("\nResults saved to:\n")
cat(normalizePath(out_dir), "\n")

cat("\nFiles generated:\n")

for (f in c(
  summary_file,
  alignment_file,
  pdf_file,
  png_file
)) {
  cat(" - ", basename(f), "\n", sep = "")
}


# ==============================================================
# OUTPUT SUMMARY
# ==============================================================
#
# OUTPUT DIRECTORY:
# R_Scripts/Results/Check_Rejected_SeaUrchin_Sequences_RESULTS/
#
# FILES GENERATED:
#
# 1. Rejected_Sequence_Alignment_Summary.csv
#    - Best retained match for each excluded sequence
#    - Excluded and retained sequence lengths
#    - Query and reference alignment coordinates
#    - Unaligned N- and C-terminal reference regions
#    - Percent amino-acid identity
#    - Query and reference coverage percentages
#    - Local alignment score
#
# 2. Rejected_Sequence_Alignments.txt
#    - Detailed local pairwise sequence alignments
#    - Excluded and retained sequence identifiers
#    - Alignment settings used for reproducibility
#
# 3. Rejected_Sequence_Coverage.pdf
#    - Vector-format alignment coverage figure
#    - Shortened labels for thesis presentation
#    - Suitable for supplementary material
#
# 4. Rejected_Sequence_Coverage.png
#    - High-resolution (300 dpi) coverage figure
#    - Shortened labels for thesis presentation
#
# INTERPRETATION:
# Reduced alignment coverage may indicate partial sequence
# representation or truncation relative to the retained
# extracellular reference sequence.
#
# These results do not independently establish which
# structural domains are missing. Domain annotations and
# full-length protein comparisons are required for that.
#
# ==============================================================

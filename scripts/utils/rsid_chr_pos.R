#!/usr/bin/env Rscript
# =============================================================================
# Script: rsid_chr_pos.R
# Purpose: Map chromosome:position identifiers to rsIDs
# Author: Rasyed
# Date: 2024
# Description:
#   This script maps SNP identifiers from format "chr:position:A1:A2" to
#   standard rsIDs using NCBI BED files. This is necessary because PRS-CS
#   requires rsIDs for matching with reference panels.
# 
# Inputs:
#   $1: File with variants in format "chr:pos:A1:A2" and "chr:pos"
#   $2: File with mapping from "chr:pos" to rsIDs
#
# Outputs:
#   tmp_chr_rsid.txt: Mapping from original variant name to rsID
#   tmp_chr_post.txt: List of variant names to extract
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# =============================================================================
# PARSE COMMAND LINE ARGUMENTS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript rsid_chr_pos.R <variant_file> <rsid_mapping_file>")
}

variant_file <- args[1]
rsid_file <- args[2]

# =============================================================================
# READ DATA
# =============================================================================

# Read variant file (format: chr:pos:A1:A2 and chr:pos)
variants <- fread(variant_file, header = FALSE, col.names = c("full_id", "chr_pos"))

# Read rsID mapping file (format: chr:pos rsID)
rsid_map <- fread(rsid_file, header = FALSE, col.names = c("chr_pos", "rsid"))

cat("Total variants:", nrow(variants), "\n")
cat("Unique chr:pos in rsID map:", length(unique(rsid_map$chr_pos)), "\n")

# =============================================================================
# HANDLE DUPLICATE MAPPINGS
# =============================================================================

# Remove duplicate chr:pos entries in rsID map (keep first occurrence)
rsid_unique <- rsid_map[!duplicated(rsid_map$chr_pos), ]

cat("Unique chr:pos after deduplication:", nrow(rsid_unique), "\n")

# =============================================================================
# MERGE AND CREATE MAPPINGS
# =============================================================================

# Merge variants with rsIDs
merged <- merge(variants, rsid_unique, by = "chr_pos", all.x = TRUE)

# Create mapping file: original variant name -> rsID
# For variants without rsID, keep original name
merged$rsid_final <- ifelse(
    is.na(merged$rsid), 
    merged$full_id, 
    merged$rsid
)

# Save mapping file
mapping_output <- merged[, .(full_id, rsid_final)]
fwrite(mapping_output, "tmp_chr_rsid.txt", col.names = FALSE, sep = "\t")

# Create extraction list (variants that successfully mapped)
extract_list <- merged[!is.na(merged$rsid), .(full_id)]
fwrite(extract_list, "tmp_chr_post.txt", col.names = FALSE)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

n_mapped <- nrow(extract_list)
n_total <- nrow(variants)
mapping_rate <- n_mapped / n_total * 100

cat("\n=== Mapping Summary ===\n")
cat(sprintf("Total variants: %d\n", n_total))
cat(sprintf("Successfully mapped to rsID: %d\n", n_mapped))
cat(sprintf("Mapping rate: %.2f%%\n", mapping_rate))
cat("========================\n")

# Optional: Save unmapped variants for review
if (n_total - n_mapped > 0) {
    unmapped <- merged[is.na(merged$rsid), .(full_id)]
    fwrite(unmapped, "tmp_chr_unmapped.txt", col.names = FALSE)
    cat(sprintf("Unmapped variants saved to: tmp_chr_unmapped.txt\n"))
}

cat("Mapping completed successfully\n")

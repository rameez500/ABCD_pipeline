#!/usr/bin/env Rscript
# =============================================================================
# Script: quality_control.R
# Purpose: Pre-imputation quality control for ABCD genotype data
# Author: Rasyed
# Date: 2024
# Description: 
#   This script performs initial QC on ABCD genotype data before imputation,
#   including missingness analysis, MAF calculation, and ancestry filtering.
# =============================================================================

# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================

library(data.table)      # Fast data manipulation for large files
library(qqman)           # QQ and Manhattan plots
library(Hmisc)           # Data analysis utilities
library(ggplot2)         # Publication-quality plots
library(ggpubr)          # Arranging multiple plots
library(tidyverse)       # Collection of data science packages
library(foreach)         # Parallel processing

# =============================================================================
# 2. SET WORKING DIRECTORY & PATHS
# =============================================================================

# Base directories (consider moving to config file)
BASE_DIR <- "/projects/bga_lab/DATA_REPOSITORIES/ABCD/dataset"
PRE_IMP_DIR <- file.path(BASE_DIR, "pre_imputation")
DATA_DIR <- file.path(BASE_DIR, "NDA/genomics_sample03/genotype_QCed")
OUTPUT_DIR <- file.path(PRE_IMP_DIR, "qc_results")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

setwd(PRE_IMP_DIR)

# =============================================================================
# 3. LOAD AND EXAMINE GENOTYPE DATA
# =============================================================================

#' Load BIM file (variant information)
#' The BIM file contains:
#'   V1: Chromosome
#'   V2: Variant ID (rsID)
#'   V3: Genetic distance (cM)
#'   V4: Base-pair position
#'   V5: Effect allele
#'   V6: Other allele
bim_file <- file.path(DATA_DIR, "ABCD_release_3.0_QCed.bim")
bim <- fread(bim_file)
cat("BIM file dimensions:", dim(bim), "\n")
head(bim)

# NOTE: Check SNP build using dbSNP
# Example: rs3131962
#   - GRCh38: 1:821224
#   - GRCh37: 1:756604
# Confirms data is in GRCh37 build

#' Load FAM file (sample information)
#' The FAM file contains:
#'   V1: Family ID (FID)
#'   V2: Individual ID (IID)
#'   V3: Paternal ID
#'   V4: Maternal ID
#'   V5: Sex
#'   V6: Phenotype
fam_file <- file.path(DATA_DIR, "ABCD_release_3.0_QCed.fam")
fam <- fread(fam_file)
cat("FAM file dimensions:", dim(fam), "\n")
head(fam)

# =============================================================================
# 4. CALCULATE MISSINGNESS STATISTICS
# =============================================================================

#' Note: Missingness calculations should be run via PLINK command line
#' PLINK command:
#' plink --bfile ${DATA_DIR}/ABCD_release_3.0_QCed \
#'   --missing \
#'   --out ${OUTPUT_DIR}/ABCD_release_3.0_QCed_missing

# Load missingness results
miss_geno <- fread(file.path(OUTPUT_DIR, "ABCD_release_3.0_QCed_missing.lmiss"))
miss_sample <- fread(file.path(OUTPUT_DIR, "ABCD_release_3.0_QCed_missing.imiss"))

cat("Variant missingness range:", range(miss_geno$F_MISS), "\n")
cat("Sample missingness range:", range(miss_sample$F_MISS), "\n")

# =============================================================================
# 5. VISUALIZE MISSINGNESS
# =============================================================================

#' Plot histogram of variant missingness
p1 <- ggplot(miss_geno, aes(x = F_MISS)) +
  geom_histogram(color = "black", fill = "white", bins = 50) +
  labs(
    title = "Variant Missingness Distribution",
    x = "Proportion Missing",
    y = "Frequency"
  ) +
  theme_minimal()

#' Plot histogram of sample missingness
p2 <- ggplot(miss_sample, aes(x = F_MISS)) +
  geom_histogram(color = "black", fill = "white", bins = 50) +
  labs(
    title = "Sample Missingness Distribution",
    x = "Proportion Missing",
    y = "Frequency"
  ) +
  theme_minimal()

# Save plots
ggsave(file.path(OUTPUT_DIR, "pre_imp_missing_geno.png"), p1, width = 5, height = 5)
ggsave(file.path(OUTPUT_DIR, "pre_imp_missing_sample.png"), p2, width = 5, height = 5)

# =============================================================================
# 6. CALCULATE ALLELE FREQUENCIES
# =============================================================================

#' PLINK command for allele frequencies:
#' plink --bfile ${DATA_DIR}/ABCD_release_3.0_QCed \
#'   --freq \
#'   --out ${OUTPUT_DIR}/ABCD_release_3.0_QCed_freq

# Load frequency results
freq <- fread(file.path(OUTPUT_DIR, "ABCD_release_3.0_QCed_freq.frq"))
cat("MAF range:", range(freq$MAF), "\n")

# Plot MAF distribution
p3 <- ggplot(freq, aes(x = MAF)) +
  geom_histogram(color = "black", fill = "white", bins = 50) +
  labs(
    title = "Minor Allele Frequency Distribution",
    x = "MAF",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave(file.path(OUTPUT_DIR, "pre_imp_MAF.png"), p3, width = 5, height = 5)

# =============================================================================
# 7. GENETIC ANCESTRY FILTERING
# =============================================================================

#' Load genetic ancestry proportions from ABCD data
#' Data structure: acspsw03.txt contains genetic ancestry information
ancestry_file <- "/projects/bga_lab/DATA_REPOSITORIES/ABCD/ABCD_3.0/acspsw03.txt"
ancestry <- fread(ancestry_file)

# Remove first row (variable descriptions)
ancestry <- ancestry[-1, ]

# Convert ancestry proportions to numeric
ancestry$genetic_af_european <- as.numeric(ancestry$genetic_af_european)

#' Create European ancestry indicator (>= 95% European)
ancestry$european_ancestry <- ifelse(
  ancestry$genetic_af_european >= 0.95, 
  1, 
  0
)

cat("European ancestry (>=95%):", sum(ancestry$european_ancestry == 1, na.rm = TRUE), "\n")
cat("Non-European ancestry:", sum(ancestry$european_ancestry == 0, na.rm = TRUE), "\n")

# Subset European ancestry samples
european_samples <- ancestry[ancestry$european_ancestry == 1, 
                              c("src_subject_id", "genetic_af_european", "european_ancestry")]

# Merge with FAM data to get FID
european_with_fid <- merge(
  european_samples, 
  fam[, .(V1, V2)], 
  by.x = "src_subject_id", 
  by.y = "V2"
)

# Create output files for PLINK filtering
# File 1: Simple keep list (FID, IID)
fwrite(european_with_fid[, .(V1, src_subject_id)], 
       file.path(OUTPUT_DIR, "eur_ancestry.txt"),
       col.names = FALSE, sep = "\t")

# File 2: For post-imputation ID update
european_with_fid$newFID <- paste(european_with_fid$V1, european_with_fid$src_subject_id, sep = "_")
european_with_fid$newID <- 0
fwrite(european_with_fid[, .(newFID, src_subject_id)], 
       file.path(OUTPUT_DIR, "eur_ancestry_imp.txt"),
       col.names = FALSE, sep = "\t")

cat("European ancestry samples saved:", nrow(european_with_fid), "\n")

# =============================================================================
# 8. OTHER ANCESTRY GROUPS (for multi-ancestry analysis)
# =============================================================================

#' Create indicators for other ancestry groups
#' This is useful for sensitivity analyses

# African ancestry
ancestry$genetic_af_african <- as.numeric(ancestry$genetic_af_african)
ancestry$african_ancestry <- ifelse(ancestry$genetic_af_african >= 0.95, 1, 0)
cat("African ancestry (>=95%):", sum(ancestry$african_ancestry == 1, na.rm = TRUE), "\n")

# East Asian ancestry
ancestry$genetic_af_east_asian <- as.numeric(ancestry$genetic_af_east_asian)
ancestry$easian_ancestry <- ifelse(ancestry$genetic_af_east_asian >= 0.95, 1, 0)
cat("East Asian ancestry (>=95%):", sum(ancestry$easian_ancestry == 1, na.rm = TRUE), "\n")

# American ancestry
ancestry$genetic_af_american <- as.numeric(ancestry$genetic_af_american)
ancestry$american_ancestry <- ifelse(ancestry$genetic_af_american >= 0.95, 1, 0)
cat("American ancestry (>=95%):", sum(ancestry$american_ancestry == 1, na.rm = TRUE), "\n")

# =============================================================================
# 9. SESSION INFORMATION
# =============================================================================

sessionInfo()
cat("Pipeline completed successfully at:", Sys.time(), "\n")

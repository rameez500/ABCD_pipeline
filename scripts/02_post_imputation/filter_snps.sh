#!/bin/bash
# =============================================================================
# Script: filter_snps.sh
# Purpose: Filter imputed SNPs based on MAF and imputation quality (Rsq)
# Author: Rasyed
# Date: 2024
# Description:
#   This script filters imputed VCF files by:
#   1. Minor Allele Frequency (MAF) > 0.01
#   2. Imputation R-squared (Rsq) > 0.8
# =============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths - UPDATE THESE FOR YOUR SYSTEM
IMPUTED_DIR="/projects/bga_lab/DATA_REPOSITORIES/ABCD/dataset/NDA/genomics_sample03/impute/imputed"
OUTPUT_DIR="/projects/bga_lab/DATA_REPOSITORIES/ABCD/dataset/post_imputation/filtered"
LOG_DIR="${OUTPUT_DIR}/logs"

# Create directories
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# Filtering thresholds
MAF_THRESHOLD=0.01
RSQ_THRESHOLD=0.8

# Chromosomes to process
CHROMOSOMES=$(seq 1 22)

# Log file
LOG_FILE="${LOG_DIR}/filter_snps_$(date +%Y%m%d_%H%M%S).log"

# =============================================================================
# FUNCTIONS
# =============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

log_message "Starting SNP filtering pipeline"
log_message "MAF threshold: ${MAF_THRESHOLD}"
log_message "Rsq threshold: ${RSQ_THRESHOLD}"

# -----------------------------------------------------------------------------
# Step 1: Extract SNPs meeting quality thresholds
# -----------------------------------------------------------------------------
# This step filters SNPs from each chromosome's imputation info file
# The awk script:
#   - Finds column indices for SNP, MAF, and Rsq from header
#   - Filters rows where MAF > threshold AND Rsq > threshold
#   - Outputs only the SNP names

log_message "Extracting high-quality SNPs from each chromosome"

for chr in ${CHROMOSOMES}; do
    log_message "Processing chromosome ${chr}"
    
    zcat "${IMPUTED_DIR}/chr${chr}.info.gz" | \
    awk -v maf_thresh="${MAF_THRESHOLD}" -v rsq_thresh="${RSQ_THRESHOLD}" '
    NR==1 {
        # Find column indices
        for (i=1; i<=NF; i++) {
            if ($i == "SNP") snp_col = i
            if ($i == "MAF") maf_col = i
            if ($i == "Rsq") rsq_col = i
        }
    }
    NR>1 && $maf_col > maf_thresh && $rsq_col > rsq_thresh {
        print $snp_col
    }' > "${OUTPUT_DIR}/info_filter_chr${chr}.txt"
    
    # Count filtered SNPs
    n_snps=$(wc -l < "${OUTPUT_DIR}/info_filter_chr${chr}.txt")
    log_message "  Chromosome ${chr}: ${n_snps} SNPs passed filters"
done

# -----------------------------------------------------------------------------
# Step 2: Combine all filtered SNPs
# -----------------------------------------------------------------------------
log_message "Combining filtered SNPs from all chromosomes"

cat ${OUTPUT_DIR}/info_filter_chr*.txt > "${OUTPUT_DIR}/maf_rsq_filter_all_chrs.txt"

total_snps=$(wc -l < "${OUTPUT_DIR}/maf_rsq_filter_all_chrs.txt")
log_message "Total SNPs passing filters: ${total_snps}"

# -----------------------------------------------------------------------------
# Step 3: Create PLINK BIM files for each chromosome
# -----------------------------------------------------------------------------
# This step extracts variant information from VCF files without loading
# all genotype data, which saves memory
# The --make-just-bim flag only creates the BIM file, not BED/BIM/FAM

log_message "Creating BIM files for each chromosome"

for chr in ${CHROMOSOMES}; do
    log_message "Processing chromosome ${chr} BIM file"
    
    plink --vcf "${IMPUTED_DIR}/chr${chr}.dose.vcf.gz" \
          --make-just-bim \
          --out "${OUTPUT_DIR}/chr_var_${chr}" \
          --extract "${OUTPUT_DIR}/maf_rsq_filter_all_chrs.txt" \
          --const-fid 0 2>&1 | tee -a "${LOG_FILE}"
done

log_message "SNP filtering pipeline completed successfully"

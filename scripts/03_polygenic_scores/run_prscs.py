#!/usr/bin/env python3
"""
Script: run_prscs.py
Purpose: Run PRS-CS (Polygenic Risk Score - Continuous Shrinkage) analysis
Author: Rasyed
Date: 2024
Description:
    This script runs PRS-CS for multiple GWAS summary statistics:
    1. Autism Spectrum Disorder (ASD)
    2. Substance Use (SU) psychopathology
    3. Non-SU Internalizing
    4. Non-SU Externalizing
    
Usage:
    python run_prscs.py --phenotype asd --phi 1e-4
"""

import os
import subprocess
import argparse
import logging
from datetime import datetime
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
PRSCS_PATH = "/home/rasyed2/Tool/PRScs/PRScs.py"
LD_REF_DIR = "/scratch/silo1/BGA_LAB/dbGaP/GTP/PRScsx/LD_ref/ldblk_1kg_eur"
BASE_DIR = "/projects/bga_lab/DATA_REPOSITORIES/ABCD/dataset/post_imputation"
PLINK_PREFIX = "ABCD_imp_hwe_chr_all"

# GWAS summary statistics files
GWAS_FILES = {
    "asd": {
        "path": "/projects/bga_lab/DATA_REPOSITORIES/ABCD/dataset/post_imputation/aut_spect_dis.txt",
        "n_gwas": 46351,
        "description": "Autism Spectrum Disorder"
    },
    "su_psych": {
        "path": "/scratch/silo1/PolySubUse/Sept_gsemGWAS/results/FS_sumstats_PGS-CS.txt",
        "n_gwas": 1734340,
        "description": "Substance Use Psychopathology"
    },
    "non_su_int": {
        "path": "/scratch/silo1/PolySubUse/Sept_gsemGWAS/results/FNSINT_sumstats_PGS-CS.txt",
        "n_gwas": 1164731,
        "description": "Non-SU Internalizing"
    },
    "non_su_ext": {
        "path": "/scratch/silo1/PolySubUse/Sept_gsemGWAS/results/FNSEXT_sumstats_PGS-CS.txt",
        "n_gwas": 730198,
        "description": "Non-SU Externalizing"
    }
}

# PRS-CS parameters
PHI_VALUES = [1e-4, 1e-2, 1]  # Shrinkage parameters to test
N_ITER = 10000
N_BURNIN = 5000

# =============================================================================
# SETUP LOGGING
# =============================================================================

def setup_logging(phenotype, phi):
    """Setup logging configuration"""
    log_dir = Path(BASE_DIR) / "prs" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"prscs_{phenotype}_phi{phi}_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

# =============================================================================
# PRS-CS FUNCTIONS
# =============================================================================

def run_prscs(phenotype, phi, bim_prefix, output_dir):
    """
    Run PRS-CS for a specific phenotype and phi value
    
    Args:
        phenotype: str, phenotype identifier (e.g., 'asd', 'su_psych')
        phi: float, shrinkage parameter
        bim_prefix: str, prefix for PLINK BIM file
        output_dir: str, output directory
    
    Returns:
        bool: True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    # Get GWAS file info
    gwas_info = GWAS_FILES[phenotype]
    sst_file = gwas_info["path"]
    n_gwas = gwas_info["n_gwas"]
    
    # Create output directory
    out_dir = Path(output_dir) / f"{phenotype}_phi{phi}"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Build PRS-CS command
    cmd = [
        "python", PRSCS_PATH,
        f"--ref_dir={LD_REF_DIR}",
        f"--bim_prefix={bim_prefix}",
        f"--sst_file={sst_file}",
        f"--n_gwas={n_gwas}",
        f"--phi={phi}",
        f"--n_iter={N_ITER}",
        f"--n_burnin={N_BURNIN}",
        f"--out_dir={out_dir}"
    ]
    
    logger.info(f"Running PRS-CS for {phenotype} (phi={phi})")
    logger.info(f"GWAS sample size: {n_gwas}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    try:
        # Run PRS-CS
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"PRS-CS completed successfully for {phenotype} (phi={phi})")
        logger.debug(f"STDOUT: {result.stdout}")
        
        # Combine chromosome files
        combine_chromosomes(out_dir, phenotype, phi)
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"PRS-CS failed for {phenotype} (phi={phi})")
        logger.error(f"STDERR: {e.stderr}")
        return False

def combine_chromosomes(out_dir, phenotype, phi):
    """
    Combine PRS-CS output from all chromosomes into a single file
    
    Args:
        out_dir: Path, output directory containing chromosome files
        phenotype: str, phenotype identifier
        phi: float, shrinkage parameter
    """
    logger = logging.getLogger(__name__)
    
    # Pattern for chromosome files
    pattern = f"{out_dir}/{phenotype}_phi{phi}_pst_eff_a1_b0.5_chr*.txt"
    combined_file = out_dir / f"{phenotype}_phi{phi}_pst_eff_all_chrs.txt"
    
    # Use shell command to concatenate files
    cmd = f"cat {pattern} > {combined_file}"
    
    logger.info(f"Combining chromosome files for {phenotype} (phi={phi})")
    subprocess.run(cmd, shell=True, check=True)
    
    # Count total variants
    n_variants = sum(1 for _ in open(combined_file))
    logger.info(f"Total variants in combined file: {n_variants}")

# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(description="Run PRS-CS analysis")
    parser.add_argument(
        "--phenotype", 
        type=str, 
        required=True,
        choices=["asd", "su_psych", "non_su_int", "non_su_ext", "all"],
        help="Phenotype to analyze (or 'all' for all phenotypes)"
    )
    parser.add_argument(
        "--phi", 
        type=float, 
        default=None,
        help="Shrinkage parameter (if not provided, runs all default values)"
    )
    parser.add_argument(
        "--bim_prefix",
        type=str,
        default=PLINK_PREFIX,
        help=f"PLINK BIM file prefix (default: {PLINK_PREFIX})"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=os.path.join(BASE_DIR, "prs"),
        help="Output directory"
    )
    
    args = parser.parse_args()
    
    # Determine which phenotypes to run
    if args.phenotype == "all":
        phenotypes = list(GWAS_FILES.keys())
    else:
        phenotypes = [args.phenotype]
    
    # Determine phi values to run
    if args.phi is not None:
        phi_values = [args.phi]
    else:
        phi_values = PHI_VALUES
    
    # Run PRS-CS for each combination
    for phenotype in phenotypes:
        for phi in phi_values:
            logger = setup_logging(phenotype, phi)
            logger.info("=" * 60)
            logger.info(f"Starting PRS-CS for: {GWAS_FILES[phenotype]['description']}")
            logger.info("=" * 60)
            
            success = run_prscs(
                phenotype=phenotype,
                phi=phi,
                bim_prefix=args.bim_prefix,
                output_dir=args.output_dir
            )
            
            if success:
                logger.info(f"✅ Successfully completed {phenotype} (phi={phi})")
            else:
                logger.error(f"❌ Failed to complete {phenotype} (phi={phi})")
    
    logger.info("PRS-CS pipeline completed")

if __name__ == "__main__":
    main()

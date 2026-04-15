# ABCD Genotype Pipeline

A comprehensive pipeline for processing ABCD (Adolescent Brain Cognitive Development) study genotype data, from pre-imputation QC to polygenic risk score (PRS) calculation.

## Overview

This pipeline processes genotype data from the ABCD study through the following stages:

1. **Pre-imputation Quality Control**: Missingness analysis, MAF calculation, ancestry filtering
2. **Post-imputation Processing**: SNP filtering (MAF > 1%, Rsq > 0.8), rsID mapping
3. **Polygenic Risk Scores**: PRS-CS analysis for multiple phenotypes

## Repository Structure
ABCD_genotype_pipeline/
├── scripts/ # Analysis scripts
├── config/ # Configuration files
├── docs/ # Documentation
├── notebooks/ # Jupyter notebooks
└── tests/ # Test scripts


## Requirements

### Software
- R (>= 4.0) with packages: data.table, ggplot2, tidyverse
- Python (>= 3.8) with packages: numpy, pandas
- PLINK (>= 1.9)
- BCFtools
- PRS-CS

### Data Access
- ABCD study data (requires data use agreement)
- GWAS summary statistics (PGC, etc.)

## Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/ABCD_genotype_pipeline.git
cd ABCD_genotype_pipeline

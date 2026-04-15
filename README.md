# Prenatal Cannabis Exposure, Genetic Predispositions, and Autism Traits in the ABCD Study

[![Paper](https://europepmc.org/article/ppr/ppr957754)](https://doi.org/10.31234/osf.io/wpng5)
[![R](https://img.shields.io/badge/Made%20with-R-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![PLINK](https://img.shields.io/badge/Tool-PLINK-orange)](https://www.cog-genomics.org/plink/)

This repository contains the complete analysis pipeline for the study:

> **Prenatal Cannabis Exposure, Genetic Predispositions, and Autism Spectrum Disorder Traits in the Adolescent Brain Cognitive Development Study**
>
> This study investigated whether prenatal cannabis exposure (PCE) predicts Autism Spectrum Disorder (ASD) traits above and beyond genetic risk (using a polygenic score, PGS) and familial confounders in children from the ABCD study.

## ABCD Genotyping & Polygenic Score Pipeline

This repository contains scripts for quality control, post-imputation processing, and polygenic score (PGS) calculation using data from the [Adolescent Brain Cognitive Development (ABCD) Study](https://abcdstudy.org/), Release 3.0.

---

## Overview

The pipeline processes ABCD genotype data through four main stages:

```
Raw QCed genotypes (PLINK BED)
        │
        ▼
1. Pre-imputation QC       ── missingness, MAF inspection, European ancestry filter
        │
        ▼
2. Post-imputation QC      ── MAF > 1%, Rsq > 0.8, HWE, rsID mapping
        │
        ▼
3. Ancestry subsets        ── European (95% threshold + FlashPCA), mixed-ancestry
        │
        ▼
4. Polygenic scores (PGS)  ── PRS-cs for ASD, SU-psychopathology, internalizing, externalizing
```

---

## Repository Structure

```
ABCD_pipeline/
│
├── README.md
│
├── scripts/                           # All pipeline scripts
│   ├── 01_pre_imputation/
│   │   ├── quality_control.R         # Pre-imputation QC
│   │   └── ancestry_filtering.R      # Genetic ancestry selection
│   │
│   ├── 02_post_imputation/
│   │   ├── filter_snps.sh            # MAF and Rsq filtering
│   │   ├── convert_to_plink.sh       # VCF to PLINK conversion
│   │   ├── map_rsids.R               # Map to rsIDs
│   │   └── merge_chromosomes.sh      # Merge chromosome files
│   │
│   ├── 03_polygenic_scores/
│   │   ├── run_prscs.py              # PRS-CS analysis
│   │   ├── calculate_prs.sh          # Calculate PRS scores
│   │   └── prs_config.yaml           # PRS-CS configuration
│   │
│   └── utils/
│       ├── rsid_chr_pos.R            # Helper function for rsID mapping
│       └── file_utils.sh             # Common file operations                   # Helper: map chr:pos to rsIDs (R)
```

---

## Scripts

| # | Folder | Language | Description |
|---|------|----------|-------------|
| 01 | `01_pre_imputation` | R | Load BIM/FAM files, compute missingness and MAF, filter European ancestry (≥95%) |
| 02 | `02_post_imputation` | Bash | Filter by MAF > 1% & Rsq > 0.8, apply HWE, map rsIDs, merge chromosomes |
| 03 | `03_polygenic_scores` | Bash | Run PRS-cs and PLINK scoring for four traits (European + mixed ancestry) |

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| R | ≥ 4.0 | Pre-imputation QC, ancestry filtering, PCA cleanup |
| PLINK | 1.9 | Genotype file manipulation, scoring |
| bcftools | ≥ 1.9 | VCF subsetting |
| FlashPCA | 2.0 | Fast principal components analysis |
| PRS-cs | [GitHub](https://github.com/getian107/PRScs) | Bayesian polygenic scoring |
| Python | ≥ 3.7 | Required by PRS-cs |

**R packages:** `data.table`, `ggplot2`, `ggpubr`, `tidyverse`, `Hmisc`, `qqman`, `foreach`, `stringr`

---

## Data Access

> ⚠️ **ABCD Study data is controlled-access.** Raw genotype files are NOT included in this repository.
> Access requires a data use agreement via the [NDA portal](https://nda.nih.gov/).

GWAS summary statistics used for PGS:
- **ASD:** PGC (Grove et al. 2019) — available at [pgc.unc.edu](https://pgc.unc.edu/for-researchers/download-results/)
- **SU-psychopathology / Internalizing / Externalizing:** Brick et al. Genomic SEM paper

---

## Pipeline Notes

- **Genome build:** GRCh37 for raw genotypes; GRCh38 for imputed VCFs (confirmed via dbSNP lookups)
- **Imputation:** Pre-imputed dataset sourced from NDA (Michigan Imputation Server)
- **rsID mapping:** NCBI BED files, build b151 GRCh38p7
- **PRS-cs reference panel:** 1000 Genomes Phase 3, European LD (`ldblk_1kg_eur`)
- **HWE filter:** Applied post-imputation within European ancestry subset (p < 1e-6)

---



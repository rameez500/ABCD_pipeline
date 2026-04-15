# Prenatal Cannabis Exposure, Genetic Predispositions, and Autism Traits in the ABCD Study

[![Paper](https://img.shields.io/badge/Paper-10.31234%2Fosf.io%2Fwpng5-blue)](https://doi.org/10.31234/osf.io/wpng5)
[![R](https://img.shields.io/badge/Made%20with-R-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![PLINK](https://img.shields.io/badge/Tool-PLINK-orange)](https://www.cog-genomics.org/plink/)

This repository contains the complete analysis pipeline for the study:

> **Prenatal Cannabis Exposure, Genetic Predispositions, and Autism Spectrum Disorder Traits in the Adolescent Brain Cognitive Development Study**
>
> This study investigated whether prenatal cannabis exposure (PCE) predicts Autism Spectrum Disorder (ASD) traits above and beyond genetic risk (using a polygenic score, PGS) and familial confounders in children from the ABCD study.

## 📚 Overview & Key Findings

This research used data from the Adolescent Brain Cognitive Development (ABCD) Study® to answer three main questions:

1.  **Genetic Prediction:** Does a polygenic score for ASD predict autism-related traits in a general population sample of 9-11 year olds?
2.  **Independent Effect of PCE:** Does prenatal cannabis exposure predict ASD traits even after accounting for the child's own genetic risk?
3.  **Role of Familial Confounders:** Does the association between PCE and ASD traits persist after controlling for broader family-level factors using a propensity score?

### Main Results from the Paper

Our analyses demonstrated that:

- **PGS is a significant predictor:** The ASD polygenic score was significantly associated with the total SRS score (β = 0.09, 95% CI [0.05, 0.13]) and most subscales.
- **PCE has an independent effect:** Prenatal cannabis exposure accounted for significant variance in ASD traits *above and beyond* the genetic risk captured by the PGS.
- **Effect remains after rigorous control:** Even in the most stringent models controlling for both the PGS *and* a propensity score (accounting for numerous family-level confounders), PCE remained positively associated with the total SRS score and several specific ASD traits.

## 🧬 Analysis Pipeline

The analysis follows a structured pipeline from raw genotype data to statistical modeling, as illustrated below:

```mermaid
graph TD
    subgraph "1. Genotype Processing (PLINK / R)"
        A[Raw ABCD Genotype Data] --> B{Pre-imputation QC};
        B --> C[Filter SNPs: MAF > 1%, Rsq > 0.8];
        C --> D[Restrict to European Ancestry];
        D --> E[Map chr:pos to rsIDs];
    end

    subgraph "2. Polygenic Score (PRS-CS / Python)"
        E --> F[ASD GWAS Summary Stats];
        F --> G[Run PRS-CS to derive posterior effect sizes];
        G --> H[Calculate ASD Polygenic Score (PGS) for each child];
    end

    subgraph "3. Statistical Modeling (R)"
        I[ABCD Phenotypic Data<br>(PCE, SRS, Covariates)] --> J[Merge with Genetic Data];
        H --> J;
        J --> K[Run Hierarchical Regression Models];
        K --> L{Key Outputs};
    end

    L --> M[Table 1: PGS association with SRS];
    L --> N[Table 2: PCE effect above PGS];
    L --> O[Table 3: PCE effect controlling for PGS + Propensity Score];

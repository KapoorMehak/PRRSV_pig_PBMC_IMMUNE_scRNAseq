# Distinct Cell-Type-Specific Gene Expression Signatures Predict Persistent Viral Infection in Pigs

---
## Overview

This repository contains the analysis code and scripts used to reproduce all figures in:

We profiled peripheral blood mononuclear cells (PBMCs) from pigs challenged with Porcine Reproductive and Respiratory Syndrome Virus type 2 (PRRSV-2) using single-cell RNA sequencing (scRNA-seq). By comparing immune cell transcriptomes across three infection outcomes — mock-infected (MI), virus extinct (VE), and persistently infected (PI) — at two time points (14 and 84 days post-infection, DPI), we identify cell-type-specific gene expression signatures that distinguish viral clearance from persistence, with a focus on monocyte-driven immune programs.

---

## Experimental Design

- **Cohort:** 36 specific-pathogen-free (SPF) 4-week-old piglets; 30 challenged with PRRSV-2, 6 mock-infected (MI)
- **Outcome classification:** RT-PCR on PBMCs and lymphoid tissues at 84 DPI stratified animals into persistently infected (PI) or virus extinct (VE)
- **scRNA-seq samples:**
  - 14 DPI: 9 samples (2 MI, 2 VE, 5 PI)
  - 84 DPI: 16 samples (4 MI, 2 VE, 10 PI)
- **Library preparation:** 10X Genomics Chromium single-cell 5′ platform

---

## Key Findings

- scRNA-seq resolved **10 major immune cell populations** in pig PBMCs (14 DPI: 30 clusters; 84 DPI: 24 clusters)
- **Monocytes** show divergent transcriptional programs between VE and PI animals at both time points: VE monocytes transition from early inflammatory activation (14 DPI) to resolution (84 DPI); PI monocytes shift from metabolic activation (14 DPI) to sustained immune engagement (84 DPI)
- Cross-timepoint analysis of **monocyte DEGs** identifies concordant and discordant gene sets that serve as candidate predictors of infection outcome
- **Flow cytometry validation** confirms single-cell-derived immune population proportions (5 major cell types; Spearman correlation, FDR < 0.05) across both time points

---

## Repository Structure

```
PRRSV_pig_PBMC_IMMUNE_scRNAseq/scripts/
├── 01_CellRanger/               # Read alignment and gene-barcode matrix generation
├── 02_ambient_RNA_removal/      # Ambient RNA decontamination (e.g., SoupX/DecontX)
├── 03_QC_filtering/             # Per-sample QC metrics and cell filtering
├── 04_doublet_detection/        # Doublet identification and removal
├── 05_downstream_analysis/      # Integration, clustering, UMAP, cell type annotation
│                                #   → Figures 2 & 3
├── 06_GO_enrichment/            # Gene Ontology enrichment (BP terms)
│                                #   → Figures 4 & 5
├── 07_concordant_discordant_analysis/  # Cross-timepoint DEG comparison
│                                #   → Figure 5
├── 08_unique_atD14_analysis/    # Unique DEGs and GO terms at 14 DPI
│                                #   → Figure 5e
├── 09_pseudobulk_analysis/      # Pseudobulk differential expression and gene detection
└── 10_Flow_cytometry/           # Flow cytometry vs. scRNA-seq proportion correlation
│                                #   → Figure 2,3,6
└── README.md
```
---

## Figure Guide

| Figure | Description | Folder |
|--------|-------------|--------|
| Fig. 1 | PRRSV challenge experimental design | *(schematic BioRender — no code)* |
| Fig. 2 | UMAP, cell type annotation, and proportion analysis at 14 DPI | `05_downstream_analysis/` → `06_downstream_analysis/`|
| Fig. 3 | UMAP, cell type annotation, and proportion analysis at 84 DPI | `05_downstream_analysis/` |
| Fig. 4 | Monocyte volcano plots and GO enrichment (VE vs PI, 14 & 84 DPI) | `11_pseudobulk_analysis/` → `08_GO_enrichment/` |
| Fig. 5 | Cross-timepoint monocyte DEG concordance, discordant sets, unique 14 DPI terms | `09_concordant_discordant_analysis/` → `10_unique_atD14_analysis/` |
| Fig. 6 | Flow cytometry validation of scRNA-seq immune cell proportions | `12_Flow_cytometry/` |

---



## Data Availability

Raw sequencing data (FASTQ files) and processed count matrices are deposited at NCBI Gene Expression Omnibus 

**Accession:** 
**URL:** 

---

# Distinct Cell-Type-Specific Gene Expression Programs associated with persistence versus clearance of Viral Infection in Pigs

---
## Overview

DOI pending: https://doi.org/10.5281/zenodo.20130635
This repository contains the analysis code and scripts used to reproduce all figures in:

We profiled peripheral blood mononuclear cells (PBMCs) from pigs challenged with Porcine Reproductive and Respiratory Syndrome Virus type 2 (PRRSV-2) using single-cell RNA sequencing (scRNA-seq). By comparing immune cell transcriptomes across three infection outcomes: mock-infected (MI), virus extinct (VE), and persistently infected (PI) at two time points (14 and 84 days post-infection, DPI). We identify cell-type-specific gene expression signatures that distinguish viral clearance from persistence, with a focus on monocyte-driven immune programs.

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
- scRNA-seq identifies monocytes as the primary cell type distinguishing PRRSV infection outcomes, with divergence driven by cell-type-specific transcriptional programs rather than differences in immune cell composition
- **Monocytes** show divergent transcriptional programs between VE and PI animals at both time points: VE monocytes transition from early inflammatory activation (14 DPI) to resolution (84 DPI); PI monocytes shift from metabolic activation (14 DPI) to sustained immune engagement (84 DPI)
- Cross-timepoint analysis of **monocyte DEGs** identifies concordant and discordant gene sets. Temporal analysis reveals that when and how immune responses are mounted which determines infection outcome rather than their magnitude 
- **Flow cytometry validation** Flow cytometry significantly correlates with scRNA-seq-derived immune cell proportions across major myeloid and lymphocyte populations at both 14 and 84 DPI, supporting the robustness of single-cell findings

---

## Repository Structure

```
PRRSV_pig_PBMC_IMMUNE_scRNAseq/
├── 01_CellRanger/
│   └── cellranger_script2.sh                       # Read alignment; gene-barcode matrix generation
├── 02_ambient_RNA_removal/
│   └── cellbender_10k_python.py                    # Ambient RNA decontamination (CellBender)
├── 03_QC_filtering/
│   ├── QC_after_cb_14dpi.ipynb                     # QC metrics and cell filtering — 14 DPI
│   └── QC_after_cb_84dpi.ipynb                     # QC metrics and cell filtering — 84 DPI
├── 04_doublet_detection/
│   └── doublet_detection_14dpi.ipynb               # Doublet identification and removal
├── 05_downstream_analysis/
│   ├── downstream_analysis_14dpi.ipynb             # Normalization, clustering, UMAP, annotation — 14 DPI  → Fig. 2
│   └── downstream_analysis_84dpi.ipynb             # Normalization, clustering, UMAP, annotation — 84 DPI  → Fig. 3
├── 06_proportion_analysis/
│   └── cell_prop_UMAP_gen.R                        # Cell type proportion plots and UMAP panels            → Figs. 2 & 3
├── 07_DEG_analysis/
│   ├── DEG_MAST_14DPI.R                            # MAST differential expression at 14 DPI
│   ├── DEG_MAST_84DPI.R                            # MAST differential expression at 84 DPI
│   └── DEG_MAST_84DPI_ABT_subset.R                 # MAST DEG — αβ T cell subset at 84 DPI
├── 08_GO_enrichment/
│   └── DEGs_GO_manuscript_mono.R                   # GO enrichment (BP/MF/CC) for monocyte DEGs           → Figs. 4 & 5
├── 09_concordant_discordant_analysis/
│   └── DE_summary_allcelltypes_allcontrast.ipynb   # Cross-timepoint DEG concordance/discordance          → Fig. 5
├── 10_unique_atD14_analysis/
│   └── unique_DEGs_GO.ipynb                        # Unique DEGs and GO terms at 14 DPI                   → Fig. 5
├── 11_Flow_cytometry/
│   └── check_flow_prop_script.R                    # Flow vs scRNA-seq proportion correlation             → Fig. 6
├── 12_pseudobulk_analysis/
│   ├── Pseudobulk_profiles.R                       # Pseudobulk expression profiles
│   └── script_for_all_pb.R                         # Pseudobulk DE across all cell types                  → Fig. 4
└── README.md
```
---

## Figure Guide

| Figure | Description | Scripts |
|--------|-------------|---------|
| Fig. 1 | PRRSV challenge experimental design | *(BioRender schematic — no code)* |
| Fig. 2 | UMAP, cell type annotation, and proportion analysis at 14 DPI | `05_downstream_analysis/downstream_analysis_14dpi.ipynb` → `06_proportion_analysis/cell_prop_UMAP_gen.R` |
| Fig. 3 | UMAP, cell type annotation, and cell type proportion analysis at 84 DPI | `05_downstream_analysis/downstream_analysis_84dpi.ipynb` → `06_proportion_analysis/cell_prop_UMAP_gen.R` |
| Fig. 4 | Monocyte volcano plots and GO enrichment (VE vs PI, 14 & 84 DPI) | `07_DEG_analysis/` → `08_GO_enrichment/DEGs_GO_manuscript_mono.R` |
| Fig. 5 | Cross-timepoint monocyte DEG concordance, discordant sets, and unique 14 DPI GO terms | `09_concordant_discordant_analysis/DE_summary_allcelltypes_allcontrast.ipynb` → `10_unique_atD14_analysis/unique_DEGs_GO.ipynb` |
| Fig. 6 | Flow cytometry validation of scRNA-seq immune cell proportions | `11_Flow_cytometry/check_flow_prop_script.R` |

---

## Data Availability

Raw sequencing data (FASTQ files) are deposited at NCBI Gene Expression Omnibus 

**Accession:** 
**URL:** 

---

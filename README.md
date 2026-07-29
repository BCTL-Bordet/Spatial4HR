# Spatial4HR+
R scripts for the spatial transcriptomics analyses presented in the Spatial4HR+ manuscript.

This repository contains the analysis code accompanying the manuscript:

** Spatial Profiling Uncovers a High-Risk Stromal-Fibrotic Subtype in HR-Positive, HER2- Negative Breast Cancer **

The workflow integrates histopathology-guided spatial annotations with spatial transcriptomics, single-cell reference mapping, cell-type deconvolution, pseudobulk transcriptomics, weighted gene co-expression network analysis (WGCNA), multimodal integration and external validation to identify the four Spatial4HR+ subtypes.


# Repository overview

```
Raw Visium spatial transcriptomics
        │
        ▼
01  Preprocessing & quality control
        │
        ▼
02  Histopathology analyses
        │
        ▼
03  Tumour architecture analysis
        │
        ▼
04  Single-cell reference assembly
        │
        ▼
05  CARD deconvolution
        │
        ▼
06  Spatial clustering
        │
        ▼
07  Pseudobulk construction
        │
        ▼
08  Molecular score calculation
        │
        ▼
09  WGCNA
        │
        ▼
10  Multimodal integration (CIMLR)
        │
        ▼
11  Spatial4HR+ genomic risk comparison
        │
        ▼
12  External validation
        │
        ▼
13  CNV & mutation analyses
```

---

# Repository structure
```
├── 00_forestplot.R
├── 01_preprocessing_import_filter_normalize_ST.R
├── 02_histopathology.R
├── 03_tumor_architecture_patterns.R
├── 04_single_cell_reference_assembly.R
├── 05_CARD.R
├── 06_unsupervised_clustering.R
├── 07_pseudobulk_construction.R
├── 08_molecular_scores.R
├── 09_WGCNA.R
├── 10_multimodal_integration.R
├── 11_spatial4hr_genomic_risk_comparison.R
├── 12_external_validation.R
└── 13_CNV_mutation.R
```


# Script descriptions

# 00_forestplot.R
Helper functions for generating forest plots used throughout the manuscript.
Used by multiple downstream analyses.


# 01_preprocessing_import_filter_normalize_ST.R
Preprocess raw Visium spatial transcriptomic data.
## Main analyses
- Import raw Visium datasets
- Import histopathology annotation fraction files
- Merge annotation fractions into Seurat metadata
- Quality control
- Spot filtering
- Data normalisation
- Histopathology-guided filtering
## Output
- Processed Seurat object
- Metadata used throughout downstream analyses


# 02_histopathology.R
Histopathological characterisation of the spatial transcriptomics cohort.
## Main analyses
- Representative H&E sections and manual annotations
- Cohort clinicopathological summary
- Histopathological compartment composition
- Histopathology–cell type associations
- Histopathology–gene expression correlations
## Figures
- Figure 1b
- Figure 1c
- Figure 1d
- Extended Data Figure 1
- Extended Data Figure 2
- Extended Data Figure 3
- Extended Data Figure 4


# 03_tumor_architecture_patterns.R
Characterise tumour architecture using tumour-content binning.
## Main analyses
- Tumour-content binning
- Hierarchical clustering of tumours
- Definition of Cell-dense, Scattered and Local-islands architecture
- Histopathological comparisons
- Cell-type comparisons
- Clinicopathological associations
- Survival analyses
## Figures
- Figure 2a–h
- Extended Data Figure 6


# 04_single_cell_reference_assembly.R
Construct the single-cell reference used for CARD deconvolution.
## Main analyses
- Wu et al. breast cancer reference
- Tang et al. adipocyte reference
- Reference integration
- Preparation of CARD reference object
## Output
Integrated single-cell reference.


# 05_CARD.R
Infer cell-type proportions using CARD.
## Main analyses
- CARD deconvolution
- Spot-level cell-type abundance estimation
- Cell-type correlation analyses
## Figures
- Extended Data Figure 5
## Output
Spot-level inferred cell-type proportions.


# 06_unsupervised_clustering.R
Identify conserved spatial transcriptional programmes.
## Main analyses
- Unsupervised clustering
- Cluster annotation
- Hallmark pathway enrichment
- Cluster abundance
- Histopathological composition
- Survival analyses
- External validation in METABRIC
- External validation in SCAN-B
## Figures
- Figure 3a–i
- Extended Data Figure 7


# 07_pseudobulk_construction.R
Generate patient-level pseudobulk transcriptomes.
## Main analyses
- Whole-sample pseudobulk
- Tumour pseudobulk
- Stromal pseudobulk
## Output
Expression matrices used by downstream analyses.


# 08_molecular_scores.R
Calculate molecular subtype and genomic risk scores.
## Main analyses
- PAM50 (AIMS)
- Genomic Grade Index (GGI)
- MammaPrint
- Oncotype DX (genefu)
## Output
Patient-level molecular subtype assignments and genomic risk scores.


# 09_WGCNA.R
Identify relapse-associated co-expression programmes.
## Main analyses
- Weighted Gene Co-expression Network Analysis
- Module detection
- Module annotation
- GSVA module scoring
- Survival analyses
- Biological process enrichment
- External validation
## Figures
- Figure 4a–f
- Extended Data Figure 8


# 10_multimodal_integration.R
Integrate histopathology, spatial clusters, WGCNA modules and cell-type composition to identify Spatial4HR+ subtypes.
## Main analyses
- CIMLR
- Multi-omic integration
- Spatial4HR+ subtype discovery
- Subtype characterisation
- Hallmark pathway analysis
- Differential expression
- EnR versus SF comparisons
## Figures
- Figure 5a–e
- Extended Data Figure 9


# 11_spatial4hr_genomic_risk_comparison.R
Compare Spatial4HR+ with established molecular prognostic assays.
## Main analyses
- PAM50 classification
- Oncotype DX
- MammaPrint
- GGI
- Survival analyses
- Stromal pseudobulk pathway analysis
- Comparison of spatial versus genomic risk
## Figures
- Figure 6a–e


# 12_external_validation.R
Validate Spatial4HR+ across independent cohorts.
## Main analyses
### METABRIC
- Data preprocessing
- Spatial4HR+ subtype prediction
- Survival analysis
- xCell cell-type inference
- Hallmark pathway analysis
### SCAN-B
- Data preprocessing
- Spatial4HR+ subtype prediction
- Survival analysis
- xCell cell-type inference
- Hallmark pathway analysis
### NeoRHEA
- Data preprocessing
- Spatial4HR+ subtype assignment
- Endocrine therapy response
- Ki67 response
## Figures
- Figure 7a–e
- Figure 7h
- Figure 7i
- Extended Data Figure 10


# 13_CNV_mutation.R
Characterise genomic alterations associated with Spatial4HR+ subtypes.
## Main analyses
- HER2 copy-number status
- Global CNV burden
- Driver mutation frequencies
- Subtype-specific genomic alterations
## Figures
- Figure 7f
- Figure 7g


---

# Main figure mapping

| Figure | Script |
|---------|--------|
| Figure 1 | 02_fig1_histopathology.R |
| Figure 2 | 03_fig2_tumor_architecture_patterns.R |
| Figure 3 | 06_unsupervised_clustering.R |
| Figure 4 | 09_WGCNA.R |
| Figure 5 | 10_multimodal_integration.R |
| Figure 6 | 11_spatial4hr_genomic_risk_comparison.R |
| Figure 7 | 12_external_validation.R, 13_CNV_mutation.R |

---

# Notes
File paths in the scripts should be updated to match local data locations before reproducing the analyses.

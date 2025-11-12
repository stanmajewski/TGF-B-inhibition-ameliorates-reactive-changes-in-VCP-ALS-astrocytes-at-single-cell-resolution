[README.md](https://github.com/user-attachments/files/23499339/README.md)
# Single-Cell RNA-seq Analysis Pipeline for VCP-ALS Astrocytes

This repository contains R scripts for analyzing single-cell RNA sequencing data from hiPSC-derived astrocytes carrying VCP mutations associated with ALS, compared to control astrocytes. The pipeline includes data loading, quality control, integration, clustering, and differential abundance analysis.

## Overview

This analysis pipeline processes Cell Ranger outputs through a complete single-cell RNA-seq workflow including:
- Quality control and filtering
- Batch correction and integration using Harmony
- Unsupervised clustering
- Cell type identification
- Differential abundance analysis using sccomp

## Requirements

### R Version
- R >= 4.2.0 (tested with R 4.3.0)

### Required R Packages

```r
# Core Seurat ecosystem
install.packages("Seurat")
install.packages("Signac")

# Tidyverse and data manipulation
install.packages("tidyverse")
install.packages("patchwork")

# Genomic annotations
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("GenomicRanges")

# Integration
install.packages("harmony")

# Visualization
install.packages("colorRamps")
install.packages("viridis")
install.packages("pheatmap")

# Statistical analysis
install.packages("lmerTest")
install.packages("sccomp")
install.packages("cmdstanr")
install.packages("DESeq2")
install.packages("clusterProfiler")
install.packages("DOSE")
BiocManager::install("org.Hs.eg.db")
install.packages("msigdbr")
install.packages("ggrepel")
install.packages("reshape2")
install.packages("gridExtra")
```

## Directory Structure

```
repository/
├── README.md                      # This file
├── scripts/                       # Analysis scripts (run in order)
│   ├── 01_load_data.R            # Load Cell Ranger outputs and QC
│   ├── 02_integrate_samples.R    # Filter, normalize, and integrate with Harmony
│   ├── 03_clustering_tests.R     # Test clustering resolutions
│   ├── 04_characterize_clusters.R # Cluster characterization and differential abundance
│   ├── 05_pseudobulk_generation.R # Generate pseudobulk samples
│   ├── 06_DESeq2_analysis.R      # Differential expression with DESeq2
│   └── 07_DEG_characterization.R # DEG enrichment and characterization
├── data/
│   ├── metadata/
│   │   └── sample_metadata.csv   # Sample information (groups, treatment, etc.)
│   └── cellranger_outputs/       # Cell Ranger output directories (user provided)
│       ├── CTRL1_U/
│       │   └── count/sample_filtered_feature_bc_matrix/
│       ├── CTRL1_SB/
│       │   └── count/sample_filtered_feature_bc_matrix/
│       └── [additional samples...]/
├── results/                       # Output directory (created automatically)
│   ├── 01_*.rda                  # Seurat objects from script 01
│   ├── 01_*.pdf                  # QC plots from script 01
│   ├── 02_*.rda                  # Integrated Seurat objects
│   ├── 02_*.csv                  # Filter statistics
│   ├── 03_*.rda                  # Clustered objects
│   ├── 03_*.pdf                  # Clustering plots
│   └── 04_*.pdf                  # Final characterization plots
└── reference_data/
    └── cell_type_markers.csv     # Cell type marker genes (user provided)
```

## Data Preparation

### 1. Cell Ranger Outputs

Place your Cell Ranger output folders in `data/cellranger_outputs/`. Each sample folder should contain:

```
data/cellranger_outputs/
├── CTRL1_U/
│   └── count/
│       └── sample_filtered_feature_bc_matrix/
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
└── [additional samples...]/
```

### 2. Sample Metadata

The file `data/metadata/sample_metadata.csv` contains information about each sample:

**Required columns:**
- `sample`: Unique sample identifier
- `geo_accession`: GEO accession number (or "N/A")
- `dataset_folder`: Path to Cell Ranger output (relative to project root)
- `Donor`: Donor identifier
- `group`: Experimental group (e.g., "ctrl", "VCP")
- `treatment`: Treatment condition (e.g., "untreated", "SB431542")
- `Sex`: Donor sex (M/F)
- `Mutation Loci`: VCP mutation information
- `Isogenic`: Isogenic line information
- `Donor Age`: Age of donor
- `Age of onset`: Age of disease onset (for patient lines)

### 3. Reference Data Files

Create the following reference files in `reference_data/`:

**Cell Type Markers** (`cell_type_markers.csv`):
- `gene`: Gene symbol
- `level`: Marker category (e.g., "cell_types", "neuronal_lineage", "NPC_patterning", etc.)

Example format:
```csv
gene,level
GFAP,cell_types
S100B,cell_types
FOXG1,NPC_patterning
```

**Transcription Factors** (`transcription_factors.csv`):
- `Symbol`: Gene symbol of transcription factor
- `Description`: Description/name of the transcription factor

Example format:
```csv
Symbol,Description
TP53,Tumor protein p53
MYC,MYC proto-oncogene
SOX2,SRY-box transcription factor 2
```

This file is used in scripts 06 and 07 to identify transcription factors among DEGs.

## Running the Analysis

### Step 1: Load and QC Data

```r
source("scripts/01_load_data.R")
```

**What it does:**
- Loads Cell Ranger outputs for all samples
- Calculates mitochondrial gene percentage
- Creates quality control plots
- Merges samples into a single Seurat object

**Outputs:**
- `results/01_cellranger_QC.csv` - Cell Ranger QC metrics
- `results/01_cellranger_QC_barplots.pdf` - QC visualizations
- `results/seur_list_patient.rda` - Individual Seurat objects
- `results/01_seur_merge_patient.rda` - Merged Seurat object
- `results/01_QC_plots.pdf` - Seurat QC plots

### Step 2: Filter, Normalize, and Integrate

```r
source("scripts/02_integrate_samples.R")
```

**What it does:**
- Filters low-quality cells (nCount_RNA 500-30000, percent.mt < 5%)
- Removes samples with >30% cells filtered or <1000 cells
- Normalizes with SCTransform
- Integrates samples using Harmony
- Performs PCA and UMAP dimensionality reduction

**Outputs:**
- `results/02_filter_stats.csv` - Cell filtering statistics
- `results/02_gr_tab_filtered.csv` - Updated sample metadata
- `results/02_seur_integr_pat.rda` - Intermediate integration
- `results/seur_integr_pat.rda` - Final integrated Seurat object

### Step 3: Clustering Tests

```r
source("scripts/03_clustering_tests.R")
```

**What it does:**
- Tests multiple clustering resolutions (0.1, 0.2, 0.25, 0.3)
- Creates UMAP visualizations for each resolution
- Generates dot plots of marker gene expression
- Creates a template for cluster annotation

**Outputs:**
- `results/03_seur_integr_filtered.rda` - Seurat object with all resolutions
- `results/03_UMAP_clustering_test_integr_dataset.pdf` - UMAP plots
- `results/03_Dotplots_clustering_test_integr_dataset_cell_type_markers.pdf`
- `results/03_Dotplots_clustering_test_integr_dataset_subtype_markers.pdf`
- `results/03_UMAP_split_by_group_vs_sample.pdf`
- `results/03_marker_expr_UMAP_integr_dataset.pdf`
- `results/03_cluster_assignment.csv` - Template for cluster annotation

**Manual step required:** Edit `results/03_cluster_assignment.csv` to assign cell types to clusters, then save as `results/03_cluster_assignment_updated.csv`

### Step 4: Characterize Clusters

**Before running:** Ensure you have:
1. Completed manual cluster annotation (saved as `results/03_cluster_assignment_updated.csv`)
2. Created the labeled Seurat object as `results/04_seur_integr_labelled.rda`

```r
source("scripts/04_characterize_clusters.R")
```

**What it does:**
- Applies cluster labels from manual annotation
- Creates comprehensive visualizations with labeled clusters
- Performs differential abundance analysis using sccomp
- Tests for significant differences in cluster proportions between groups
- Generates publication-quality figures

**Outputs:**
- `results/04_seur_integr_filtered.rda` - Labeled Seurat object
- `results/04_UMAP_clusters_labelled.pdf` - Labeled UMAP plots
- `results/04_Cell_markers_dotplot_clusters_labelled.pdf` - Marker expression
- `results/04_cell_abundance_by_sample_cluster.csv` - Cell counts and proportions
- `results/04_cell_abundance_by_sample_cluster.pdf` - Bar plots
- `results/04_cell_abundance_by_sample_cluster_heatmap.pdf` - Heatmaps
- `results/04_sccomp_*` - Comprehensive differential abundance results and plots

### Step 5: Generate Pseudobulk Samples

```r
source("scripts/05_pseudobulk_generation.R")
```

**What it does:**
- Aggregates single-cell counts to pseudobulk by cluster and sample
- Filters pseudobulks with <10 cells
- Creates metadata for pseudobulk samples
- Prepares data for differential expression analysis

**Outputs:**
- `results/05_bulk_meta.csv` - Pseudobulk sample metadata
- `results/05_bulk_data.rda` - Pseudobulk count matrix and metadata

### Step 6: Differential Expression Analysis with DESeq2

```r
source("scripts/06_DESeq2_analysis.R")
```

**What it does:**
- Performs differential expression analysis using DESeq2
- Tests multiple comparisons (VCP vs control, treated vs untreated)
- Filters genes by expression level
- Identifies significantly differentially expressed genes (DEGs)
- Separates DEGs by cluster and direction (up/down)

**Outputs:**
- `results/06_bulk_data_with_DESeq_results.rda` - Full DESeq2 results
- `results/06_DESeq2_up_down_by_cluster_genes.csv` - List of DEGs by cluster
- `results/06_DESeq2_up_down_by_cluster_N_genes.csv` - DEG counts summary

### Step 7: DEG Characterization and Enrichment

```r
source("scripts/07_DEG_characterization.R")
```

**What it does:**
- Performs GO (Gene Ontology) enrichment analysis
- Tests pathway enrichment using MSigDB
- Analyzes transcription factor enrichment
- Creates heatmaps of DEG expression
- Quantifies treatment effects on VCP-dysregulated genes
- Generates volcano plots and enrichment visualizations

**Outputs:**
- `results/07_DEGs_numbers_by_comp_and_cluster_deg_characterization.pdf` - DEG counts
- `results/07_*_GO_enrichment.pdf` - GO enrichment plots
- `results/07_*_pathway_enrichment.pdf` - Pathway analysis plots
- `results/07_*_heatmap.pdf` - Expression heatmaps
- `results/07_treatment_effect_*.pdf` - Treatment effect analyses
- `results/07_bulk_data_with_treatment_analysis.rda` - Complete results object

## Analysis Parameters

### Quality Control Thresholds
- **nCount_RNA**: 500 - 30,000 UMIs per cell
- **percent.mt**: < 5% mitochondrial content
- **Sample filtering**: Retain samples with ≥1000 cells and ≤30% cells filtered

### Normalization
- **Method**: SCTransform
- **Variable features**: 2000 genes
- **ncells**: 3000 cells for parameter estimation

### Integration
- **Method**: Harmony
- **Batch variable**: Sample ID
- **Dimensions**: 1-30 PCs

### Clustering
- **Algorithm**: Louvain (algorithm 1)
- **Tested resolutions**: 0.1, 0.2, 0.25, 0.3
- **Default resolution**: 0.3 (user selectable)

### Dimensionality Reduction
- **PCA dimensions**: 1-30
- **UMAP dimensions**: 2D (from Harmony-corrected PCA)

## Expected Runtime

Runtime varies by dataset size and computational resources:
- **Script 01** (Load): ~10-30 minutes
- **Script 02** (Integration): ~30-120 minutes (most computationally intensive)
- **Script 03** (Clustering): ~20-60 minutes
- **Script 04** (Characterization): ~30-90 minutes
- **Script 05** (Pseudobulk): ~10-30 minutes
- **Script 06** (DESeq2): ~20-60 minutes
- **Script 07** (DEG characterization): ~30-90 minutes

**Total**: 2.5-8 hours for complete pipeline

## Memory Requirements

- **Minimum RAM**: 32 GB
- **Recommended RAM**: 64 GB or more for large datasets (>100,000 cells)

Scripts include memory management options:
- `options(future.globals.maxSize = ...)` controls SCTransform memory
- Subsampling for plotting reduces memory usage

## Troubleshooting

### Common Issues

**1. "Cannot find Cell Ranger outputs"**
- Verify paths in `data/metadata/sample_metadata.csv` are correct
- Ensure folder structure matches: `data/cellranger_outputs/SAMPLE/count/sample_filtered_feature_bc_matrix/`

**2. "Memory allocation error"**
- Increase available RAM or use a computing cluster
- Reduce number of samples in initial test run

**3. "Package not found"**
- Install missing packages (see Requirements section)
- Update Bioconductor packages if using older versions

**4. "Harmony integration fails"**
- Ensure sufficient variation between samples
- Check that sample IDs are correctly specified

**5. "sccomp analysis fails"**
- Verify cmdstanr is properly installed
- Check that cluster assignments are valid factors

## Citation

If you use this pipeline, please cite:
- Seurat: Hao et al., Cell 2021
- Harmony: Korsunsky et al., Nature Methods 2019
- sccomp: Mangiola et al., PNAS 2023

## Contact

For questions about this analysis pipeline, please open an issue in this repository.

## License

This analysis code is provided for research purposes. Please refer to individual package licenses for usage restrictions.

# C9 snRNAseq Analysis Pipeline

Comprehensive single-nucleus RNA sequencing analysis pipeline for C9-ALS/FTD cortex data processing, quality control, integration, clustering, and cell type annotation.

## Overview

This pipeline processes Cell Ranger outputs through to annotated, integrated single-nucleus RNA-seq datasets with differential abundance analysis. The workflow implements Harmony batch correction, multi-resolution clustering, and automated cell type annotation using scType.

## Pipeline Structure

### C01_load_cellranger_data.R
**Purpose:** Load and quality control Cell Ranger outputs

**Functions:**
- Loads Cell Ranger metrics and count matrices
- Creates individual Seurat objects per sample
- Calculates mitochondrial percentages
- Merges samples into unified dataset
- Generates comprehensive QC visualizations

**Inputs:**
- `data/metadata/sample_metadata.csv`: Sample information and file paths
- Cell Ranger H5 files (paths specified in metadata)

**Outputs:**
- `results/01_data_loading/cellranger_qc.csv`
- `results/01_data_loading/cellranger_qc_barplots.pdf`
- `results/01_data_loading/seurat_list.rda`
- `results/01_data_loading/merged_seurat.rda`
- `results/01_data_loading/qc_plots.pdf`

**Key Parameters:**
- No filtering applied at this stage
- Mitochondrial genes identified by "^MT-" pattern

---

### C02_harmony_integration.R
**Purpose:** Filter low-quality cells and integrate samples with Harmony

**Functions:**
- Applies quality filters (UMI counts, mitochondrial %)
- Removes samples with excessive cell loss
- Performs SCTransform normalization
- Runs Harmony batch correction
- Generates UMAP embedding and nearest neighbor graph

**Inputs:**
- `results/01_data_loading/merged_seurat.rda`
- `data/metadata/sample_metadata.csv`

**Outputs:**
- `results/02_integration/filter_statistics.csv`
- `results/02_integration/filtered_metadata.csv`
- `results/02_integration/integrated_seurat.rda`

**Key Parameters:**
- UMI count range: 500 - 30,000
- Mitochondrial percentage: < 10%
- Sample retention criteria: <50% cells removed, ≥1,000 cells remaining
- SCTransform: 3,000 cells, 2,000 variable features
- Harmony dimensions: 1-30

---

### C03_clustering_analysis.R
**Purpose:** Multi-resolution clustering and initial characterization

**Functions:**
- Tests clustering resolutions (0.3 - 1.5)
- Generates UMAP visualizations for each resolution
- Creates marker gene dot plots
- Checks sample integration quality
- Produces cluster assignment template

**Inputs:**
- `results/02_integration/integrated_seurat.rda`
- `results/02_integration/filtered_metadata.csv`
- `data/reference/cell_type_markers.csv`

**Outputs:**
- `results/03_clustering/clustered_seurat.rda`
- `results/03_clustering/umap_clustering_test.pdf`
- `results/03_clustering/dotplot_cell_type_markers.pdf`
- `results/03_clustering/dotplot_subtype_markers.pdf`
- `results/03_clustering/umap_by_sample.pdf`
- `results/03_clustering/umap_marker_expression.pdf`
- `results/03_clustering/cluster_assignment_template.csv`

**Key Parameters:**
- Test resolutions: 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.5
- Clustering algorithm: Louvain (algorithm = 1)

**Manual Step Required:**
After running this script, manually annotate clusters in `cluster_assignment_template.csv` based on marker expression patterns.

---

### C04_cluster_characterization.R
**Purpose:** Apply cluster labels and perform differential abundance analysis

**Functions:**
- Applies user-annotated cluster labels
- Generates labeled visualizations
- Quantifies cell abundance across samples/clusters
- Performs statistical testing (t-tests)
- Runs sccomp compositional analysis

**Inputs:**
- `results/03_clustering/clustered_seurat.rda`
- `results/02_integration/filtered_metadata.csv`
- `results/03_clustering/cluster_assignment_template.csv` (user-annotated)
- `data/reference/cell_type_markers.csv`

**Outputs:**
- `results/04_characterization/labeled_seurat.rda`
- `results/04_characterization/umap_labeled_clusters.pdf`
- `results/04_characterization/dotplot_labeled_clusters.pdf`
- `results/04_characterization/umap_by_sample_labeled.pdf`
- `results/04_characterization/abundance_by_sample_cluster.csv`
- `results/04_characterization/abundance_barplot.pdf`
- `results/04_characterization/abundance_heatmap.pdf`
- `results/04_characterization/abundance_ttest.csv`
- `results/04_characterization/sccomp_results.csv`
- `results/04_characterization/sccomp_plots.pdf`

**Key Parameters:**
- Default clustering resolution: 0.1 (adjust in script if needed)
- Statistical tests: Pairwise t-tests with Benjamini-Hochberg correction
- sccomp: Compositional analysis with group as factor

---

### C05_sctype_annotation.R
**Purpose:** Automated cell type annotation using scType

**Functions:**
- Loads scType framework from GitHub
- Prepares brain-specific marker gene sets
- Performs automated cell type classification
- Generates comprehensive visualization suite
- Extracts astrocyte subset for downstream analysis

**Inputs:**
- `results/04_characterization/labeled_seurat.rda`
- `results/02_integration/filtered_metadata.csv`
- scType database (loaded from GitHub)

**Outputs:**
- `results/05_sctype/sctype_annotated_seurat.rda`
- `results/05_sctype/astrocytes_subset.rda`
- `results/05_sctype/sctype_umap.pdf`
- `results/05_sctype/cell_type_distribution.pdf`
- `results/05_sctype/marker_expression_heatmap.pdf`
- `results/05_sctype/violin_plots.pdf`
- `results/05_sctype/feature_plots.pdf`
- `results/05_sctype/cluster_composition.pdf`
- `results/05_sctype/sctype_complete_report.pdf`
- `results/05_sctype/cell_type_summary.csv`

**Key Parameters:**
- Uses resolution 0.1 clustering for classification
- Brain-specific marker sets for 8 major cell types
- Astrocyte markers: GFAP, AQP4, ALDH1L1, S100B, GLUL, SLC1A3

---

## Required Input Files

### Sample Metadata
Create `data/metadata/sample_metadata.csv` with columns:
- `sample`: Unique sample identifier
- `group`: Experimental group (e.g., "Control", "C9-ALS")
- `cellranger_path`: Path to Cell Ranger output directory
- `h5_path`: Path to filtered_feature_bc_matrix.h5 file
- Additional metadata columns as needed

Example:
```csv
sample,group,cellranger_path,h5_path
Sample1,Control,/path/to/sample1/outs,/path/to/sample1/outs/filtered_feature_bc_matrix.h5
Sample2,C9-ALS,/path/to/sample2/outs,/path/to/sample2/outs/filtered_feature_bc_matrix.h5
```

### Cell Type Markers (Optional)
Create `data/reference/cell_type_markers.csv` with columns:
- `gene`: Gene symbol
- `level`: Marker category (e.g., "cell_types", "neuronal_lineage", "cortical_layers")

If not provided, scripts will use default marker sets.

---

## Running the Pipeline

### Step 1: Prepare Environment
```bash
# Create directory structure
mkdir -p data/metadata
mkdir -p data/reference
mkdir -p results

# Prepare sample metadata file
# Edit data/metadata/sample_metadata.csv with your sample information
```

### Step 2: Load and QC Data
```bash
Rscript C01_load_cellranger_data.R
```
Review QC plots to assess data quality before proceeding.

### Step 3: Integration
```bash
Rscript C02_harmony_integration.R
```
Check filter statistics and verify sample retention.

### Step 4: Clustering
```bash
Rscript C03_clustering_analysis.R
```
**Manual intervention required:** Review clustering plots and annotate `cluster_assignment_template.csv`

### Step 5: Characterization
```bash
Rscript C04_cluster_characterization.R
```
Requires completed cluster annotations from Step 4.

### Step 6: Automated Annotation
```bash
Rscript C05_sctype_annotation.R
```

---

## Dependencies

### R Packages
- **Seurat** (≥4.0): Core single-cell analysis
- **tidyverse**: Data manipulation and visualization
- **harmony**: Batch correction integration
- **sccomp**: Compositional analysis
- **viridis**: Color palettes
- **pheatmap**: Heatmap generation
- **patchwork**: Plot composition
- **colorRamps**: Additional color palettes

### External Resources
- scType framework (auto-loaded from GitHub)
- scType brain tissue database (auto-loaded from GitHub)

### Installation
```r
# Install from CRAN
install.packages(c("tidyverse", "viridis", "pheatmap", "patchwork", 
                   "colorRamps", "HGNChelper", "openxlsx"))

# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Seurat", "harmony"))

# Install sccomp
devtools::install_github("MangiolaLaboratory/sccomp")
```

---

## Key Features

### Quality Control
- Cell Ranger metrics tracking
- Comprehensive filtering (UMI, gene counts, mitochondrial %)
- Sample-level quality assessment
- Automatic removal of low-quality samples

### Integration
- Harmony batch correction
- SCTransform normalization
- Layer-wise processing for memory efficiency
- UMAP embedding with Harmony-corrected dimensions

### Clustering
- Multi-resolution testing (9 resolutions)
- Visual assessment tools (UMAP, dot plots)
- Sample integration verification
- Template-based cluster annotation

### Analysis
- Cell abundance quantification
- Differential abundance testing (t-tests, sccomp)
- Comprehensive visualizations (UMAP, heatmaps, bar plots)
- Automated cell type annotation with scType

---

## Output Organization

```
results/
├── 01_data_loading/          # Initial QC and merged data
├── 02_integration/           # Filtered and integrated data
├── 03_clustering/            # Multi-resolution clustering
├── 04_characterization/      # Labeled clusters and abundance
└── 05_sctype/               # Automated annotations
```

---

## Computational Requirements

### Memory
- Minimum: 32 GB RAM
- Recommended: 64+ GB RAM for large datasets (>100,000 cells)

### Storage
- ~10-20 GB per sample for intermediate files
- Final integrated objects: 5-10 GB depending on dataset size

### Runtime
- C01: 30-60 minutes
- C02: 1-2 hours
- C03: 1-2 hours
- C04: 30-60 minutes
- C05: 30-60 minutes

---

## Troubleshooting

### Memory Issues
- Reduce `ncells` parameter in SCTransform
- Process samples in batches
- Increase system swap space

### Integration Quality
- Check filter statistics - ensure sufficient cells per sample
- Review UMAP split by sample
- Adjust Harmony parameters if needed

### Clustering
- Test multiple resolutions
- Verify marker gene expression patterns
- Check for over/under-clustering

---

## Citation

If using this pipeline, please cite:
- Seurat: Hao et al., Cell 2021
- Harmony: Korsunsky et al., Nat Methods 2019
- sccomp: Mangiola et al., PNAS 2023
- scType: Ianevski et al., Nat Commun 2022

---

## Contact

For questions or issues with this pipeline, please refer to:
- Seurat documentation: https://satijalab.org/seurat/
- Harmony documentation: https://github.com/immunogenomics/harmony
- scType repository: https://github.com/IanevskiAleksandr/sc-type

---

## License

This pipeline is provided for academic research use. Please ensure compliance with individual package licenses.

---

**Last Updated:** November 2025
**Pipeline Version:** 1.0

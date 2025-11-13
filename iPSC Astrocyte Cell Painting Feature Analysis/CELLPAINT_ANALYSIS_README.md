# Cell Painting Analysis Pipeline

Comprehensive analysis pipeline for Cell Painting high-content imaging data comparing VCP-ALS patient-derived astrocytes (treated and untreated) versus healthy control astrocytes.

## Overview

This pipeline analyzes morphological and intensity features extracted from Cell Painting experiments to:
1. Identify cellular features that define the VCP disease signature
2. Evaluate treatment effects on disease-associated features
3. Characterize responses across different feature categories (AreaShape, Intensity, Texture, etc.)

## Scripts

### 01. Comprehensive Analysis (`cellpaint_01_comprehensive_analysis.R`)

**Purpose:** Complete statistical analysis of differential Cell Painting features with treatment effect classification.

**Key Features:**
- Volcano plots of differential features
- Top up/down-regulated feature identification
- Treatment impact classification (Corrective/Exacerbating/Neutral)
- Effect size and p-value change analysis
- Category-level summaries

**Usage:**
```bash
Rscript cellpaint_01_comprehensive_analysis.R \
  --untreated ALS_Untreated_vs_Control_statistical_results.csv \
  --treated ALS_Treated_vs_Control_statistical_results.csv \
  --output cellpaint_analysis_results \
  --pvalue 0.05
```

**Outputs:**
- `01_volcano_plot_enhanced.png` - Enhanced volcano plot with significance annotations
- `02a_top_upregulated_features.png` - Top 15 increased features
- `02b_top_downregulated_features.png` - Top 15 decreased features
- `03_category_distribution.png` - Feature distribution by category and direction
- `05_treatment_effect_scatter.png` - Effect size comparison plot
- `06_significance_change_scatter.png` - P-value change plot
- `09_treatment_impact_by_category.png` - Treatment effects by category
- `treatment_comparison_results.csv` - Complete comparison data
- `analysis_summary_statistics.csv` - Summary metrics
- `top_corrective_features.csv` - Features improved by treatment
- `top_exacerbating_features.csv` - Features worsened by treatment

---

### 02. Treatment Signature Analysis (`cellpaint_02_treatment_signature_analysis.R`)

**Purpose:** Focused analysis of treatment effects on features defining the VCP disease signature.

**Key Features:**
- VCP signature definition (significantly differential features)
- Refined treatment impact classification with effect size thresholds
- Corrective vs exacerbating vs unchanged vs mixed categorization
- Category-specific treatment responses

**Usage:**
```bash
Rscript cellpaint_02_treatment_signature_analysis.R \
  --untreated ALS_Untreated_vs_Control_statistical_results.csv \
  --treated ALS_Treated_vs_Control_statistical_results.csv \
  --biological biological_replicates_data.csv \
  --output cellpaint_treatment_analysis \
  --pvalue 0.05 \
  --effect-threshold 0.1
```

**Classification Criteria:**
- **Corrective:** P-value increases AND effect size decreases (with meaningful change ≥ threshold)
  - OR direction reversal with meaningful change
- **Exacerbating:** P-value decreases AND effect size increases (with meaningful change)
- **Unchanged:** Effect size change below threshold
- **Mixed:** Other patterns

**Outputs:**
- `01_treatment_impact_overview.png` - Summary of impact categories
- `02_effect_size_changes.png` - Effect size scatter plot
- `03_pvalue_changes.png` - Statistical significance changes
- `04_impact_by_category.png` - Treatment effects by feature category
- `05_top_corrective_features.png` - Top 20 corrective features
- `06_top_exacerbating_features.png` - Top 20 exacerbating features
- `vcp_disease_signature.csv` - Features defining VCP signature
- `treatment_effects_on_signature.csv` - Complete treatment analysis
- `corrective_features.csv` - All corrective features
- `exacerbating_features.csv` - All exacerbating features
- `unchanged_features.csv` - Features unaffected by treatment
- `mixed_features.csv` - Features with mixed responses

---

### 03. Category-Specific Analysis (`cellpaint_03_category_analysis.R`)

**Purpose:** Analysis organized by CellProfiler measurement categories to understand compartment and feature-type specific responses.

**Key Features:**
- CellProfiler measurement type classification:
  - AreaShape (morphology)
  - Intensity (fluorescence)
  - Texture (spatial patterns)
  - Granularity (texture at multiple scales)
  - RadialDistribution (spatial organization)
  - Correlation (channel co-localization)
  - ObjectSkeleton (shape skeleton features)
- Category-specific volcano plots
- Treatment impact stratified by category

**Usage:**
```bash
Rscript cellpaint_03_category_analysis.R \
  --untreated ALS_Untreated_vs_Control_statistical_results.csv \
  --treated ALS_Treated_vs_Control_statistical_results.csv \
  --output cellpaint_category_analysis \
  --pvalue 0.05 \
  --top-features 3
```

**Outputs:**
- `01_volcano_plot_by_category.png` - Volcano plot colored by category
- `02_category_distribution.png` - Significant features by category
- `03_top_features_per_category.png` - Top features in each category
- `04_treatment_impact_by_category.png` - Treatment effects by category
- `05_effect_size_by_category.png` - Effect size distributions
- `06_corrective_features_by_category.png` - Category-specific heatmap
- `07_volcano_plots_faceted.png` - Faceted volcano plots
- `category_summary_statistics.csv` - Per-category statistics
- `treatment_effects_by_category.csv` - Treatment impact by category
- `treatment_comparison_with_categories.csv` - Full annotated results

## Input Data Format

### Statistical Results Files

Required columns:
- `Feature` - Feature name (e.g., "Cells_AreaShape_Area", "Nuclei_Intensity_IntegratedIntensity_DAPI")
- `P_value` - P-value from statistical test
- `P_value_FDR` - FDR-corrected p-value
- `Mean_difference` - Mean difference between conditions
- `Effect_size` - Standardized effect size (e.g., Cohen's d)
- `Significant_FDR` - Logical/character indicating FDR significance

### Optional: Biological Replicates Data

For box plot generation (Script 02):
- `Metadata_Patient_Line` - Patient identifier
- `Condition` - Experimental condition ("Control Untreated", "VCP Untreated", "VCP Treated", etc.)
- Feature columns with numeric values

## Parameters

### Common Parameters
- `--untreated` / `-u` - Path to untreated vs control results (default: `ALS_Untreated_vs_Control_statistical_results.csv`)
- `--treated` / `-t` - Path to treated vs control results (default: `ALS_Treated_vs_Control_statistical_results.csv`)
- `--output` / `-o` - Output directory name
- `--pvalue` / `-p` - P-value threshold for significance (default: 0.05)

### Script-Specific Parameters
- Script 02:
  - `--biological` / `-b` - Path to biological replicates data (optional)
  - `--effect-threshold` / `-e` - Effect size change threshold (default: 0.1)
- Script 03:
  - `--top-features` / `-n` - Number of top features per category (default: 3)

## Output Organization

Each script creates a separate output directory containing:
- **PNG files** - High-resolution plots (300 DPI) for publication
- **CSV files** - Detailed results tables for further analysis
- **Log output** - Timestamped progress messages and statistics

## Feature Categories

The pipeline automatically classifies features into categories based on CellProfiler naming:

| Category | Description | Example Features |
|----------|-------------|------------------|
| **AreaShape** | Morphological measurements | Area, Perimeter, Compactness, Eccentricity |
| **Intensity** | Fluorescence intensity | IntegratedIntensity, MeanIntensity, MedianIntensity |
| **Texture** | Spatial texture patterns | Haralick features (contrast, correlation, entropy) |
| **Granularity** | Multi-scale texture | Granularity at different pixel scales |
| **RadialDistribution** | Spatial distribution | Radial intensity profiles, zernike moments |
| **Correlation** | Channel co-localization | Pearson/Spearman correlation between channels |
| **ObjectSkeleton** | Skeleton-based features | Branch lengths, endpoints |

## Dependencies

Required R packages:
```r
tidyverse      # Data manipulation and visualization
dplyr          # Data manipulation
ggplot2        # Plotting
gridExtra      # Multiple plot layouts
pheatmap       # Heatmaps
RColorBrewer   # Color palettes
corrplot       # Correlation plots
VennDiagram    # Venn diagrams
ggrepel        # Label repulsion
scales         # Scale functions
optparse       # Command-line argument parsing
```

Install with:
```r
install.packages(c("tidyverse", "gridExtra", "pheatmap", "RColorBrewer", 
                   "corrplot", "VennDiagram", "ggrepel", "scales", "optparse"))
```

## Treatment Impact Classification

The pipeline uses consistent criteria across scripts to classify treatment effects:

### Corrective Features
Treatment normalizes the feature toward control values:
- P-value increases (less significant in treated)
- AND effect size decreases
- OR direction reversal (sign change)
- Requires meaningful change (≥ effect size threshold in Script 02)

### Exacerbating Features  
Treatment worsens the disease signature:
- P-value decreases (more significant in treated)
- AND effect size increases in same direction
- Requires meaningful change (in Script 02)

### Neutral/Unchanged Features
- Minimal change between conditions (Script 02 uses effect size threshold)
- Changes don't meet corrective or exacerbating criteria

### Mixed Features
- Complex patterns not fitting other categories
- Only in Script 02 with refined classification

## Example Workflow

```bash
# 1. Run comprehensive analysis
Rscript cellpaint_01_comprehensive_analysis.R \
  --untreated data/ALS_Untreated_vs_Control_statistical_results.csv \
  --treated data/ALS_Treated_vs_Control_statistical_results.csv \
  --output results/comprehensive \
  --pvalue 0.05

# 2. Analyze treatment effects on disease signature
Rscript cellpaint_02_treatment_signature_analysis.R \
  --untreated data/ALS_Untreated_vs_Control_statistical_results.csv \
  --treated data/ALS_Treated_vs_Control_statistical_results.csv \
  --output results/treatment_signature \
  --pvalue 0.05 \
  --effect-threshold 0.1

# 3. Examine category-specific patterns
Rscript cellpaint_03_category_analysis.R \
  --untreated data/ALS_Untreated_vs_Control_statistical_results.csv \
  --treated data/ALS_Treated_vs_Control_statistical_results.csv \
  --output results/categories \
  --pvalue 0.05 \
  --top-features 5
```

## Interpretation Guidelines

### Effect Size Interpretation
- **Small effect:** |d| < 0.2
- **Medium effect:** 0.2 ≤ |d| < 0.5
- **Large effect:** |d| ≥ 0.5
- **Very large effect:** |d| ≥ 0.8

### P-value Thresholds
- Standard significance: p < 0.05
- Stringent threshold: p < 0.01
- High confidence: p < 0.001
- FDR-corrected: Controls for multiple testing

### Category Insights
- **AreaShape changes** suggest morphological alterations
- **Intensity changes** indicate expression level differences
- **Texture patterns** reflect spatial organization
- **Correlation changes** show altered co-localization

## Troubleshooting

### Common Issues

**Missing features between datasets:**
- Scripts automatically handle with `complete.cases()`
- Check for consistent feature naming

**Memory issues with large datasets:**
- Filter to significant features first
- Reduce plot resolution (change `dpi` parameter)

**Category classification errors:**
- Verify CellProfiler naming conventions
- Check for custom feature names

**No significant features:**
- Adjust p-value threshold
- Check input data quality
- Verify statistical power

## Citation

If you use this pipeline, please cite the original publication and acknowledge the analysis methods used.

## Contact

For questions or issues with the analysis pipeline, please open an issue in the repository.

## License

This analysis pipeline is provided for research use. See LICENSE file for details.

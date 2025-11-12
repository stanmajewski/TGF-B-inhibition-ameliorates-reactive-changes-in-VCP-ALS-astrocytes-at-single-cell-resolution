# Metadata Conversion Guide

## Original to Anonymized Metadata Mapping

This document explains how the original metadata was converted to the anonymized, public-release version.

---

## Column Name Changes

### Original → Anonymized

| Original Column | Anonymized Column | Change Reason |
|----------------|-------------------|---------------|
| `sample` | `sample` | Converted from codes (A1, C1) to descriptive names |
| `geo_accession` | `geo_accession` | **Kept** - Public identifiers |
| `dataset_folder` | `h5_path` | Separated into two columns with relative paths |
| (new) | `cellranger_path` | Added for compatibility with standard workflows |
| `Donor` | `donor_id` | Renamed for clarity |
| `group` | `group` | Changed values: CTRL → Control |
| `sex` | `sex` | **Unchanged** |
| `age` | `age` | **Unchanged** |
| `sequencing_batch` | `sequencing_batch` | Renamed values: batch_1 → Batch_1 |
| (new) | `tissue` | Added to document tissue type |

---

## Sample Name Mapping

### C9-ALS Patients

| Original | Anonymized | Rationale |
|----------|------------|-----------|
| A1 | C9ALS_Patient1 | More descriptive, indicates condition |
| A2 | C9ALS_Patient2 | Consistent naming scheme |
| A3 | C9ALS_Patient3 | Clear patient identification |
| A4 | C9ALS_Patient4 | Self-documenting |
| A5 | C9ALS_Patient5 | Professional naming |
| A6 | C9ALS_Patient6 | Explicit grouping |

### Control Samples

| Original | Anonymized | Rationale |
|----------|------------|-----------|
| C1 | Control_1 | Clear control designation |
| C2 | Control_2 | Avoids confusion with C9 |
| C3 | Control_3 | Consistent format |
| C4 | Control_4 | Professional standard |
| C5 | Control_5 | Self-explanatory |
| C6 | Control_6 | Clear identification |

---

## Group Label Changes

### Original → Anonymized

| Original | Anonymized | Reason |
|----------|------------|--------|
| `C9-ALS` | `C9-ALS` | **Unchanged** - Standard terminology |
| `CTRL` | `Control` | Expanded abbreviation for clarity |

Both formats are valid, but "Control" is more explicit and follows best practices for documentation.

---

## Batch Label Formatting

### Original → Anonymized

| Original | Anonymized | Change |
|----------|------------|--------|
| `batch_1` | `Batch_1` | Capitalized for consistency |
| `batch_2` | `Batch_2` | Professional formatting |
| `batch_3` | `Batch_3` | Standard naming |
| `batch_4` | `Batch_4` | Clear hierarchy |
| `batch_5` | `Batch_5` | Improved readability |
| `batch_6` | `Batch_6` | Consistent style |

---

## File Path Anonymization

### Original Paths (Institutional)

```
/nemo/project/proj-luscombn-patani/working/public-data/scrnaseq-motor-cortex-c9-ALS/processed_files/GSM6781907_snRNA_C9ALS_A1_MotorCortex_cellBender_corrected_filtered.h5
```

**Issues:**
- Contains institutional server path (`/nemo/`)
- Includes project identifier (`proj-luscombn-patani`)
- References specific directory structure
- Not portable to other systems

### Anonymized Paths (Portable)

```
data/h5_files/GSM6781907_snRNA_C9ALS_A1_MotorCortex_cellBender_corrected_filtered.h5
```

**Improvements:**
- Relative to project root
- Generic directory structure
- Works on any system
- No institutional information
- Preserves GEO accession and file details

---

## Path Structure Breakdown

### H5 Path Components

```
data/h5_files/[GEO_ACCESSION]_snRNA_[GROUP]_[SAMPLE]_MotorCortex_cellBender_corrected_filtered.h5
```

Example:
```
data/h5_files/GSM6781907_snRNA_C9ALS_A1_MotorCortex_cellBender_corrected_filtered.h5
│              │          │       │    │            │                       │
│              │          │       │    │            │                       └─ Format
│              │          │       │    │            └───────────────────────── Processing
│              │          │       │    └──────────────────────────────────── Tissue
│              │          │       └───────────────────────────────────────── Original Sample
│              │          └───────────────────────────────────────────────── Group
│              └──────────────────────────────────────────────────────────── GEO ID
└─────────────────────────────────────────────────────────────────────────── Directory
```

### Cell Ranger Path

```
data/cellranger_outputs/[ANONYMIZED_SAMPLE_NAME]/
```

This structure assumes Cell Ranger was run with anonymized sample names. If using original Cell Ranger outputs, update paths accordingly.

---

## Metadata File Comparison

### Original File Structure
```csv
sample,geo_accession,dataset_folder,Donor,group,sex,age,sequencing_batch
A1,https://...,/nemo/project/.../GSM6781907_...h5,A1,C9-ALS,M,72,batch_5
```

### Anonymized File Structure
```csv
sample,geo_accession,h5_path,cellranger_path,donor_id,group,sex,age,sequencing_batch,tissue
C9ALS_Patient1,GSM6781907,data/h5_files/GSM6781907_...h5,data/cellranger_outputs/C9ALS_Patient1,Donor_A1,C9-ALS,M,72,Batch_5,Motor_Cortex
```

---

## Column Purpose and Usage

### Required Columns (Used by Pipeline)

1. **sample** - Primary identifier used throughout analysis
   - Must be unique
   - Used for file naming
   - Used in visualizations

2. **h5_path** - Location of count matrix
   - Must point to valid H5 file
   - Can be absolute or relative path
   - Used by C01 script

3. **group** - Experimental condition
   - Used for differential analysis
   - Used in visualizations
   - Batch correction factor

4. **sequencing_batch** - Technical batch identifier
   - Used by Harmony for batch correction
   - Important for data quality

### Optional Columns (Metadata only)

5. **geo_accession** - Public data identifier
   - For reproducibility
   - Data provenance
   - Not used in analysis code

6. **cellranger_path** - Original Cell Ranger output
   - For accessing additional QC files
   - Metrics summary location
   - Optional if only using H5 files

7. **donor_id** - Biological specimen identifier
   - For record keeping
   - Linking samples across studies
   - Not used in core analysis

8. **sex**, **age** - Demographic information
   - For population description
   - Potential covariates
   - Quality control factors

9. **tissue** - Anatomical source
   - Documentation
   - Study description
   - Not used in single-tissue analyses

---

## Usage Notes

### In Analysis Scripts

The scripts primarily use these columns:
- **sample**: Sample identification and labeling
- **h5_path**: Loading count data
- **group**: Grouping for statistics and visualization
- **sequencing_batch**: Batch correction (optional)

### Customization

You can add additional columns for:
- Disease subtype
- Treatment history
- Post-mortem interval
- RNA integrity number (RIN)
- Additional technical metadata

Just ensure required columns remain intact.

---

## Validation Checklist

Before using the metadata file:

- [ ] All sample names are unique
- [ ] H5 paths point to existing files
- [ ] Group labels are consistent
- [ ] No missing values in required columns
- [ ] File is proper CSV format
- [ ] Column names match expected format
- [ ] Paths are appropriate for your system

---

## Converting Your Own Metadata

If you have data in a different format, use this template:

```r
# R code to convert metadata
original <- read.csv("original_metadata.csv")

converted <- data.frame(
  sample = paste0("Sample_", 1:nrow(original)),  # Create new names
  geo_accession = original$GEO_ID,                # Keep public IDs
  h5_path = file.path("data/h5_files", basename(original$file_path)),  # Relative paths
  cellranger_path = file.path("data/cellranger_outputs", paste0("Sample_", 1:nrow(original))),
  donor_id = paste0("Donor_", original$subject_id),
  group = original$condition,
  sex = original$sex,
  age = original$age_years,
  sequencing_batch = paste0("Batch_", original$batch),
  tissue = "Motor_Cortex"  # Or your tissue type
)

write.csv(converted, "sample_metadata.csv", row.names = FALSE)
```

---

## Common Issues and Solutions

### Issue 1: Sample Names with Spaces
**Problem:** `"C9-ALS Patient 1"` causes parsing errors  
**Solution:** Use underscores: `"C9ALS_Patient_1"`

### Issue 2: Absolute Paths
**Problem:** `/home/user/data/file.h5` not portable  
**Solution:** Use relative: `data/h5_files/file.h5`

### Issue 3: Missing H5 Files
**Problem:** Path in metadata doesn't match actual files  
**Solution:** Verify with `file.exists()` in R or check paths

### Issue 4: Group Label Typos
**Problem:** Mix of "Control", "control", "CTRL"  
**Solution:** Standardize to one format, ensure consistency

### Issue 5: Missing Required Columns
**Problem:** Pipeline expects columns that don't exist  
**Solution:** Add columns with default values if not available

---

## Best Practices

1. **Use Descriptive Names**: `C9ALS_Patient1` > `A1`
2. **Relative Paths**: Portable across systems
3. **Consistent Formatting**: Same style throughout
4. **Document Changes**: Keep this guide updated
5. **Version Control**: Track metadata changes
6. **Validate Early**: Check format before analysis
7. **Backup Original**: Keep original metadata file

---

**Last Updated:** November 2025  
**Version:** 1.0  
**Compatible with:** C9 snRNAseq Analysis Pipeline v1.0

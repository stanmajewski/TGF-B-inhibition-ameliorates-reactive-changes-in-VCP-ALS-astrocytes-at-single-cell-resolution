# C9-ALS Motor Cortex snRNAseq Dataset - Metadata Guide

## Dataset Information

This pipeline analyzes publicly available single-nucleus RNA sequencing data from C9orf72 ALS/FTD motor cortex samples.

**Original Publication:**  
Divergent impacts of C9orf72 repeat expansion on neurons and glia in ALS and FTD

**GEO Accession:** GSE219281 (parent series)

**Tissue:** Motor Cortex  
**Technology:** 10x Genomics single-nucleus RNA-seq  
**Processing:** CellBender corrected, filtered count matrices

---

## Sample Metadata Structure

### Metadata File: `data/metadata/sample_metadata.csv`

This file contains information about all samples in the analysis.

### Column Descriptions

| Column | Type | Description |
|--------|------|-------------|
| `sample` | string | Unique sample identifier (e.g., "C9ALS_Patient1", "Control_1") |
| `geo_accession` | string | GEO accession number for public data retrieval |
| `h5_path` | string | Path to CellBender-corrected H5 file |
| `cellranger_path` | string | Path to Cell Ranger output directory |
| `donor_id` | string | Anonymized donor identifier |
| `group` | string | Experimental group: "C9-ALS" or "Control" |
| `sex` | string | Biological sex: "M" (Male) or "F" (Female) |
| `age` | integer | Age at time of death (years) |
| `sequencing_batch` | string | Sequencing batch identifier (Batch_1 through Batch_6) |
| `tissue` | string | Tissue type: "Motor_Cortex" |

---

## Sample Cohorts

### C9-ALS Patients (n=6)
- **C9ALS_Patient1**: Male, 72 years, Batch_5
- **C9ALS_Patient2**: Female, 66 years, Batch_4
- **C9ALS_Patient3**: Female, 71 years, Batch_2
- **C9ALS_Patient4**: Female, 58 years, Batch_1
- **C9ALS_Patient5**: Male, 63 years, Batch_6
- **C9ALS_Patient6**: Male, 61 years, Batch_3

**Demographics:**
- Sex: 3 Male, 3 Female
- Age range: 58-72 years
- Mean age: 65.2 years

### Control Samples (n=6)
- **Control_1**: Female, 67 years, Batch_1
- **Control_2**: Male, 64 years, Batch_3
- **Control_3**: Male, 70 years, Batch_5
- **Control_4**: Male, 60 years, Batch_6
- **Control_5**: Female, 74 years, Batch_4
- **Control_6**: Female, 54 years, Batch_2

**Demographics:**
- Sex: 3 Male, 3 Female
- Age range: 54-74 years
- Mean age: 64.8 years

---

## Data Organization

### Directory Structure

```
project_root/
├── data/
│   ├── metadata/
│   │   └── sample_metadata.csv          # This file
│   ├── h5_files/                        # CellBender-corrected H5 files
│   │   ├── GSM6781907_snRNA_C9ALS_A1_MotorCortex_cellBender_corrected_filtered.h5
│   │   ├── GSM6781909_snRNA_C9ALS_A2_MotorCortex_cellBender_corrected_filtered.h5
│   │   └── ... (all sample H5 files)
│   └── cellranger_outputs/             # Cell Ranger output directories
│       ├── C9ALS_Patient1/
│       ├── C9ALS_Patient2/
│       └── ... (all samples)
└── results/                             # Analysis outputs
```

---

## Data Acquisition

### Option 1: Download from GEO

The original data files are publicly available from GEO:

```bash
# Example for downloading one sample
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6781907&format=file"

# Or use SRA toolkit for fastq files
prefetch SRR_number
fasterq-dump SRR_number
```

### Option 2: Use Preprocessed H5 Files

If you have access to the CellBender-corrected H5 files:

1. Place all H5 files in `data/h5_files/`
2. Ensure filenames match those in the metadata file
3. Update paths in `sample_metadata.csv` if needed

---

## Important Notes

### CellBender Processing

These data have been processed with CellBender to remove ambient RNA contamination. This is why the files are labeled as:
- `cellBender_corrected_filtered.h5`

The H5 files contain:
- Filtered count matrices
- Gene names (ENSEMBL IDs)
- Barcode information

### File Naming Convention

GEO filenames follow this pattern:
```
GSM[accession]_snRNA_[Group]_[Sample]_MotorCortex_cellBender_corrected_filtered.h5
```

Example:
- `GSM6781907`: GEO accession
- `C9ALS`: Patient group
- `A1`: Sample identifier
- `MotorCortex`: Tissue type
- `cellBender_corrected_filtered`: Processing status

---

## Batch Information

### Sequencing Batches

Samples were sequenced across 6 batches. Batch effects will be corrected using Harmony integration in the analysis pipeline.

| Batch | C9-ALS Samples | Control Samples |
|-------|----------------|-----------------|
| Batch_1 | C9ALS_Patient4 | Control_1 |
| Batch_2 | C9ALS_Patient3 | Control_6 |
| Batch_3 | C9ALS_Patient6 | Control_2 |
| Batch_4 | C9ALS_Patient2 | Control_5 |
| Batch_5 | C9ALS_Patient1 | Control_3 |
| Batch_6 | C9ALS_Patient5 | Control_4 |

**Note:** Each batch contains one C9-ALS and one Control sample, providing balanced experimental design.

---

## Usage in Analysis Pipeline

### Step 1: Prepare Metadata

Ensure your `sample_metadata.csv` is located at:
```
data/metadata/sample_metadata.csv
```

### Step 2: Update File Paths

If your H5 files are in a different location, update the `h5_path` column:

```bash
# Example: If files are in /path/to/data/
sed -i 's|data/h5_files/|/path/to/data/|g' data/metadata/sample_metadata.csv
```

### Step 3: Run Pipeline

The pipeline scripts will automatically:
1. Read sample information from metadata file
2. Load H5 files based on `h5_path` column
3. Group samples by `group` column (C9-ALS vs Control)
4. Account for sequencing batch in Harmony integration

---

## Quality Expectations

Based on public data standards, expect:

- **Cells per sample:** 3,000 - 15,000
- **Median genes per cell:** 1,000 - 3,000
- **Median UMIs per cell:** 2,000 - 6,000
- **Mitochondrial %:** 2-8% (after CellBender correction)

---

## Citation

If using this dataset, please cite the original publication:

**Divergent impacts of C9orf72 repeat expansion on neurons and glia in ALS and FTD**

And the data repository:
- GEO Series: GSE219281

---

## Metadata Verification Checklist

Before running the pipeline, verify:

- [ ] All H5 files are downloaded and accessible
- [ ] File paths in `h5_path` column are correct
- [ ] Metadata file is properly formatted (CSV)
- [ ] No missing values in required columns
- [ ] Sample names are unique
- [ ] Group labels are consistent (C9-ALS vs Control)
- [ ] Sex values are valid (M or F)
- [ ] Ages are numeric
- [ ] Batch information is present

---

## Troubleshooting

### Issue: H5 files not found

**Solution:**
1. Check file paths in metadata
2. Verify files are downloaded
3. Ensure correct file permissions
4. Use absolute paths if needed

### Issue: Sample names mismatch

**Solution:**
1. Verify sample column matches across all pipeline steps
2. Ensure no special characters in sample names
3. Check for trailing spaces or line endings

### Issue: Group labels inconsistent

**Solution:**
1. Use exactly "C9-ALS" and "Control" (case-sensitive)
2. Check for typos or extra spaces
3. Validate with: `cut -d',' -f6 sample_metadata.csv | sort | uniq`

---

## Additional Resources

- **GEO Series Page:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219281
- **Original Publication:** Divergent impacts of C9orf72 repeat expansion on neurons and glia in ALS and FTD
- **CellBender Documentation:** https://github.com/broadinstitute/CellBender
- **10x Genomics File Formats:** https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview

---

**Last Updated:** November 2025  
**Dataset Version:** CellBender-corrected filtered matrices  
**Compatible with:** C9 snRNAseq Analysis Pipeline v1.0

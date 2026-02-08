# Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This folder contains scripts and documentation for running **single-cell RNA-seq preprocessing, QC, integration, and downstream cell-type–specific analyses** on the Minerva HPC cluster.

The pipeline follows a strict, staged order from raw FASTQs to integrated objects and targeted analyses (e.g., microglia).

---

## Pipeline overview (exact order)

1. **Cell Ranger** — alignment and gene-level quantification  
2. **CellBender** — ambient/background RNA removal  
3. **QC and filtering** (per-sample)  
4. **Seurat + Harmony integration, clustering, and annotation**  
5. **Cell-type–specific downstream analyses** (e.g., microglia)

---

Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh
│   ├── run_cellbender_all.sh
│   ├── miBrain_QC_filtering.R
│   ├── QC_filtering_integration_clustering.R
│   ├── Seurat_integration_harmony.R
│   ├── microglia_subset_signatures.R
│   ├── Pseudobulk_Dream_Pathway_analysis.R
│   ├── Pseudobulk_edgeR.R
│   └── Propeller_Camera.R
└── docs/
    └── (analysis notes, troubleshooting, version info)

---
## 1. Cell Ranger

**Script:**  
Single_Cell/scripts/run_cellranger_all.sh

Runs 10x Genomics Cell Ranger for all samples. This step:

- Aligns reads  
- Performs basic QC  
- Produces gene-level count matrices (one per sample)

Example usage (on Minerva):

    ./Single_Cell/scripts/run_cellranger_all.sh

See comments in the script for required input paths and sample lists.

---

## 2. CellBender

**Script:**  
Single_Cell/scripts/run_cellbender_all.sh

Runs CellBender on Cell Ranger outputs to remove ambient/background RNA.

Example usage (on Minerva):

    ./Single_Cell/scripts/run_cellbender_all.sh

This produces cleaned HDF5 matrices (one per sample), which can then be loaded directly into R / Seurat.

---

## 3. QC and filtering (per-sample)

**Script:**  
Single_Cell/scripts/miBrain_QC_filtering.R

This script is run **independently on each sample** after CellBender. It performs:

- Cell-level QC filtering (features, counts, mitochondrial content)
- Removal of low-quality cells
- Creation of QC’d per-sample Seurat objects

This step **does not perform cell-type annotation** and **must be run before integration**.

---

## 4. Seurat + Harmony integration, clustering, and annotation

**Script:**  
Single_Cell/scripts/Seurat_integration_harmony.R

Template script for integrating multiple QC’d single-cell samples using Seurat and Harmony. It:

- Loads a list of QC’d Seurat objects (one per sample)
- Runs SCTransform normalization
- Performs SCT-based integration across samples
- Applies Harmony batch correction on PCA embeddings
- Builds neighbors, clusters, and UMAP embeddings
- Computes cluster markers with FindAllMarkers
- Performs cell-type annotation

Example usage (on Minerva):

    module load R/4.2.0
    Rscript Single_Cell/scripts/Seurat_integration_harmony.R

Before running, edit the header of the script to set:

- qc_list_path — path to the saved list of QC’d per-sample Seurat objects  
- output_dir — directory where integrated objects, plots, and marker tables are written  

---

## 5. Microglia subset and gene signature analysis

**Script:**  
Single_Cell/scripts/microglia_subset_signatures.R

This script is run **after integration** and operates on an integrated Seurat object. It performs:

- Identification of microglia using module scores
- Subsetting and reprocessing of microglia
- Microglial subclustering
- Module scoring and cluster-level summaries
- Genotype enrichment analyses
- Generation of publication-ready plots and heatmaps
- Saving of a microglia-specific Seurat object

---

## 6. Microglia pseudobulk analyses

All scripts in this section operate on **pseudobulk microglia data** derived from the integrated Seurat object. These analyses are performed after microglia identification, subsetting, and reclustering.

---

## 6.1 Pseudobulk differential expression and pathway analysis (DREAM)

**Script:**  
Single_Cell/scripts/Pseudobulk_Dream_Pathway_analysis.R

This script performs microglia pseudobulk differential expression using **DREAM**, accounting for sample-level random effects.

Downstream pathway and gene-set enrichment analyses are performed on DREAM results.

---

## 7 Pseudobulk differential expression (edgeR)

**Script:**  
Single_Cell/scripts/Pseudobulk_edgeR.R

This script constructs microglia pseudobulk expression profiles per sample and performs differential expression using **edgeR**.

It serves as an alternative DEG framework to DREAM and is useful for method comparison and robustness checks.

---

## 8 Abundance testing and gene-set analysis (Propeller + CAMERA)

**Script:**  
Single_Cell/scripts/Propeller_Camera.R

This script performs complementary pseudobulk-based statistical analyses on microglia:

**Abundance testing (Propeller-style limma):**
- Tests cluster-level abundance differences between genotypes
- Tests state-level abundance differences (e.g., Homeostatic vs DAM)
- Uses asin(sqrt(proportion)) transformation
- Applies cell-count–based precision weights

**Gene-set testing (CAMERA):**
- Tests curated microglial gene modules
- Accounts for inter-gene correlation
- Operates on voom-transformed pseudobulk expression data
- Reports directionality and FDR-adjusted significance

All analyses in this script are restricted to microglia and assume pseudobulk inputs.

---








## Docs and notes

Use the Single_Cell/docs/ directory for:

- Version notes (R, Seurat, Harmony, Cell Ranger, CellBender)
- Common error messages and fixes
- Run logs and Minerva-specific configuration details

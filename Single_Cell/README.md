# RNA-seq Pipelines (MiBrain / Minerva)

This repository contains pipelines and documentation for running both **single-cell** and **bulk RNA-seq** analyses on the Minerva HPC cluster.

- `Single_Cell/` — 10x Genomics single-cell RNA-seq preprocessing, QC, integration, and downstream cell-type–specific analyses  
- `Bulk/` — bulk RNA-seq analysis using RAPiD, DESeq2, and pathway analysis (GSEA + ORA)

---

## 1. Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This pipeline covers preprocessing, QC, integration, and downstream analyses for MiBrain single-cell RNA-seq data on Minerva.

### Overview of steps (exact order)

1. **Cell Ranger** — alignment and gene-level quantification  
2. **CellBender** — ambient/background RNA removal  
3. **QC and filtering** (per-sample)  
4. **Seurat + Harmony integration, clustering, and annotation**  
5. **Cell-type–specific downstream analyses** (e.g., microglia)

---

### Folder structure

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh
│   ├── run_cellbender_all.sh
│   ├── miBrain_QC_filtering.R
│   ├── Seurat_integration_harmony.R
│   └── microglia_subset_signatures.R
└── docs/
    └── (analysis notes, troubleshooting, etc.)


1.1 Cell Ranger
Script:
Single_Cell/scripts/run_cellranger_all.sh
Runs 10x Genomics Cell Ranger for all samples in the project. This step:
Aligns reads
Performs basic QC
Generates gene-level count matrices (per sample)
Example usage (on Minerva):
./Single_Cell/scripts/run_cellranger_all.sh
(See comments inside the script for required input paths and sample lists.)
1.2 CellBender
Script:
Single_Cell/scripts/run_cellbender_all.sh
Runs CellBender on the Cell Ranger outputs to remove ambient/background RNA.
Example usage (on Minerva):
./Single_Cell/scripts/run_cellbender_all.sh
This produces cleaned HDF5 matrices (one per sample), which can then be loaded directly into R / Seurat.
1.3 QC and filtering (per-sample)
Script:
Single_Cell/scripts/miBrain_QC_filtering.R
This script is run on individual samples after CellBender and performs:
Cell-level QC filtering (features, counts, mitochondrial content)
Removal of low-quality cells
Creation of QC’d per-sample Seurat objects
This step does not perform cell-type annotation and must be run before integration.
1.4 Seurat + Harmony integration, clustering, and annotation
Script:
Single_Cell/scripts/Seurat_integration_harmony.R
This script integrates multiple QC’d single-cell samples using Seurat and Harmony. It:
Loads a list of QC’d Seurat objects (one per sample)
Runs SCTransform normalization
Performs SCT-based integration across samples
Applies Harmony batch correction on PCA embeddings
Computes neighbors, clusters, and UMAP
Finds cluster markers with FindAllMarkers
Performs cell-type annotation
Example usage (on Minerva):
module load R/4.2.0
Rscript Single_Cell/scripts/Seurat_integration_harmony.R
Before running, edit the top of Seurat_integration_harmony.R and set:
qc_list_path – path to the saved list of QC’d per-sample Seurat objects
output_dir – directory where integrated objects, plots, and marker tables will be written
1.5 Microglia subset and gene signature analysis
Script:
Single_Cell/scripts/microglia_subset_signatures.R
This script is run after integration and operates on an integrated Seurat object. It performs:
Identification of microglia using module scores
Subsetting and reprocessing of microglia
Microglial subclustering
Module scoring and cluster-level summaries
Generation of figures and heatmaps
Saving of a microglia-specific Seurat object
1.6 Docs and notes (Single-cell)
Additional notes, version info, and troubleshooting tips for the single-cell pipeline live under:
Single_Cell/docs/
Suggested files:
Seurat_integration_harmony_notes.md
cellranger_cellbender_pipeline_notes.md
1.7 Example MiBrain project layout on Minerva
For reference, a typical project directory for MiBrain P2RY12 single-cell data:
/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
├── 00_fastq/
├── 01_pilot_fastq/
├── 02_cellranger/
├── 03_cellbender/
├── 04_seurat_objects/
├── logs/
└── metadata/
Convenience commands:
./Single_Cell/scripts/run_cellranger_all.sh
./Single_Cell/scripts/run_cellbender_all.sh
2. Bulk RNA-seq Pipeline Overview
The bulk RNA-seq pipeline is organized into four main stages:
RAPiD — preprocessing, QC, alignment, and quantification
Collate samples (coming soon) — combine featureCounts outputs
DEG analysis — differential expression with DESeq2
Pathway analysis — GSEA + ORA
Detailed documentation for each stage lives in Bulk/docs/.
2.1 Repository structure (Bulk)
Bulk/
├── scripts/
│   ├── 00_run_RAPiD.sh
│   ├── 01_collate_samples.R
│   ├── BULK_rnaseq_DEG.R
│   └── BULK_rnaseq_pathway_1.R
└── docs/
    ├── RAPiD_details.md
    ├── RAPiD_tips_jen.md
    └── DEG_pathway_details.md
2.2 RAPiD — preprocessing / alignment / QC
RAPiD performs:
Adapter trimming (Trimmomatic)
QC (FastQC, Picard)
Alignment (STAR)
Gene quantification (featureCounts, Salmon, Kallisto, RSEM)
Splicing analysis (LeafCutter)
Aggregated QC reports (MultiQC)
See:
Bulk/docs/RAPiD_details.md
Bulk/docs/RAPiD_tips_jen.md
2.3 Collate samples (coming soon)
Planned stage to:
Read RAPiD featureCounts outputs
Combine them into a unified gene × sample matrix
Prepare inputs for DESeq2
Planned locations:
Bulk/scripts/01_collate_samples.R
Bulk/docs/collate_samples_details.md
2.4 DEG analysis — DESeq2
Script:
Bulk/scripts/BULK_rnaseq_DEG.R
Takes as input:
featureCounts matrix
Sample metadata
Produces:
results_raw.rds
results_shrunken_apeglm.rds
QC plots and result tables
See Bulk/docs/DEG_pathway_details.md.
2.5 Pathway analysis — GSEA + ORA
Script:
Bulk/scripts/BULK_rnaseq_pathway_1.R
Uses DESeq2 output to perform:
GSEA (fgseaMultilevel)
ORA (KEGG, GO BP, Hallmark) via clusterProfiler
Standard enrichment plots
2.6 Bulk documentation summary
All bulk RNA-seq documentation lives under Bulk/docs/.

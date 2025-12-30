# RNA-seq Pipelines (MiBrain / Minerva)

This repository contains pipelines and documentation for running **single-cell** and **bulk RNA-seq** analyses for MiBrain projects on the Minerva HPC cluster.

The repository is organized into two primary workflows:

- **Single-cell RNA-seq** — 10x Genomics data preprocessing, QC, integration, and downstream cell-type–specific analyses  
- **Bulk RNA-seq** — alignment, quantification, differential expression, and pathway analysis  

---

# 1. Single-cell RNA-seq Pipeline (MiBrain / Minerva)

The single-cell pipeline covers preprocessing, QC, integration, clustering, annotation, and downstream analyses for MiBrain scRNA-seq data.

## 1.1 Pipeline overview (exact order)

1. **Cell Ranger** — alignment and gene-level quantification  
2. **CellBender** — ambient/background RNA removal  
3. **QC and filtering** (per-sample)  
4. **Seurat + Harmony integration, clustering, and annotation**  
5. **Cell-type–specific downstream analyses** (e.g., microglia)

---

## 1.2 Folder structure (Single-cell)

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh
│   ├── run_cellbender_all.sh
│   ├── miBrain_QC_filtering.R
│   ├── Seurat_integration_harmony.R
│   └── microglia_subset_signatures.R
└── docs/
    └── (analysis notes, troubleshooting, version info)



1.3 Key scripts (Single-cell)
Cell Ranger
Script: Single_Cell/scripts/run_cellranger_all.sh
Runs 10x Genomics Cell Ranger for all samples:
Alignment
Basic QC
Gene-level count matrices (per sample)
CellBender
Script: Single_Cell/scripts/run_cellbender_all.sh
Removes ambient/background RNA from Cell Ranger outputs, producing cleaned HDF5 matrices.
QC and filtering (per-sample)
Script: Single_Cell/scripts/miBrain_QC_filtering.R
Cell-level QC filtering (features, counts, mitochondrial content)
Removal of low-quality cells
Outputs QC’d per-sample Seurat objects
No cell-type annotation performed at this stage
Seurat + Harmony integration and annotation
Script: Single_Cell/scripts/Seurat_integration_harmony.R
SCTransform normalization
Sample integration
Harmony batch correction
Clustering and UMAP
Marker detection
Cell-type annotation
Microglia subset and gene signature analysis
Script: Single_Cell/scripts/microglia_subset_signatures.R
Downstream analysis performed after integration, including:
Identification and subsetting of microglia
Microglial reprocessing and subclustering
Module scoring and cluster-level summaries
Genotype enrichment analyses
Publication-ready plots and heatmaps
Saving of a microglia-specific Seurat object

1.4 Example MiBrain single-cell project layout (Minerva)
/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
├── 00_fastq/          # raw FASTQs
├── 01_pilot_fastq/    # pilot FASTQs / symlinks
├── 02_cellranger/     # Cell Ranger outputs
├── 03_cellbender/     # CellBender outputs
├── 04_seurat_objects/ # downstream R analysis
├── logs/              # LSF logs
└── metadata/          # sample sheets, etc.








# 2. **Bulk RNA-seq Pipeline (MiBrain / Minerva)**

The bulk RNA-seq pipeline supports alignment, quantification, differential expression, and pathway analysis for bulk RNA-seq datasets.

---

## 2.1 **Pipeline overview**

1. **RAPiD** — preprocessing, QC, alignment, and quantification  
2. **Collate samples** (planned) — combine count matrices  
3. **Differential expression analysis** — DESeq2  
4. **Pathway analysis** — GSEA and ORA  

---

## 2.2 **Folder structure (Bulk)**

```text
Bulk/
├── scripts/
│   ├── 00_run_RAPiD.sh
│   ├── 01_collate_samples.R         # coming soon
│   ├── BULK_rnaseq_DEG.R
│   └── BULK_rnaseq_pathway_1.R
└── docs/
    ├── RAPiD_details.md
    ├── RAPiD_tips_jen.md
    └── DEG_pathway_details.md



2.3 RAPiD — preprocessing / alignment / QC
RAPiD performs:
Adapter trimming (Trimmomatic)
QC (FastQC, Picard)
Alignment (STAR)
Gene quantification (featureCounts, Salmon, Kallisto, RSEM)
Splicing analysis (LeafCutter)
Aggregated QC reports (MultiQC)
Documentation:
Bulk/docs/RAPiD_details.md
Bulk/docs/RAPiD_tips_jen.md
2.4 Differential expression analysis (DESeq2)
Script: Bulk/scripts/BULK_rnaseq_DEG.R
Inputs:
featureCounts matrix
Sample metadata
Outputs:
Raw and shrunken DEG results
QC plots and result tables
2.5 Pathway analysis (GSEA + ORA)
Script: Bulk/scripts/BULK_rnaseq_pathway_1.R
GSEA (fgseaMultilevel)
Over-representation analysis (KEGG, GO, Hallmark)
Standard enrichment plots
2.6 Bulk documentation summary
All bulk RNA-seq documentation lives under:
Bulk/docs/
Notes
This repository is designed for use on the Minerva HPC cluster
Script headers contain required input paths and configuration details
Documentation is intentionally modular so pipelines can be adapted to new projects

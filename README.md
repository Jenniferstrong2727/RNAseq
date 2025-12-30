# RNA-seq Pipelines (MiBrain / Minerva)

This repository contains pipelines and documentation for running **bulk** and **single-cell RNA-seq** analyses for MiBrain projects on the Minerva HPC cluster.

The repository is organized into two primary workflows aligned with project aims:

- **Aim 1 — Bulk RNA-seq**: alignment, quantification, differential expression, and pathway analysis  
- **Aim 2 — Single-cell RNA-seq**: preprocessing, QC, integration, and downstream cell-type–specific analyses  

---

# 1. Bulk RNA-seq Pipeline (MiBrain / Minerva)

The bulk RNA-seq pipeline supports alignment, quantification, differential expression, and pathway analysis for bulk RNA-seq datasets.

---

## 1.1 Pipeline overview

1. **RAPiD** — preprocessing, QC, alignment, and quantification  
2. **Collate samples** (planned) — combine count matrices  
3. **Differential expression analysis** — DESeq2  
4. **Pathway analysis** — GSEA and ORA  

---

## 1.2 Folder structure (Bulk)

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

---

## 1.3 RAPiD — preprocessing / alignment / QC

RAPiD performs:
- Adapter trimming (Trimmomatic)
- QC (FastQC, Picard)
- Alignment (STAR)
- Gene quantification (featureCounts, Salmon, Kallisto, RSEM)
- Splicing analysis (LeafCutter)
- Aggregated QC reports (MultiQC)

Documentation:
- Bulk/docs/RAPiD_details.md  
- Bulk/docs/RAPiD_tips_jen.md  

---

## 1.4 Differential expression analysis (DESeq2)

**Script:**  
Bulk/scripts/BULK_rnaseq_DEG.R

Inputs:
- featureCounts matrix  
- Sample metadata  

Outputs:
- Raw DEG results  
- Shrunk logFC results  
- QC plots and result tables  

---

## 1.5 Pathway analysis (GSEA + ORA)

**Script:**  
Bulk/scripts/BULK_rnaseq_pathway_1.R

Performs:
- GSEA (fgseaMultilevel)
- Over-representation analysis (KEGG, GO, Hallmark)
- Standard enrichment plots

---

# 2. Single-cell RNA-seq Pipeline (MiBrain / Minerva)

**Preprocessing, QC, integration, clustering, annotation, and downstream analyses for MiBrain scRNA-seq data**

---

## 2.1 Pipeline overview (exact order)

1. **Cell Ranger** — alignment and gene-level quantification  
2. **CellBender** — ambient/background RNA removal  
3. **QC and filtering** (per-sample)  
4. **Seurat + Harmony integration, clustering, and annotation**  
5. **Cell-type–specific downstream analyses** (e.g., microglia)

---

## 2.2 Folder structure (Single-cell)

    Single_Cell/
    ├── scripts/
    │   ├── run_cellranger_all.sh
    │   ├── run_cellbender_all.sh
    │   ├── miBrain_QC_filtering.R
    │   ├── Seurat_integration_harmony.R
    │   └── microglia_subset_signatures.R
    └── docs/
        └── (analysis notes, troubleshooting, version info)

---

## 2.3 Key scripts (Single-cell)

### Cell Ranger

**Script:**  
Single_Cell/scripts/run_cellranger_all.sh

- Alignment  
- Basic QC  
- Gene-level count matrices (per sample)

---

### CellBender

**Script:**  
Single_Cell/scripts/run_cellbender_all.sh

Removes ambient/background RNA from Cell Ranger outputs.

---

### QC and filtering (per-sample)

**Script:**  
Single_Cell/scripts/miBrain_QC_filtering.R

- Cell-level QC filtering (features, counts, mitochondrial content)  
- Removal of low-quality cells  
- Outputs QC’d per-sample Seurat objects  
- **No cell-type annotation performed at this stage**

---

### Seurat + Harmony integration and annotation

**Script:**  
Single_Cell/scripts/Seurat_integration_harmony.R

- SCTransform normalization  
- Sample integration  
- Harmony batch correction  
- Clustering and UMAP  
- Marker detection  
- Cell-type annotation  

---

### Microglia subset and gene signature analysis

**Script:**  
Single_Cell/scripts/microglia_subset_signatures.R

Downstream analysis performed **after integration**, including:
- Identification and subsetting of microglia  
- Microglial reprocessing and subclustering  
- Module scoring and cluster-level summaries  
- Genotype enrichment analyses  
- Publication-ready plots and heatmaps  
- Saving of a microglia-specific Seurat object  

---

## 2.4 Example MiBrain single-cell project layout (Minerva)

    /sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
    ├── 00_fastq/
    ├── 01_pilot_fastq/
    ├── 02_cellranger/
    ├── 03_cellbender/
    ├── 04_seurat_objects/
    ├── logs/
    └── metadata/

---

## Notes

- This repository is designed for use on the **Minerva HPC cluster**
- Script headers contain required input paths and configuration details
- Documentation is intentionally modular so pipelines can be adapted to new projects

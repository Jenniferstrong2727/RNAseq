# RNA-seq Pipelines (MiBrain / Minerva)

This repository contains end-to-end pipelines for **single-cell** and **bulk RNA-seq** analyses performed on the **Minerva HPC cluster**, with a focus on the **P2RY12 MiBrain project**.

---

## Repository Overview

* `Single_Cell/` — MiBrain single-cell RNA-seq preprocessing, integration, and downstream analyses
* `Bulk/` — bulk RNA-seq processing, differential expression, and pathway analysis

---

# 1. Single-cell RNA-seq Pipeline 

This pipeline processes raw 10x Genomics data through integration and **cell-type–specific analyses**, with modular scripts for reproducibility and scalability.

---

## Pipeline Overview 

1. **Cell Ranger** — alignment and gene quantification
2. **CellBender** — ambient RNA correction
3. **QC & filtering** — per-sample Seurat objects
4. **Integration (Seurat + Harmony)** — clustering and annotation
5. **Downstream analyses** — cell-type–specific (e.g., microglia)
6. **Pseudobulk + pathway analysis**
7. **Pilot vs postmortem comparison**

---

## Folder Structure

```text
Single_Cell/
└── P2RY12/
    ├── New_miBrain/
    │   └── SCRIPTS/
    │       ├── run_cellranger_all.sh
    │       ├── run_cellbender_all.sh
    │       ├── miBrain_QC_filtering.R
    │       ├── Seurat_integration_harmony.R
    │       ├── microglia_subset_signatures.R
    │       ├── Pseudobulk_edgeR.R
    │       ├── Pseudobulk_Dream_Pathway_analysis.R
    │       └── Propeller_Camera.R
    │
    └── pilot_postmortem_comparison/
        ├── Merged_pilot_analysis.R
        ├── PM_abundance.R
        ├── PM_pilot_abundance_comparison.R
        ├── PM_pilot_DEG_correlation_concordance.R
        └── state_specific_DEG.R
```

---

## Scripts (New_miBrain)

### 1. Cell Ranger

**Script:** `run_cellranger_all.sh`

* Aligns FASTQ files
* Generates gene × cell count matrices
* Performs initial QC

**Run:**

```bash
./run_cellranger_all.sh
```

---

### 2. CellBender

**Script:** `run_cellbender_all.sh`

* Removes ambient/background RNA contamination
* Produces corrected `.h5` matrices

```bash
./run_cellbender_all.sh
```

---

### 3. QC and Filtering

**Script:** `miBrain_QC_filtering.R`

* Filters low-quality cells
* Applies thresholds (features, counts, mitochondrial %)
* Outputs per-sample Seurat objects

---

### 4. Integration, Clustering, Annotation

**Script:** `Seurat_integration_harmony.R`

* SCTransform normalization
* Integration across samples
* Harmony batch correction
* UMAP + clustering
* Marker detection (`FindAllMarkers`)
* Cell-type annotation

**Run:**

```bash
module load R/4.2.0
Rscript Seurat_integration_harmony.R
```

**Edit inside script:**

* `qc_list_path`
* `output_dir`

---

### 5. Microglia Subclustering & Signatures

**Script:** `microglia_subset_signatures.R`

* Identifies microglia via module scoring
* Subsets and reclusters microglia
* Computes functional gene modules
* Generates figures and summary tables

---

### 6. Pseudobulk Differential Expression

**Script:** `Pseudobulk_edgeR.R`

* Aggregates counts per sample/cell type
* Performs DEG analysis using **edgeR**
* Enables genotype comparisons (e.g., G/G vs A/A)

---

### 7. Pathway Analysis (Dream / GSEA)

**Script:** `Pseudobulk_Dream_Pathway_analysis.R`

* Linear modeling with **dream (variancePartition)**
* Gene set enrichment (GSEA/GO/KEGG)
* Identifies pathway-level differences

---

### 8. Cell Type Composition Analysis

**Script:** `Propeller_Camera.R`

* Differential abundance testing (propeller)
* Gene set testing (camera)
* Evaluates compositional changes across conditions

---

## Pilot vs Postmortem Comparison

**Location:**
`Single_Cell/P2RY12/pilot_postmortem_comparison/`

### Key Analyses

**Data Integration**

* `Merged_pilot_analysis.R`

  * Merges pilot + postmortem datasets
  * Harmonizes metadata and annotations

**Cell Type Abundance**

* `PM_abundance.R`
* `PM_pilot_abundance_comparison.R`

  * Compares cell-type proportions
  * Identifies shifts across datasets

**DEG Concordance**

* `PM_pilot_DEG_correlation_concordance.R`

  * Correlates DEG signatures between MiBrain and postmortem datasets

**State-Specific DEG Analysis**

* `state_specific_DEG.R`

  * Identifies DEGs within specific cell states
  * Enables cross-system validation

---

## Example Project Layout (Minerva)

```text
/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
├── 00_fastq/
├── 01_pilot_fastq/
├── 02_cellranger/
├── 03_cellbender/
├── 04_seurat_objects/
├── logs/
└── metadata/
```

---

## Bulk RNA-seq Pipeline

### Workflow

1. RAPiD — preprocessing, QC, alignment
2. Collation — combine count matrices
3. DEG analysis — DESeq2
4. Pathway analysis — GSEA + ORA

---

### Folder Structure

```text
Bulk/
├── scripts/
│   ├── 00_run_RAPiD.sh
│   ├── 01_collate_samples.R
│   ├── BULK_rnaseq_DEG.R
│   └── BULK_rnaseq_pathway_1.R
└── docs/
```

---

### Key Tools

* Alignment: STAR
* Quantification: featureCounts, Salmon, Kallisto
* DEG: DESeq2
* Pathways: fgsea, clusterProfiler

---

## Documentation

* `Bulk/docs/`
* Recommended additions:

  * `integration_notes.md`
  * `cellbender_qc_notes.md`

---

## Notes

* Designed for **Minerva HPC execution**
* Paths should be customized per project
* Pipeline is modular and stepwise
* Optimized for **MiBrain P2RY12 variant analysis and cross-system validation**


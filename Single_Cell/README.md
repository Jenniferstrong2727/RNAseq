# Single-Cell RNA-seq Pipeline (MiBrain / Minerva)

This directory contains the **single-cell RNA-seq analysis pipeline** for the MiBrain P2RY12 project, including preprocessing, integration, and downstream cell-type–specific analyses.

---

## Pipeline Overview

This workflow processes raw 10x Genomics data through integration and biological interpretation:

1. **Cell Ranger** — alignment and gene quantification
2. **CellBender** — ambient RNA correction
3. **QC & filtering** — per-sample Seurat objects
4. **Integration (Seurat + Harmony)** — clustering and annotation
5. **Cell-type–specific analyses** (e.g., microglia)
6. **Pseudobulk differential expression and pathway analysis**
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

## Core Pipeline (New_miBrain)

### Preprocessing

**Cell Ranger** (`run_cellranger_all.sh`)

* Alignment and count matrix generation

**CellBender** (`run_cellbender_all.sh`)

* Ambient RNA removal

---

### QC and Integration

**QC Filtering** (`miBrain_QC_filtering.R`)

* Filters low-quality cells
* Outputs per-sample Seurat objects

**Integration** (`Seurat_integration_harmony.R`)

* SCTransform normalization
* Harmony batch correction
* Clustering, UMAP, and annotation

**Run:**

```bash
module load R/4.2.0
Rscript Seurat_integration_harmony.R
```

---

### Cell-type–specific Analysis

**Microglia Analysis** (`microglia_subset_signatures.R`)

* Subclustering and module scoring
* Functional annotation

---

### Pseudobulk & Pathway Analysis

**DEG Analysis** (`Pseudobulk_edgeR.R`)

* Pseudobulk aggregation
* Differential expression (edgeR)

**Pathway Analysis** (`Pseudobulk_Dream_Pathway_analysis.R`)

* Linear modeling (dream)
* GSEA / GO / KEGG

**Composition Analysis** (`Propeller_Camera.R`)

* Cell-type proportion testing
* Gene set testing

---

## Pilot vs Postmortem Comparison

Location:
`Single_Cell/P2RY12/pilot_postmortem_comparison/`

### Analyses

* **Merged_pilot_analysis.R** — dataset integration
* **PM_abundance.R** — cell-type proportions in postmortem data
* **PM_pilot_abundance_comparison.R** — comparison of cell-type proportions between MiBrain and postmortem datasets
* **PM_pilot_DEG_correlation_concordance.R** — DEG concordance between MiBrain and postmortem datasets
* **state_specific_DEG.R** — state-level differential expression across datasets

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

## Notes

* Designed for **Minerva HPC execution**
* Paths should be customized per project
* Pipeline is modular and stepwise
* Optimized for **MiBrain P2RY12 variant analysis and cross-system validation**

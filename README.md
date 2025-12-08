# RNA-seq Pipelines (MiBrain / Minerva)

This repository contains pipelines and documentation for running both **single-cell** and **bulk** RNA-seq analysis on the Minerva HPC cluster.

- `Single_Cell/` — 10x Genomics single-cell RNA-seq preprocessing + integration (Cell Ranger → CellBender → Seurat + Harmony).
- `Bulk/` — bulk RNA-seq analysis using RAPiD, DESeq2, and pathway analysis (GSEA + ORA).

---

## 1. Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This pipeline covers preprocessing and integration for MiBrain single-cell RNA-seq data on Minerva.

Currently includes:

- **Cell Ranger** (10x Genomics) for alignment and gene-level quantification  
- **CellBender** (Broad) for removal of ambient/background RNA  
- **Seurat + Harmony** for downstream integration, clustering, marker detection, and simple module scores

### Folder structure

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh        # submit Cell Ranger jobs for all samples
│   ├── run_cellbender_all.sh        # submit CellBender jobs for all samples
│   └── Seurat_integration_harmony.R # Seurat + Harmony integration, clustering, markers
└── docs/
    └── (analysis notes, troubleshooting, etc.)
1.1 Cell Ranger
Script: Single_Cell/scripts/run_cellranger_all.sh
Runs 10x Genomics Cell Ranger for all samples in the project. This step:

Aligns reads
Performs basic QC
Generates gene-level count matrices (per sample)
Example usage (on Minerva):
# From your project directory on Minerva
./Single_Cell/scripts/run_cellranger_all.sh
(See comments inside the script for required input paths and sample lists.)
1.2 CellBender
Script: Single_Cell/scripts/run_cellbender_all.sh
Runs CellBender on the Cell Ranger outputs to remove ambient/background RNA.

Example usage (on Minerva):

# From your project directory on Minerva
./Single_Cell/scripts/run_cellbender_all.sh
This produces cleaned HDF5 matrices (one per sample), which can then be loaded into R / Seurat.
1.3 Seurat + Harmony Integration
Script: Single_Cell/scripts/Seurat_integration_harmony.R
This script is a general template for integrating multiple scRNA-seq samples after CellBender using Seurat and Harmony. It:

Loads a list of per-sample QC’d Seurat objects (one per sample)
Runs SCTransform normalization
Performs SCT-based integration across samples
Applies Harmony batch correction on the PCA embeddings
Computes neighbors, clusters, and UMAP
Finds cluster markers with FindAllMarkers
Optionally computes simple module scores for broad cell types
(e.g., neurons, astrocytes, microglia)
Example usage (on Minerva):
module load R/4.2.0

Rscript Single_Cell/scripts/Seurat_integration_harmony.R
Before running, edit the top of Seurat_integration_harmony.R and set:
qc_list_path – path to your saved Seurat list of QC’d samples
output_dir – directory where integrated objects, plots, markers, and module scores will be saved
1.4 Docs and Notes (Single-cell)
Additional notes, version info, and troubleshooting tips for the single-cell pipeline live under:
Single_Cell/docs/
Suggested files:
Seurat_integration_harmony_notes.md – R / Seurat / Harmony versions, common errors, and fixes
cellranger_cellbender_pipeline_notes.md – run logs, flags, and cluster config details
1.5 Example MiBrain project layout on Minerva
For reference, a typical project directory for the MiBrain P2RY12 single-cell data:
/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
├── 00_fastq/          # new MiBrain FASTQs (RunMiBrain-1..8)
├── 01_pilot_fastq/    # symlinks to pilot FASTQs (miBrain-1, miBrain-2)
├── 02_cellranger/     # Cell Ranger outputs (one folder per sample)
├── 03_cellbender/     # CellBender outputs (one folder per sample)
├── 04_seurat_objects/ # downstream analysis in R
├── logs/              # LSF logs for Cell Ranger / CellBender
└── metadata/          # sample sheets, etc.
Convenience commands:
# Submit all Cell Ranger jobs
./Single_Cell/scripts/run_cellranger_all.sh

# After Cell Ranger finishes (or partially finishes):
./Single_Cell/scripts/run_cellbender_all.sh
2. Bulk RNA-seq Pipeline Overview
The bulk RNA-seq pipeline is organized into four main stages:
RAPiD — preprocessing, QC, alignment, and quantification
Collate Samples (coming soon) — combine featureCounts outputs into a unified matrix
DEG Analysis — differential expression with DESeq2
Pathway Analysis — GSEA + ORA using fgsea and clusterProfiler
Detailed documentation for each stage lives in the Bulk/docs/ directory.
2.1 Repository structure (Bulk)
Bulk/
├── scripts/
│   ├── 00_run_RAPiD.sh              # run RAPiD preprocessing
│   ├── 01_collate_samples.R         # (coming soon) create combined count matrix
│   ├── BULK_rnaseq_DEG.R            # differential expression (DESeq2)
│   └── BULK_rnaseq_pathway_1.R      # GSEA + ORA pathway analysis
│
└── docs/
    ├── RAPiD_details.md             # what RAPiD does and key outputs
    ├── RAPiD_tips_jen.md            # Jen’s practical tips for running RAPiD on Minerva
    └── DEG_pathway_details.md       # combined DEG + pathway documentation
2.2 RAPiD — Preprocessing / Alignment / QC
RAPiD performs:
Adapter trimming (Trimmomatic)
QC (FastQC, Picard)
Alignment (STAR)
Gene quantification (featureCounts, Salmon, Kallisto, RSEM)
Splicing analysis (LeafCutter)
Aggregated QC reports (MultiQC HTML)
Further details:
Bulk/docs/RAPiD_details.md — full description of each step
Bulk/docs/RAPiD_tips_jen.md — practical tips for running RAPiD on Minerva
2.3 Collate Samples — (Coming soon)
This stage will:
Read RAPiD featureCounts outputs
Combine them into a unified gene × sample matrix (.tsv)
Standardize column names
Prepare the counts matrix for DESeq2
Planned locations:
Script: Bulk/scripts/01_collate_samples.R
Docs: Bulk/docs/collate_samples_details.md
2.4 DEG Analysis — Differential Expression (DESeq2)
Script: Bulk/scripts/BULK_rnaseq_DEG.R
Takes as input:

featureCounts matrix
Sample metadata (sample, genotype, handler, etc.)
Produces:
results_raw.rds
results_shrunken_apeglm.rds
QC plots and top tables
Details and examples: Bulk/docs/DEG_pathway_details.md.
2.5 Pathway Analysis — GSEA + ORA
Script: Bulk/scripts/BULK_rnaseq_pathway_1.R
Uses the DESeq2 output to perform:

GSEA using fgseaMultilevel
Over-representation analysis (ORA) for KEGG, GO BP, Hallmark, etc., via clusterProfiler
Standard dot plots and enrichment plots
More documentation: Bulk/docs/DEG_pathway_details.md.
2.6 Bulk documentation summary
All bulk RNA-seq docs:
Bulk/docs/RAPiD_details.md — full RAPiD step descriptions
Bulk/docs/RAPiD_tips_jen.md — Jen’s RAPiD tips on Minerva
Bulk/docs/DEG_pathway_details.md — DESeq2 + GSEA/ORA workflow
(planned) Bulk/docs/collate_samples_details.md — collating featureCounts outputs

---

### B. New `Single_Cell/README.md`

This one should basically be the **single-cell section only**, so users who click into that folder still see the full description.

In `Single_Cell/README.md`, replace everything with:

```markdown
# Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This folder contains scripts and documentation for running the **single-cell RNA-seq preprocessing and integration pipeline** on the Minerva HPC cluster.

Currently includes:

- **Cell Ranger** (10x Genomics) for alignment and gene quantification  
- **CellBender** (Broad) for ambient/background RNA removal  
- **Seurat + Harmony** for downstream integration, clustering, marker detection, and simple module scores

---

## Folder structure

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh        # submit Cell Ranger jobs for all samples
│   ├── run_cellbender_all.sh        # submit CellBender jobs for all samples
│   └── Seurat_integration_harmony.R # Seurat + Harmony integration, clustering, markers
└── docs/
    └── (analysis notes, troubleshooting, etc.)
1. Cell Ranger
Script: Single_Cell/scripts/run_cellranger_all.sh
Runs 10x Genomics Cell Ranger for all samples. This step:

Aligns reads
Performs QC
Produces gene-level count matrices for each sample
Example usage (on Minerva):
# From your project directory on Minerva
./Single_Cell/scripts/run_cellranger_all.sh
See comments in the script for required input paths and sample lists.
2. CellBender
Script: Single_Cell/scripts/run_cellbender_all.sh
Runs CellBender on Cell Ranger outputs to remove ambient/background RNA.

Example usage (on Minerva):

# From your project directory on Minerva
./Single_Cell/scripts/run_cellbender_all.sh
Outputs cleaned HDF5 matrices (one per sample) that can be loaded directly into R / Seurat.
3. Seurat + Harmony Integration
Script: Single_Cell/scripts/Seurat_integration_harmony.R
Template script for integrating multiple scRNA-seq samples after CellBender using Seurat and Harmony. It:

Loads a list of QC’d Seurat objects (one per sample)
Runs SCTransform normalization
Performs SCT-based integration across samples
Applies Harmony to correct batch effects on the PCA space
Builds neighbors, clusters, and UMAP embeddings
Computes cluster markers with FindAllMarkers
Optionally computes simple module scores for broad cell types
(e.g., neurons, astrocytes, microglia)
Example usage:
module load R/4.2.0

Rscript Single_Cell/scripts/Seurat_integration_harmony.R
Before running, edit the header of the script to set:
qc_list_path – path to the saved list of per-sample QC’d Seurat objects
output_dir – where integrated objects, plots, markers, and module score outputs are written
4. Docs and notes
Use Single_Cell/docs/ for:
Version notes (R, Seurat, Harmony, CellBender, Cell Ranger)
Common error messages and fixes
Run logs and LSF config details for Minerva


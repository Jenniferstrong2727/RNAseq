# Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This folder contains scripts and docs for running the single-cell RNA-seq preprocessing and integration pipeline on the Minerva HPC cluster.

Currently includes:

- **Cell Ranger** (10x Genomics) for alignment + quantification  
- **CellBender** (Broad) for ambient RNA/background removal  
- **Seurat + Harmony** for downstream integration, clustering, and marker detection  

---

## Folder structure

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh        # Submit Cell Ranger jobs for all samples
│   ├── run_cellbender_all.sh        # Submit CellBender jobs for all samples
│   └── Seurat_integration_harmony.R # Seurat + Harmony integration, clustering, markers
└── docs/
    └── (analysis notes, troubleshooting, etc.)





1. Cell Ranger
Script:
Single_Cell/scripts/run_cellranger_all.sh
Runs 10x Genomics Cell Ranger for all samples in the project.
This step performs alignment and gene-level quantification.
Example usage:
# From your project directory on Minerva
./Single_Cell/scripts/run_cellranger_all.sh
(See comments inside the script for required input paths and sample lists.)
2. CellBender
Script:
Single_Cell/scripts/run_cellbender_all.sh
Runs CellBender on the Cell Ranger outputs to remove ambient RNA and background.
Example usage:
# From your project directory on Minerva
./Single_Cell/scripts/run_cellbender_all.sh
This produces cleaned HDF5 matrices (one per sample), which can then be loaded into R / Seurat.
3. Seurat + Harmony Integration
Script:
Single_Cell/scripts/Seurat_integration_harmony.R
This script is a generic template for integrating multiple scRNA-seq samples after CellBender using Seurat and Harmony. It:
Loads a list of per-sample QC’d Seurat objects (one per sample)
Runs SCTransform normalization
Performs SCT-based integration across samples
Applies Harmony batch correction on the PCA embeddings
Computes neighbors, clusters, and UMAP
Finds cluster markers with FindAllMarkers
Optionally computes simple module scores for broad cell types
(e.g., neurons, astrocytes, microglia)
Usage on Minerva (example):
module load R/4.2.0

Rscript Single_Cell/scripts/Seurat_integration_harmony.R
Before running, edit the top of Seurat_integration_harmony.R and set:
qc_list_path – path to your saved Seurat list of QC’d samples
output_dir – where to save integrated objects, plots, markers, and module scores
4. Docs and Notes
Any additional notes, version info, and troubleshooting tips for this pipeline can go under:
Single_Cell/docs/
For example, you might add:
Seurat_integration_harmony_notes.md – R / Seurat / Harmony versions, common errors, and fixes
cellranger_cellbender_pipeline_notes.md – run logs and config details

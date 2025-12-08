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



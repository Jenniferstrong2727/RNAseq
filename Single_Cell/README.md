# Single-cell RNA-seq Pipeline (MiBrain / Minerva)

This folder contains scripts and docs for running the single-cell RNA-seq preprocessing pipeline on Minerva:

- **Cell Ranger** (10x Genomics) for alignment + quantification  
- **CellBender** (Broad) for ambient RNA/background removal

## Folder structure

```text
Single_Cell/
├── scripts/
│   ├── run_cellranger_all.sh     # Submit Cell Ranger jobs for all samples
│   └── run_cellbender_all.sh     # Submit CellBender jobs for all samples
└── docs/
    └── (tips / notes)



/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12
├── 00_fastq/          # new MiBrain FASTQs (RunMiBrain-1..8)
├── 01_pilot_fastq/    # symlinks to pilot FASTQs (miBrain-1, miBrain-2)
├── 02_cellranger/     # Cell Ranger outputs (one folder per sample)
├── 03_cellbender/     # CellBender outputs (one folder per sample)
├── 04_seurat_objects/ # downstream analysis in R
├── logs/              # LSF logs for Cell Ranger / CellBender
└── metadata/          # sample sheets, etc.



# Submit all Cell Ranger jobs
./Single_Cell/scripts/run_cellranger_all.sh

# After Cell Ranger finishes (or partially finishes):
./Single_Cell/scripts/run_cellbender_all.sh








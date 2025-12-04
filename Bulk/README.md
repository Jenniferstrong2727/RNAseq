# Bulk RNA-seq Pipeline Overview  
A modular pipeline for bulk RNA-seq analysis, organized into four main stages:

1. **RAPiD** — preprocessing, QC, alignment, and quantification  
2. **Collate Samples** *(coming soon)* — combine featureCounts outputs into a unified matrix  
3. **DEG Analysis** — differential expression using DESeq2  
4. **Pathway Analysis** — GSEA + ORA using FGSEA and clusterProfiler  

Detailed documentation for each stage lives in the `docs/` directory.

---

## 📁 Repository Structure

```text
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
    └── DEG_pathway_details.md       # combined DEG + Pathway documentation







1. RAPiD — Preprocessing / Alignment / QC
RAPiD runs:
Adapter trimming (Trimmomatic)
QC (FastQC, Picard)
Alignment (STAR)
Gene quantification (featureCounts, Salmon, Kallisto, RSEM)
Splicing (LeafCutter)
MultiQC (summary HTML report)
Full details: docs/RAPiD_details.md
Cluster tips: docs/RAPiD_tips_jen.md
2. Collate Samples — (Coming Soon)
This step will:
Read RAPiD’s featureCounts outputs
Combine them into a unified gene × sample matrix (.tsv)
Standardize column names
Prepare data for DESeq2
Script will live in: scripts/01_collate_samples.R
Docs will live in: docs/collate_samples_details.md
3. DEG Analysis — Differential Expression (DESeq2)
Script: scripts/BULK_rnaseq_DEG.R
Takes:
featureCounts matrix
metadata (sample, genotype, handler)
Produces:
results_raw.rds
results_shrunken_apeglm.rds
QC plots & top tables
Details: docs/DEG_pathway_details.md
4. Pathway Analysis — GSEA + ORA
Script: scripts/BULK_rnaseq_pathway_1.R
Takes the DESeq2 outputs and performs:
GSEA via fgseaMultilevel
ORA (KEGG, GO BP, Hallmark) via clusterProfiler
Dotplots + enrichment plots
Details: docs/DEG_pathway_details.md
Documentation
All documentation for this pipeline:
docs/RAPiD_details.md — full RAPiD step descriptions
docs/RAPiD_tips_jen.md — Jen’s tips for RAPiD on Minerva
docs/DEG_pathway_details.md — DESeq2 + GSEA/ORA guide
(coming soon) docs/collate_samples_details.md

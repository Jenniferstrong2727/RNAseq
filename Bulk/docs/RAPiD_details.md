# RAPiD — Bulk RNA-seq Preprocessing Pipeline

This document describes what the **RAPiD** pipeline does *after you have FASTQ files*, and where to find the outputs it creates.

High-level step order:

no_merge → trimmomatic → fastqc → star → markdup & filter_sort → QC/quant
(picard, featureCounts, rsem, kallisto, salmon, leafcutter) → multiqc



---

## 1. no_merge (internal setup)

**What it is**  
Internal RAPiD step that:

- reads the FASTQ manifest  
- organizes FASTQ pairs (R1/R2)  
- prepares jobs for downstream tools  

**Output**  
- Ready-to-process FASTQ pairs passed to trimming/alignment steps.

---

## 2. trimmomatic

**What it is**  
Adapter and quality trimming step (e.g. `--trimAdapter NexteraPE-PE`), which:

- removes adapters
- trims low-quality bases from read ends

**Output**  
- Cleaned, trimmed `.fastq.gz` files for each sample.

---

## 3. fastqc

**What it is**  
Per-sample quality assessment on the (usually trimmed) FASTQs:

- GC content  
- per-base quality  
- adapter content, etc.

**Output**  
- One `FastQC` HTML report per sample (`*.html`)  
- Summary text files per run  

These are later summarized by `MultiQC`.

---

## 4. star (alignment)

**What it is**  
Aligns reads to a reference genome (e.g. **GRCh38.Gencode.v30**):

- spliced alignment of reads to genome  
- handles introns/exons  

**Output**  
- A coordinate-sorted (or sortable) BAM file per sample (e.g. `sample.Aligned.out.bam` or similar)  

These BAMs then go into duplication marking and filtering.

---

## 5. markdup & filter_sort (Picard + Samtools)

**What it is**  

- **Picard** is used to mark PCR duplicates.  
- **Samtools** or similar is used to sort and filter the BAM.

This step produces BAM files that are ready for:

- QC metrics (Picard)  
- gene quantification (featureCounts, RSEM, etc.)

**Output**  

- Final processed BAM per sample, e.g. `sample.rmdup.sorted.bam`.

---

## 6. QC & Quantification

These steps often run in parallel. The key tools:

### Picard metrics

Calculates quality and alignment metrics, e.g.:

- insert size  
- duplication rate  
- RNA-seq specific metrics  

**Output**  

- Various `*.metrics` files under `QC/` (exact paths depend on config).

---

### featureCounts

**What it is**  

- Gene-level counting: how many reads overlap each gene.

**Output**  

- A gene count matrix (e.g. `counts.txt`) with:
  - rows = genes  
  - columns = samples  

> 💡 This is the table you feed into your **DESeq2** / `BULK_rnaseq_DEG.R` step.

---

### RSEM, kallisto, salmon

**What they do**  

- These tools perform transcript- and/or gene-level quantification using different statistical models.

**Outputs**  

- Tool-specific quantification files (TPM, FPKM, counts) in subfolders:
  - `Processed/rsem/`
  - `Processed/kallisto/`
  - `Processed/salmon/`

You can use these for alternative downstream analyses if you don’t want to rely solely on featureCounts.

---

### leafcutter (splicing)

**What it is**  

- Focused on **alternative splicing**, using intron clusters and junction reads.  
- Detects differential splicing patterns between groups.

**Output**  

- Files under `Processed/leafcutter/`  
- Junction counts and splicing cluster data for further analysis.

---

## 7. multiqc (final aggregated QC)

**What it is**  

- Scans all log files and QC outputs from:
  - FastQC
  - STAR
  - Picard
  - Trimmomatic
  - etc.

- Combines everything into **one** interactive HTML report.

**Output**  

- `QC/multiqc_report.html` — the first file you should open to assess run quality.

---

# 📁 Where to find the important outputs

Assuming you are inside your `RAPiD-nf` project directory, expect something like:

```text
RAPiD-nf/
├── QC/
│   ├── fastqc/                 # per-sample FastQC reports
│   ├── picard_metrics/         # picard outputs (names vary)
│   └── multiqc_report.html     # combined QC dashboard
└── Processed/
    ├── RAPiD/bams/             # final BAMs (rmdup, sorted)
    ├── featureCounts/          # gene-level counts (e.g. counts.txt)
    ├── rsem/                   # RSEM quant outputs
    ├── kallisto/               # Kallisto quant outputs
    ├── salmon/                 # Salmon quant outputs
    └── leafcutter/             # splicing analysis outputs


The most important files for typical downstream DE analysis are:
Processed/featureCounts/ (gene count matrix)
QC/multiqc_report.html (QC overview)
Processed/RAPiD/bams/ (aligned, cleaned BAMs if needed later)



From there, you continue with:
RAPiD → generate counts + QC
BULK_rnaseq_DEG.R → DESeq2
BULK_rnaseq_pathway_1.R → GSEA + ORA



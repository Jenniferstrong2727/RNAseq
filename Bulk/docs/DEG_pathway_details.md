# BULK_rnaseq_DEG.R — Differential Expression (DESeq2)

This document describes the inputs, outputs, model design, usage, and troubleshooting for the script **`BULK_rnaseq_DEG.R`**, which performs differential expression (DESeq2) for bulk RNA-seq based on a simple **case vs control** genotype comparison.

---

# 📥 Inputs

## 1. featureCounts table (`--counts`)

A TSV file with:

- one row per gene  
- one column per sample  
- a gene ID column named **Geneid** or **geneid**

Example header:

Geneid SampleA SampleB SampleC
ENSG0001 123 456 789


Sample column names MUST match the `sample` column in metadata.

---

## 2. metadata table (`--metadata`)

An `.xlsx` file with REQUIRED columns:

| column   | description |
|----------|-------------|
| `sample` | must match count matrix column names (case-insensitive) |
| `genotype` | used to assign case vs control groups |
| `handler` | optional batch or processing covariate |

Example:

| sample | genotype | handler |
|--------|----------|---------|
| PL-1   | G/G      | H1      |
| PL-2   | A/A      | H1      |
| PL-3   | G/G      | H2      |

---

# 🧩 Name Alignment Rules (Important!)

If you get an error like:

Cannot align metadata to counts header


Check:

- No typos  
- No `.bam`, `_Aligned.out`, or other suffixes in column names  
- Hyphens and underscores match exactly (`PL-1` ≠ `PL_1`)  
- Extra spaces are removed  

---

# ▶️ Typical Usage

Rscript scripts/BULK_rnaseq_DEG.R
--counts /path/to/featureCounts.tsv
--metadata /path/to/metadata.xlsx
--outdir /path/to/OutputRoot
--project MyProject
--case-geno "A/A"
--ctrl-geno "G/G"


### Different genotypes?
Example:


### Different genotypes?
Example:

--case-geno "T/C" --ctrl-geno "T/T"


---

# 📦 Outputs

DEG results are written to:

<outdir>/<project>/
└── DEG_handler_adjusted/
├── results_raw.rds
├── results_shrunken_apeglm.rds
├── Top100_.csv
├── QC_.png
└── ...


### What each file means

| file | meaning |
|------|---------|
| **results_raw.rds** | raw `DESeq2::results()` output (stat, pvalue, padj) |
| **results_shrunken_apeglm.rds** | lfcShrink with apeglm → stable effect sizes |
| **Top100_*.csv** | optional filtered tables |
| **QC_*.png** | quality plots such as p-value histograms |

---

# 🧠 DE Model Design

Internally the script models:

~ handler + group


Where:

- `group` = case vs control
- control group becomes reference level
- contrast usually = `group_case_vs_control`

This is printed in logs.

---

# ⚙️ Command-line Options (DEG)

| option | description |
|--------|-------------|
| `--counts` | path to featureCounts .tsv |
| `--metadata` | path to metadata .xlsx |
| `--outdir` | output root directory |
| `--project` | project name (subfolder name) |
| `--case-geno` | genotype string representing CASE samples |
| `--ctrl-geno` | genotype string representing CONTROL samples |

Outputs generated:

- DESeq2 results (raw + shrunk)
- QC plots
- cleaned metadata used for modeling

---

# ❗ Troubleshooting

## Problem: *“Cannot align metadata to counts header”*

Fix:

- ensure exact sample name match  
- remove `.bam`, `.fq`, `_Aligned.out`  
- convert spaces → underscores  
- check for Excel auto-formatting issues

---

## Problem: apeglm `coef not found`

Run inside R:


Where:

- `group` = case vs control
- control group becomes reference level
- contrast usually = `group_case_vs_control`

This is printed in logs.

---

# ⚙️ Command-line Options (DEG)

| option | description |
|--------|-------------|
| `--counts` | path to featureCounts .tsv |
| `--metadata` | path to metadata .xlsx |
| `--outdir` | output root directory |
| `--project` | project name (subfolder name) |
| `--case-geno` | genotype string representing CASE samples |
| `--ctrl-geno` | genotype string representing CONTROL samples |

Outputs generated:

- DESeq2 results (raw + shrunk)
- QC plots
- cleaned metadata used for modeling

---

# ❗ Troubleshooting

## Problem: *“Cannot align metadata to counts header”*

Fix:

- ensure exact sample name match  
- remove `.bam`, `.fq`, `_Aligned.out`  
- convert spaces → underscores  
- check for Excel auto-formatting issues

---

## Problem: apeglm `coef not found`

Run inside R:


Where:

- `group` = case vs control
- control group becomes reference level
- contrast usually = `group_case_vs_control`

This is printed in logs.

---

# ⚙️ Command-line Options (DEG)

| option | description |
|--------|-------------|
| `--counts` | path to featureCounts .tsv |
| `--metadata` | path to metadata .xlsx |
| `--outdir` | output root directory |
| `--project` | project name (subfolder name) |
| `--case-geno` | genotype string representing CASE samples |
| `--ctrl-geno` | genotype string representing CONTROL samples |

Outputs generated:

- DESeq2 results (raw + shrunk)
- QC plots
- cleaned metadata used for modeling

---

# ❗ Troubleshooting

## Problem: *“Cannot align metadata to counts header”*

Fix:

- ensure exact sample name match  
- remove `.bam`, `.fq`, `_Aligned.out`  
- convert spaces → underscores  
- check for Excel auto-formatting issues

---

## Problem: apeglm `coef not found`

Run inside R:

resultsNames(dds)


If you do NOT see `group_case_vs_control`, fix by:

- releveling metadata factor in the sheet  
- swapping `--case-geno` and `--ctrl-geno`  

---

## Problem: Very few or zero significant genes

- check DESeq2 size factors  
- inspect QC: p-value histogram should NOT be flat  
- consider loosening padj threshold  

---

# 📜 Reproducibility

Store session info:

```r
writeLines(
  capture.output(sessionInfo()),
  file.path(<outdir>, <project>, "sessionInfo_DEG.txt")
)

*.rds
*.csv
*.tsv
*.xlsx
*.png
*.pdf
Output/
*/DEG_handler_adjusted/
*/PATHWAYS/
*/QC/


If using this workflow in published work, consider citing:
DESeq2
apeglm
clusterProfiler (if performing ORA later)
fgsea (if pathway analysis later)
org.Hs.eg.db or relevant organism database

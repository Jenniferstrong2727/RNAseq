# Bulk RNA-seq (DEG + Pathways)

Two scripts:
- `BULK_rnaseq_DEG.R` – runs DESeq2, saves raw and shrunken results, plus basic QC.
- `BULK_rnaseq_pathway_1.R` – runs GSEA (fgseaMultilevel) and ORA (GO BP, KEGG, Hallmark) from those results; saves tables and plots.

---

## Requirements

R ≥ 4.2 and the following packages:

```r
install.packages(c("optparse","readr","readxl","dplyr","tidyr","stringr","janitor",
                   "ggplot2","ggrepel"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2","AnnotationDbi","org.Hs.eg.db",
                       "msigdbr","fgsea","clusterProfiler","enrichplot"))



Rscript BULK_rnaseq_DEG.R \
  --counts   /path/to/featureCounts.tsv \
  --metadata /path/to/metadata.xlsx \
  --outdir   /path/to/OutputRoot \
  --project  MyProject \
  --case-geno "A/A" \
  --ctrl-geno "G/G"


Rscript BULK_rnaseq_pathway_1.R \
  --res_raw   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_raw.rds \
  --res_shr   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_shrunken_apeglm.rds \
  --outdir    /path/to/OutputRoot \
  --project   MyProject \
  --msig_collection H \
  --padj_cutoff 0.05 \
  --lfc_cutoff  0.5

# Bulk RNA-seq (DEG + Pathways)

Two scripts:
- `BULK_rnaseq_DEG.R` – runs DESeq2, saves raw and shrunken results, plus basic QC.
- `BULK_rnaseq_pathway_1.R` – runs GSEA (fgseaMultilevel) and ORA (GO BP, KEGG, Hallmark) from those results; saves tables and plots.

---

## Requirements

R ≥ 4.2 and the following packages:

```r
install.packages(c("optparse","readr","readxl","dplyr","tidyr","stringr","janitor",
                   "ggplot2","ggrepel"))
install.packages("BiocManager")
BiocManager::install(c("DESeq2","AnnotationDbi","org.Hs.eg.db",
                       "msigdbr","fgsea","clusterProfiler","enrichplot"))



Rscript BULK_rnaseq_DEG.R \
  --counts   /path/to/featureCounts.tsv \
  --metadata /path/to/metadata.xlsx \
  --outdir   /path/to/OutputRoot \
  --project  MyProject \
  --case-geno "A/A" \
  --ctrl-geno "G/G"


Rscript BULK_rnaseq_pathway_1.R \
  --res_raw   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_raw.rds \
  --res_shr   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_shrunken_apeglm.rds \
  --outdir    /path/to/OutputRoot \
  --project   MyProject \
  --msig_collection H \
  --padj_cutoff 0.05 \
  --lfc_cutoff  0.5

##############
1) featureCounts table (--counts)
TSV with one row per gene and one column per sample.
A gene id column named Geneid (case-insensitive) or geneid.
Sample column headers should match the sample column in metadata (after simple cleaning like trimming spaces).
Example (header only):
Geneid    SampleA    SampleB    SampleC
ENSG...   123        456        789
2) metadata (--metadata)
Excel (.xlsx) with at least these columns (case-insensitive; the script cleans names):
sample – must match the counts column headers (see “Name alignment” below)
genotype – used to map case vs control (e.g., A/A, G/G, or T/C, T/T)
handler – optional batch/handler covariate (string or factor)
Minimal example:
| sample | genotype | handler |
|-------|----------|---------|
| PL-1 | G/G | H1 |
| PL-2 | A/A | H1 |
| PL-3 | G/G | H2 |
Name alignment tips
The DEG script does light sanitization (trims spaces and normalizes obvious suffixes); however, keep names consistent where possible.
If you see an error like “Cannot align metadata to counts header”, double-check:
typos in sample
extra file suffixes in counts headers (e.g., remove .bam, _Aligned.out etc.)
consistent use of hyphens/underscores
Typical run commands
DEG (DESeq2)
Rscript BULK_rnaseq_DEG.R \
  --counts   /path/to/featureCounts.tsv \
  --metadata /path/to/metadata.xlsx \
  --outdir   /path/to/OutputRoot \
  --project  MyProject \
  --case-geno "A/A" \
  --ctrl-geno "G/G"
Pathways (GSEA + ORA)
Rscript BULK_rnaseq_pathway_1.R \
  --res_raw   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_raw.rds \
  --res_shr   /path/to/OutputRoot/MyProject/DEG_handler_adjusted/results_shrunken_apeglm.rds \
  --outdir    /path/to/OutputRoot \
  --project   MyProject \
  --msig_collection H \
  --padj_cutoff 0.05 \
  --lfc_cutoff  0.5
🧪 Different genotypes? Just change --case-geno and --ctrl-geno, e.g. for LRRK2:
--case-geno "T/C" --ctrl-geno "T/T"
Output structure
<outdir>/<project>/
├─ DEG_handler_adjusted/
│  ├─ results_raw.rds                      # DESeq2::results()
│  ├─ results_shrunken_apeglm.rds          # lfcShrink() with apeglm
│  ├─ Top100 CSVs, QC histograms, etc.
│  └─ ...
├─ PATHWAYS/
│  ├─ <project>_fgsea_hallmark_rankByWaldStat.(csv|rds)
│  ├─ <project>_fgsea_hallmark_rankByShrunkenLFC.(csv|rds)
│  ├─ <project>_ORA_GO_BP.(csv|rds)
│  ├─ <project>_ORA_KEGG.(csv|rds)
│  ├─ <project>_ORA_Hallmark.(csv|rds)     # Hallmark via msigdbr ENTREZ
│  └─ plots/
│     ├─ GSEA_dotplot_*.png                # top 20 FDR
│     ├─ GSEA_ES_*.png                     # enrichment curves
│     ├─ ORA_GO_BP_dotplot_top20.png
│     ├─ ORA_KEGG_dotplot_top20.png
│     └─ ORA_Hallmark_dotplot_top20.png
What each file is for (quick glossary)
results_raw.rds — raw Wald stats, p-values, padj; best for ranking by stat (GSEA) and for provenance.
results_shrunken_apeglm.rds — effect sizes shrunk by apeglm; best for ranking by LFC (GSEA) and filtering significant genes for ORA (more stable LFCs).
GSEA outputs — pathway-level enrichment using a continuous ranked list (two flavors: by raw stat and by shrunken LFC).
ORA outputs — pathway over-representation using significant gene set (padj + |LFC| thresholds).
plots/ — communication-ready figures: dotplots (overview) and ES curves (examples).
Options reference
BULK_rnaseq_DEG.R
--counts (required) – featureCounts .tsv
--metadata (required) – metadata .xlsx
--outdir (required) – output root
--project (required) – project name (folder tag)
--case-geno / --ctrl-geno (required) – map metadata genotype → case/control
Design: ~ handler + group, with control as the reference level
Exports: DESeq2 objects (RDS), top tables (CSV), QC p-value histograms (PNG)
BULK_rnaseq_pathway_1.R
--res_raw (required) – results_raw.rds from DEG
--res_shr (required) – results_shrunken_apeglm.rds
--outdir, --project (required)
--msig_collection – MSigDB collection (e.g., H, C2, C5; default H)
--padj_cutoff – FDR for ORA sig set (default 0.05)
--lfc_cutoff – |LFC| threshold for ORA sig set (default 0.5)
--make_cnet – optional cnetplot (requires enrichplot, ggraph, ggtangle)
Notes & recommendations
GSEA method: script uses fgseaMultilevel automatically (no nperm argument).
Hallmark ORA: uses msigdbr ENTREZ IDs directly (avoids SYMBOL↔ENTREZ ambiguity).
Gene ID mapping rates: it’s normal to see 75–85% SYMBOL→ENTREZ mapping depending on gene set and annotation version.
KEGG: clusterProfiler::enrichKEGG() expects ENTREZ IDs and organism "hsa".
Too many plots? You can delete subsets you don’t need from <project>/PATHWAYS/plots/—the tables (CSV/RDS) are the authoritative outputs.
Troubleshooting
“Cannot align metadata to counts header.”
Ensure metadata$sample exactly matches counts column headers (case/spacing/hyphenation). Remove file extensions like .bam from counts headers if present.
apeglm/coef error (e.g., coef not found).
Check resultsNames(dds) and verify the script detected group_case_vs_control. If levels are inverted, re-level group or flip --case-geno/--ctrl-geno.
msigdbr warnings about category deprecated.
Use collection = "H" (the scripts already do this).
cnetplot errors about ggtangle.rdb corruption.
Reinstall:
install.packages("BiocManager")
BiocManager::install(c("enrichplot","ggraph","ggtangle"), ask = FALSE)
Or skip --make_cnet.
Low or zero ORA genes.
Loosen thresholds (e.g., --padj_cutoff 0.1 and/or --lfc_cutoff 0.2) and re-run.
Reproducibility
Save your session info alongside outputs:
writeLines(capture.output(sessionInfo()),
           file.path(<outdir>, <project>, "sessionInfo.txt"))
Consider using a lockfile/renv for strict package versions.
.gitignore (suggested)
Create a .gitignore in this folder:
# never commit raw data or large outputs
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
License & citation (optional)
Add a LICENSE 
If you publish, cite: DESeq2, apeglm, fgsea, msigdbr, clusterProfiler, enrichplot, org.Hs.eg.db.

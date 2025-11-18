# Bulk RNA-seq DEG

Generic DESeq2 pipeline (single-project per run).

**Run example:**
./bulk_deg.R \
  --counts /path/to/featureCounts.tsv \
  --metadata /path/to/metadata.xlsx \
  --outdir /path/to/Output \
  --project MyProject \
  --case-geno "A/A" --ctrl-geno "G/G"

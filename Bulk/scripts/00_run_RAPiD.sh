#!/usr/bin/env bash
set -euo pipefail

# ---- USER CONFIG ----
PROJECT_DIR=/sc/arion/projects/ad-omics/Jennifer/RAPiD_P2RY12_bulk
FASTQ_DIR=/sc/arion/projects/ad-omics/Jennifer/bulk_rna_seq/00_fastq
RAPID_SOURCE=/sc/arion/projects/ad-omics/data/software/RAPiD-nf
ACCOUNT=acc_ad-omics
GENOME=GRCh38.Gencode.v30
STRANDED=none
# ---------------------

# 1. Create and enter project directory
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR"

# 2. Copy lab's RAPiD pipeline
cp -R "$RAPID_SOURCE" ./RAPiD-nf
cd RAPiD-nf

# 3. Create sample_fastq_table.tsv from FASTQ filenames
python3 - <<PY
import os, re, sys
base = "$FASTQ_DIR"
files = [f for f in os.listdir(base) if f.endswith(".fastq.gz")]
pairs = {}
pat = re.compile(r"(.+)_L\\d{3}_R([12])_001\\.fastq\\.gz$")
for f in files:
    m = pat.match(f)
    if not m:
        print(f"# WARNING: unexpected filename: {f}", file=sys.stderr)
        continue
    sample, read = m.group(1), m.group(2)
    pairs.setdefault(sample, {})[read] = os.path.join(base, f)
samples = sorted(pairs.keys())
with open("sample_fastq_table.tsv", "w") as w:
    w.write("sample\tf1\tf2\n")
    for s in samples:
        w.write(f"{s}\t{pairs[s]['1']}\t{pairs[s]['2']}\n")
print("Wrote sample_fastq_table.tsv with", len(samples), "samples")
PY

# 4. Set up Nextflow cache and framework
mkdir -p "$PROJECT_DIR/.nextflow"
export NXF_HOME="$PROJECT_DIR/.nextflow"
export NXF_WORK=$(pwd)/work

mkdir -p ../.nextflow/framework/19.10.0
cd ../.nextflow/framework/19.10.0
curl -L https://www.nextflow.io/releases/v19.10.0/nextflow-19.10.0-one.jar -o nextflow-19.10.0-one.jar

# return to RAPiD directory
cd "$PROJECT_DIR/RAPiD-nf"

# 5. Create per-sample directory structure and symlink FASTQs
tail -n +2 sample_fastq_table.tsv | while IFS=$'\t' read -r sample f1 f2; do
  mkdir -p "$sample/Raw/Illumina"
  ln -sf "$f1" "$sample/Raw/Illumina/$(basename "$f1")"
  ln -sf "$f2" "$sample/Raw/Illumina/$(basename "$f2")"
done

# 6. Meta + logs
mkdir -p Meta Logs
cp sample_fastq_table.tsv Meta/
sed -i 's/acc_apollo/'"$ACCOUNT"'/g' config/lsf.config

# 7. Submit RAPiD job
bsub -J rapid_p2ry12 -n 1 -W 48:00 -q premium -P "$ACCOUNT" \
  -o "$PROJECT_DIR/RAPiD-nf/Logs/rapid_%J.out" \
  -e "$PROJECT_DIR/RAPiD-nf/Logs/rapid_%J.err" \
  -R "rusage[mem=3750] span[hosts=1]" \
  "/bin/bash -lc '
     set -euo pipefail
     export NXF_OFFLINE=true
     export NXF_HOME=$PROJECT_DIR/.nextflow
     export NXF_WORK=\$(pwd)/work

     source /etc/profile.d/modules.sh || true
     module load java || true

     /hpc/packages/minerva-common/nextflow/19.10.0/bin/nextflow run \$(pwd)/RAPiD.nf \
       --run \$(pwd) \
       --genome $GENOME \
       --stranded $STRANDED \
       -profile chimera \
       --qc --fastqc --leafcutter --featureCounts --rsem --kallisto --salmon \
       --trimAdapter NexteraPE-PE \
       -resume
  '"


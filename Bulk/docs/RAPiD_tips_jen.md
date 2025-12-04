# RAPiD Tips
These are practical notes from running **RAPiD** on the cluster (Minerva-style environment), including manifest preparation, `screen`, and `bsub` usage.

---

## 1. Getting a working copy of RAPiD

Example:

```bash
# Copy lab version of RAPiD into your project directory
cp -R /sc/arion/projects/ad-omics/data/software/RAPiD-nf $PWD

# Optionally clone additional helper scripts (Raj lab RNA-pipelines)
git clone git@github.com:RajLabMSSM/RNA-pipelines.git


Sample manifest: full paths for FASTQ files
RAPiD sometimes expects full paths in the f1/f2 columns of your sample manifest.
Example R code used to prepend full paths:

d <- read.table("MyND_samples.tsv", header = TRUE)  # original manifest
head(d)  # check

d$f1 <- paste0("/sc/arion/projects/pd-omics/30-1165110063/fastq/", d$f1)
d$f2 <- paste0("/sc/arion/projects/pd-omics/30-1165110063/fastq/", d$f2)

file.exists(d$f1)  # sanity check paths
file.exists(d$f2)

library(readr)
write_tsv(d, "MyND4_samples.tsv")  # save updated manifest

d <- read.table("MyND_samples.tsv", header = TRUE)  # original manifest
head(d)  # check

d$f1 <- paste0("/sc/arion/projects/pd-omics/30-1165110063/fastq/", d$f1)
d$f2 <- paste0("/sc/arion/projects/pd-omics/30-1165110063/fastq/", d$f2)

file.exists(d$f1)  # sanity check paths
file.exists(d$f2)

library(readr)
write_tsv(d, "MyND4_samples.tsv")  # save updated manifest


Use this pattern with your own base FASTQ directory.

cd RAPiD-nf/
python  # ensure Python is available (or use python3)

# Example using helper from RajLab RNA-pipelines repo
python ../RNA-pipelines/RAPiD-tips/create_rapid_structure.py ../MyND4_samples.tsv

cd RAPiD-nf/
python  # ensure Python is available (or use python3)

# Example using helper from RajLab RNA-pipelines repo
python ../RNA-pipelines/RAPiD-tips/create_rapid_structure.py ../MyND4_samples.tsv

This script sets up metadata + directory structure for RAPiD (e.g., Meta/, Raw/Illumina/, etc.).

# Start a screen named RAPiD
screen -S RAPiD

# Start a screen named RAPiD
screen -S RAPiD

Inside screen, you can run bsub to submit jobs. Useful keybinds:
Detach: Ctrl + a, then d
Reattach: screen -r RAPiD
Terminate: Ctrl + d from within that screen
You can have multiple screens for different runs.

bsub -I -n 2 -W 144:00 -q long -P acc_pd-omics \
  -R "rusage[mem=3750]" -R "span[hosts=1]" \
  "/sc/arion/projects/H_PBG/nextflow/bin/nextflow run /path/to/RAPiD-nf/RAPiD.nf \
     --run `pwd` \
     --genome GRCh38.Gencode.v30 \
     --stranded none \
     -profile chimera \
     --qc --fastqc --leafcutter --featureCounts --rsem --kallisto --salmon \
     --trimAdapter NexteraPE-PE \
     --rawPath Raw/Illumina \
     -resume"

bsub -I -n 2 -W 144:00 -q long -P acc_pd-omics \
  -R "rusage[mem=3750]" -R "span[hosts=1]" \
  "/sc/arion/projects/H_PBG/nextflow/bin/nextflow run /path/to/RAPiD-nf/RAPiD.nf \
     --run `pwd` \
     --genome GRCh38.Gencode.v30 \
     --stranded none \
     -profile chimera \
     --qc --fastqc --leafcutter --featureCounts --rsem --kallisto --salmon \
     --trimAdapter NexteraPE-PE \
     --rawPath Raw/Illumina \
     -resume"

Notes:
--rawPath Raw/Illumina should match where your FASTQs are linked.
--genome and --stranded must match your experiment and RAPiD config.
-resume lets Nextflow continue where it left off if partially completed.

bjobs -w         # list jobs
bpeek <JOBID>    # peek at job output
bpeek -f <JOBID> # follow job output
ls -l Logs/      # check log files if you've defined a Logs/ directory


bjobs -w         # list jobs
bpeek <JOBID>    # peek at job output
bpeek -f <JOBID> # follow job output
ls -l Logs/      # check log files if you've defined a Logs/ directory


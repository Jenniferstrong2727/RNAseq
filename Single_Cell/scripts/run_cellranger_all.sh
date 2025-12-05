#!/bin/bash
set -euo pipefail

#############################################
# Cell Ranger multi-sample submission script
# - Designed for LSF (bsub) on Minerva
# - One job per sample
#############################################

########## USER CONFIG (EDIT THESE) ##########

# Base project directory
PROJECT_DIR="/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12"

# Where FASTQs live
FASTQ_NEW="${PROJECT_DIR}/00_fastq"              # new 10x runs
FASTQ_PILOT_BASE="/sc/arion/projects/ad-omics/miBrain/data"  # pilot symlinks

# Cell Ranger reference
TRANSCRIPTOME="/sc/arion/projects/ad-omics/sc_mic_rawdata/cellranger_files/refdata-cellranger-GRCh38-3.0.0"

# Output + logs
OUT_BASE="${PROJECT_DIR}/02_cellranger"
LOG_DIR="${PROJECT_DIR}/logs"

# LSF settings
QUEUE="premium"
PROJECT_CODE="acc_ad-omics"
N_CORES=8
MEM_MB=64000
WALL_HOURS=72

#############################################

mkdir -p "${OUT_BASE}" "${LOG_DIR}"

# Helper to submit one Cell Ranger job
submit_cr_job () {
  local SAMPLE_ID="$1"      # e.g. PL1_new
  local FASTQ_DIR="$2"      # e.g. /.../00_fastq or pilot dir
  local SAMPLE_PREFIX="$3"  # e.g. RunMiBrain-1 or miBrain-1

  local OUTDIR="${OUT_BASE}/${SAMPLE_ID}"

  echo "Submitting job for ${SAMPLE_ID} with prefix ${SAMPLE_PREFIX}"

  bsub <<EOF
#BSUB -J CR_${SAMPLE_ID}
#BSUB -n ${N_CORES}
#BSUB -M ${MEM_MB}
#BSUB -W ${WALL_HOURS}:00
#BSUB -P ${PROJECT_CODE}
#BSUB -q ${QUEUE}
#BSUB -o ${LOG_DIR}/CR_${SAMPLE_ID}_%J.out
#BSUB -e ${LOG_DIR}/CR_${SAMPLE_ID}_%J.err

module load cellranger/9.0.1
export MRO_DISK_SPACE_CHECK=disable   # disable disk-space preflight if needed

cellranger count \\
  --id=${SAMPLE_ID} \\
  --transcriptome=${TRANSCRIPTOME} \\
  --fastqs=${FASTQ_DIR} \\
  --sample=${SAMPLE_PREFIX} \\
  --include-introns=true \\
  --create-bam=false \\
  --output-dir=${OUTDIR}
EOF
}

#############################################
# PROJECT-SPECIFIC EXAMPLE (MiBrain P2RY12)
#############################################

echo "Submitting Cell Ranger jobs for PILOT samples..."
# Pilot samples (FASTQs already demultiplexed with sample IDs in read names)
submit_cr_job "PL1_pilot" "${FASTQ_PILOT_BASE}/miBrain-1" "miBrain-1"
submit_cr_job "B3_1_pilot" "${FASTQ_PILOT_BASE}/miBrain-2" "miBrain-2"

echo "Submitting Cell Ranger jobs for NEW samples..."
# New samples (RunMiBrain-X naming scheme)
submit_cr_job "PL1_new" "${FASTQ_NEW}" "RunMiBrain-1"
submit_cr_job "PL2_new" "${FASTQ_NEW}" "RunMiBrain-2"
submit_cr_job "PL3_new" "${FASTQ_NEW}" "RunMiBrain-3"
submit_cr_job "B3_1_new" "${FASTQ_NEW}" "RunMiBrain-4"
submit_cr_job "B3_2_new" "${FASTQ_NEW}" "RunMiBrain-5"
submit_cr_job "G1_1_new" "${FASTQ_NEW}" "RunMiBrain-7"
submit_cr_job "G1_2_new" "${FASTQ_NEW}" "RunMiBrain-8"

echo "All Cell Ranger jobs submitted."

#!/bin/bash
set -euo pipefail

######################################################
# CellBender remove-background multi-sample launcher
# - GPU version (cellbender/0.2.1-gpu on Minerva)
# - Submits one job per sample via LSF (bsub)
######################################################

########## USER CONFIG (EDIT THESE) ##########

PROJECT_DIR="/sc/arion/projects/ad-omics/Jennifer/scRNA_mibrain_P2RY12"

CR_OUT_BASE="${PROJECT_DIR}/02_cellranger"
CB_OUT_BASE="${PROJECT_DIR}/03_cellbender"
LOG_DIR="${PROJECT_DIR}/logs"

# Samples expected to have Cell Ranger outputs in:
#   ${CR_OUT_BASE}/${SAMPLE_ID}/outs/raw_feature_bc_matrix.h5
SAMPLES=(
  PL1_new
  PL2_new
  PL3_new
  B3_1_new
  B3_2_new
  G1_1_new
  G1_2_new
)

# Global hyperparameters (tune per dataset/project)
EXPECTED_CELLS=205
TOTAL_DROPLETS_INCLUDED=5000

# LSF settings
CB_QUEUE="gpu"
CB_PROJECT="acc_ad-omics"
CB_NCORES=4
CB_MEM_MB=64000
CB_WALL_HOURS=24

#############################################

mkdir -p "${CB_OUT_BASE}" "${LOG_DIR}"

submit_cb_job () {
  local SAMPLE_ID="$1"

  local CR_OUT="${CR_OUT_BASE}/${SAMPLE_ID}/outs"
  local INPUT_H5="${CR_OUT}/raw_feature_bc_matrix.h5"
  local CB_OUTDIR="${CB_OUT_BASE}/${SAMPLE_ID}"
  local OUTPUT_PREFIX="${CB_OUTDIR}/${SAMPLE_ID}_cellbender"

  # Skip if Cell Ranger not done yet
  if [[ ! -f "${INPUT_H5}" ]]; then
    echo "Skipping ${SAMPLE_ID}: raw_feature_bc_matrix.h5 not found."
    return
  fi

  # Skip if CellBender already completed
  if [[ -f "${OUTPUT_PREFIX}.h5" ]]; then
    echo "Skipping ${SAMPLE_ID}: CellBender output already exists."
    return
  fi

  mkdir -p "${CB_OUTDIR}"

  echo "Submitting CellBender job for ${SAMPLE_ID}..."

  bsub <<EOF
#BSUB -J CB_${SAMPLE_ID}
#BSUB -q ${CB_QUEUE}
#BSUB -gpu "num=1"
#BSUB -n ${CB_NCORES}
#BSUB -M ${CB_MEM_MB}
#BSUB -W ${CB_WALL_HOURS}:00
#BSUB -P ${CB_PROJECT}
#BSUB -o ${LOG_DIR}/CB_${SAMPLE_ID}_%J.out
#BSUB -e ${LOG_DIR}/CB_${SAMPLE_ID}_%J.err

module load cellbender/0.2.1-gpu

cellbender remove-background \\
  --input "${INPUT_H5}" \\
  --output "${OUTPUT_PREFIX}.h5" \\
  --expected-cells ${EXPECTED_CELLS} \\
  --total-droplets-included ${TOTAL_DROPLETS_INCLUDED} \\
  --cuda
EOF
}

echo "Submitting CellBender jobs..."
for SAMPLE_ID in "${SAMPLES[@]}"; do
  submit_cb_job "${SAMPLE_ID}"
done
echo "All CellBender jobs submitted."

#!/usr/bin/bash -l

module load singularity


# MEGAHIT metagenomic assembly script. So arguably you could
# get away with just running exec_assembly.sh, but if you
# want to be thorough and guage the quality of resulting assemblies,
# its best to run exec_megahit.sh and exec_

# Usage: ./exec_megahit.sh <input_dir> <output_dir>
# Input: interleaved filtered fastq files
#        (*_interleaved_filtered.fastq.gz)


set -Eeuo pipefail
shopt -s nullglob

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    echo "[INFO] SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

log() {
    echo "[$(date +'%F %T')] $*"
}

trap 'rc=$?; log "[ERROR] line ${LINENO}: command failed: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR

if [[ $# -ne 2 ]]; then
  log "[ERROR] Usage: $0 <input_dir> <output_dir>"
  exit 2
fi

INPUT_DIR=$1
OUT_DIR=$2

MEGAHIT_THREADS=${SLURM_CPUS_PER_TASK:-48}
MEGAHIT_MEM=${MEGAHIT_MEM:-0.9}
MEGAHIT_IMAGE="megahit.sif"

if [[ ! -f "${MEGAHIT_IMAGE}" ]]; then
    log "[ERROR] MEGAHIT singularity image not found at ${MEGAHIT_IMAGE}."
    exit 1
fi

FILTERED_FASTQS=( ${INPUT_DIR}/*_interleaved_filtered.fastq.gz )
mkdir -p "${OUT_DIR}"

if [[ ${#FILTERED_FASTQS[@]} -eq 0 ]]; then
  log "[ERROR] No interleaved filtered files found in ${INPUT_DIR}/"
  exit 1
fi

for i in "${FILTERED_FASTQS[@]}"; do
    BASE="$(basename "${i}" _interleaved_filtered.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}/megahit_output"
    SAMPLE_OUT_SENS="${SAMPLE_OUT_DIR}/sensitive"
    SAMPLE_OUT_LARGE="${SAMPLE_OUT_DIR}/large"

    if [[ -f "${SAMPLE_OUT_DIR}/final.contigs.fa" ]]; then
        log "[SKIP] MEGAHIT output already exists for ${BASE}, skipping."
        continue
    fi

    rm -rf "${SAMPLE_OUT_DIR}"

    log "[RUN] megahit: ${BASE}"

    # meta-sensitive
    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${MEGAHIT_IMAGE}" \
      megahit \
      --12 "${i}" \
      -t "${MEGAHIT_THREADS}" \
      -m "${MEGAHIT_MEM}" \
      --min-count 2 \
      --k-list 1,29,39,49,59,69,79,89,99,129,141 \
      -o "${SAMPLE_OUT_SENS}"

      # meta-large
    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${MEGAHIT_IMAGE}" \
      megahit \
      --12 "${i}" \
      -t "${MEGAHIT_THREADS}" \
      -m "${MEGAHIT_MEM}" \
      --k-min 27 \
      --k-max 127 \
      --k-step 10 \
      -o "${SAMPLE_OUT_LARGE}"
done

log "MEGAHIT assembly completed."

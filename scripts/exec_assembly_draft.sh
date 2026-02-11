#!/usr/bin/bash -l

module load singularity
module load openjdk/17.0.5_8


######################################################
# three phases of execution:
# 1. Bloom Filtering with bbtools
# 2. SPADES
# 3. fungalrelease.sh
######################################################
# Usage: ./exec_assembly.sh <input_dir> <output_dir>
######################################################


# Phase 0: Set Up
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

# Configuration, be sure you know your paths!:
INPUT_DIR=$1 # expects interleaved filter files
OUT_DIR=$2

# java virt max should be set to 80% of available memory
MEM_FOR_JAVA_TOOLS=-Xmx160g
BBTOOLS_IMAGE="bbtools_38.86.sif"
SPADES_IMAGE="spades_3.15.2.sif"



FILTERED_FASTQS=( ${INPUT_DIR}/*_interleaved_filtered.fastq.gz )
mkdir -p ${OUT_DIR}


if [[ ${#FILTERED_FASTQS[@]} -eq 0 ]]; then
  log "[ERROR] No interleaved files were produced in ${INPUT_DIR}/"
  exit 1
fi

touch readlen.txt

# Phase 1: Bloom filtering with BBTools

BLOOM_FILTERED=()

# Dependencies: BBTools, Java 17

for i in "${FILTERED_FASTQS[@]}"; do
    log "Running BBCMS & Assembly on ${i} ..."

    BASE="$(basename "${i}" _interleaved_filtered.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}"
    SAMPLE_COUNTS="${SAMPLE_OUT_DIR}/counts.metadata.json"

    SPADES_IN1="${SAMPLE_OUT_DIR}/bbcms_output.FWD.fastq.gz"
    SPADES_IN2="${SAMPLE_OUT_DIR}/bbcms_output.REV.fastq.gz"

    mkdir -p "${SAMPLE_OUT_DIR}"
    
    log "[RUN] bbcms.sh: $BASE"

    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      readlength.sh -Xmx3g \
      in="${i}" \
      out="${SAMPLE_OUT_DIR}/${BASE}_readlength.txt" \
      overwrite

    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      bbcms.sh ${MEM_FOR_JAVA_TOOLS} \
        mincount=2 \
        highcountfraction=0.6 \
        in="${i}" \
        out1="${SPADES_IN1}" \
        out2="${SPADES_IN2}" \
        1> "${SAMPLE_OUT_DIR}/bbcms.${BASE}.log" \
        2> "${SAMPLE_OUT_DIR}/bbcms.${BASE}.err"
    
    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      reformat.sh -Xmx5g \
        in1="${SPADES_IN1}" \
        in2="${SPADES_IN2}" \
        out="${SAMPLE_OUT_DIR}/${BASE}_bbcms_interleaved.fastq.gz" \
        overwrite=t
    BLOOM_FILTERED_INT+=( "${SAMPLE_OUT_DIR}/${BASE}_bbcms_interleaved.fastq.gz" )
    grep Uniq ${SAMPLE_OUT_DIR}/bbcms.${BASE}.err | awk '{print $NF}' > "${SAMPLE_OUT_DIR}/${BASE}_unique31mer.txt"
    
done

log "Bloom filtering completed. Processed files: ${#BLOOM_FILTERED[@]}"

# conditional if spades flag in position 3 is "spades"
if [[ "$3" == "spades" || $# -lt 3 ]]; then
    log "Starting SPADES Assembly Phase ..."
else
    log "Skipping SPADES Assembly Phase as per user request ..."
    exit 0
fi

for i in "${BLOOM_FILTERED_INT[@]}"; do
    log "Running SPADES Assembly on ${i} ..."

    BASE="$(basename "${i}" _bbcms_interleaved.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}" # should be same as bloom sample dir
    SAMPLE_TMP_DIR="${OUT_DIR}/${BASE}/tmp"

    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${SPADES_IMAGE}" \
      spades.py \
      -m ${SPADES_MEM} \
      -t ${SLURM_NTASKS_PER_NODE:-8} \
      --tmp-dir ${SAMPLE_TMP_DIR} \
      -o "${SAMPLE_OUT_DIR}/spades_output" \
      --only-assembler \
      -k 33,55,77,99,127 \
      --meta \
      -1 "${SAMPLE_OUT_DIR}/bbcms_output.FWD.fastq.gz" \
      -2 "${SAMPLE_OUT_DIR}/bbcms_output.REV.fastq.gz"
done

for i in "${BLOOM_FILTERED_INT[@]}"; do
    BASE="$(basename "${i}" _bbcms_interleaved.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}"
    SPADES_OUT="${SAMPLE_OUT_DIR}/spades_output/scaffolds.fasta"
    
    FUNGAL_RELEASE_SCAFFOLDS="${SAMPLE_OUT_DIR}/${BASE}_scaffolds.fasta"
    FUNGAL_RELEASE_CONTIGS="${SAMPLE_OUT_DIR}/${BASE}_contigs.fasta"
    AGP_OUT="${SAMPLE_OUT_DIR}/${BASE}.agp"
    LEGEND_OUT="${SAMPLE_OUT_DIR}/${BASE}scaffold.legend"


    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      fungalrelease.sh -Xmx40g \
        in="${SPADES_OUT}" \
        out="${FUNGAL_RELEASE_SCAFFOLDS}" \
        outc="${FUNGAL_RELEASE_CONTIGS}" \
        agp="${AGP_OUT}" \
        legend="${LEGEND_OUT}" \
        mincontig=200 \
        minscaffold=200 \
        sortscaffolds=t \
        sortcontigs=t \
        overwrite=t
done
#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mem=500G   # maximum memory per node
#SBATCH --mail-user=  # email address
#SBATCH --mail-type=END
#SBATCH --output="s-%j-count.out" # job standard output file (%j replaced by job id)

# Base directory where samples are located
BASE_DIR="./PRRSV"

# Path to the cellranger executable
CELLRANGER="./cellranger/apps/cellranger-7.2.0/cellranger"

# Transcriptome reference
TRANSCRIPTOME="./reference_files/PRRSV"

# Iterate over each sample directory
for SAMPLE_DIR in ${BASE_DIR}/Sample_*; do
  # Extract the sample name from the directory path
  SAMPLE_NAME=$(basename ${SAMPLE_DIR})

  # Define the FASTQ directory for this sample
  FASTQ_DIR="${BASE_DIR}/${SAMPLE_NAME}"

  # Extract the sample identifier from the first FASTQ file name
  SAMPLE_ID=$(ls ${FASTQ_DIR} | grep -oP '^\d+' | head -1)

  # output directory
  OUTPUT_DIR="${SAMPLE_DIR}/${SAMPLE_NAME}_output"
  if [ -d "${OUTPUT_DIR}" ]; then
      echo "Removing existing output directory ${OUTPUT_DIR} to ensure a fresh start."
      rm -rf "${OUTPUT_DIR}"
  fi

  # cellranger count
  echo "Processing ${SAMPLE_NAME} with sample ID ${SAMPLE_ID}..."
  ${CELLRANGER} count --id=${SAMPLE_NAME} --transcriptome=${TRANSCRIPTOME} --fastqs=${FASTQ_DIR} --sample=${SAMPLE_ID} --output-dir=${OUTPUT_DIR}
  echo "Finished processing ${SAMPLE_NAME}"
done


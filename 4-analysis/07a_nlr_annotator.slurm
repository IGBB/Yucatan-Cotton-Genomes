#!/bin/bash

#SBATCH --job-name=nlr_annotator
#SBATCH --output=joblog/%x_%A_%a.log
#SBATCH --error=joblog/%x_%A_%a.log
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH -A jfw-lab

# ══════════════════════════════════════════════════════════════════════════════
# 16_nlr_annotator.slurm
# Runs NLR-Annotator on each of the five core genomes to identify NLR loci
# directly from the genome sequence (no gene annotation required).
#
# PREREQUISITES:
#   git clone https://github.com/steuernb/NLR-Annotator  (already done)
#     NLR-Annotator/NLR-Annotator-v2.1b.jar
#     NLR-Annotator/src/mot.txt
#     NLR-Annotator/src/store.txt
#
# INPUTS:
#   genomes/{genome}.fa
#
# OUTPUTS:
#   nlr_annotator/{genome}/{genome}.nlr.txt   — tabular NLR loci
#   nlr_annotator/{genome}/{genome}.nlr.gff3  — GFF3 NLR loci
#   nlr_annotator/{genome}/{genome}.nlr.bed   — BED NLR loci
# ══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

GENOMES=(TX2094 Chorlitos Consorcio Preseco Pretropical)
GENOME="${GENOMES[$SLURM_ARRAY_TASK_ID]}"

JAR="NLR-Annotator/NLR-Annotator-v2.1b.jar"
MOT="NLR-Annotator/src/mot.txt"
STORE="NLR-Annotator/src/store.txt"

OUTDIR="nlr_annotator/${GENOME}"
mkdir -p "${OUTDIR}" joblog

echo "Date: $(date --iso-8601=seconds)"
echo "Job: $SLURM_JOB_NAME  $SLURM_JOB_ID  array task $SLURM_ARRAY_TASK_ID"
echo "Host: $HOSTNAME"
echo "Genome: ${GENOME}"
echo

ml openjdk/21.0.3_9-vngib7s

java -Xmx245G \
    -jar "${JAR}" \
    -i  "genomes/${GENOME}.fa" \
    -x  "${MOT}" \
    -y  "${STORE}" \
    -t  "${SLURM_CPUS_PER_TASK}" \
    -o  "${OUTDIR}/${GENOME}.nlr.txt" \
    -g  "${OUTDIR}/${GENOME}.nlr.gff3" \
    -b  "${OUTDIR}/${GENOME}.nlr.bed"

echo
echo "$(date --iso-8601=seconds) NLR-Annotator complete for ${GENOME}"
echo "Output:"
ls -lh "${OUTDIR}/"

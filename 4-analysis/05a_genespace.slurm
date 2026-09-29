#!/bin/bash

#SBATCH --job-name=genespace
#SBATCH --output=joblog/%x_%j.log
#SBATCH --error=joblog/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH -A jfw-lab

set -euo pipefail

echo "Date: $(date --iso-8601=seconds)"
echo "Job: $SLURM_JOB_NAME  ID: $SLURM_JOB_ID"
echo "Host: $HOSTNAME"
echo

ml apptainer 
source /lustre/hdd/LAS/mhufford-lab/corrinne/app/env/genome_assembly

mkdir -p joblog genespace

Rscript 05b_genespace.R

echo "$(date --iso-8601=seconds) genespace complete"

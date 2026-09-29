#!/bin/bash

#SBATCH --job-name=resistify
#SBATCH --output=joblog/%x_%j.log
#SBATCH --error=joblog/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH --array=0-4
#SBATCH -A jfw-lab


thr=$SLURM_CPUS_PER_TASK

ml apptainer 
#apptainer exec --env HF_HUB_DISABLE_XET=1  --bind /lustre/hdd/LAS/jfw-lab/corrinne/FourGenomes/fromTony:/root/.cache/huggingface resistify.sif resistify download

samples=(Chorlitos Consorcio Preseco Pretropical TX2094)
name=${samples[$SLURM_ARRAY_TASK_ID]}

apptainer exec \
  --env HF_HUB_OFFLINE=1 \
  --env HF_HUB_DISABLE_XET=1 \
  --env HF_HUB_CACHE=/lustre/hdd/LAS/jfw-lab/corrinne/FourGenomes/fromTony \
  --bind /lustre/hdd/LAS/jfw-lab/corrinne/FourGenomes/fromTony:/root/.cache/huggingface \
  resistify.sif resistify nlr proteins/${name}.proteins.fa -o resistify/${name}_resistify
  

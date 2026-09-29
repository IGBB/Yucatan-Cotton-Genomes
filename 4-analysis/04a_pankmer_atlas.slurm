#!/bin/bash

#SBATCH --job-name=pankmer
#SBATCH --output=joblog/%x_%j.log
#SBATCH --error=joblog/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH -A jfw-lab

# Pankmer builds a k-mer index then produces a pairwise similarity matrix
# and a clustered heatmap SVG.

set -euo pipefail

INDEX=manyGenomes/fromAtlas/genomes_index.tar
MATRIX=manyGenomes/fromAtlas/pankmer_adj_matrix.csv
HEATMAP=manyGenomes/fromAtlas/pankmer_clustermap.svg

echo "Date: $(date --iso-8601=seconds)"
echo "Job: $SLURM_JOB_NAME  ID: $SLURM_JOB_ID"
echo "Host: $HOSTNAME"
echo

ml apptainer

thr=$SLURM_CPUS_PER_TASK

mkdir -p joblog


# Step 1: Build k-mer index
echo "=== Step 1: pankmer index ==="
apptainer exec pankmer.sif pankmer index -g manyGenomes/fromAtlas/*.fa.gz -o "$INDEX" -t $thr --time
echo "Index complete: $INDEX"
echo

# Step 2: Pairwise adjacency matrix
echo "=== Step 2: pankmer adj-matrix ==="
apptainer exec pankmer.sif pankmer adj-matrix -i "$INDEX" -o "$MATRIX" --time
echo "Matrix complete: $MATRIX"
echo

# Step 3: Clustered heatmap
echo "=== Step 3: pankmer clustermap ==="
apptainer exec pankmer.sif pankmer clustermap -i "$MATRIX" -o "$HEATMAP"
echo "Heatmap complete: $HEATMAP"
echo

echo "$(date --iso-8601=seconds) pankmer complete"

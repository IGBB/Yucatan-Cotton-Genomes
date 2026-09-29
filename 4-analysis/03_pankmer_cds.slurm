#!/bin/bash

#SBATCH --job-name=pankmer_cds
#SBATCH --output=joblog/%x_%j.log
#SBATCH --error=joblog/%x_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH -A jfw-lab

set -euo pipefail

GENOMES=(
    cds/Chorlitos.cds.fa
    cds/Consorcio.cds.fa
    cds/Preseco.cds.fa
    cds/Pretropical.cds.fa
    cds/TX2094.cds.fa
)

INDEX=pankmer/cds_index.tar
MATRIX=pankmer/pankmer_cds_adj_matrix.csv
HEATMAP=pankmer/pankmer_cds_clustermap.svg

echo "Date: $(date --iso-8601=seconds)"
echo "Job: $SLURM_JOB_NAME  ID: $SLURM_JOB_ID"
echo "Host: $HOSTNAME"
echo

ml apptainer

thr=$SLURM_CPUS_PER_TASK

mkdir -p pankmer joblog

for fa in "${GENOMES[@]}"; do
    [[ -f "$fa" ]] || { echo "ERROR: missing CDS FASTA: $fa"; exit 1; }
done

# Step 1: Build k-mer index
echo "=== Step 1: pankmer index ==="
apptainer exec pankmer.sif pankmer index -g "${GENOMES[@]}" -o "$INDEX" -t $thr --time
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

echo "$(date --iso-8601=seconds) pankmer_cds complete"

#!/bin/bash

#SBATCH --job-name=synteny_pairwise
#SBATCH --output=joblog/%x_%A_%a.log
#SBATCH --error=joblog/%x_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=18:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH -A jfw-lab
#SBATCH --array=0-23,25-29

# Pairwise synteny dotplots: all Yuc × Ref genome combinations.
# Pre-filter each genome to chromosomes-only with filter_chromosomes.sh before submitting.
#
# Outputs (per pair):
#   $out_dir/{yuc}_vs_{ref}.paf   PAF alignment
#   $out_dir/{yuc}_vs_{ref}.html  dotplot

set -euo pipefail

# ── genome sets ───────────────────────────────────────────────────────────────
# Edit paths to point at your chromosome-filtered FASTAs (output of filter_chromosomes.sh).

YUC_NAMES=(Chorlitos     Consorcio     Preseco     Pretropical     TX2094)
YUC_FA=(
    genomes/Chorlitos.fa
    genomes/Consorcio.fa
    genomes/Preseco.fa
    genomes/Pretropical.fa
    genomes/TX2094.fa
)

REF_NAMES=(TX2094         Latifolium         Jin668         TM1         YZ1         ZM113)
REF_FA=(
    genomes/TX2094.fa
    genomes/Latifolium.26.fa
    genomes/Jin668.T2T.fasta
    genomes/TM1.T2T.26.fa
    genomes/YZ1.T2T.fasta
    genomes/ZM113.T2T.fa
)

# ── parameters ────────────────────────────────────────────────────────────────
asm=asm5          # minimap2 preset: asm5 / asm10 / asm20
min=100000          # dotplotly minimum alignment length (-q)
out_dir=dotplots   # output directory for PAFs and HTMLs

# ── pair mapping ──────────────────────────────────────────────────────────────
# Task ID → (yuc_idx, ref_idx) via row-major order over 5×6 matrix.
n_ref=${#REF_NAMES[@]}
yuc_idx=$(( SLURM_ARRAY_TASK_ID / n_ref ))
ref_idx=$(( SLURM_ARRAY_TASK_ID % n_ref ))

yuc_name="${YUC_NAMES[$yuc_idx]}"
ref_name="${REF_NAMES[$ref_idx]}"
yuc_fa="${YUC_FA[$yuc_idx]}"
ref_fa="${REF_FA[$ref_idx]}"

pair="${yuc_name}_vs_${ref_name}"

# ── environment ───────────────────────────────────────────────────────────────
echo "Date: $(date --iso-8601=seconds)"
echo "Job: $SLURM_JOB_NAME  Array: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Host: $HOSTNAME"
echo "Pair: $pair  (task $SLURM_ARRAY_TASK_ID)"
echo

ml apptainer samtools r/4.5.1-py311-ggqalmr
source /lustre/hdd/LAS/mhufford-lab/corrinne/app/env/genome_assembly

thr=$SLURM_CPUS_PER_TASK
echo "Threads: $thr  asm: $asm  min: $min"
echo "Yuc: $yuc_fa"
echo "Ref: $ref_fa"
echo

# ── checks ────────────────────────────────────────────────────────────────────
[[ -f "$yuc_fa" ]] || { echo "ERROR: missing Yuc FASTA: $yuc_fa"; exit 1; }
[[ -f "$ref_fa" ]] || { echo "ERROR: missing Ref FASTA: $ref_fa"; exit 1; }

mkdir -p "$out_dir" joblog

paf="${out_dir}/${pair}.paf"

# ── minimap2 ──────────────────────────────────────────────────────────────────
echo "Start minimap2: $pair"
minimap2 -t "$thr" -x "$asm" -c "$ref_fa" "$yuc_fa" > "$paf"
echo "minimap2 complete"

# Dotplots are rendered in batch by 02b_dotplot_pdf.R after all array tasks finish.

echo "$(date --iso-8601=seconds) task ${SLURM_ARRAY_TASK_ID} ($pair) complete"

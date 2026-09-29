#!/bin/bash
# filter_chromosomes.sh — Extract the first N sequences from a FASTA (chromosomes only).
# Usage: bash filter_chromosomes.sh input.fa output.fa [n=26]
#
# Run this once per genome before submitting synteny_pairwise.slurm.
# Requires: samtools

set -euo pipefail

input="${1:?Usage: $0 input.fa output.fa [n_chromosomes]}"
output="${2:?Usage: $0 input.fa output.fa [n_chromosomes]}"
n_chr="${3:-26}"

[[ -f "$input" ]] || { echo "ERROR: missing input FASTA: $input"; exit 1; }

echo "Indexing $input"
samtools faidx "$input"

mapfile -t chr_names < <(awk -v n="$n_chr" 'NR <= n { print $1 }' "${input}.fai")

echo "Extracting ${#chr_names[@]} sequences -> $output"
samtools faidx "$input" "${chr_names[@]}" > "$output"
samtools faidx "$output"

echo "Done: $output"

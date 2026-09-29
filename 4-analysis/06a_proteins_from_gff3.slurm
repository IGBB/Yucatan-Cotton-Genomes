#!/bin/bash

#SBATCH --job-name=proteins_from_gff3
#SBATCH --output=joblog/%x_%A_%a.log
#SBATCH --error=joblog/%x_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH --array=0-4
#SBATCH -A jfw-lab

ml gffread

samples=(Chorlitos Consorcio Preseco Pretropical TX2094)
name=${samples[$SLURM_ARRAY_TASK_ID]}

# Step 1: GFF3 → fixed GTF  (-T writes GTF format)
# Uncomment if starting from GFF3; skip if .fixed.gtf already exists in gtf/
# gffread gff/${name}.gff3 -T -o gtf/${name}.gtf
# awk '$3 == "transcript" {
#     $9 = sprintf("transcript_id \"%s\"", $9)} 1' \
#     FS="\t" OFS="\t" gtf/${name}.gtf > gtf/${name}.fixed.gtf

# Step 2: fixed GTF + genome FASTA → protein FASTA
echo "Processing: $name"
gffread gtf/${name}.fixed.gtf \
    -g genomes/${name}.fa \
    -y proteins/${name}.proteins.fa

# Step 3: fixed GTF + genome FASTA → CDS FASTA
gffread gtf/${name}.fixed.gtf \
    -g genomes/${name}.fa \
    -x cds/${name}.cds.fa

echo "Done: $name"

# Analysis Code and Output

Scripts and selected output files supporting the analysis of four new *Gossypium hirsutum* genome assemblies (Chorlitos, Consorcio, Preseco, Pretropical) and the reference accession TX2094.

---

## scripts

### Pairwise synteny

| File | Description |
|------|-------------|
| `01_filter_chromosomes.sh` | Subsets each genome FASTA to the 26 named chromosomes (Ah_01–Ah_13, Dh_01–Dh_13); run once per genome before submitting `02a` |
| `02a_synteny_pairwise.slurm` | Array job running minimap2 pairwise whole-genome alignments for all genome pairs |
| `02b_dotplot_pdf.R` | Renders PAF alignments into dotplots using ggplot2 and patchwork; excludes alignments shorter than 100 kb |

### Genome uniqueness and annotation similarity (PanKmer)

| File | Description |
|------|-------------|
| `04a_pankmer_atlas.slurm` | Builds a k-mer index, pairwise adjacency matrix, and clustered heatmap SVG for all 167 *G. hirsutum* genomes (4 new + TX2094 + 162 from CottonGen); reads from and writes to `manyGenomes/fromAtlas/` |
| `04b_pankmer_heatmap_atlas.R` | Produces main (~60 genome) and supplemental (~167 genome) Jaccard similarity heatmaps from the atlas adjacency matrix using pheatmap |
| `03_pankmer_cds.slurm` | Runs PanKmer on CDS sequences to estimate annotation similarity among the five focal genomes; writes index, adjacency matrix, and SVG heatmap to `pankmer/` |

### Orthology and synteny (GENESPACE)

| File | Description |
|------|-------------|
| `05a_genespace.slurm` | Runs GENESPACE v1.2.3 on the five focal genomes with A and D subgenomes treated independently; calls `05b_genespace.R` |
| `05b_genespace.R` | GENESPACE configuration and execution: `query_pangenes` writes `genespace/pangenes/pangene_db.tsv`; `plot_riparian` produces per-genome riparian plots |
| `05c_reformatMatrix.R` | Builds the presence/absence (PAV) matrix from `pangene_db.tsv` and saves the family-count and private-gene summary tables |

### NLR immune receptor annotation (Resistify)

| File | Description |
|------|-------------|
| `06a_proteins_from_gff3.slurm` | Array job extracting protein and CDS sequences from each BRAKER3 GTF + genome FASTA |
| `06b_resistify.slurm` | Array job running Resistify on predicted protein sequences from each genome; models cached locally via `--bind` mount |
| `06c_resistify_pav.R` | Joins Resistify NLR predictions with the GENESPACE pangene table to compute presence/absence variation (PAV) across the five genomes; saves `resistify/nlr_pg_braker.rds` |
| `06d_resistify_braker_plots.R` | Publication figures: NLR counts by PAV class, NLR classification by genome, A vs. D subgenome breakdown, MADA/CJID motif prevalence |

### NLR immune receptor annotation (NLR-Annotator) and integration

| File | Description |
|------|-------------|
| `07a_nlr_annotator.slurm` | Array job running NLR-Annotator v2.1b on each genome FASTA to identify NLR loci from sequence motifs; requires the NLR-Annotator JAR and motif files cloned to `NLR-Annotator/` |
| `07b_nlr_annotator_compare.R` | Intersects BRAKER3 gene models with NLR-Annotator loci; classifies each gene as identified by both methods, Resistify only, or NLR-Annotator only; integrates results with the GENESPACE pangene table; identifies NLR clusters (consecutive NLRs on the same chromosome separated by ≤100 kb); produces summary tables and figures |

---

## output/

### dotplots/

| File | Description |
|------|-------------|
| `FourNewCottonGenomes_SuppFig3.jpg` | Tiled overview of pairwise whole-genome synteny for the four new genomes vs. TX2094 and five reference accessions |
| `all_comparisons.pdf` | Full-resolution per-comparison dotplots (one page per genome pair) |

### pankmer/

| File | Description |
|------|-------------|
| `FourNewCottonGenomes_Figure2.jpg` | Main-text hierarchically clustered k-mer similarity heatmap (~60 genomes; Jing and Yang107 sets excluded) |
| `pankmer_heatmap_atlas_main.pdf` | PDF version of the main heatmap |
| `pankmer_heatmap_atlas_supp.pdf` | Supplemental heatmap including all ~167 genomes |
| `FourNewCottonGenomes_SuppFig4.jpg` | CDS-based k-mer similarity heatmap for the five focal genomes |
| `pankmer_cds_clustermap.svg` | SVG version of the CDS similarity heatmap |
| `pankmer_cds_adj_matrix.csv` | Raw pairwise CDS k-mer adjacency matrix for the five focal genomes |

### genespace_riparian/

Riparian plots from GENESPACE showing syntenic block structure for each genome relative to TX2094. 

### genespace_pav/

Pangene extraction from `genespace/pangenes/pangene_db.tsv`, produced by `05c_reformatMatrix.R`.

| File | Description |
|------|-------------|
| `pav_summary_table.tsv` | Gene family counts by number of genomes present (one–five), plus `SCfive` (single-copy in all five genomes) |
| `unique_genes_table.tsv` | Genes private to each genome — total families, and Gohir.A / Gohir.D breakdown |
| `unique_genes_pct_table.tsv` | Private genes as a percentage of each genome's total GTF gene count |
| `unique_genes.tsv` | Gene-level IDs of every private gene, by genome |

### resistify/

| File | Description |
|------|-------------|
| `resistifyGenespacePAV.txt` | Summary of NLR PAV classes (core / near-core / variable / private) from Resistify × GENESPACE integration |
| `nlr_gene_pgrow.tsv` | Per-gene NLR PAV table with pangene family assignments |
| `panelA_pav_per_genome.jpg` | NLR gene counts per genome, stacked by PAV class (core / variable / private) |
| `panelB_class_per_genome.jpg` | NLR classification (CNL, TNL, etc.) per genome, shaded by PAV class |
| `panelC_subgenome_class.jpg` | A vs. D subgenome NLR classification breakdown per genome |
| `panelD_mada_cjid.jpg` | MADA and CJID motif prevalence by NLR classification |

### nlr_annotator/

Summary tables and figures from `07b_nlr_annotator_compare.R`.

| File | Description |
|------|-------------|
| `nlr_compare.tsv` | Per-gene classification table: genome, gene ID, coordinates, Resistify class, NLR-Annotator locus, and category (both / resistify_only / annotator_gene) |
| `nlr_compare_summary.tsv` | Per-genome counts of each category plus total Resistify NLRs, total NLR-Annotator loci, and unannotated loci |
| `nlr_pangene_pav.tsv` | Per-pangene-family NLR statistics: PAV class, number of genomes with double-confirmed NLRs, number with any NLR evidence |
| `nlr_clusters.tsv` | Per-gene cluster assignments (cluster ID, size, category) |
| `nlr_cluster_summary.tsv` | Per-genome cluster summary: total NLRs, number of clusters, singleton count, cluster size breakdown |
| `nlr_cluster_coords.tsv` | Per-cluster coordinates: chromosome, span (kb), size, method breakdown |
| `nlr_chr_distribution.jpg` | Bar chart of NLR pangene families per chromosome, faceted by PAV class (A and D subgenomes side by side) |
| `nlr_chr_locations.jpg` | Dot plot of NLR family positions along all 26 chromosomes, coloured by PAV class |
| `nlr_tracks_chr11.jpg` | Linear genome track plot for chromosome 11 (A and D), one track per genome, showing NLR positions |
| `nlr_unannotated.tsv` | NLR-Annotator loci with no overlapping BRAKER3 gene model |

# Software Used

## Sequence alignment & synteny
| Program | Version | Role | Used in |
|---|---|---|---|
| minimap2 | v2.28 | Pairwise whole-genome alignment (`-x asm5`, `-c`) | `02a_synteny_pairwise.slurm` |
| samtools | v1.23.1 | Indexing / general sequence handling | `02a_synteny_pairwise.slurm` |

## Orthology / pangenes
| Program | Version | Role | Used in |
|---|---|---|---|
| GENESPACE | v1.2.3 | Synteny-anchored orthology and pangene table construction | `05a_genespace.slurm`, `05b_genespace.R` |
| gffread | v0.12.7 | Protein/CDS extraction | `06a_proteins_from_gff3.slurm` |

## K-mer similarity
| Program | Version | Role | Used in |
|---|---|---|---|
| PanKmer | v0.20.4 | k-mer index, adjacency matrix, heatmap | `03_pankmer_cds.slurm`, `04a_pankmer_atlas.slurm` |

## NLR / immune receptor annotation
| Program | Version | Role | Used in |
|---|---|---|---|
| Resistify | v2.1.0  | Protein-based NLR classification | `06b_resistify.slurm`, `06c_resistify_pav.R` |
| NLR-Annotator | v2.1b | Motif-based NLR locus detection | `07a_nlr_annotator.slurm` |
| OpenJDK | 21.0.3_9-vngib7s | Java runtime for NLR-Annotator | `07a_nlr_annotator.slurm` |

## R environment
| Program | Version |
|---|---|
| R | 4.5.1  |
| R | 4.4.3  |

### R packages (local R 4.4.3)
| Package | Version |
|---|---|
| data.table | 1.17.0 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.0 |
| pheatmap | 1.0.12 |
| RColorBrewer | 1.1.3 |



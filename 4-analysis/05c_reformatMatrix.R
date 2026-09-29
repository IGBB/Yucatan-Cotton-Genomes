# 05c_reformatMatrix.R
#
# Pangene extraction from the GENESPACE pangene table: builds the
# presence/absence (PAV) matrix across the five focal genomes and saves
# the summary tables that back the paper's PAV / core-gene numbers.
#
# Adapted from the root-level 05_reformatMatrix.R, which only printed
# summary_table / unique_table to the console — nothing was ever saved
# to a file. This version writes everything out.
#
# Input:  genespace/pangenes/pangene_db.tsv  (from 05a/05b_genespace)
# Output (FourNewCottonGenomes/output/genespace_pav/):
#   pav_summary_table.tsv       — family counts by # genomes present (1-5) + SCfive
#   unique_genes_table.tsv      — genes private to each genome (families / Gohir.A / Gohir.D)
#   unique_genes_pct_table.tsv  — private genes as % of each genome's GTF gene count
#   unique_genes.tsv            — gene-level IDs of all private genes, by genome
#
# Run from the fromTony working directory:
#   Rscript FourNewCottonGenomes/scripts/05c_reformatMatrix.R

library(data.table)

PANGENE_DB <- "genespace/pangenes/pangene_db.tsv"
OUT_DIR    <- "FourNewCottonGenomes/output/genespace_pav"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pg <- fread(PANGENE_DB)

# Identify subgenome columns
at_cols <- grep("\\.At$", names(pg), value = TRUE)
dt_cols <- grep("\\.Dt$", names(pg), value = TRUE)
genomes <- sub("\\.At$", "", at_cols)   # Chorlitos, Consorcio, etc.

# Stack A and D subgenomes into one set of genome columns
pgAh <- pg[, c("interpChr", at_cols, dt_cols), with = FALSE]
pgDh <- copy(pgAh)
pgDh[, interpChr := gsub("Ah", "Dh", interpChr)]
pgAh[, (dt_cols) := NULL]
pgDh[, (at_cols) := NULL]
setnames(pgAh, at_cols, genomes)
setnames(pgDh, dt_cols, genomes)
pg_combined <- rbind(pgAh, pgDh)

# PAV matrix — presence = any Gohir gene in that genome's cell
pav <- pg_combined[, lapply(.SD, function(x) as.integer(grepl("Gohir", x))), .SDcols = genomes]
keep <- rowSums(pav) > 0
pav     <- pav[keep]
pg_filt <- pg_combined[keep]

# Single-copy flag: any cell in the row contains | (multi-gene)
has_pipe <- rowSums(pg_filt[, lapply(.SD, function(x) grepl("\\|", x)), .SDcols = genomes]) > 0

# SCfive: all 5 genomes present, no multi-copy gene in any column
sc_five <- sum(rowSums(pav) == 5 & !has_pipe)

# ── PAV summary table ─────────────────────────────────────────────────────────
hist_counts <- tabulate(rowSums(pav), nbins = length(genomes))
summary_table <- data.frame(
  one    = hist_counts[1],
  two    = hist_counts[2],
  three  = hist_counts[3],
  four   = hist_counts[4],
  five   = hist_counts[5],
  SCfive = sc_five
)
fwrite(summary_table, file.path(OUT_DIR, "pav_summary_table.tsv"), sep = "\t")
cat("PAV summary table:\n"); print(summary_table)

# ── Unique (private) genes per genome ─────────────────────────────────────────
unique_rows        <- rowSums(pav) == 1
unique_per_genome  <- colSums(pav[unique_rows])

count_gohir_type <- function(genome, pattern) {
  cells <- pg_filt[[genome]][unique_rows & pav[[genome]] == 1L]
  sum(grepl(pattern, unlist(strsplit(as.character(cells), "\\|"))))
}
unique_A <- sapply(genomes, count_gohir_type, pattern = "Gohir\\.A")
unique_D <- sapply(genomes, count_gohir_type, pattern = "Gohir\\.D")

unique_table <- rbind(families = unique_per_genome, GohirA = unique_A, GohirD = unique_D)
fwrite(as.data.table(unique_table, keep.rownames = "metric"),
       file.path(OUT_DIR, "unique_genes_table.tsv"), sep = "\t")
cat("\nUnique (private) genes per genome:\n"); print(unique_table)

# Gene counts from GTF files (ground truth — CLAUDE.md, alphabetical genome order)
gtf_table <- rbind(
  total  = c(77610, 77232, 80273, 77620, 85371),
  GohirA = c(38216, 38146, 40166, 38130, 41720),
  GohirD = c(39394, 39086, 40107, 39490, 43651)
)
colnames(gtf_table) <- c("Chorlitos", "Consorcio", "Preseco", "Pretropical", "TX2094")

pct_table <- unique_table / gtf_table[c("total", "GohirA", "GohirD"), genomes] * 100
rownames(pct_table) <- c("families_pct", "GohirA_pct", "GohirD_pct")
fwrite(as.data.table(pct_table, keep.rownames = "metric"),
       file.path(OUT_DIR, "unique_genes_pct_table.tsv"), sep = "\t")
cat("\nUnique genes as % of each genome's total GTF gene count:\n"); print(pct_table)

# ── Gene-level IDs of all private genes, by genome ────────────────────────────
unique_genes <- rbindlist(lapply(genomes, function(g) {
  rows  <- unique_rows & pav[[g]] == 1L
  cells <- pg_filt[[g]][rows]
  ids   <- unlist(strsplit(as.character(cells), "\\|"))
  ids   <- ids[grepl("Gohir", ids)]
  data.table(genome = g, gene_id = gsub("\\.\\d+$", "", ids))
}))
fwrite(unique_genes, file.path(OUT_DIR, "unique_genes.tsv"), sep = "\t")
cat(sprintf("\nPrivate gene IDs saved: %s (%d genes)\n",
            file.path(OUT_DIR, "unique_genes.tsv"), nrow(unique_genes)))

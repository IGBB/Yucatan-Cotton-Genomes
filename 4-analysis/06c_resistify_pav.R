# 08_resistify_pav.R
#
# Overview:
#   Integrates Resistify NLR annotations with the GENESPACE pangene PAV matrix
#   to characterize presence/absence variation of NLR immune receptor genes
#   across the five Gossypium hirsutum genomes.
#
# Inputs:
#   genespace/pangenes/pangene_db.tsv        — GENESPACE pangene table (all 5 genomes)
#   resistify/{genome}_resistify/results.tsv — Resistify per-gene NLR classification
#                                              for each of the 5 genomes
#
# Outputs:
#   resistify/resistifyGenespacePAV.txt

library(data.table)

# ── Rebuild PAV matrix (mirrors 4_reformatMatrix.R) ──────────────────────────
pg <- fread("genespace/pangenes/pangene_db.tsv")
at_cols <- grep("\\.At$", names(pg), value = TRUE)
dt_cols <- grep("\\.Dt$", names(pg), value = TRUE)
genomes <- sub("\\.At$", "", at_cols)

pgAh <- pg[, c("interpChr", at_cols, dt_cols), with = FALSE]
pgDh <- copy(pgAh)
pgDh[, interpChr := gsub("Ah", "Dh", interpChr)]
pgAh[, (dt_cols) := NULL]
pgDh[, (at_cols) := NULL]
setnames(pgAh, at_cols, genomes)
setnames(pgDh, dt_cols, genomes)
pg_combined <- rbind(pgAh, pgDh)

pav <- pg_combined[, lapply(.SD, function(x) as.integer(grepl("Gohir", x))),
                   .SDcols = genomes]
keep <- rowSums(pav) > 0
pav     <- pav[keep]
pg_filt <- pg_combined[keep]

# ── Load Resistify results ────────────────────────────────────────────────────
nlr_raw <- rbindlist(lapply(genomes, function(g) {
  dt <- fread(file.path("resistify", paste0(g, "_resistify"), "results.tsv"))
  dt[, genome := g]
}))

# Strip genome prefix and isoform suffix to get bare gene IDs
# e.g. Chorlitos_Gohir.D11G380000.2 → Gohir.D11G380000
nlr_raw[, gene_id := sub("\\.[0-9]+$", "",
                         sub("^[^_]+_", "", Sequence))]

# One entry per gene: keep the isoform with the most domains
nlr_raw[, domain_len := nchar(Domains)]
setorder(nlr_raw, gene_id, -domain_len)
nlr <- unique(nlr_raw, by = c("gene_id", "genome"))
nlr_raw[, domain_len := NULL]

# ── Build gene → pg_filt row index map ───────────────────────────────────────
gene_map <- rbindlist(lapply(genomes, function(g) {
  vals    <- as.character(pg_filt[[g]])
  split_v <- strsplit(vals, "\\|")
  idx     <- rep(seq_along(split_v), lengths(split_v))
  genes   <- unlist(split_v)
  ok      <- grepl("Gohir", genes)
  data.table(pg_row  = idx[ok],
             gene_id = sub("\\.[0-9]+$", "", sub("^[^_]+_", "", genes[ok])),
             genome  = g)
}))

# ── Join NLRs to pangene rows ─────────────────────────────────────────────────
nlr_pg <- merge(nlr, gene_map, by = c("gene_id", "genome"), all.x = TRUE)

# Attach number of genomes sharing each pangene family
nlr_pg[!is.na(pg_row), n_genomes := rowSums(pav)[pg_row]]

# Subgenome: Ah → A, Dh → D (interpChr contains "Ah" or "Dh")
nlr_pg[!is.na(pg_row),
  subgenome := ifelse(grepl("Ah", pg_filt$interpChr[pg_row]), "A", "D")]

# ── Early diagnostics: gene lists and pangene table ───────────────────────────
# One row per NLR gene (isoforms already combined); pg_row is NA if unplaced
fwrite(nlr_pg[order(genome, gene_id), .(genome, gene_id, pg_row)],
  "resistify/nlr_gene_pgrow.tsv", sep = "\t", na = "NA")

# Full filtered pangene table with row index for cross-referencing
pg_out <- copy(pg_filt)
pg_out[, pg_row := .I]
fwrite(pg_out, "resistify/pg_filt.tsv", sep = "\t", na = "NA")

# PAV class (variable = 2–4 genomes, combining dispensable)
nlr_pg[, pav_class := fcase(
  n_genomes == 5L,    "core",
  n_genomes %in% 2:4, "variable",
  n_genomes == 1L,    "private",
  is.na(pg_row),      "unplaced"
)]
nlr_pg[, pav_class := factor(pav_class,
  levels = c("core", "variable", "private", "unplaced"))]

# ── Summaries ─────────────────────────────────────────────────────────────────
out <- file("resistify/resistifyGenespacePAV.txt", open = "wt")

cat("=== Total unique NLR genes per genome (all, regardless of pangene matrix) ===\n",
    file = out)
write.table(nlr[, .N, by = genome][order(genome)],
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)
cat("\n", file = out)

cat("=== NLR classification per genome (all) ===\n", file = out)
write.table(
  dcast(nlr, genome ~ Classification,
        fun.aggregate = length, value.var = "gene_id", fill = 0L),
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)
cat("\n", file = out)

cat("=== Pangene matrix coverage ===\n", file = out)
cat("Matched to pangene matrix:", sum(!is.na(nlr_pg$pg_row)), "\n", file = out)
cat("Unplaced (not in matrix):", sum(is.na(nlr_pg$pg_row)), "\n\n", file = out)

cat("NLR Classification x PAV class\n", file = out)
write.table(dcast(nlr_pg[!is.na(pg_row)],
  Classification ~ pav_class, fun.aggregate = length, value.var = "gene_id"),
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)

cat("\nPer-genome NLR counts by PAV class\n", file = out)
write.table(dcast(nlr_pg[!is.na(pg_row)],
  genome ~ pav_class, fun.aggregate = length, value.var = "gene_id"),
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)

cat("\nSubgenome x PAV class\n", file = out)
write.table(dcast(nlr_pg[!is.na(pg_row)],
  subgenome ~ pav_class, fun.aggregate = length, value.var = "gene_id"),
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)

# ── Family-level NLR PAV ──────────────────────────────────────────────────────
# A pangene row is an "NLR row" if any gene in it was classified as NLR by
# Resistify. Representation = how many genomes have *any* gene in that row
# (from the full pav matrix, not just NLR genes), so paralogs don't inflate
# the count. Core/variable/private are then counted at the family level.

nlr_row_ids <- sort(unique(nlr_pg[!is.na(pg_row), pg_row]))
pav_nlr     <- pav[nlr_row_ids]
pav_nlr[, pg_row    := nlr_row_ids]
pav_nlr[, n_rep     := rowSums(.SD), .SDcols = genomes]
pav_nlr[, row_pav   := fcase(
  n_rep == 5L,    "core",
  n_rep %in% 2:4, "variable",
  n_rep == 1L,    "private"
)]

cat("\n=== Family-level NLR PAV (one count per pangene row) ===\n", file = out)
cat("Total NLR pangene families:", nrow(pav_nlr), "\n", file = out)
write.table(
  pav_nlr[, .N, by = row_pav][order(factor(row_pav,
    levels = c("core", "variable", "private")))],
  file = out, quote = FALSE, sep = "\t", row.names = FALSE)
cat("\n", file = out)

cat("Family-level NLR counts per genome\n", file = out)
fam_counts <- rbindlist(lapply(genomes, function(g) {
  data.table(
    genome   = g,
    core     = sum(pav_nlr$n_rep == 5L),
    variable = sum(pav_nlr$n_rep %in% 2:4 & pav_nlr[[g]] == 1L),
    private  = sum(pav_nlr$n_rep == 1L    & pav_nlr[[g]] == 1L)
  )
}))
write.table(fam_counts, file = out, quote = FALSE, sep = "\t", row.names = FALSE)

close(out)
saveRDS(nlr_pg, "resistify/nlr_pg_braker.rds")
message("Report written to resistify/resistifyGenespacePAV.txt")

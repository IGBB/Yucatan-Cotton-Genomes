# 17_nlr_annotator_compare.R
#
# Cross-reference NLR-Annotator loci with BRAKER gene models and Resistify NLR
# predictions across the five core genomes.
#
# Inputs:
#   gtf/{genome}.fixed.gtf                       — BRAKER gene models
#   nlr_annotator/{genome}/{genome}.nlr.txt       — NLR-Annotator loci (1-based coords)
#   resistify/{genome}_resistify/results.tsv      — Resistify NLR classifications
#
# Outputs:
#   nlr_annotator/nlr_compare.tsv         — per-gene classification table
#   nlr_annotator/nlr_compare_summary.tsv — per-genome category counts (wide)
#   nlr_annotator/nlr_unannotated.tsv     — NLR-Annotator loci with no overlapping gene

library(data.table)
library(ggplot2)

genomes <- c("TX2094", "Chorlitos", "Consorcio", "Preseco", "Pretropical")

# ── 1. Gene models from GTF ────────────────────────────────────────────────────
# Gene features in .fixed.gtf have a bare gene ID in col 9 (no key="value" quoting).
# Strip the genome-specific prefix (everything up to and including the first _) to
# get canonical Gohir.XNNGNNNNNNN IDs that match Resistify.
load_genes <- function(genome) {
  f <- file.path("gtf", paste0(genome, ".fixed.gtf"))
  gtf <- fread(f, sep = "\t", header = FALSE, quote = "",
               col.names = c("chr","source","feature","start","end",
                             "score","strand","frame","attr"),
               showProgress = FALSE)
  gtf <- gtf[!startsWith(chr, "#") & feature == "gene"]
  gtf[, prefixed_id := attr]
  gtf[, gene_id := sub("^[^_]+_", "", attr)]   # Gohir.A01G000100
  gtf[, .(genome = genome, chr,
          start   = as.integer(start),
          end     = as.integer(end),
          strand, gene_id, prefixed_id)]
}

# ── 2. NLR-Annotator loci ─────────────────────────────────────────────────────
# TXT columns: chr  nlr_id  nlr_class  start  end  strand  motifs
# Coordinates appear to use the same scale as the BED output; treat as 1-based
# (a 1-bp edge difference is immaterial for gene-level overlaps).
load_nlr_annotator <- function(genome) {
  f <- file.path("nlr_annotator", genome, paste0(genome, ".nlr.txt"))
  loci <- fread(f, sep = "\t", header = FALSE,
                col.names = c("chr","nlr_id","nlr_class","start","end","strand","motifs"),
                showProgress = FALSE)
  loci[, `:=`(start = as.integer(start), end = as.integer(end), genome = genome)]
  loci[, .(genome, chr, start, end, strand, nlr_id, nlr_class)]
}

# ── 3. Resistify NLR genes ────────────────────────────────────────────────────
# results.tsv columns: Sequence  Length  LRR_Length  Motifs  Domains
#                      Classification  NBARC_motifs  MADA  CJID
# Sequence format: {prefix}_Gohir.XNNGNNNNNNN.isoform
# Non-NLR rows have empty Classification.
load_resistify <- function(genome) {
  f <- file.path("resistify", paste0(genome, "_resistify"), "results.tsv")
  res <- fread(f, showProgress = FALSE)
  res[, gene_id := sub("\\.[0-9]+$", "", sub("^[^_]+_", "", Sequence))]
  res[, genome  := genome]
  res <- res[!is.na(Classification) & Classification != ""]
  res[, .(genome, gene_id, resistify_class = Classification, MADA, CJID)]
}

# ── 4. Load all data ───────────────────────────────────────────────────────────
cat("Loading gene models...\n")
genes <- rbindlist(lapply(genomes, load_genes))
cat(sprintf("  %d genes across %d genomes\n", nrow(genes), length(genomes)))

cat("Loading NLR-Annotator loci...\n")
nlr_loci <- rbindlist(lapply(genomes, load_nlr_annotator))
cat("  Loci per genome:\n")
print(nlr_loci[, .N, by = genome])

cat("Loading Resistify NLR genes...\n")
resistify_nlr <- rbindlist(lapply(genomes, load_resistify))
cat("  Resistify NLR genes per genome:\n")
print(resistify_nlr[, .N, by = genome])

# ── 5. Overlap: which genes fall inside an NLR-Annotator locus? ───────────────
# Process per genome to avoid chromosome name collisions across assemblies.
# A gene is considered overlapping if any part of [gene.start, gene.end]
# intersects [locus.start, locus.end] (type = "any").

gene_with_locus <- rbindlist(lapply(genomes, function(g) {
  g_genes <- genes[genome == g, .(chr, start, end, strand, gene_id)]
  g_loci  <- nlr_loci[genome == g, .(chr, start, end, nlr_id, nlr_class)]

  if (nrow(g_loci) == 0) {
    return(g_genes[, .(chr, start, end, strand, gene_id,
                       nlr_id = NA_character_, nlr_class = NA_character_)])
  }

  setkey(g_loci, chr, start, end)
  ov <- foverlaps(g_genes, g_loci,
                  by.x = c("chr", "start", "end"),
                  by.y = c("chr", "start", "end"),
                  type = "any", nomatch = NA, mult = "first")
  ov[, .(chr, start, end, strand, gene_id, nlr_id, nlr_class)]
}), idcol = "genome")
gene_with_locus[, genome := genomes[genome]]

# ── 6. Overlap: which NLR-Annotator loci have no gene? ────────────────────────
loci_covered <- rbindlist(lapply(genomes, function(g) {
  g_genes <- genes[genome == g, .(chr, start, end, gene_id)]
  g_loci  <- nlr_loci[genome == g, .(chr, start, end, nlr_id, nlr_class)]

  if (nrow(g_genes) == 0 || nrow(g_loci) == 0) return(data.table(nlr_id = character()))

  setkey(g_genes, chr, start, end)
  ov <- foverlaps(g_loci, g_genes,
                  by.x = c("chr", "start", "end"),
                  by.y = c("chr", "start", "end"),
                  type = "any", nomatch = NA, mult = "first")
  ov[!is.na(gene_id), .(nlr_id)]
}))

unannotated <- nlr_loci[!nlr_id %in% loci_covered$nlr_id]

# ── 7. Classify each gene ─────────────────────────────────────────────────────
result <- merge(gene_with_locus, resistify_nlr, by = c("genome", "gene_id"), all.x = TRUE)

result[, category := fcase(
  !is.na(resistify_class) & !is.na(nlr_id), "both",
  !is.na(resistify_class) &  is.na(nlr_id), "resistify_only",
  is.na(resistify_class)  & !is.na(nlr_id), "annotator_gene",
  default = NA_character_
)]

# Keep only genes relevant to NLR analysis
nlr_result <- result[!is.na(category)]

# ── 8. Summary ────────────────────────────────────────────────────────────────
cat("\n=== Classification summary ===\n")
summary_long <- nlr_result[, .N, by = .(genome, category)]
summary_wide <- dcast(summary_long, genome ~ category, value.var = "N", fill = 0L)

# Add totals
summary_wide <- merge(summary_wide,
  resistify_nlr[, .(resistify_total = .N), by = genome],
  by = "genome", all.x = TRUE)
summary_wide <- merge(summary_wide,
  nlr_loci[, .(nlr_annotator_loci = .N), by = genome],
  by = "genome", all.x = TRUE)
summary_wide <- merge(summary_wide,
  unannotated[, .(unannotated_loci = .N), by = genome],
  by = "genome", all.x = TRUE)

print(summary_wide)

cat("\nUnannotated NLR-Annotator loci (no overlapping gene):",
    nrow(unannotated), "\n")
cat("  By genome:\n")
print(unannotated[, .N, by = genome])

# ── 9. NLR-Annotator class breakdown for 'both' genes ────────────────────────
cat("\n=== NLR-Annotator class × Resistify class (genes in 'both') ===\n")
print(nlr_result[category == "both",
                 .N, by = .(nlr_class, resistify_class)][order(-N)])

# ── 10. Pangene PAV for "both" confirmed NLR families ────────────────────────
# For every pangene family containing at least one "both" gene, compute:
#   n_genomes_pav  — how many of the 5 genomes have any gene in this family
#   n_genomes_both — how many of those genomes also have a "both"-confirmed NLR
#
# n_genomes_pav tells you the PAV class of the NLR family (core / near-core / variable).
# n_genomes_both tells you in how many genomes the NLR identity was independently
# confirmed by both protein-based (Resistify) and genome-scan (NLR-Annotator) approaches.
# A family where n_genomes_pav=5 but n_genomes_both=3 means 2 genomes have the gene
# but it was not doubly-confirmed there — potentially diverged or pseudogenized orthologs.

cat("\n=== Pangene PAV for confirmed NLR families ===\n")

# Rebuild pg_filt — same logic as 08_reformatMatrix.R / 08e_nlr_cds_families.R
pg      <- fread("genespace/pangenes/pangene_db.tsv", showProgress = FALSE)
at_cols <- grep("\\.At$", names(pg), value = TRUE)
dt_cols <- grep("\\.Dt$", names(pg), value = TRUE)
genomes_pg <- sub("\\.At$", "", at_cols)

pgAh <- pg[, c("interpChr", "interpOrd", at_cols, dt_cols), with = FALSE]
pgDh <- copy(pgAh)
pgDh[, interpChr := gsub("Ah", "Dh", interpChr)]
pgAh[, (dt_cols) := NULL]
pgDh[, (at_cols) := NULL]
setnames(pgAh, at_cols, genomes_pg)
setnames(pgDh, dt_cols, genomes_pg)
pg_combined <- rbind(pgAh, pgDh)

keep    <- rowSums(pg_combined[, lapply(.SD, function(x) as.integer(grepl("Gohir", x))),
                               .SDcols = genomes_pg]) > 0
pg_filt <- pg_combined[keep]
pg_filt[, pg_row := .I]

# Build gene_id → pg_row lookup (strips prefix + isoform suffix, same as other scripts)
gene_pg <- rbindlist(lapply(genomes_pg, function(g) {
  vals  <- as.character(pg_filt[[g]])
  split <- strsplit(vals, "\\|", fixed = FALSE)
  raw   <- unlist(split)
  idx   <- rep(seq_len(nrow(pg_filt)), lengths(split))
  ok    <- grepl("Gohir", raw)
  if (!any(ok)) return(NULL)
  data.table(
    genome  = g,
    gene_id = sub("\\.[0-9]+$", "", sub("^[^_]+_", "", raw[ok])),
    pg_row  = idx[ok]
  )
}))
gene_pg <- unique(gene_pg)

# Join "both" genes onto their pangene family
both_pg <- merge(
  nlr_result[category == "both", .(genome, gene_id)],
  gene_pg,
  by = c("genome", "gene_id"),
  all.x = TRUE
)
cat(sprintf("  'both' genes matched to a pangene family: %d / %d\n",
            sum(!is.na(both_pg$pg_row)), nrow(both_pg)))
if (anyNA(both_pg$pg_row))
  cat(sprintf("  'both' genes with no pangene match: %d (may be in filtered-out families)\n",
              sum(is.na(both_pg$pg_row))))

# Per-family statistics
focus_rows <- unique(both_pg$pg_row[!is.na(both_pg$pg_row)])

pg_sub <- pg_filt[pg_row %in% focus_rows]

pav_mat <- pg_sub[, lapply(.SD, function(x) as.integer(grepl("Gohir", x))),
                  .SDcols = genomes_pg]
pg_sub[, n_genomes_pav := rowSums(pav_mat)]

# Count distinct genomes with a "both" NLR per family
both_per_fam <- both_pg[!is.na(pg_row), .(n_genomes_both = uniqueN(genome)), by = pg_row]

# Count distinct genomes where the family has ANY NLR-Annotator-overlapping gene
# (category "both" OR "annotator_gene") — NLR-Annotator structural evidence regardless
# of whether Resistify also confirmed it
annotator_pg <- merge(
  nlr_result[category %in% c("both", "annotator_gene"), .(genome, gene_id)],
  gene_pg,
  by = c("genome", "gene_id"),
  all.x = TRUE
)
annotator_per_fam <- annotator_pg[!is.na(pg_row) & pg_row %in% focus_rows,
                                   .(n_genomes_nlr_annotator = uniqueN(genome)),
                                   by = pg_row]

# Count distinct genomes where the family has a Resistify-only NLR
# (protein evidence but no NLR-Annotator locus overlap)
resistify_pg <- merge(
  nlr_result[category == "resistify_only", .(genome, gene_id)],
  gene_pg,
  by = c("genome", "gene_id"),
  all.x = TRUE
)
resistify_per_fam <- resistify_pg[!is.na(pg_row) & pg_row %in% focus_rows,
                                   .(n_genomes_resistify_only = uniqueN(genome)),
                                   by = pg_row]

# Count distinct genomes where the family has any NLR evidence from either method
# (union of all three categories — every row in nlr_result)
any_nlr_pg <- merge(
  nlr_result[, .(genome, gene_id)],
  gene_pg,
  by = c("genome", "gene_id"),
  all.x = TRUE
)
any_nlr_per_fam <- any_nlr_pg[!is.na(pg_row) & pg_row %in% focus_rows,
                                .(n_genomes_any_nlr = uniqueN(genome)),
                                by = pg_row]

pg_stats <- merge(
  pg_sub[, .(pg_row, interpChr, interpOrd, n_genomes_pav)],
  both_per_fam,
  by = "pg_row", all.x = TRUE
)
pg_stats <- merge(pg_stats, annotator_per_fam,  by = "pg_row", all.x = TRUE)
pg_stats <- merge(pg_stats, resistify_per_fam,  by = "pg_row", all.x = TRUE)
pg_stats <- merge(pg_stats, any_nlr_per_fam,    by = "pg_row", all.x = TRUE)
pg_stats[is.na(n_genomes_both),           n_genomes_both           := 0L]
pg_stats[is.na(n_genomes_nlr_annotator),  n_genomes_nlr_annotator  := 0L]
pg_stats[is.na(n_genomes_resistify_only), n_genomes_resistify_only := 0L]
pg_stats[is.na(n_genomes_any_nlr),        n_genomes_any_nlr        := 0L]

# One column per genome:
#   "NLR"     — gene present AND predicted NLR by either method (any nlr_result category)
#   "present" — gene present but no NLR evidence from either method
#   NA        — gene absent
nlr_genome_pg <- unique(any_nlr_pg[!is.na(pg_row), .(genome, pg_row)])

presence <- copy(pg_sub[, .(pg_row)])
for (g in genomes_pg) {
  is_present <- grepl("Gohir", pg_sub[[g]])
  is_nlr     <- pg_sub$pg_row %in% nlr_genome_pg[genome == g, pg_row]
  presence[, (g) := fifelse(!is_present, NA_character_,
                             fifelse(is_nlr, "NLR", "present"))]
}
pg_stats <- merge(pg_stats, presence, by = "pg_row", all.x = TRUE)

# Assign PAV class labels
pg_stats[, pav_class := fcase(
  n_genomes_pav == 5, "core",
  n_genomes_pav == 4, "near_core",
  n_genomes_pav == 3, "variable_3",
  n_genomes_pav == 2, "variable_2",
  n_genomes_pav == 1, "private"
)]

cat(sprintf("\n  Pangene families with ≥1 'both' gene: %d\n", nrow(pg_stats)))
cat(sprintf("  of which single-genome (private):     %d\n",
            pg_stats[n_genomes_pav == 1, .N]))

# Cross-tab: n_genomes_pav × n_genomes_both
cat("\n  Cross-tab  n_genomes_pav × n_genomes_both\n")
cat("  (rows = PAV; cols = genomes with Resistify + NLR-Annotator)\n\n")
xtab <- dcast(pg_stats[, .N, by = .(n_genomes_pav, n_genomes_both)],
              n_genomes_pav ~ n_genomes_both,
              value.var = "N", fill = 0L)
print(xtab)

# Cross-tab: n_genomes_pav × n_genomes_nlr_annotator
cat("\n  Cross-tab  n_genomes_pav × n_genomes_nlr_annotator\n")
cat("  (rows = PAV; cols = genomes with any NLR-Annotator overlap)\n\n")
xtab2 <- dcast(pg_stats[, .N, by = .(n_genomes_pav, n_genomes_nlr_annotator)],
               n_genomes_pav ~ n_genomes_nlr_annotator,
               value.var = "N", fill = 0L)
print(xtab2)

# Summarise by PAV class
cat("\n  Summary by PAV class:\n")
pav_summary <- pg_stats[, .(
  n_families           = .N,
  mean_nlr_annotator   = round(mean(n_genomes_nlr_annotator), 2),
  mean_both            = round(mean(n_genomes_both), 2),
  pct_all_nlr_annotator = round(100 * mean(n_genomes_nlr_annotator == n_genomes_pav), 1),
  pct_all_both         = round(100 * mean(n_genomes_both == n_genomes_pav), 1)
), by = pav_class][order(-n_families)]
print(pav_summary)

# ── 11. Write outputs ─────────────────────────────────────────────────────────
fwrite(nlr_result,    "nlr_annotator/nlr_compare.tsv",         sep = "\t")
fwrite(summary_wide,  "nlr_annotator/nlr_compare_summary.tsv", sep = "\t")
fwrite(unannotated,   "nlr_annotator/nlr_unannotated.tsv",     sep = "\t")
fwrite(pg_stats,      "nlr_annotator/nlr_pangene_pav.tsv",     sep = "\t")

cat("\nOutputs written:\n")
cat(sprintf("  nlr_annotator/nlr_compare.tsv         — %d genes\n",          nrow(nlr_result)))
cat(sprintf("  nlr_annotator/nlr_compare_summary.tsv — %d genomes\n",        length(genomes)))
cat(sprintf("  nlr_annotator/nlr_unannotated.tsv     — %d unannotated loci\n", nrow(unannotated)))
cat(sprintf("  nlr_annotator/nlr_pangene_pav.tsv     — %d NLR pangene families\n", nrow(pg_stats)))

# ── 12. Plot: NLR pangene families per chromosome ────────────────────────────
# Ah and Dh subgenomes are plotted side-by-side for each chromosome number.
# Faceted into four panels: all families, core (5 genomes), variable (2–4), private (1).

pg_plot <- copy(pg_stats)
pg_plot[, subgenome := sub("_.*", "", interpChr)]      # "Ah" or "Dh"
pg_plot[, chr_num   := sub("^[AD]h_", "", interpChr)]  # "01" – "13"
pg_plot[, partition := fcase(
  n_genomes_pav == 5, "core",
  n_genomes_pav  > 1, "variable",
  default = "private"
)]

# Duplicate rows for the "all" panel
plot_data <- rbind(
  pg_plot[, .(chr_num, subgenome, partition)],
  pg_plot[, .(chr_num, subgenome, partition = "all")]
)
plot_data[, partition := factor(partition, levels = c("all", "core", "variable", "private"))]
plot_data[, chr_num   := factor(chr_num,   levels = sprintf("%02d", 1:13))]

chr_counts <- plot_data[, .N, by = .(partition, chr_num, subgenome)]

ggplot(chr_counts, aes(x = chr_num, y = N, fill = subgenome)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.7) +
  facet_wrap(~ partition, ncol = 2, scales = "free_y",
             labeller = labeller(partition = c(
               all      = "All NLR families",
               core     = "Core (5 genomes)",
               variable = "Variable (2–4 genomes)",
               private  = "Private (1 genome)"
             ))) +
  scale_fill_manual(values = c(Ah = "#4472C4", Dh = "#ED7D31"),
                    labels = c(Ah = "A subgenome", Dh = "D subgenome")) +
  labs(x = "Chromosome", y = "NLR pangene families", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))

ggsave("nlr_annotator/nlr_chr_distribution.jpg", width = 9, height = 7, dpi = 300)
cat("  nlr_annotator/nlr_chr_distribution.jpg\n")

# ── 13. Plot: NLR family positions along each chromosome ─────────────────────
# One facet per chromosome (2 rows × 13 cols: Ah on top, Dh on bottom).
# X-axis is interpOrd (gene order along chromosome); each point is one NLR
# pangene family coloured by PAV class.

chr_levels <- c(paste0("Ah_", sprintf("%02d", 1:13)),
                paste0("Dh_", sprintf("%02d", 1:13)))

pg_loc <- pg_stats[, .(interpChr, interpOrd, pav_class)]
pg_loc[, interpChr := factor(interpChr, levels = chr_levels)]
pg_loc[, pav_class := factor(pav_class, levels = c("core", "near_core",
                                                     "variable_3", "variable_2",
                                                     "private"))]

pav_colours <- c(core       = "#2166AC",
                 near_core  = "#74ADD1",
                 variable_3 = "#F4A582",
                 variable_2 = "#D6604D",
                 private    = "#B2182B")
pav_labels  <- c(core       = "Core (5)",
                 near_core  = "Near-core (4)",
                 variable_3 = "Variable (3)",
                 variable_2 = "Variable (2)",
                 private    = "Private (1)")

ggplot(pg_loc, aes(x = interpOrd, y = pav_class, colour = pav_class)) +
  geom_point(size = 1, alpha = 0.6, shape = 16) +
  facet_wrap(~ interpChr, ncol = 13, nrow = 2) +
  scale_colour_manual(values = pav_colours, labels = pav_labels) +
  scale_y_discrete(labels = pav_labels) +
  labs(x = "Gene order along chromosome", y = NULL, colour = "PAV class") +
  theme_bw(base_size = 9) +
  theme(axis.text.x    = element_blank(),
        axis.ticks.x   = element_blank(),
        panel.grid     = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text     = element_text(size = 7))

ggsave("nlr_annotator/nlr_chr_locations.jpg", width = 14, height = 6, dpi = 300)
cat("  nlr_annotator/nlr_chr_locations.jpg\n")

# ── 14. Track plot: Ah_11 and Dh_11 ──────────────────────────────────────────
# Linear circos-style stacked tracks, one per genome. Tick marks are placed at
# actual genomic coordinates (base pairs) from the GTF, not gene-order indices.
# NLR = predicted by either method; present = in an NLR pangene family but no
# NLR call for this genome's gene at this locus.

chr11        <- c("Ah_11", "Dh_11")
genome_order <- rev(genomes)   # rev so TX2094 draws at top (ggplot y = 1 is bottom)

# NLR genes: use actual coordinates already in nlr_result
nlr_chr11 <- nlr_result[chr %in% chr11,
                         .(genome, gene_id, chr, start, status = "NLR")]

# "Present, not NLR" genes: in an NLR pangene family on chr 11 but not in nlr_result.
# Join gene models → gene_pg (limited to NLR families) then anti-join nlr_result.
present_chr11 <- merge(
  genes[chr %in% chr11, .(genome, gene_id, chr, start)],
  gene_pg[pg_row %in% focus_rows, .(genome, gene_id)],
  by = c("genome", "gene_id")
)[!nlr_chr11[, .(genome, gene_id)], on = c("genome", "gene_id")]
present_chr11[, status := "present"]

track_data <- rbind(nlr_chr11, present_chr11)
track_data[, genome  := factor(genome, levels = genome_order)]
track_data[, chr     := factor(chr,    levels = chr11)]
track_data[, status  := factor(status, levels = c("NLR", "present"))]

n_g      <- length(genomes)
bg_bands <- data.table(y    = seq_len(n_g),
                        fill = rep(c("white", "grey93"), length.out = n_g))

ggplot(track_data) +
  geom_rect(data = bg_bands,
            aes(xmin = -Inf, xmax = Inf,
                ymin = y - 0.5, ymax = y + 0.5, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity(guide = "none") +
  geom_segment(aes(x      = start,
                   xend   = start,
                   y      = as.integer(genome) - 0.42,
                   yend   = as.integer(genome) + 0.42,
                   colour = status),
               linewidth = 0.45) +
  facet_wrap(~ chr, ncol = 1, scales = "free_x",
             labeller = labeller(chr = c(
               Ah_11 = "A subgenome  —  Chr 11",
               Dh_11 = "D subgenome  —  Chr 11"
             ))) +
  scale_colour_manual(values = c(NLR = "#B2182B", present = "grey65"),
                      labels = c(NLR = "Predicted NLR", present = "Present, not NLR")) +
  scale_x_continuous(labels = scales::label_comma(scale = 1e-6, suffix = " Mb")) +
  scale_y_continuous(breaks = seq_len(n_g), labels = genome_order,
                     expand = c(0, 0)) +
  labs(x = "Genomic position", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "grey60"),
        legend.position  = "bottom",
        strip.background = element_rect(fill = "grey85"),
        strip.text       = element_text(face = "bold", size = 10))

ggsave("nlr_annotator/nlr_tracks_chr11.jpg", width = 12, height = 5, dpi = 300)
cat("  nlr_annotator/nlr_tracks_chr11.jpg\n")

# ── 15. NLR clustering ────────────────────────────────────────────────────────
# Two NLRs on the same chromosome are in the same cluster when the gap between
# the end of one and the start of the next is ≤ MAX_GAP bp.  Cluster IDs are
# assigned by sorting genes on (genome, chr, start) and incrementing a counter
# whenever a new chromosome starts or the gap exceeds the threshold.
# Includes all categories (both / resistify_only / annotator_gene).

MAX_GAP <- 100000L   # 100 kb

nlr_clust <- nlr_result[!is.na(start) & !is.na(end) & category != "annotator_gene",
                         .(genome, gene_id, chr, start, end, category)
                        ][order(genome, chr, start)]

nlr_clust[, cluster_id := {
  prev_end <- shift(end,  fill = NA_integer_)
  prev_chr <- shift(chr,  fill = "")
  is_new   <- chr != prev_chr | is.na(prev_end) | (start - prev_end) > MAX_GAP
  cumsum(is_new)
}, by = genome]

clust_sz <- nlr_clust[, .(cluster_size = .N), by = .(genome, cluster_id)]
nlr_clust <- merge(nlr_clust, clust_sz, by = c("genome", "cluster_id"))

# Per-genome summary
clust_summary <- nlr_clust[, .(
  n_nlr              = .N,
  n_clusters         = uniqueN(cluster_id),
  n_singleton        = sum(cluster_size == 1  & !duplicated(cluster_id)),
  n_cluster_2_5      = sum(cluster_size %in% 2:5   & !duplicated(cluster_id)),
  n_cluster_6_20     = sum(cluster_size %in% 6:20  & !duplicated(cluster_id)),
  n_supercluster_gt20 = sum(cluster_size > 20 & !duplicated(cluster_id)),
  pct_in_gt5         = round(100 * sum(cluster_size >  5) / .N, 1),
  pct_in_gt20        = round(100 * sum(cluster_size > 20) / .N, 1),
  max_cluster_size   = max(cluster_size)
), by = genome]

cat(sprintf("\n=== NLR clustering (gap threshold = %s bp) ===\n",
            format(MAX_GAP, big.mark = ",")))
print(clust_summary)

# Overall (across all genomes / pangene families are not deduplicated here —
# these are raw per-genome gene counts)
cat(sprintf(
  "\nOverall:  %.1f%% of NLR gene calls occur in clusters of >5 NLRs\n",
  100 * sum(nlr_clust$cluster_size > 5) / nrow(nlr_clust)
))
cat(sprintf(
  "          %.1f%% occur in super-clusters of >20 NLRs\n",
  100 * sum(nlr_clust$cluster_size > 20) / nrow(nlr_clust)
))
cat(sprintf(
  "          %d super-clusters (>20 NLRs) total across all genomes\n",
  sum(clust_summary$n_supercluster_gt20)
))

# Cluster size distribution table (clusters, not genes)
size_cat <- clust_sz[, size_cat := fcase(
  cluster_size == 1,  "singleton",
  cluster_size <= 5,  "small_2_5",
  cluster_size <= 20, "large_6_20",
  default             = "super_gt20"
)]
cat("\n  Cluster counts by size category (per genome):\n")
print(dcast(size_cat[, .N, by = .(genome, size_cat)],
            genome ~ size_cat, value.var = "N", fill = 0L))

# Per-cluster coordinate summary
clust_coords <- nlr_clust[, .(
  chr         = chr[1],
  cluster_size = cluster_size[1],
  start_bp    = min(start),
  end_bp      = max(end),
  span_kb     = round((max(end) - min(start)) / 1e3, 1),
  n_both      = sum(category == "both"),
  n_resistify = sum(category == "resistify_only"),
  n_annotator = sum(category == "annotator_gene")
), by = .(genome, cluster_id)]

clust_coords[, subgenome  := sub("_.*", "", chr)]
clust_coords[, start_mb   := round(start_bp / 1e6, 2)]
clust_coords[, end_mb     := round(end_bp   / 1e6, 2)]

top_n <- 10L
for (sg in c("Ah", "Dh")) {
  top <- clust_coords[subgenome == sg][order(-cluster_size)][seq_len(min(.N, top_n))]
  cat(sprintf("\n  Top %d clusters — %s subgenome:\n", nrow(top),
              if (sg == "Ah") "A" else "D"))
  print(top[, .(genome, chr, start_mb, end_mb, span_kb,
                cluster_size, n_both, n_resistify, n_annotator)])
}

# Write outputs
fwrite(nlr_clust[, .(genome, gene_id, chr, start, end, category,
                     cluster_id, cluster_size)],
       "nlr_annotator/nlr_clusters.tsv",         sep = "\t")
fwrite(clust_summary,
       "nlr_annotator/nlr_cluster_summary.tsv",  sep = "\t")
fwrite(clust_coords[order(subgenome, -cluster_size),
                    .(genome, chr, subgenome, start_mb, end_mb, span_kb,
                      cluster_size, n_both, n_resistify, n_annotator)],
       "nlr_annotator/nlr_cluster_coords.tsv",   sep = "\t")
cat("\n  nlr_annotator/nlr_clusters.tsv\n")
cat("  nlr_annotator/nlr_cluster_summary.tsv\n")
cat("  nlr_annotator/nlr_cluster_coords.tsv\n")

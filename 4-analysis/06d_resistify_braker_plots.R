# 08d_resistify_braker_plots.R
#
# NLR PAV visualization from BRAKER3-based resistify results.
# Input: resistify/nlr_pg_braker.rds  (saved by 08_resistify_pav.R)
#
# Outputs (all to resistify/):
#   panelA_pav_per_genome.jpg       — NLR counts per genome, stacked by PAV class
#   panelB_class_per_genome.jpg     — NLR classification per genome, stacked by PAV class
#   panelC_subgenome_class.jpg      — A vs D subgenome NLR type breakdown per genome
#   panelD_mada_cjid.jpg            — MADA / CJID motif prevalence by classification

library(data.table)
library(ggplot2)

nlr <- readRDS("resistify/nlr_pg_braker.rds")

# Collapse dispensable into variable for 3-level PAV display
nlr[, pav3 := fcase(
  pav_class == "core",                         "core",
  pav_class %in% c("variable", "dispensable"), "variable",
  pav_class == "private",                      "private",
  default = NA_character_
)]
nlr[, pav3 := factor(pav3, levels = c("core", "variable", "private"))]

# Ordered classification factor (most-complex → least)
class_levels <- c("CNL", "CN", "TNL", "TN", "NL", "RNL", "N")
nlr[, Classification := factor(Classification,
                               levels = intersect(class_levels, unique(Classification)))]

# Exclude unplaced genes from all plots
nlr_p <- nlr[!is.na(pav3)]

# ── Colour palettes ───────────────────────────────────────────────────────────
pal_pav  <- c(core = "#1b7837", variable = "#762a83", private = "#d6604d")
pal_class <- c(CNL = "#4E79A7", CN = "#86BCd9",
               TNL = "#E15759", TN = "#F9A8A9",
               NL  = "#F28E2B", RNL = "#59A14F", N = "#B07AA1")

genome_order <- c("TX2094", "Chorlitos", "Consorcio", "Preseco", "Pretropical")
nlr_p[, genome := factor(genome, levels = intersect(genome_order, unique(genome)))]

# ── Panel A: NLR count per genome, stacked by PAV class ──────────────────────
pa_dt <- nlr_p[, .N, by = .(genome, pav3)]

pA <- ggplot(pa_dt, aes(x = genome, y = N, fill = pav3)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = pal_pav, name = NULL) +
  labs(x = NULL, y = "NLR gene count",
       title = "A  Resistance gene presence across genomes") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "right")

ggsave("resistify/panelA_pav_per_genome.jpg", pA,
       width = 6, height = 5, dpi = 150)
message("Saved panelA_pav_per_genome.jpg")

# ── Panel B: NLR classification per genome, stacked by PAV class ─────────────
# Facets = genome; x = Classification; fill = pav3
pb_dt <- nlr_p[, .N, by = .(genome, Classification, pav3)]

pB <- ggplot(pb_dt, aes(x = Classification, y = N, fill = pav3)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = pal_pav, name = NULL) +
  facet_wrap(~ genome, nrow = 1) +
  labs(x = "NLR classification", y = "Count",
       title = "B  NLR classification per genome, shaded by PAV class") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"))

ggsave("resistify/panelB_class_per_genome.jpg", pB,
       width = 13, height = 5, dpi = 150)
message("Saved panelB_class_per_genome.jpg")

# ── Panel C: Subgenome (A vs D) NLR type breakdown per genome ────────────────
# Requires subgenome assignment; excludes any NA subgenome
pc_dt <- nlr_p[!is.na(subgenome), .N, by = .(genome, subgenome, Classification)]
pc_dt[, subgenome := factor(subgenome, levels = c("A", "D"))]

pC <- ggplot(pc_dt, aes(x = subgenome, y = N, fill = Classification)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = pal_class, name = "NLR class") +
  facet_wrap(~ genome, nrow = 1) +
  labs(x = "Subgenome", y = "NLR gene count",
       title = "C  NLR classification by subgenome") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "grey90"))

ggsave("resistify/panelC_subgenome_class.jpg", pC,
       width = 11, height = 5, dpi = 150)
message("Saved panelC_subgenome_class.jpg")

# ── Panel D: MADA / CJID motif prevalence by classification ──────────────────
# MADA and CJID mark functionally distinct CC-NBS subclades; stored as
# Python-style True/False strings — convert before counting.
if (is.character(nlr_p$MADA)) {
  nlr_p[, MADA := MADA == "True"]
  nlr_p[, CJID := CJID == "True"]
}

motif_dt <- nlr_p[, .(
  MADA_pos = sum(MADA,  na.rm = TRUE),
  CJID_pos = sum(CJID,  na.rm = TRUE),
  total    = .N
), by = Classification]
motif_dt[, MADA_pct := MADA_pos / total * 100]
motif_dt[, CJID_pct := CJID_pos / total * 100]

motif_long <- melt(motif_dt,
  id.vars      = "Classification",
  measure.vars = c("MADA_pct", "CJID_pct"),
  variable.name = "motif",
  value.name    = "pct"
)
motif_long[, motif := sub("_pct$", "", motif)]

pD <- ggplot(motif_long, aes(x = Classification, y = pct, fill = motif)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.65),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(MADA = "#4E79A7", CJID = "#F28E2B"),
                    name = "Motif") +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = "NLR classification", y = "% of genes with motif",
       title = "D  MADA and CJID motif prevalence by NLR class") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

ggsave("resistify/panelD_mada_cjid.jpg", pD,
       width = 7, height = 5, dpi = 150)
message("Saved panelD_mada_cjid.jpg")

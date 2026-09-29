# 03b_pankmer_heatmap_atlas.R
#
# Whole-genome k-mer similarity heatmap for the Atlas genome set.
# Produces two figures from the same adjacency matrix:
#
#   MAIN (60 genomes): Jing_ and Yang107_ excluded
#     pankmer_heatmap_main.pdf / .jpg
#
#   SUPPLEMENTAL (~167 genomes): Jing_ excluded only
#     pankmer_heatmap_supp.pdf / .jpg
#
# Run from the fromTony working directory:
#   Rscript 03b_pankmer_heatmap_atlas.R

library(pheatmap)
library(RColorBrewer)

MATRIX_FILE <- "manyGenomes/fromAtlas/pankmer_adj_matrix.csv"
OUT_DIR     <- "manyGenomes/fromAtlas"

# ── load full matrix ──────────────────────────────────────────────────────────

raw_all <- as.matrix(read.csv(MATRIX_FILE, row.names = 1, check.names = FALSE))
storage.mode(raw_all) <- "double"
cat(sprintf("Full matrix: %d genomes\n", nrow(raw_all)))

# ── helpers ───────────────────────────────────────────────────────────────────

jacc_from_raw <- function(raw) {
  total <- diag(raw)
  denom <- outer(total, total, "+") - raw
  j     <- raw / denom
  j     <- (j + t(j)) / 2
  diag(j) <- 1
  j
}

get_source <- function(nm) {
  core <- c("Chorlitos", "Consorcio", "Preseco", "Pretropical")
  ifelse(nm %in% core,              "present study",
  ifelse(grepl("^HAU35_",  nm),     "HAU35",
  ifelse(grepl("^Yang107_", nm),    "Yang107",
  ifelse(grepl("^HAU_",    nm),     "HAU",
  ifelse(grepl("^CRI_",    nm),     "CRI",
  ifelse(grepl("^HGS_",    nm),     "HGS",
  ifelse(grepl("^HEAU_",   nm),     "HEAU",
  ifelse(grepl("^NSF_",    nm),     "NSF",
  ifelse(grepl("^TM1_",    nm),     "TM1",
  ifelse(grepl("^USDA_",   nm),     "USDA",
  ifelse(grepl("^XAAS_",   nm),     "XAAS",
  ifelse(grepl("^ZJU_",    nm),     "ZJU",
         "other"))))))))))))
}

make_annot <- function(nms) {
  sources <- get_source(nms)
  lvls    <- c("present study", sort(setdiff(unique(sources), "present study")))
  annot   <- data.frame(Source = factor(sources, levels = lvls), row.names = nms)
  levels  <- levels(annot$Source)
  palette <- setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(length(levels)),
    levels
  )
  palette["present study"] <- "black"
  list(df = annot, colors = list(Source = palette))
}

# ── wild / semi-wild / domesticated status ────────────────────────────────────
# HAU35 and Yang107 status come from their sample tables; the remaining
# genomes in the main figure are assigned here. Genomes with no assignment
# are shown as "Unclassified".

HAU35_FILE <- "manyGenomes/fromAtlas/HAU35_samples.tsv"
YANG_FILE  <- "manyGenomes/fromAtlas/Yang107_samples.tsv"

hau35 <- read.delim(HAU35_FILE, stringsAsFactors = FALSE)
yang  <- read.delim(YANG_FILE,  stringsAsFactors = FALSE)

status_of <- function(status, ...) {
  nms <- c(...)
  setNames(rep(status, length(nms)), nms)
}

STATUS <- c(
  setNames(hau35$Status,   hau35$MatrixName),
  setNames(yang$Category,  yang$MatrixName),
  status_of("Wild", "Chorlitos", "Consorcio", "Preseco", "Pretropical", "USDA_TX2094"),
  status_of("Semi-wild", "CRI_TX1000", "HAU_TX31"),
  status_of("Cultivar", "CRI_ZM113", "CRI_ZM24", "HAU_Jin668", "HAU_YZ1", "HAU_EM22",
            "HEAU_NDM8", "HGS_CSX8308", "HGS_UA48", "HGS_UGA230", "USDA_PSC355",
            "XAAS_YM11", "NSF_B713", "NSF_Bar32", "ZJU_HPF17",
            "TM1_CRI", "TM1_HAU", "TM1_UTXv31", "TM1_ZJU_T2T")
)
stopifnot(all(STATUS %in% c("Wild", "Semi-wild", "Cultivar")))

make_annot_status <- function(nms) {
  st <- unname(STATUS[nms])
  st[is.na(st)] <- "Unclassified"
  if (any(st == "Unclassified"))
    cat("  Unclassified:", paste(nms[st == "Unclassified"], collapse = ", "), "\n")
  st <- droplevels(factor(st,
    levels = c("Wild", "Semi-wild", "Cultivar", "Unclassified"),
    labels = c("Wild", "Semi-wild", "Domesticated (cultivar)", "Unclassified")))
  pal <- c("Wild"                    = "#1B9E77",
           "Semi-wild"               = "#D95F02",
           "Domesticated (cultivar)" = "#7570B3",
           "Unclassified"            = "#BDBDBD")
  list(df     = data.frame(Status = st, row.names = nms),
       colors = list(Status = pal[levels(st)]))
}

# Highlight only the four new genomes, in forest green; every other genome
# gets no annotation colour (white, blends into the background).
make_annot_highlight <- function(nms, core = c("Chorlitos", "Consorcio", "Preseco", "Pretropical")) {
  mark  <- factor(ifelse(nms %in% core, "New genome", ""), levels = c("New genome", ""))
  annot <- data.frame(" " = mark, row.names = nms, check.names = FALSE)
  list(df = annot, colors = list(" " = setNames(c("#0B5D1E", "white"), c("New genome", ""))))
}

draw_heatmap <- function(jacc, stem, fs_general, fs_labels, w, h, color_by = "source",
                         legends = TRUE) {
  annot  <- switch(color_by,
                    status    = make_annot_status(rownames(jacc)),
                    highlight = make_annot_highlight(rownames(jacc)),
                    make_annot(rownames(jacc)))
  lo     <- floor(min(jacc[lower.tri(jacc)]) * 100) / 100
  breaks <- seq(lo, 1.0, length.out = 101)
  colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
  dist_j <- as.dist(1 - jacc)

  args <- list(
    mat                      = jacc,
    color                    = colors,
    breaks                   = breaks,
    annotation_row           = annot$df,
    annotation_col           = annot$df,
    annotation_colors        = annot$colors,
    legend                   = legends,
    annotation_legend        = legends,
    clustering_distance_rows = dist_j,
    clustering_distance_cols = dist_j,
    clustering_method        = "average",
    display_numbers          = FALSE,
    fontsize                 = fs_general,
    fontsize_row             = fs_labels,
    fontsize_col             = fs_labels
  )

  cat(sprintf("  Writing %s.pdf (%d × %d in)...\n", stem, w, h))
  do.call(pheatmap, c(args, list(filename = file.path(OUT_DIR, paste0(stem, ".pdf")),
                                 width = w, height = h)))

  cat(sprintf("  Writing %s.jpg...\n", stem))
  do.call(pheatmap, c(args, list(filename = file.path(OUT_DIR, paste0(stem, ".jpg")),
                                 width = w, height = h)))
}

# ── MAIN figure: exclude Jing_ and Yang107_ (~60 genomes) ────────────────────

keep_main <- !grepl("^Jing_|^Yang107_", rownames(raw_all))
raw_main  <- raw_all[keep_main, keep_main]
cat(sprintf("\nMain figure: %d genomes (Jing_ + Yang107_ excluded)\n", nrow(raw_main)))

draw_heatmap(
  jacc       = jacc_from_raw(raw_main),
  stem       = "pankmer_heatmap_main",
  fs_general = 24,
  fs_labels  = 20,
  w          = 38,
  h          = 32,
  color_by   = "highlight",
  legends    = FALSE
)

# ── SUPPLEMENTAL figure: exclude Jing_ only (~167 genomes) ───────────────────

keep_supp <- !grepl("^Jing_", rownames(raw_all))
raw_supp  <- raw_all[keep_supp, keep_supp]
cat(sprintf("\nSupplemental figure: %d genomes (Jing_ excluded)\n", nrow(raw_supp)))

draw_heatmap(
  jacc       = jacc_from_raw(raw_supp),
  stem       = "pankmer_heatmap_supp",
  fs_general = 24,
  fs_labels  = 18,
  w          = 74,
  h          = 68
)

message("\nDone.")


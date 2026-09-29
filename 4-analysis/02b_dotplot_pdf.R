# 02b_dotplot_pdf.R
#
# Render pairwise dotplots from PAF files into:
#   1. dotplots/all_comparisons_tiled.jpg  — 6-row × 5-col overview image
#   2. dotplots/all_comparisons.pdf        — one full page per comparison
#
# Run after all minimap2 array tasks are complete:
#   Rscript 02b_dotplot_pdf.R

library(data.table)
library(ggplot2)
library(patchwork)

PAF_DIR    <- "dotplots"
OUT_TILED  <- "dotplots/all_comparisons_tiled.jpg"
OUT_PDF    <- "dotplots/all_comparisons.pdf"
MIN_SEQLEN <- 100000L   # minimum sequence length — mirrors dotplotly -q

# Grid definition: rows = reference genomes, columns = query genomes
REF_NAMES <- c("TX2094", "Latifolium", "Jin668", "TM1", "YZ1", "ZM113")
YUC_NAMES <- c("Chorlitos", "Consorcio", "Preseco", "Pretropical", "TX2094")

STRAND_COLS <- c("+" = "#2166AC", "-" = "#D73027")

# Display labels for the figure (file names keep the assembly names)
DISPLAY <- c(Chorlitos   = "Río Lagartos (RiCh)",
             Consorcio   = "Centro (CeCo)",
             Preseco     = "Sisal (SiPr)",
             Pretropical = "Celestún (CePr)")
disp <- function(x) ifelse(x %in% names(DISPLAY), DISPLAY[x], x)

# ── helpers ───────────────────────────────────────────────────────────────────

read_paf <- function(f) {
  paf <- fread(f, sep = "\t", header = FALSE, quote = "",
               select = 1:12, showProgress = FALSE, fill = TRUE,
               col.names = c("qname", "qlen", "qstart", "qend", "strand",
                             "tname", "tlen", "tstart", "tend",
                             "nmatch", "alen", "mapq"))
  paf[qlen >= MIN_SEQLEN & tlen >= MIN_SEQLEN]
}

# Cumulative axis offsets so sequences sit end-to-end
seq_offsets <- function(names, lens) {
  dt <- unique(data.table(name = names, len = lens))[order(name)]
  dt[, offset := c(0L, cumsum(as.numeric(len))[-.N])]
  dt[, mid    := offset + len / 2]
  dt[, label  := sub("^([AD])h_?0*", "\\1", name)]   # Ah_01 → A1
  dt
}

# Full-page dotplot (for PDF)
make_dotplot_full <- function(paf, title) {
  q_off <- seq_offsets(paf$qname, paf$qlen)
  t_off <- seq_offsets(paf$tname, paf$tlen)
  paf   <- merge(paf, q_off[, .(qname = name, q_off = offset)], by = "qname")
  paf   <- merge(paf, t_off[, .(tname = name, t_off = offset)], by = "tname")
  paf[, `:=`(x1 = qstart + q_off, x2 = qend + q_off,
              y1 = tstart + t_off, y2 = tend  + t_off)]
  paf[strand == "-", `:=`(y1 = tend + t_off, y2 = tstart + t_off)]

  ggplot(paf) +
    geom_vline(xintercept = q_off$offset[-1], colour = "grey78", linewidth = 0.2) +
    geom_hline(yintercept = t_off$offset[-1], colour = "grey78", linewidth = 0.2) +
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, colour = strand),
                 linewidth = 0.3, alpha = 0.85) +
    scale_colour_manual(values = STRAND_COLS,
                        labels = c("+" = "Forward", "-" = "Reverse")) +
    scale_x_continuous(breaks = q_off$mid, labels = q_off$label,
                       limits = c(0, sum(as.numeric(q_off$len))), expand = c(0, 0)) +
    scale_y_continuous(breaks = t_off$mid, labels = t_off$label,
                       limits = c(0, sum(as.numeric(t_off$len))), expand = c(0, 0)) +
    labs(title = title, x = "Query", y = "Reference", colour = "Strand") +
    theme_bw(base_size = 12) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y     = element_text(size = 9),
          plot.title      = element_text(face = "bold", size = 13),
          panel.grid      = element_blank(),
          legend.position = "bottom")
}

# Small tile dotplot (for tiled image) — no axis text, minimal chrome
make_dotplot_tile <- function(paf) {
  q_off <- seq_offsets(paf$qname, paf$qlen)
  t_off <- seq_offsets(paf$tname, paf$tlen)
  paf   <- merge(paf, q_off[, .(qname = name, q_off = offset)], by = "qname")
  paf   <- merge(paf, t_off[, .(tname = name, t_off = offset)], by = "tname")
  paf[, `:=`(x1 = qstart + q_off, x2 = qend + q_off,
              y1 = tstart + t_off, y2 = tend  + t_off)]
  paf[strand == "-", `:=`(y1 = tend + t_off, y2 = tstart + t_off)]

  ggplot(paf) +
    geom_vline(xintercept = q_off$offset[-1], colour = "grey82", linewidth = 0.1) +
    geom_hline(yintercept = t_off$offset[-1], colour = "grey82", linewidth = 0.1) +
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, colour = strand),
                 linewidth = 0.2, alpha = 0.85) +
    scale_colour_manual(values = STRAND_COLS, guide = "none") +
    scale_x_continuous(limits = c(0, sum(as.numeric(q_off$len))), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, sum(as.numeric(t_off$len))), expand = c(0, 0)) +
    theme_bw(base_size = 7) +
    theme(axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          axis.title  = element_blank(),
          panel.grid  = element_blank(),
          plot.margin = margin(1, 1, 1, 1))
}

blank_tile <- function(label = "") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             size = 3.5, colour = "grey60", fontface = "italic") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw(base_size = 7) +
    theme(axis.text       = element_blank(),
          axis.ticks      = element_blank(),
          axis.title      = element_blank(),
          panel.grid      = element_blank(),
          panel.background = element_rect(fill = "grey96", colour = NA),
          plot.margin     = margin(1, 1, 1, 1))
}

# ── build tiled image ─────────────────────────────────────────────────────────
cat("Building tiled image...\n")

tile_list <- list()

for (ref in REF_NAMES) {
  for (yuc in YUC_NAMES) {
    key      <- sprintf("%s_vs_%s", yuc, ref)
    paf_file <- file.path(PAF_DIR, paste0(key, ".paf"))
    is_self  <- yuc == ref

    if (is_self) {
      p <- blank_tile("(self)")
    } else if (!file.exists(paf_file)) {
      p <- blank_tile("no PAF")
    } else {
      paf <- tryCatch(read_paf(paf_file), error = function(e) NULL)
      if (is.null(paf) || nrow(paf) == 0) {
        p <- blank_tile("no aln")
      } else {
        p <- make_dotplot_tile(paf)
      }
    }
    tile_list[[length(tile_list) + 1]] <- p
    cat(sprintf("  [%s vs %s] %s\n", yuc, ref,
                if (is_self) "self" else sprintf("%d alns", if (!is.null(paf)) nrow(paf) else 0)))
  }
}

# Column headers (query names) across top; row headers (ref names) down left side.
# Achieved by labelling the top row plots and prepending a narrow label column.

# Add column headers to the top row (first 5 tiles = first ref row)
for (i in seq_along(YUC_NAMES)) {
  tile_list[[i]] <- tile_list[[i]] +
    labs(title = disp(YUC_NAMES[i])) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                    margin = margin(b = 2)))
}

# Build row-label plots (narrow, one per reference)
row_labels <- lapply(REF_NAMES, function(ref) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = disp(ref),
             angle = 90, size = 5.5, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
})

# Interleave: for each row, prepend its label then the 5 tile plots
ordered_plots <- list()
for (r in seq_along(REF_NAMES)) {
  ordered_plots <- c(ordered_plots,
                     list(row_labels[[r]]),
                     tile_list[((r - 1) * 5 + 1):(r * 5)])
}

tiled <- wrap_plots(ordered_plots, ncol = 6, nrow = 6,
                    widths  = c(0.14, 1, 1, 1, 1, 1)) +
  plot_annotation(
    title   = "Pairwise whole-genome synteny",
    caption = sprintf("Sequences ≥ %d kb shown. Blue = forward, red = reverse strand.",
                      MIN_SEQLEN / 1000L),
    theme   = theme(plot.title   = element_text(size = 18, face = "bold"),
                    plot.caption = element_text(size = 11))
  )

ggsave(OUT_TILED, tiled, width = 22, height = 18, dpi = 200)
cat(sprintf("Tiled image → %s\n\n", OUT_TILED))

# ── multi-page PDF (full detail per comparison) ───────────────────────────────
cat("Building per-comparison PDF...\n")
pdf(OUT_PDF, width = 12, height = 10)
n_ok <- 0L

for (ref in REF_NAMES) {
  for (yuc in YUC_NAMES) {
    if (yuc == ref) next
    paf_file <- file.path(PAF_DIR, sprintf("%s_vs_%s.paf", yuc, ref))
    if (!file.exists(paf_file)) next
    paf <- tryCatch(read_paf(paf_file), error = function(e) NULL)
    if (is.null(paf) || nrow(paf) == 0) next
    print(make_dotplot_full(paf, sprintf("%s  vs  %s", disp(yuc), disp(ref))))
    n_ok <- n_ok + 1L
  }
}

dev.off()
cat(sprintf("PDF (%d pages) → %s\n", n_ok, OUT_PDF))

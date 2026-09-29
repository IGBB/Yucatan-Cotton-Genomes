# Drop host R library injected via ~/.Renviron — use only container libraries
.libPaths(.libPaths()[.libPaths() != "/lustre/hdd/LAS/mhufford-lab/corrinne/app/R/library"])

library(GENESPACE)
library(data.table)

load("genespace/results/gsParams.rda")

# Get rid of Tony's paths and substitute my own
old_wd <- gsParam$paths$wd
new_wd <- normalizePath("genespace")

gsParam$paths <- lapply(gsParam$paths, function(x) gsub(old_wd, new_wd, x, fixed = TRUE))

# ── Pangene database ───────────────────────────────────────────────────────────
pg <- query_pangenes(
  gsParam = gsParam,
  refGenome = "TX2094.At",
  transform = TRUE,
  showArrayMem = FALSE,
  showNSOrtho = FALSE,
  maxMem2Show = Inf,
  showUnPlacedPgs = FALSE
)

fwrite(pg, file = "genespace/pangenes/pangene_db.tsv", sep = "\t")

# ── Riparian plots ─────────────────────────────────────────────────────────────
pdf("genespace/riparian_TX2094At.pdf", width = 10, height = 8)
plot_riparian(
  gsParam   = gsParam,
  refGenome = "TX2094.At",
  useRegions = FALSE
)
dev.off()

pdf("genespace/riparian_TX2094At_ordered.pdf", width = 13, height = 6)
plot_riparian(
  gsParam           = gsParam,
  refGenome         = "TX2094.At",
  genomeIDs         = c("TX2094.At", "TX2094.Dt",
                        "Chorlitos.At", "Chorlitos.Dt",
                        "Consorcio.At", "Consorcio.Dt",
                        "Preseco.At",   "Preseco.Dt",
                        "Pretropical.At", "Pretropical.Dt"),
  reorderBySynteny  = FALSE,
  forceRecalcBlocks = FALSE,
  useOrder          = FALSE,
  minChrLen2plot    = 10,
  chrLabFontSize    = 6
)
dev.off()

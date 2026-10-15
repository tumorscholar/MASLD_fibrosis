###############################################################################
# 41 - External validation: T-cell CLONAL EXPANSION increases with MASLD stage
#
#      Independent cohort: Nat. Genet. 2025, s41588-025-02407-8
#      ("Spatially resolved multi-omics of human MASLD"; 61 donors, liver scRNA +
#       TCR). Source Data Fig. 2f (Supplementary Table / MOESM6), sheet "Fig. 2f":
#       per-cell  sample | disease_status | cell_type_lvl2 | clone_id | clone_id_size
#       where clone_id = TCRa/b CDR3-similarity clonotype CLUSTER and
#       clone_id_size = no. cells in that cluster (their Fig. 2f selects size >= 20).
#
# This is a SEPARATE validation plot (our data is NOT merged in). It shows, in an
# external cohort, the same direction we report: healthy livers carry very few /
# no (hyper)expanded clones vs MASLD.
#
# Stage mapping to our axis:  Healthy = CTRL, No_fibrosis = MASL, Fibrosis = MASH.
#
# Two panels only:
#   (a) % of T cells in expanded clones (clonal size >= 20)   <- their exact cutoff
#   (b) largest clonal size per sample (cells)                <- (hyper)expansion
# dot = sample (size ~ no. cells) · bar = median · Wilcoxon (Mann-Whitney) pairwise.
#
# NOTE on seq depth: the >=20 cutoff is absolute and CTRL libraries are smaller,
# so it is biased against detecting expansion in controls; panel (b) (max clone
# size) is reported alongside precisely because it still separates CTRL cleanly.
###############################################################################

suppressMessages({
  library(readxl); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork)
})
source("00_nature_style.R")   # theme_nature(), save_nature(), stage_cols, npg_pal

## ---- paths -----------------------------------------------------------------
## Point XLSX at the downloaded Source Data file for Fig. 2 of the paper.
XLSX    <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData/41588_2025_2407_MOESM6_ESM.xlsx"
SHEET   <- "Fig. 2f"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
OUTNAME <- "ExtDataFig_ExternalValidation_TCR_clonalExpansion_NatGenet2025"

## ---- load & per-sample clonal-expansion metrics ----------------------------
d <- as.data.table(read_excel(XLSX, sheet = SHEET))
setnames(d, 1, "barcode")                       # first (unnamed) col = cell barcode
d[, clone_id := as.character(clone_id)]

per <- d[, .(stage     = disease_status[1],
             n_cells   = .N,
             pct_exp20 = 100 * sum(clone_id_size >= 20) / .N,  # % cells in clusters (global) >=20
             max_clone = max(clone_id_size)),                  # largest clone touching this sample
         by = sample]

## stage factor: CTRL -> MASL -> MASH, labelled with our equivalent axis
lev <- c("CTRL", "MASL", "MASH")
lab <- c("CTRL\n(Healthy)", "MASL\n(No_fibrosis)", "MASH\n(Fibrosis)")
per[, stage := factor(stage, levels = lev, labels = lab)]
setorder(per, stage, sample)
fwrite(per, file.path(tab_dir, paste0(OUTNAME, "_perSample.csv")))

## ---- colours (our stage palette, re-keyed to the labelled stages) ----------
stage_fill <- setNames(unname(stage_cols[c("Healthy","No_fibrosis","Fibrosis")]), lab)
cmp <- list(lab[1:2], c(lab[1], lab[3]), lab[2:3])   # CTRL-MASL, CTRL-MASH, MASL-MASH

## ---- panel builder ---------------------------------------------------------
base_panel <- function(dt, yvar, ylab, logy = FALSE) {
  p <- ggplot(dt, aes(stage, .data[[yvar]])) +
    geom_jitter(aes(size = n_cells, fill = stage),
                position = position_jitter(width = 0.14, height = 0, seed = 1),
                shape = 21, colour = "black", stroke = 0.25, alpha = 0.85) +
    stat_summary(fun = median, fun.min = median, fun.max = median,   # thin median bar ON TOP of dots
                 geom = "crossbar", width = 0.42, linewidth = 0.25, fatten = 1, colour = "black") +
    scale_fill_manual(values = stage_fill, guide = "none") +
    scale_size_area("No. cells", max_size = 4,
                    breaks = c(100, 500, 1000, 2000)) +
    stat_compare_means(comparisons = cmp, method = "wilcox.test",
                       label = "p.signif", size = 2.4, tip.length = 0.01,
                       bracket.size = 0.3) +
    labs(x = NULL, y = ylab) +
    theme_nature(7) +
    theme(legend.position = "right",
          legend.key.size = unit(3, "mm"),
          axis.text.x = element_text(size = 6))
  if (logy) p <- p + scale_y_log10()
  p
}

pa <- base_panel(per, "pct_exp20", "% of T cells in expanded\nclones (clonal size ≥ 20)") +
  theme(legend.position = "none")
pb <- base_panel(per, "max_clone", "Largest clonal size\nper sample (cells)", logy = TRUE)

p <- (pa | pb) +
  plot_annotation(
    title = "T-cell clonal expansion increases with MASLD stage (external cohort, Li et al, Nat. Genet. 2025)",
    subtitle = sprintf("Fig. 2f source data · dot = sample (size ∝ no. cells) · bar = median · Wilcoxon · CTRL n=%d, MASL n=%d, MASH n=%d",
                       sum(grepl("CTRL", per$stage)), sum(grepl("MASL", per$stage)),
                       sum(grepl("MASH", per$stage))),
    theme = theme(plot.title    = element_text(size = 7, hjust = 0.5),
                  plot.subtitle = element_text(size = 6, hjust = 0.5)))

save_nature(p, file.path(fig_dir, OUTNAME), w_cm = 14, h_cm = 7)
message("Wrote ", file.path(fig_dir, OUTNAME), ".pdf/.tiff and _perSample.csv in ", tab_dir)

## ---- console: medians + pairwise Wilcoxon ----------------------------------
cat("\nMedian by stage:\n")
print(per[, .(n_cells   = as.numeric(median(n_cells)),
              pct_exp20 = round(median(pct_exp20), 2),
              max_clone = as.numeric(median(max_clone))), by = stage])

for (v in c("pct_exp20", "max_clone")) {
  cat("\n", v, " (Wilcoxon):\n", sep = "")
  for (pair in cmp) {
    pv <- wilcox.test(per[[v]][per$stage == pair[1]],
                      per[[v]][per$stage == pair[2]])$p.value
    cat(sprintf("  %-18s vs %-18s  p = %.4g\n",
                gsub("\n", " ", pair[1]), gsub("\n", " ", pair[2]), pv))
  }
}

# End of the script

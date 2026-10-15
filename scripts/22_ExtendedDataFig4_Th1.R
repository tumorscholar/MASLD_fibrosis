###############################################################################
# 22 - Extended Data Figure 4 - known MASLD-associated T-cell subsets (Th1 / auto-aggressive)
#
#  a : PBMC Th1 cells (CXCR3+ TBX21+), % of PBMC cells, per patient, by stage
#  b : liver auto-aggressive CD8 (CD8A+ CD69+ CXCR6+), % of liver, by stage
#
# V7: percentages are COMPUTED from SeuObjV7 (the old script hard-coded the Prism
# values); positive = RNA log-norm > 0.25. Bar = mean, error = SEM, points = per
# patient.
###############################################################################
library(Seurat); library(tidyverse)
source("00_nature_style.R")

seu_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE); dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
stage_lev <- c("Healthy", "No_fibrosis", "Fibrosis")

SeuObj <- readRDS(seu_rds); DefaultAssay(SeuObj) <- "RNA"

# per-patient % positive within a tissue, then bar+SEM+jitter by stage
pct_per_patient <- function(obj, tissue, genes, thr = 0.25) {
  o <- subset(obj, subset = Tissue == tissue)
  d <- FetchData(o, vars = c(genes, "Patient_ID", "Stage"), layer = "data")
  d$pos <- Reduce(`&`, lapply(genes, function(g) d[[g]] > thr))
  d %>% group_by(Patient_ID, Stage) %>%
    summarise(pos = sum(pos), total = n(), pct = 100 * pos / total, .groups = "drop") %>%
    mutate(Stage = factor(Stage, levels = stage_lev))
}
bar_plot <- function(df, ylab, ttl) ggplot(df, aes(Stage, pct, fill = Stage)) +
  stat_summary(fun = mean, geom = "col", width = 0.6, colour = "black", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, colour = "black", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 1.1, colour = "black", show.legend = FALSE) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  labs(x = NULL, y = ylab, title = ttl) +
  theme_nature(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- ED Fig 3a: Th1 in PBMC ------------------------------------------------
a <- pct_per_patient(SeuObj, "PBMC", c("CXCR3", "TBX21"))
write.csv(a, file.path(tab_dir, "EDFig3a_Th1_PBMC_perPatient_v7.csv"), row.names = FALSE)
save_nature(bar_plot(a, "CXCR3+ TBX21+", "% Th1 cells in PBMC"),
            file.path(fig_dir, "ExtendedDataFig3a_Th1_PBMC_v7"), 4.64, 6.95)

# ---- ED Fig 3b: auto-aggressive CD8 in liver -------------------------------
b <- pct_per_patient(SeuObj, "LIVER", c("CD8A", "CD69", "CXCR6"))
write.csv(b, file.path(tab_dir, "EDFig3b_AutoAggr_liver_perPatient_v7.csv"), row.names = FALSE)
save_nature(bar_plot(b, "CD8+ CD69+ CXCR6+", "% Auto-aggressive CD8 T cells\nin liver"),
            file.path(fig_dir, "ExtendedDataFig3b_AutoAggr_liver_v7"), 5.52, 6.97)
cat("\nDone Extended Data Fig 3.\n")
# End of the script

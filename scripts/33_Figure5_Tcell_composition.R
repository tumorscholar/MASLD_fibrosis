###############################################################################
# 33 - Figure 5 - T-cell composition shift, baseline vs follow-up (liver + PBMC)
#
# The raw FNA composition is confounded: follow-up FNA barely captures stroma /
# parenchyma, so every T-cell % is inflated by sampling. This panel therefore
# normalises WITHIN T cells (each T-cell subtype as a % of all T cells captured
# in that visit), which removes the stromal-dropout confound and shows the true
# shift in T-cell composition. Recomputed directly from the FNA objects.
###############################################################################
suppressMessages({library(tidyverse)})
source("00_nature_style.R")

# read the SAME saved composition tables that 30_Figure5_liverFNA.R / 31_Figure5_pbmcFNA.R write
# (one source of truth - the tables behind the raw Fig 5 composition bar plot).
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
comp_files <- c(LIVER = file.path(fig_dir, "Fig5_liver_composition_v7.csv"),
                PBMC  = file.path(fig_dir, "Fig5_pbmc_composition_v7.csv"))
stopifnot(all(file.exists(comp_files)))   # run 18_/19_ first if this fails

t_types <- c("Naive/TCM T cells","Memory T cells","Activated T cells",
             "Activated CD8 T cells","Cytotoxic CD8 T cells")

comp <- imap_dfr(comp_files, function(f, tis)
  readr::read_csv(f, show_col_types = FALSE) %>% mutate(Tissue = tis)) %>%
  filter(cell_type %in% t_types) %>%
  group_by(Tissue, Visit) %>% mutate(pct_T = 100 * n / sum(n)) %>% ungroup() %>%  # normalise WITHIN T cells
  mutate(Visit = factor(Visit, levels = c("Baseline","Followup")),
         cell_type = factor(cell_type, levels = t_types),
         Tissue = factor(Tissue, levels = c("LIVER","PBMC")))

write.csv(comp, file.path(tab_dir, "Fig5_Tcell_composition_normalised_v7.csv"), row.names = FALSE)
cat("T-cell-normalised composition (% of captured T cells):\n")
print(comp %>% arrange(Tissue, Visit, cell_type), n = 40)

# NOTE: the plotted panel (within-T normalised composition) is SUPERSEDED by
# 33b_Figure5_immune_composition.R (all immune cell types, cluster-labelled, which
# also resolves the five T-cell clusters). This script now only writes its data
# table; the figure is not plotted, so re-running 30-40 produces no stale panel.
# To restore the within-T figure, uncomment the block below.
#
# t_cols <- c("Naive/TCM T cells"="#9ecae1","Memory T cells"="#6baed6",
#             "Activated T cells"="#fdae6b","Activated CD8 T cells"="#e6550d",
#             "Cytotoxic CD8 T cells"="#a50f15")
# p <- ggplot(comp, aes(Visit, pct_T, fill = cell_type)) +
#   geom_col(width = 0.7) + facet_wrap(~ Tissue) +
#   scale_fill_manual(values = t_cols, name = NULL) +
#   labs(x = NULL, y = "% of captured T cells",
#        title = "T-cell composition, baseline vs follow-up (normalised within T cells)") +
#   theme_nature(7) +
#   theme(legend.position = "right", legend.key.size = unit(0.35,"cm"),
#         strip.background = element_rect(fill = "grey92", colour = NA))
# save_nature(p, file.path(fig_dir, "Fig5_Tcell_composition_normalised_v7"), 11, 7)
cat("\nDone Fig 5 T-cell composition table (figure superseded by 33b).\n")
# End of the script

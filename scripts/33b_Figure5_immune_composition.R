###############################################################################
# 33b - Figure 5c - IMMUNE-cell composition, baseline vs follow-up
#        (ALL immune cell types, each as a % of all immune / CD45+ cells)
#        LIVER  -> main Fig 5c        PBMC -> supplementary
#
# The raw FNA composition is confounded by stromal/parenchymal drop-out at
# follow-up. This panel instead asks, WITHIN the captured immune compartment,
# how the make-up of the leukocytes shifts after surgery: every immune cell type
# is a % of all CD45+ (immune) cells captured in that visit. Non-immune lineages
# (endothelium, fibroblasts, hepatocytes, cholangiocytes, adipocytes, mesothelium,
# erythrocytes) are excluded from the denominator.
#
# STACKED, cluster-wise: each bar is stacked by cell type and every segment is
# labelled with its Fig-1 cluster, e.g. "Activated T cells (C20)", so the five
# T-cell clusters (C0/C1/C2/C6/C20) are individually distinguishable.
#
# Reads the composition tables written by 30_/31_ (one source of truth):
#   Fig5_liver_composition_v7.csv , Fig5_pbmc_composition_v7.csv
# Figures (only the two we use): Fig5c_immune_composition_LIVER_v7 (main),
#                                Fig5c_immune_composition_PBMC_v7  (supp)
# Tables : Fig5c_immune_composition_fine_v7.csv (Tissue,Visit,cell_type,cluster,n,pct_immune)
#          Fig5c_immune_composition_lineage_v7.csv (+ baseline->follow-up delta)
###############################################################################
suppressMessages({ library(tidyverse) })
source("00_nature_style.R")

fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
ann_csv <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1/cluster_annotation_V7.csv"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

comp_files <- c(LIVER = file.path(fig_dir, "Fig5_liver_composition_v7.csv"),
                PBMC  = file.path(fig_dir, "Fig5_pbmc_composition_v7.csv"))
stopifnot(all(file.exists(comp_files)))   # run 30_/31_ first if this fails

# ---- cell_type -> Fig-1 cluster label (e.g. "Activated T cells (C20)") ------
ann <- read.csv(ann_csv, stringsAsFactors = FALSE)          # columns: cluster (C#), cell_type
clabel <- ann %>%
  mutate(cid_num = suppressWarnings(as.integer(sub("^C", "", cluster)))) %>%
  arrange(cid_num) %>%
  group_by(cell_type) %>%
  summarise(cid = paste(cluster, collapse = "/"), .groups = "drop") %>%   # C19/C29 for multi-cluster
  mutate(label = paste0(cell_type, " (", cid, ")"))
lab_of <- setNames(clabel$label, clabel$cell_type)
labelise <- function(ct) ifelse(ct %in% names(lab_of), unname(lab_of[ct]), ct)

# ---- immune vs non-immune (CD45+) -------------------------------------------
nonimmune_rx <- paste(c("Endothelial","Fibroblast","Hepatocyte","Cholangiocyte",
                        "Adipocyte","Mesothelial","Erythrocyte"), collapse = "|")
is_immune <- function(ct) !grepl(nonimmune_rx, ct, ignore.case = TRUE)

lineage_of <- function(ct) dplyr::case_when(
  grepl("T cell",      ct)                     ~ "T cells",
  grepl("NK",          ct)                     ~ "NK cells",
  grepl("Naive B|B cell", ct)                  ~ "B cells",
  grepl("Plasma",      ct)                     ~ "Plasma cells",
  grepl("Monocyte",    ct)                     ~ "Monocytes",
  grepl("Macrophage",  ct)                     ~ "Macrophages",
  grepl("cDC|pDC|Dendritic", ct)               ~ "Dendritic cells",
  grepl("Neutrophil",  ct)                     ~ "Neutrophils",
  grepl("Mast",        ct)                     ~ "Mast cells",
  grepl("Cycling|Proliferat", ct)              ~ "Cycling cells",
  TRUE                                         ~ "Other immune")

# stable cell-type order (lineage-grouped) + palette; keyed by BASE cell_type,
# then re-keyed to the "(C#)" labels once the label map is known.
fine_levels <- c(
  "Naive/TCM T cells","Memory T cells","Activated T cells",
  "Activated CD8 T cells","Cytotoxic CD8 T cells",              # T (reds)
  "Cytotoxic NK cells",                                          # NK
  "Naive B cells","Plasma cells",                               # B / plasma
  "Classical Monocytes","Non-classical Monocytes","Inflammatory Monocytes",
  "SPP1+ Macrophages","FOLR2+ Resident Macrophages","TREM2+ Macrophages",
  "cDC2","pDC","Neutrophils","Mast cells","Cycling cells")
fine_cols <- c(
  "Naive/TCM T cells"="#7f0000","Memory T cells"="#b30000","Activated T cells"="#d7301f",
  "Activated CD8 T cells"="#ef6548","Cytotoxic CD8 T cells"="#fc8d59",
  "Cytotoxic NK cells"="#6a51a3",
  "Naive B cells"="#3182bd","Plasma cells"="#6baed6",
  "Classical Monocytes"="#00441b","Non-classical Monocytes"="#238b45","Inflammatory Monocytes"="#66c2a4",
  "SPP1+ Macrophages"="#005824","FOLR2+ Resident Macrophages"="#41ab5d","TREM2+ Macrophages"="#a1d99b",
  "cDC2"="#fec44f","pDC"="#fe9929","Neutrophils"="#737373","Mast cells"="#bdbdbd",
  "Cycling cells"="#dd3497")
# re-key order + palette to the cluster-labelled names
fine_levels_lab <- labelise(fine_levels)
names(fine_cols) <- labelise(names(fine_cols))

# ---- load, keep immune, renormalise WITHIN immune (CD45+) -------------------
raw <- imap_dfr(comp_files, function(f, tis)
  readr::read_csv(f, show_col_types = FALSE) %>% mutate(Tissue = tis))

imm <- raw %>%
  filter(is_immune(cell_type)) %>%
  group_by(Tissue, Visit) %>% mutate(pct_immune = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(lineage = lineage_of(cell_type),
         cluster = clabel$cid[match(cell_type, clabel$cell_type)],
         ct_lab  = labelise(cell_type),
         Visit   = factor(Visit,  levels = c("Baseline","Followup")),
         Tissue  = factor(Tissue, levels = c("LIVER","PBMC")))

# fine table (each immune cell type as % of CD45+, with its cluster) ----------
fine <- imm %>%
  mutate(ct_lab = factor(ct_lab, levels = intersect(fine_levels_lab, unique(ct_lab)))) %>%
  arrange(Tissue, Visit, ct_lab)
write.csv(fine %>% select(Tissue, Visit, cell_type, cluster, n, pct_immune),
          file.path(tab_dir, "Fig5c_immune_composition_fine_v7.csv"), row.names = FALSE)

# lineage table + baseline->follow-up delta (kept as a data table only) -------
lin <- imm %>% group_by(Tissue, Visit, lineage) %>%
  summarise(pct_immune = sum(pct_immune), .groups = "drop")
lin_delta <- lin %>%
  pivot_wider(names_from = Visit, values_from = pct_immune, values_fill = 0) %>%
  mutate(delta_FU_minus_BL = round(Followup - Baseline, 2),
         Baseline = round(Baseline, 2), Followup = round(Followup, 2)) %>%
  arrange(Tissue, desc(abs(delta_FU_minus_BL)))
write.csv(lin_delta, file.path(tab_dir, "Fig5c_immune_composition_lineage_v7.csv"), row.names = FALSE)
cat("\n== immune-lineage composition (% of CD45+), baseline vs follow-up ==\n")
print(as.data.frame(lin_delta), row.names = FALSE)

# ---- stacked, cluster-wise bar (segment labels show % for readable segments)-
stack_fine <- function(TIS) {
  d <- droplevels(filter(fine, Tissue == TIS))
  ggplot(d, aes(Visit, pct_immune, fill = ct_lab)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    geom_text(aes(label = ifelse(pct_immune >= 4, sprintf("%.0f", pct_immune), "")),
              position = position_stack(vjust = 0.5), size = 1.9, colour = "white") +
    scale_fill_manual(values = fine_cols, name = "Cell type (cluster)", drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "% of immune (CD45+) cells",
         title = sprintf("%s FNA - immune composition", TIS)) +
    theme_nature(7) +
    theme(legend.position = "right", legend.key.size = unit(0.30, "cm"),
          legend.text = element_text(size = 5.2))
}

save_nature(stack_fine("LIVER"), file.path(fig_dir, "Fig5c_immune_composition_LIVER_v7"), 9, 8)  # main 5c
save_nature(stack_fine("PBMC"),  file.path(fig_dir, "Fig5c_immune_composition_PBMC_v7"),  9, 8)  # supp

cat("\nDone Fig 5c immune composition (LIVER main / PBMC supp; stacked, cluster-labelled).\n")
# End of the script

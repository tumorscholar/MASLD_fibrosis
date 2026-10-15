###############################################################################
# 24 - Extended Data Figure 7 - non-expanded T cells (CD4/CD8 + Trm markers)
#
#  a : CD4 vs CD8 proportion of NON-EXPANDED T cells, per tissue x stage
#  b : Tissue-residency (Trm) marker DotPlot2 across tissue x stage
#
# (The original script also had an fgsea translation-pathway panel (6c); that
#  needs the MSigDB gene sets and is left as a separate analysis - not ported
#  here. Ask if you want it.)
#
# Input : TcellObjV7.rds
###############################################################################
library(Seurat); library(SeuratExtend); library(tidyverse)
source("00_nature_style.R")

tcell_rds <- "/data/Blizard-AlazawiLab/rk/seurat/TcellObjV7.rds"
fig_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_Fibrosis","9991"="Fibrosis",
  "9999"="No_Fibrosis","10113-1"="Fibrosis","10113-2"="No_Fibrosis","10202"="Fibrosis",
  "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
  "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
  "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")
stage_lev <- c("Healthy","No_Fibrosis","Fibrosis")

SeuObj <- readRDS(tcell_rds)
if (is.null(SeuObj$Stage)) SeuObj$Stage <- unname(stage_map[as.character(SeuObj$Patient_ID)])
SeuObj$Stage <- gsub("No_fibrosis", "No_Fibrosis", SeuObj$Stage)
DefaultAssay(SeuObj) <- "RNA"

# CD4 / CD8 gate (log-norm > 0.25)
d4 <- FetchData(SeuObj, vars = "CD4", layer = "data")
d8 <- FetchData(SeuObj, vars = c("CD8A","CD8B"), layer = "data")
SeuObj$T_cell_type <- dplyr::case_when(
  d4$CD4 > 0.25 ~ "CD4",
  (d8$CD8A > 0.25 | d8$CD8B > 0.25) ~ "CD8", TRUE ~ "Other")

# ---- ED Fig 6a: CD4 vs CD8 of NON-expanded T cells --------------------------
pd <- SeuObj@meta.data %>%
  filter(Clonality == "Non-expanded", T_cell_type %in% c("CD4","CD8")) %>%
  count(Tissue, Stage, T_cell_type) %>%
  group_by(Tissue, Stage) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(Stage = factor(Stage, levels = stage_lev),
         Tissue = factor(Tissue, levels = c("LIVER","PBMC","SAT","VAT")))
p6a <- ggplot(pd, aes(Tissue, pct, fill = T_cell_type)) +
  geom_col(position = "dodge") + facet_grid(Stage ~ ., switch = "y") +
  scale_fill_manual(values = c("CD4" = npg_pal[1], "CD8" = npg_pal[2]), name = "T cell") +
  labs(x = NULL, y = "Percentage (%)") + theme_nature(7) +
  theme(strip.placement = "outside", strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
save_nature(p6a, file.path(fig_dir, "ExtendedDataFig6a_NonExpanded_CD4CD8_v7"), 8, 8)

# ---- ED Fig 6b: Trm marker dotplot ------------------------------------------
trm_genes <- c("NR4A1","BHLHE40","RUNX3","KLF2","S1PR1","SELL","CD44","ITGB1",
               "CD101","CXCR6","ITGA1","ITGAE","CD69")
SeuObj$Tissue.Stage <- factor(paste(SeuObj$Tissue, SeuObj$Stage, sep = "."),
  levels = as.vector(t(outer(c("LIVER","PBMC","SAT","VAT"), stage_lev, paste, sep = "."))))
p6b <- DotPlot2(SeuObj, features = trm_genes, group.by = "Tissue.Stage") +
  labs(x = NULL, y = NULL) + theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_nature(p6b, file.path(fig_dir, "ExtendedDataFig6b_NonExpanded_Trm_v7"), 10, 10)
cat("\nDone Extended Data Fig 6 (a,b). fgsea panel 6c not ported.\n")
# End of the script

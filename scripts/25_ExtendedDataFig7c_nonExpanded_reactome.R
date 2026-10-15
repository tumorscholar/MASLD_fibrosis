###############################################################################
# 25 - Extended Data Figure 7c - Reactome FGSEA of non-expanded T cells
#
# Ported from the original 32_ExtendedDataFigure6.R (panel 6c). Panels 6a/6b are
# already plotted elsewhere - this script builds 6c ONLY.
#
# Non-expanded T cells from No_Fibrosis + Fibrosis are ranked PER TISSUE by mean
# expression and run through fgsea against Reactome (msigdbr C2 CP:REACTOME).
# The pathways shown are then PICKED FROM THE FGSEA RESULT (top significant
# pathways), not a hard-coded list.
#
# HONEST METHOD NOTE: ranking by mean expression (not a differential contrast)
# is dominated by the most abundant transcripts, so ribosomal / translation
# pathways score highest in EVERY tissue essentially by construction. This panel
# therefore describes the dominant expressed programme of non-expanded T cells,
# NOT a difference between stages/tissues. (A stage-contrast fgsea would be the
# way to make a differential claim - ask if you want that version.)
#
# Input : TcellObjV7.rds
###############################################################################
library(Seurat); library(tidyverse); library(msigdbr); library(fgsea)
library(viridis); library(Matrix)
source("00_nature_style.R")

tcell_rds <- "/data/Blizard-AlazawiLab/rk/seurat/TcellObjV7.rds"
fig_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
tab_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

n_show <- 20            # how many pathways to display (picked from the fgsea result)

SeuObj <- readRDS(tcell_rds)
SeuObj$Stage <- gsub("No_fibrosis", "No_Fibrosis", SeuObj$Stage)
DefaultAssay(SeuObj) <- "RNA"
SeuObj <- JoinLayers(SeuObj)          # v5: collapse split layers so counts is one matrix

# ---- non-expanded T cells, No_Fibrosis + Fibrosis only ----------------------
nonexp <- subset(SeuObj, subset = Stage %in% c("No_Fibrosis", "Fibrosis") &
                  Clonality == "Non-expanded")
cat("Non-expanded cells (No_Fib+Fib) per tissue:\n"); print(table(nonexp$Tissue))

# ---- Reactome gene sets (defensive across msigdbr API versions) -------------
reactome_df <- tryCatch(
 msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
 error = function(e)
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
cat("Reactome gene sets:", length(reactome_sets), "\n")

# ---- per-tissue fgsea on mean-expression ranks ------------------------------
tissues <- c("LIVER", "PBMC", "SAT", "VAT")
fgsea_all <- list()
for (i in seq_along(tissues)) {
 tt <- tissues[i]
 set.seed(123 + i)                                   # reproducible tie-jitter/permutation
 cells_tt <- colnames(nonexp)[nonexp$Tissue == tt]
 if (length(cells_tt) < 50) { cat("  <50 cells for", tt, "- skipped.\n"); next }
 obj_t <- subset(nonexp, cells = cells_tt)
 ranks <- Matrix::rowMeans(GetAssayData(obj_t, assay = "RNA", layer = "counts"))
 ranks <- ranks + runif(length(ranks), -1e-6, 1e-6)  # break ties reproducibly
 ranks <- sort(ranks, decreasing = TRUE)
 fg <- fgsea(pathways = reactome_sets, stats = ranks,
             minSize = 15, maxSize = 500, nPermSimple = 10000)
 fg$Tissue <- tt
 fgsea_all[[tt]] <- fg
 cat("  ", tt, ": ", nrow(obj_t@meta.data), " cells, ",
     sum(fg$padj < 0.05, na.rm = TRUE), " pathways padj<0.05\n", sep = "")
}
fgsea_all <- bind_rows(fgsea_all) %>% filter(!is.na(padj))

# full results table + per-tissue top-10 (for the record)
write.csv(fgsea_all %>% select(Tissue, pathway, NES, padj, pval, size),
          file.path(tab_dir, "EDFig6c_fgsea_reactome_nonexpanded_ALL_v7.csv"), row.names = FALSE)
top10 <- fgsea_all %>% group_by(Tissue) %>% slice_min(padj, n = 10, with_ties = FALSE) %>% ungroup()
write.csv(top10 %>% select(Tissue, pathway, NES, padj, pval, size),
          file.path(tab_dir, "EDFig6c_fgsea_reactome_top10_perTissue_v7.csv"), row.names = FALSE)

# ---- PICK pathways FROM the fgsea result ------------------------------------
# significant (padj<0.05) in >=1 tissue; take the n_show with the smallest best
# padj across tissues, then order rows by mean NES.
picked <- fgsea_all %>% group_by(pathway) %>%
 summarise(best_padj = min(padj), mean_NES = mean(NES), .groups = "drop") %>%
 filter(best_padj < 0.05) %>% slice_min(best_padj, n = n_show, with_ties = FALSE) %>%
 arrange(mean_NES)
cat("Displaying", nrow(picked), "pathways picked from the fgsea result.\n")

plot_df <- fgsea_all %>% filter(pathway %in% picked$pathway) %>%
 mutate(pathway = factor(gsub("REACTOME_", "", pathway),
                         levels = gsub("REACTOME_", "", picked$pathway)),
        Tissue  = factor(Tissue, levels = tissues))

p <- ggplot(plot_df, aes(Tissue, pathway, size = -log10(padj), color = NES)) +
 geom_point(alpha = 0.9) +
 scale_color_viridis(option = "viridis", direction = 1) +
 scale_size(range = c(2, 6)) +
 labs(x = "Tissue", y = "Reactome pathway", color = "NES",
      size = expression(-log[10]~adj.~p),
      title = "Non-expanded T cells - Reactome enrichment") +
 theme_nature(7) +
 theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 8),
       panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 8))
save_nature(p, file.path(fig_dir, "ExtendedDataFig6c_NonExpanded_Reactome_fgsea_v7"), 20.3, 12.7)

cat("\nDone Extended Data Fig 6c (pathways picked from fgsea result).\n")
# End of the script
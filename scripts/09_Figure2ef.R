###############################################################################
# Figure 2e,f - markers of clonal expansion in T cells
#
# Fig 2e : RNA markers,            expanded vs non-expanded T cells (by tissue)
# Fig 2f : surface protein (ADT) markers, expanded vs non-expanded (by tissue)
#
# Clonal-expansion status (Expanded = CDR3aa observed in more than one cell) is
# read from the T-cell object. ADT uses the cleaned, CLR-normalised panel
# (isotype controls and hashtag oligos removed).
#
# Input : T-cell Seurat object with clonal-expansion annotation (RNA + ADT)
# Output: differential-expression tables and Figures 2e, 2f (PDF + SVG + TIFF)
###############################################################################

library(Seurat)
library(tidyverse)
library(SeuratExtend)          # DotPlot2
library(patchwork)
source("00_nature_style.R")    # theme_nature(), npg_pal

# ---- CONFIG (v7) ------------------------------------------------------------
tcell_rds <- "/data/Blizard-AlazawiLab/rk/seurat/TcellObjV7.rds"
out_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"   # figure panels
tab_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"    # data tables
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# panel sizes (cm). Your original was 5x5 cm; with 22 RNA features that is very
# cramped, so defaults are a little larger - change here if you want the 5x5.
e_w <- 16;  e_h <- 9     # Fig 2e (RNA, 22 genes) - wider so dots aren't squeezed
f_w <- 11;  f_h <- 7     # Fig 2f (ADT, 10 markers)

obj <- readRDS(tcell_rds)
obj$Clonality <- factor(obj$Clonality, levels = c("Non-expanded", "Expanded"))

# saver: PDF + SVG + TIFF (same three formats as your original) ---------------
save3 <- function(p, name, w_cm, h_cm) {
  grDevices::cairo_pdf(file.path(out_dir, paste0(name, ".pdf")), width = w_cm/2.54, height = h_cm/2.54)
  print(p); dev.off()
  svglite::svglite(file.path(out_dir, paste0(name, ".svg")), width = w_cm/2.54, height = h_cm/2.54)
  print(p); dev.off()
  tiff(file.path(out_dir, paste0(name, ".tiff")), width = w_cm, height = h_cm,
       units = "cm", res = 600, compression = "lzw")
  print(p); dev.off()
  message("Saved: ", name)
}

# common styling -------------------------------------------------------------
# - drop the DotPlot2 default axis titles ("Var1" / "group"); the tick labels
#   already name the genes/ADT (y) and tissues (x).
# - cap the z-score colour scale symmetrically (limits + squish) so the ADT
#   panel isn't over-saturated (it otherwise runs to ~2 vs ~1.5 for RNA).
zlim <- c(-1.5, 1.5)
style_dot <- function(p) {
  p +
    scale_fill_gradientn(colours = c("white", "#c6dbef", "#4292c6", "#08519c", "#3f007d"),
                         limits = zlim, oob = scales::squish, name = "z-score") +
    scale_colour_manual(values = c("Expanded" = npg_pal[1], "Non-expanded" = npg_pal[2]),
                        name = "Clonality") +
    labs(x = NULL, y = NULL) +                        # remove "group" / "Var1"
    theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Expanded vs Non-expanded T cells")
}

# =============================== Fig 2e (RNA) ================================
DefaultAssay(obj) <- "RNA"

expanded_genes <- c("FGFBP2","GZMH","GNLY","ADGRG1","FCRL6","GZMB","CX3CR1","S1PR5",
                    "NKG7","KLRD1","PRF1","TBX21","PLEK","FCGR3A","ZEB2","CST7",
                    "CCL5","EFHD2","CTSW","SPON2","CD8A","CD4")
miss_rna <- setdiff(expanded_genes, rownames(obj[["RNA"]]))
if (length(miss_rna)) message("RNA features not found (skipped): ", paste(miss_rna, collapse = ", "))
genes_use <- expanded_genes[expanded_genes %in% rownames(obj[["RNA"]])]   # keep your order

# differential expression (reference table, same thresholds as your script)
de_rna <- FindMarkers(obj, ident.1 = "Expanded", ident.2 = "Non-expanded",
                      group.by = "Clonality", logfc.threshold = 0.25, min.pct = 0.25)
de_rna <- de_rna[order(de_rna$avg_log2FC, decreasing = TRUE), ]
write.csv(de_rna, file.path(tab_dir, "DE_expanded_vs_nonexpanded_RNA_V7.csv"))
cat("Top RNA (Expanded vs Non-expanded):\n"); print(head(de_rna, 10))

p_e <- style_dot(DotPlot2(obj, features = genes_use, group.by = "Tissue",
                          split.by = "Clonality", legend_order = c("fill","color","size")))
save3(p_e, "Fig2e_Expanded_vs_Nonexpanded_RNA_V7", e_w, e_h)

# =============================== Fig 2f (ADT) ================================

if (!"ADTonly" %in% names(obj@assays))
  stop("No 'ADTonly' assay in TcellObjV7. Assays present: ",
       paste(names(obj@assays), collapse = ", "))
DefaultAssay(obj) <- "ADTonly"

adt_features <- c("Hu.CD57","Hu.CD8","Hu.KLRG1","Hu.CD11a","Hu.CD18",
                  "Hu.CD4-RPA.T4","Hu.CD7","Hu.CD27","Hu.CD62L")
miss_adt <- setdiff(adt_features, rownames(obj[["ADTonly"]]))
if (length(miss_adt)) message("ADT features not found (skipped): ", paste(miss_adt, collapse = ", "))
adt_use <- adt_features[adt_features %in% rownames(obj[["ADTonly"]])]

de_adt <- FindMarkers(obj, ident.1 = "Expanded", ident.2 = "Non-expanded",
                      group.by = "Clonality", logfc.threshold = 0.25, min.pct = 0.25)
de_adt <- de_adt[order(de_adt$avg_log2FC, decreasing = TRUE), ]
write.csv(de_adt, file.path(tab_dir, "DE_expanded_vs_nonexpanded_ADT_V7.csv"))
cat("\nTop ADT (Expanded vs Non-expanded):\n"); print(head(de_adt, 10))

p_f <- style_dot(DotPlot2(obj, features = adt_use, group.by = "Tissue",
                          split.by = "Clonality", legend_order = c("fill","color","size")))
save3(p_f, "Fig2f_Expanded_vs_Nonexpanded_ADT_V7", f_w, f_h)

# End of the script
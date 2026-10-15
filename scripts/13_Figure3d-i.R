###############################################################################
# Figure 3d-i - transcriptional / surface remodelling of shared-tissue clones
#
# Reads the shared-class-annotated T-cell object built by 11_Figure3a-c.R
# (expandedTcellFinalV7.rds), which labels every expanded clone as
# "Shared tissue" / "One tissue" / "PBMC".
#
#   d : migration/adhesion GENES   across PBMC / One tissue / Shared tissue
#   e : migration/adhesion ADTs    across PBMC / One tissue / Shared tissue
#   f : effector-memory GENES      One tissue vs Shared tissue (solid tissues)
#   g : innate / interleukin GENES One tissue vs Shared tissue (solid tissues)
#   h : Fibrosis vs No Fibrosis GENES  within shared-tissue cells (by tissue)
#   i : Fibrosis vs No Fibrosis ADTs   within shared-tissue cells (by tissue)
#
# DotPlot2 default axis names ("Var1"/"group") are dropped; fill = z-score
# (white -> blue), size = % expressed. Curated marker panels kept as published.
###############################################################################

library(Seurat)
library(SeuratExtend)
library(tidyverse)
source("00_nature_style.R")            # theme_nature(), npg_pal, save_nature()

# ---- v7 data locations ------------------------------------------------------
obj_rds  <- "/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinalV7.rds"
fig_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure3"
tab_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure3"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# PDF (vector) + SVG (editable) + 600-dpi TIFF
save3 <- function(p, file, w_cm, h_cm) {
  grDevices::cairo_pdf(file.path(fig_dir, paste0(file, ".pdf")), width = w_cm/2.54, height = h_cm/2.54)
  print(p); grDevices::dev.off()
  svglite::svglite(file.path(fig_dir, paste0(file, ".svg")), width = w_cm/2.54, height = h_cm/2.54)
  print(p); grDevices::dev.off()
  grDevices::tiff(file.path(fig_dir, paste0(file, ".tiff")), width = w_cm, height = h_cm,
                  units = "cm", res = 600, compression = "lzw")
  print(p); grDevices::dev.off()
  message("Saved: ", file)
}

# common DotPlot2 styling: drop "Var1"/"group" titles; z-score capped to a
# symmetric [-1.5, 1.5] scale (squish) so ADT panels don't over-saturate -
# matches the Fig 2e/f convention.
zscore_fill <- scale_fill_gradientn(
  colours = c("white", "#c6dbef", "#4292c6", "#08519c", "#3f007d"),
  limits = c(-1.5, 1.5), oob = scales::squish, name = "z-score")
style_dot <- function(p, flip = TRUE) {
  p <- p +
    zscore_fill +
    labs(x = NULL, y = NULL) +
    theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (flip) p <- p + coord_flip()
  p
}

obj <- readRDS(obj_rds)
# object is saved with original labels (a_/b_/c_); recode to friendly names here
obj$expanded_shared_class <- factor(dplyr::recode(as.character(obj$expanded_shared_class),
    "a_PBMC" = "PBMC", "b_One_solid_tissue" = "One tissue",
    "c_Shared_between_solid_tissues" = "Shared tissue"),
  levels = c("Shared tissue", "One tissue", "PBMC"))
cat("shared-class cells:\n"); print(table(obj$expanded_shared_class))

# ============================ Fig 3d: migration genes =======================
migration_genes <- c("CCR7","SELL","S1PR1",
                     "ITGAL","ITGB2","ITGA4","ITGB7","CD44","CXCR3",
                     "CXCR6","CCR5","ITGAE","CX3CR1",
                     "ICAM1","VCAM1","ITGB1","ADGRE5")
DefaultAssay(obj) <- "RNA"
p_d <- style_dot(DotPlot2(obj, features = migration_genes, group.by = "expanded_shared_class") +
                 ggtitle("Migration and adhesion associated genes"))
save3(p_d, "Fig3d_migration_genes_v7", w_cm = 12, h_cm = 8)

# ============================ Fig 3e: migration ADTs ========================
migration_ADTs <- c("Hu.CD11a","Hu.CD11b","Hu.CD11c","Hu.CD18","Hu.CD29",
                    "Hu.CD49a","Hu.CD49b","Hu.CD49d","Hu.CX3CR1","Hu.CD62L")
DefaultAssay(obj) <- "ADTonly"
p_e <- style_dot(DotPlot2(obj, features = migration_ADTs, group.by = "expanded_shared_class",
                          legend_order = c("fill","size")) +
                 ggtitle("Migration and adhesion associated ADTs"))
save3(p_e, "Fig3e_migration_ADTs_v7", w_cm = 12, h_cm = 8)

# ==================== Fig 3f/g: solid-tissue One vs Shared ===================
obj_solid <- subset(obj, cells = colnames(obj)[obj$Tissue %in% c("LIVER","SAT","VAT") &
                                               obj$expanded_shared_class %in% c("Shared tissue","One tissue")])
obj_solid$expanded_shared_class <- factor(obj_solid$expanded_shared_class,
                                          levels = c("Shared tissue","One tissue"))
DefaultAssay(obj_solid) <- "RNA"
cat("\nsolid-tissue One vs Shared cells:\n"); print(table(obj_solid$expanded_shared_class))

effector_memory_genes <- c("GZMH","GNLY","FGFBP2","GZMB","PRF1","CCL5","NKG7","CTSW","CST7",
                           "TBX21","KLF2","HCST","HOPX","EFHD2","RASSF1","PLEK")
p_f <- style_dot(DotPlot2(obj_solid, features = effector_memory_genes, group.by = "expanded_shared_class",
                          legend_order = c("fill","size")) +
                 ggtitle("Effector memory activation genes"))
save3(p_f, "Fig3f_effector_memory_genes_v7", w_cm = 12, h_cm = 8)

innate_il_genes <- c("TLR2","TLR4","NLRP3","IL1B","IL6","IL10","IL12A","IL12B","IL18",
                     "IFNB1","CXCL8","CCL2","CCL3","STAT1","IRF7")
p_g <- style_dot(DotPlot2(obj_solid, features = innate_il_genes, group.by = "expanded_shared_class",
                          legend_order = c("fill","size")) +
                 ggtitle("Innate and interleukin genes"))
save3(p_g, "Fig3g_innate_il_genes_v7", w_cm = 12, h_cm = 8)

# ==================== Fig 3h/i: Fibrosis vs No Fibrosis (shared) =============
# shared-tissue cells in solid tissues, MASLD only (drop Healthy).
obj_shared <- subset(obj_solid, cells = colnames(obj_solid)[
  obj_solid$expanded_shared_class == "Shared tissue" & obj_solid$Stage != "Healthy"])
obj_shared$Stage  <- factor(as.character(obj_shared$Stage), levels = c("Fibrosis","No_fibrosis"))
obj_shared$Tissue <- factor(obj_shared$Tissue, levels = c("LIVER","SAT","VAT"))
cat("\nshared-tissue MASLD cells (Fib/No-fib):\n"); print(table(obj_shared$Stage))

stage_outline <- scale_colour_manual(values = c("Fibrosis" = "#084594", "No_fibrosis" = "#9ecae1"),
                                     name = "Stage")

# ---- Fig 3h: genes, Fibrosis vs No Fibrosis --------------------------------
DefaultAssay(obj_shared) <- "RNA"
mk_gene <- FindMarkers(obj_shared, ident.1 = "Fibrosis", ident.2 = "No_fibrosis",
                       group.by = "Stage", test.use = "wilcox")
write.csv(mk_gene %>% rownames_to_column("gene") %>% filter(p_val_adj < 0.05),
          file.path(tab_dir, "Fig3h_shared_Fib_vs_NoFib_genes_v7.csv"), row.names = FALSE)

fibrosis_genes <- c("TRAV38-1",
  "IFITM3",
  "ZEB2",
  "XCL2",
  "CCL3L1",
  "CCL4L2",
  "JUN",
  "FOS",
  "EGR1",
  "CD27",
  "TMIGD2"
)


p_h <- DotPlot2(obj_shared, features = fibrosis_genes, group.by = "Tissue", split.by = "Stage",
                legend_order = c("fill","color","size")) +
  zscore_fill + stage_outline +
  labs(x = NULL, y = NULL) + theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fibrosis vs No Fibrosis - shared T cells")
save3(p_h, "Fig3h_fibrosis_shared_genes_v7", w_cm = 12, h_cm = 8)

# ---- Fig 3i: ADTs, Fibrosis vs No Fibrosis ---------------------------------
DefaultAssay(obj_shared) <- "ADTonly"
adt_cells  <- colnames(obj_shared[["ADTonly"]])
obj_adt    <- subset(obj_shared, cells = adt_cells)
obj_adt$Tissue <- factor(obj_adt$Tissue, levels = c("LIVER","SAT","VAT"))
mk_adt <- FindMarkers(obj_adt, ident.1 = "Fibrosis", ident.2 = "No_fibrosis",
                      group.by = "Stage", test.use = "wilcox")
write.csv(mk_adt %>% rownames_to_column("adt") %>% filter(p_val_adj < 0.05),
          file.path(tab_dir, "Fig3i_shared_Fib_vs_NoFib_ADTs_v7.csv"), row.names = FALSE)

fibrosis_adts <- c("Hu.CD7","Hu.CD45RA",
  "Hu.CD36",
  "Hu.CD40",
  "Hu.CD95",
  "Hu.CD69",
  "Hu.HLA.ABC",
  "Hu.KLRG1",
  "Hu.CD25",
  "Hu.CD39"
)

p_i <- DotPlot2(obj_adt, features = fibrosis_adts, group.by = "Tissue", split.by = "Stage",
                legend_order = c("fill","color","size")) +
  zscore_fill + stage_outline +
  labs(x = NULL, y = NULL) + theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Fibrosis vs No Fibrosis - shared T cells")
save3(p_i, "Fig3i_fibrosis_shared_ADTs_v7", w_cm = 12, h_cm = 8)

# End of the script

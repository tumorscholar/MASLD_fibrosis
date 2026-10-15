###############################################################################
# 18 - Figure 4a,f  (CellChat visualisation, all tissues in one run)
#
# Loops LIVER -> SAT -> VAT (order fixed) and makes:
#   Fig 4a : 2-D signaling-role scatter (outgoing vs incoming), one per stage,
#            per tissue. Shared, widened x/y limits so panels are comparable.
#   Fig 4f : differential ligand-receptor bubbles, Endothelial <-> T cell,
#            No Fibrosis vs Fibrosis - TWO per tissue (Endo->T, T->Endo) = 6,
#            saved individually AND as one ordered 2 x 3 panel
#            (row 1 = Endo->T, row 2 = T->Endo; cols = LIVER, SAT, VAT).
#
# Sources and targets are specified by cell-type name (endothelial cell types
# and "T cell"), not by cluster index.
#
# Input : cellchat_<tissue>_<stage>_V7.rds   (from script 17)
# Output: Fig4a per tissue + Fig4f 6 bubbles (+ combined) under Figure4/
###############################################################################
library(Seurat); library(CellChat); library(ggplot2); library(patchwork)
source("00_nature_style.R")
options(stringsAsFactors = FALSE)

# ---- Fig-1a cell-type colours (name -> colour), + pooled "T cell" -----------
# so each cell type keeps ONE colour across every stage and tissue.
fig1_celltype_cols <- function() {
  ann  <- read.csv("/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1/cluster_annotation_V7.csv",
                   stringsAsFactors = FALSE)
  cl2t <- setNames(ann$cell_type, sub("^C", "", ann$cluster))
  lineage_order_ids <- c(0,2,1,6,20, 4, 8,19,29, 3,11,23, 7,9,26, 14,25, 13,18,
                         12,27,24,28, 5,10,15, 17,21, 16,22)
  my_colors <- c("#7f0000","#b30000","#d7301f","#ef6548","#fc8d59","#6a51a3",
                 "#08519c","#3182bd","#6baed6","#00441b","#238b45","#66c2a4",
                 "#41ab5d","#74c476","#a1d99b","#addd8e","#78c679","#fec44f",
                 "#fe9929","#8c510a","#bf812d","#dfc27d","#f6e8c3","#01665e",
                 "#35978f","#80cdc1","#d94801","#fd8d3c","#737373","#bdbdbd")
  nm <- unname(cl2t[as.character(lineage_order_ids)]); keep <- !duplicated(nm)
  c(setNames(my_colors[keep], nm[keep]), "T cell" = "#b30000")   # pooled T cells -> one red
}
ct_cols <- fig1_celltype_cols()

tissues    <- c("LIVER", "SAT", "VAT")                       # fixed panel order
cc_root    <- "/data/Blizard-AlazawiLab/rk/cellchat_v7"
fig_dir    <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure4"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
endo_types <- c("Capillary Endothelial cells", "Venous Endothelial cells",
                "TAGLN Endothelial cells")

rd <- function(tissue, st) {
  f <- file.path(cc_root, tolower(tissue), "sharedTcells",
                 paste0("cellchat_", tolower(tissue), "_", st, "_V7.rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

# coord_flip() puts the L-R pairs on the horizontal axis (matches the original
# 22_Figure4a,f script). Combined with a single fixed canvas size for all six
# plots (see the save loop below), every Fig 4f panel comes out the SAME size.
style_bubble <- function(p, ttl) p +
  coord_flip() +
  scale_color_gradient(low = "mistyrose", high = "maroon", name = "Comm. prob.") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) + ggtitle(ttl) +
  theme(legend.position = "right", legend.key.size = unit(0.8, "lines"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 9, hjust = 0.5))
# panel height grows with the number of L-R pairs so dots never squish
n_pairs <- function(p) { d <- p$data
nm <- intersect(c("interaction_name","interaction_name_2","interaction"), colnames(d))
if (length(nm)) length(unique(d[[nm[1]]])) else 25 }
bub_h <- function(p) max(8, 0.42 * n_pairs(p) + 3)          # cm
save3 <- function(p, file, w = 16, h = 6) {
  grDevices::cairo_pdf(file.path(fig_dir, paste0(file, ".pdf")), width = w/2.54, height = h/2.54); print(p); dev.off()
  svglite::svglite(file.path(fig_dir, paste0(file, ".svg")), width = w/2.54, height = h/2.54); print(p); dev.off()
  ggsave(file.path(fig_dir, paste0(file, ".tiff")), p, width = w/2.54, height = h/2.54,
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
}

# save a bubble with a FIXED panel (bubble grid) size so all six Fig 4f dot plots
# have IDENTICAL panel width regardless of how long the y-axis labels are or how
# many L-R pairs a tissue has. Uses base ggplotGrob + gtable (ships with ggplot2,
# no extra package) to pin the panel column/row to a fixed cm size; the panel
# HEIGHT still grows with the number of pairs so dots stay legible.
PANEL_W_CM <- 6            # <- common bubble-grid width for every Fig 4f panel
save_fixed <- function(p, file) {
  np <- n_pairs(p)
  panel_h <- max(4, 0.32 * np)                       # cm, grows with number of pairs
  g <- ggplotGrob(p)
  pcells <- which(g$layout$name == "panel")
  g$widths[unique(g$layout$l[pcells])]  <- grid::unit(PANEL_W_CM, "cm")
  g$heights[unique(g$layout$t[pcells])] <- grid::unit(panel_h,    "cm")
  W <- PANEL_W_CM + 9        # + room for y-labels and the legend (constant)
  H <- panel_h + 5           # + room for x-labels and title
  grDevices::cairo_pdf(file.path(fig_dir, paste0(file, ".pdf")), width = W/2.54, height = H/2.54)
  grid::grid.newpage(); grid::grid.draw(g); dev.off()
  ggsave(file.path(fig_dir, paste0(file, ".tiff")), g, width = W/2.54, height = H/2.54,
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
}

bub_endoT <- list(); bub_Tendo <- list()                    # collect for the combined panel

for (tissue in tissues) {
  ccH <- rd(tissue, "Healthy"); ccN <- rd(tissue, "No_fibrosis"); ccF <- rd(tissue, "Fibrosis")
  
  # ---- Fig 4a: role scatter (widened, shared limits across the tissue's stages)
  obj.list <- Filter(Negate(is.null), list(Healthy = ccH, No_fibrosis = ccN, Fibrosis = ccF))
  if (length(obj.list) == 0) { message("no CellChat objects for ", tissue); next }
  num.link <- unlist(sapply(obj.list, function(x)
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)))
  wMinMax  <- c(min(num.link, na.rm = TRUE), max(num.link, na.rm = TRUE))
  # Fig-1a colours per cell type (fixed across stage/tissue), fixed 0-40 axes,
  # bigger axis titles (outgoing / incoming interaction strength).
  gg <- lapply(names(obj.list), function(nm) {
    ob  <- obj.list[[nm]]
    grp <- levels(ob@idents)                                 # groups present in this stage
    cu  <- unname(ct_cols[grp]); cu[is.na(cu)] <- "grey70"   # Fig-1a colour per group, in order
    netAnalysis_signalingRole_scatter(ob, title = nm, weight.MinMax = wMinMax, color.use = cu,
                                      label.size = 3, font.size = 8, font.size.title = 9) +
      scale_x_continuous(breaks = seq(0, 25, 5)) + scale_y_continuous(breaks = seq(0, 25, 5)) +
      coord_cartesian(xlim = c(0, 25), ylim = c(0, 25)) +    # same scale everywhere (no data dropped)
      theme_nature(7) +
      theme(axis.title = element_text(size = 10, face = "bold"),
            axis.text  = element_text(size = 9),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.line = element_line(colour = "black"),
            legend.position = "right", legend.key.size = unit(0.9, "lines"))
  })
  fig4a <- wrap_plots(gg, nrow = 1)
  ggsave(file.path(fig_dir, paste0("Fig4a_signalingRole_", tolower(tissue), "_v7.png")),
         fig4a, dpi = 300, width = 18, height = 6.5, bg = "white")
  save_nature(fig4a, file.path(fig_dir, paste0("Fig4a_signalingRole_", tolower(tissue), "_v7")), 45.72, 16.5)
  
  dis <- Filter(Negate(is.null), list(No_fibrosis = ccN, Fibrosis = ccF))
  if (length(dis) != 2) { message("comparison skipped for ", tissue, " - need both disease stages"); next }
  # lift both stages to the UNION of cell groups so differential matrices are
  # conformable (SAT/VAT can have a group present in one stage but not the other).
  grp.union <- union(levels(dis[[1]]@idents), levels(dis[[2]]@idents))
  dis <- setNames(lapply(names(dis), function(n)
    if (setequal(levels(dis[[n]]@idents), grp.union)) dis[[n]] else liftCellChat(dis[[n]], grp.union)),
    names(dis))
  cellchat <- mergeCellChat(dis, add.names = names(dis))
  
  # ---- Extended Data Fig 8: differential interaction heatmap (red up / blue down)
  # red  = increased in Fibrosis vs No Fibrosis; blue = decreased. Count & weight.
  for (meas in c("count", "weight")) {
    ht <- tryCatch(
      netVisual_heatmap(cellchat, comparison = c(1, 2), measure = meas,
                        title.name = paste0("Differential ", meas, " - ", tissue)),
      error = function(e) { message("ED8 ", meas, " ", tissue, " failed: ", e$message); NULL })
    if (is.null(ht)) next
    for (dev_fun in list(
      function() grDevices::cairo_pdf(file.path(fig_dir, paste0("EDFig8_diff_", meas, "_", tolower(tissue), "_v7.pdf")),
                                      width = 13/2.54, height = 11/2.54),
      function() png(file.path(fig_dir, paste0("EDFig8_diff_", meas, "_", tolower(tissue), "_v7.png")),
                     width = 13, height = 11, units = "cm", res = 300))) {
      dev_fun(); ComplexHeatmap::draw(ht); dev.off()
    }
  }
  
  # ---- Fig 4f: Endo <-> T cell bubbles (No Fibrosis vs Fibrosis) --------------
  endo <- intersect(endo_types, levels(cellchat@idents$joint))
  cat(tissue, "endothelial partners:", paste(endo, collapse = ", "), "\n")
  
  # Exclude the CD4-receptor artefacts from BOTH directions. CellChat scores
  # MHC-II (HLA-DR/DP/DQ) and IL16 exclusively onto the CD4 receptor; neither the
  # endothelial clusters nor the cytotoxic-CD8 T cells here are genuinely CD4+, and
  # the endothelial-T axis of interest is MHC-I -> CD8. Dropping these two pathways
  # removes every spurious "... -> CD4" row from all Fig 4f bubbles (Endo->T, T->Endo
  # and the combined all-tissue panel).
  drop_paths <- c("MHC-II", "IL16")
  keepP <- tryCatch(setdiff(unique(unlist(lapply(cellchat@netP, function(x) x$pathways))), drop_paths),
                    error = function(e) NULL)
  # Fig 4f shows the relevant interactions driving the paper: every pathway CellChat
  # detected between endothelial and T cells, MINUS the CD4-receptor artefacts
  # (MHC-II, IL16) already dropped in keepP. A hardcoded pathway allow-list was
  # tried but the adhesion/chemokine pathway names in this dataset did not match it
  # and only MHC-I survived; keepP keeps the full antigen-presentation / adhesion /
  # chemokine / integrin picture the legend describes.
  focus_paths <- keepP
  
  g_endoT <- netVisual_bubble(cellchat, sources.use = endo, targets.use = "T cell",
                              signaling = focus_paths,
                              comparison = c(1, 2), max.dataset = 2, angle.x = 45, remove.isolate = TRUE,
                              dot.size.min = 3, dot.size.max = 7,
                              title.name = paste0("Endothelial -> T cell (", tissue, ")")) |>
    style_bubble(paste0("Endothelial -> T cell (", tissue, ")"))
  g_Tendo <- netVisual_bubble(cellchat, sources.use = "T cell", targets.use = endo,
                              signaling = focus_paths,
                              comparison = c(1, 2), max.dataset = 2, angle.x = 45, remove.isolate = TRUE,
                              dot.size.min = 3, dot.size.max = 7,
                              title.name = paste0("T cell -> Endothelial (", tissue, ")")) |>
    style_bubble(paste0("T cell -> Endothelial (", tissue, ")"))
  
  # collect; all six are saved together below at ONE identical canvas size.
  bub_endoT[[tissue]] <- g_endoT; bub_Tendo[[tissue]] <- g_Tendo
}

# ---- save all six Fig 4f dot plots at ONE identical canvas size --------------
# Matches the original 22_Figure4a,f script: every bubble plot is written to the
# SAME width x height, so all six panels are exactly the same size when placed
# side-by-side in the figure. The common width is set from the busiest plot so
# no plot is cramped; every plot then uses that same width (and a fixed height).
# max.dataset = 2 already restricts each bubble to interactions ENRICHED in
# Fibrosis (dataset 2), which is the paper's point.
all_plots <- c(bub_endoT, bub_Tendo)
all_plots <- all_plots[!vapply(all_plots, is.null, logical(1))]
maxN <- max(vapply(all_plots, n_pairs, numeric(1)), na.rm = TRUE)
Wc   <- max(16, 0.42 * maxN + 6)     # common width  (cm) - same for all six
Hc   <- 6                            # common height (cm) - same for all six
for (tissue in tissues) {
  if (!is.null(bub_endoT[[tissue]]))
    save3(bub_endoT[[tissue]], paste0("Fig4f_EndoToT_", tolower(tissue), "_v7"), w = Wc, h = Hc)
  if (!is.null(bub_Tendo[[tissue]]))
    save3(bub_Tendo[[tissue]], paste0("Fig4f_TToEndo_", tolower(tissue), "_v7"), w = Wc, h = Hc)
}
cat(sprintf("\nDone Fig 4a/f - six individual Fig 4f dot plots, all %.1f x %.1f cm (Endo->T, T->Endo x liver/SAT/VAT).\n", Wc, Hc))
# End of the script
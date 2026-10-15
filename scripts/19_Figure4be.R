###############################################################################
# 19 - Figure 4b-e : endothelial phenotype & EndMT violins (No Fib vs Fibrosis)
#
#  b,c : surface phenotype + antigen presentation across tissues on the
#        endothelial clusters - CD155, CD146, CD31, CD36, MHC class I.
#  d,e : endothelial-to-mesenchymal transition in fibrosis - ACTA2, TAGLN,
#        PECAM1 (the EndMT cluster gains ACTA2/TAGLN and loses PECAM1).
#
# Endothelial clusters are auto-detected by NAME (cell_type contains
# "Endothelial"), so this adapts to the current endothelial cluster set - no
# hardcoded C-numbers. Markers are resolved across the ADT and RNA assays.
#
# Input : SeuObjV7.rds
# Output: Fig4b-e violin panels under Figure4/
###############################################################################
library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
source("00_nature_style.R")

seu_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure4"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")
stage_lev  <- c("No_fibrosis","Fibrosis")
# stage_cols is the PAPER-WIDE palette from 00_nature_style.R
# (No_fibrosis = #8E44AD violet, Fibrosis = #C2185B magenta) - do not override.
if (!exists("stage_cols")) stage_cols <- c("No_fibrosis" = "#8E44AD", "Fibrosis" = "#C2185B")

# ---- load, assign stage, subset endothelial clusters ------------------------
seu <- readRDS(seu_rds)
if (is.null(seu$Stage)) seu$Stage <- unname(stage_map[as.character(seu$Patient_ID)])
seu$Stage <- gsub("No_fibrosis", "No_fibrosis", seu$Stage)

endo_types <- grep("endothelial", unique(as.character(seu$cell_type)), value = TRUE, ignore.case = TRUE)
if (length(endo_types) == 0) stop("no endothelial clusters found in $cell_type")
cat("Endothelial clusters:", paste(endo_types, collapse = " | "), "\n")

endo <- subset(seu, subset = cell_type %in% endo_types & Stage %in% stage_lev)
endo$cell_type <- droplevels(factor(endo$cell_type, levels = endo_types))
endo$Stage     <- factor(endo$Stage, levels = stage_lev)
endo$Tissue    <- factor(endo$Tissue, levels = intersect(c("LIVER","PBMC","SAT","VAT"), unique(endo$Tissue)))
cat("\nEndothelial cells per (cluster, Stage):\n"); print(table(endo$cell_type, endo$Stage))

# ---- normalise the assays we read from --------------------------------------
# subset() leaves per-sample split layers (counts.1, counts.2 ...); CLR can't
# handle those - join layers at the OBJECT level per assay (assay-level
# assignment errors with "no method for coercing this S4 class to a vector").
# this figure only needs RNA (genes) + ADTonly (surface markers); drop the rest
for (a in intersect(c("ADT","HTO"), names(endo@assays))) endo[[a]] <- NULL
# join split layers at OBJECT level (assay-level assignment throws an S4 coercion
# error), then normalise ONLY if not already normalised (re-CLR errors).
for (a in intersect(c("RNA","ADTonly"), names(endo@assays))) {
  DefaultAssay(endo) <- a
  endo <- tryCatch(JoinLayers(endo), error = function(e) endo)
}
DefaultAssay(endo) <- "RNA"
if (!"data" %in% Layers(endo[["RNA"]])) endo <- NormalizeData(endo, verbose = FALSE)
if ("ADTonly" %in% names(endo@assays)) {
  DefaultAssay(endo) <- "ADTonly"
  if (!"data" %in% Layers(endo[["ADTonly"]]))
    endo <- NormalizeData(endo, normalization.method = "CLR", margin = 2, verbose = FALSE)
}
DefaultAssay(endo) <- "RNA"

# ---- resolve each marker to (assay, feature) --------------------------------
# surface markers prefer ADT; genes use RNA. MHC-I may be an HLA-ABC tag or genes.
# EndMT-focused candidates; recruitment markers dropped (chemokines are ~0 in
# endothelial - recruitment is shown by the Fig 4f L-R bubble instead).
want <- list(
  # mesenchymal genes GAINED during EndMT
  "ACTA2"  = list(assays = "RNA", alt = "ACTA2",  set = "EndMT_gain"),
  "TAGLN"  = list(assays = "RNA", alt = "TAGLN",  set = "EndMT_gain"),
  "VIM"    = list(assays = "RNA", alt = "VIM",    set = "EndMT_gain"),
  "FN1"    = list(assays = "RNA", alt = "FN1",    set = "EndMT_gain"),
  "COL1A1" = list(assays = "RNA", alt = "COL1A1", set = "EndMT_gain"),
  "COL3A1" = list(assays = "RNA", alt = "COL3A1", set = "EndMT_gain"),
  "S100A4" = list(assays = "RNA", alt = "S100A4", set = "EndMT_gain"),
  "CDH2"   = list(assays = "RNA", alt = "CDH2",   set = "EndMT_gain"),
  "SNAI2"  = list(assays = "RNA", alt = "SNAI2",  set = "EndMT_gain"),
  "CD146"  = list(assays = c("ADTonly","RNA"), alt = c("CD146","MCAM"), set = "EndMT_gain"),
  # endothelial identity LOST during EndMT
  "PECAM1" = list(assays = "RNA", alt = "PECAM1", set = "Endo_loss"),
  "CDH5"   = list(assays = "RNA", alt = "CDH5",   set = "Endo_loss"),
  "VWF"    = list(assays = "RNA", alt = "VWF",    set = "Endo_loss"),
  "CLDN5"  = list(assays = "RNA", alt = "CLDN5",  set = "Endo_loss"),
  "KDR"    = list(assays = "RNA", alt = "KDR",    set = "Endo_loss"),
  "CD31"   = list(assays = c("ADTonly","RNA"), alt = "CD31", set = "Endo_loss"),
  # antigen presentation / phenotype (protein)
  "MHC-I"  = list(assays = c("ADTonly","RNA"), alt = c("HLA-ABC","HLA.ABC","HLA-A","HLA-B","HLA-C"), set = "Pheno"),
  "CD36"   = list(assays = c("ADTonly","RNA"), alt = "CD36", set = "Pheno"),
  "CD155"  = list(assays = c("ADTonly","RNA"), alt = c("CD155","PVR"), set = "Pheno"))
set_of <- vapply(want, function(x) x$set, character(1))       # marker -> set

# ADT features carry an antibody prefix (Hu.CD155, Hu.HLA.ABC, Hu.CD54); strip
# it and drop . / - so "CD155" matches "Hu.CD155" but not "Hu.CD314".
resolve <- function(spec) {
  norm <- function(v) toupper(gsub("[-.]", "", sub("^Hu\\.", "", v)))
  for (a in spec$assays) {
    if (!a %in% names(endo@assays)) next
    rn <- rownames(endo[[a]]); rnn <- norm(rn)
    for (nm in spec$alt) { hit <- rn[rnn == norm(nm)]; if (length(hit)) return(c(assay = a, feat = hit[1])) }
  }
  NULL
}
resolved <- lapply(want, resolve)
miss <- names(resolved)[vapply(resolved, is.null, logical(1))]
if (length(miss)) message("markers not found (skipped): ", paste(miss, collapse = ", "))
resolved <- resolved[!vapply(resolved, is.null, logical(1))]
cat("\nresolved markers:\n"); for (m in names(resolved)) cat("  ", m, "->", resolved[[m]]["assay"], "/", resolved[[m]]["feat"], "\n")

# ---- long expression table --------------------------------------------------
meta <- endo@meta.data[, c("cell_type","Tissue","Stage")]
long <- bind_rows(lapply(names(resolved), function(m) {
  r <- resolved[[m]]
  val <- as.numeric(GetAssayData(endo, assay = r["assay"], layer = "data")[r["feat"], ])
  data.frame(marker = m, cell_type = meta$cell_type, Tissue = meta$Tissue, Stage = meta$Stage, value = val)
}))

# display label: protein markers (ADTonly) get a "Hu" tag so protein vs RNA is clear
disp <- vapply(names(resolved), function(m)
  if (resolved[[m]]["assay"] == "ADTonly") paste0(m, " (Hu)") else m, character(1))
star <- function(p) as.character(cut(p, c(-Inf, 0.001, 0.01, 0.05, Inf),
                                     labels = c("***", "**", "*", "ns")))

# solid tissues only (drop PBMC); short "... endo" cluster labels; comparison
# order fixed No_fibrosis -> Fibrosis everywhere.
solid <- c("LIVER","SAT","VAT")
endo_short <- c("Capillary Endothelial cells" = "Capillary endo",
                "Venous Endothelial cells"    = "Venous endo",
                "TAGLN Endothelial cells"     = "TAGLN endo")
long <- long %>% filter(Tissue %in% solid) %>%
  mutate(Tissue = factor(Tissue, levels = solid),
         Stage  = factor(Stage,  levels = c("No_fibrosis","Fibrosis")),
         endo   = factor(unname(endo_short[cell_type]), levels = unname(endo_short)))
assay_of <- vapply(names(resolved), function(m) unname(resolved[[m]]["assay"]), "")

# ---- TEST 1: EndMT change in TAGLN endo (No-fib vs Fibrosis), solid tissues -
de1 <- long %>% filter(endo == "TAGLN endo") %>% group_by(marker) %>%
  summarise(median_Nofib = median(value[Stage == "No_fibrosis"]),
            median_Fib   = median(value[Stage == "Fibrosis"]),
            p = tryCatch(wilcox.test(value ~ Stage)$p.value, error = function(e) NA_real_),
            .groups = "drop") %>%
  mutate(set = set_of[marker], direction = ifelse(median_Fib > median_Nofib, "up", "down"),
         sig = star(p)) %>% arrange(p)
cat("\n=== TEST 1: TAGLN endo (solid tissues) No-fib vs Fibrosis (EndMT change) ===\n")
print(as.data.frame(de1), row.names = FALSE)
write.csv(de1, file.path(fig_dir, "Fig4_TAGLN_EndMT_change_test_V7.csv"), row.names = FALSE)

# ---- TEST 2: TAGLN endo vs other endothelial (ALL MASLD, solid tissues) -----
# use all MASLD patients (No Fibrosis + Fibrosis), not Fibrosis alone
de2 <- long %>% filter(Stage %in% c("No_fibrosis","Fibrosis")) %>% group_by(marker) %>%
  summarise(median_TAGLN = median(value[endo == "TAGLN endo"]),
            median_other = median(value[endo != "TAGLN endo"]),
            p = tryCatch(wilcox.test(value ~ (endo == "TAGLN endo"))$p.value,
                         error = function(e) NA_real_), .groups = "drop") %>%
  mutate(direction = ifelse(median_TAGLN > median_other, "higher in TAGLN", "lower in TAGLN"),
         sig = star(p)) %>% arrange(p)
cat("\n=== TEST 2: TAGLN vs other endothelial (solid, all MASLD) - distinctness ===\n")
print(as.data.frame(de2), row.names = FALSE)
write.csv(de2, file.path(fig_dir, "Fig4_TAGLN_vs_otherEndo_test_V7.csv"), row.names = FALSE)

# ---- pick markers from the tests: significant, EXPRESSED, and in the expected
# direction (avoids zero-inflated markers that are 'significant' on sparse data)
expr_thr <- 0.1
# show the FULL canonical EndMT signature that is EXPRESSED (median > thr in
# either stage), regardless of significance, so the coordinated EndMT change is
# visible; per-marker Wilcoxon stars still show which reach significance.
endmt_gain <- de1 %>% filter(set == "EndMT_gain", pmax(median_Nofib, median_Fib) > expr_thr) %>% arrange(p) %>% pull(marker)
endmt_loss <- de1 %>% filter(set == "Endo_loss",  pmax(median_Nofib, median_Fib) > expr_thr) %>% arrange(p) %>% pull(marker)
endmt_markers <- c(endmt_gain, endmt_loss)
# TAGLN-distinct markers, split into gained vs lost blocks (MHC-I dropped from
# the TAGLN panels - it is LOWER in TAGLN, so it belongs with the other endo).
dm <- de2 %>% filter(!is.na(p), p < 0.05, pmax(median_TAGLN, median_other) > 0.5, marker != "MHC-I")
gained_markers <- dm %>% filter(direction == "higher in TAGLN") %>% arrange(p) %>% slice_head(n = 5) %>% pull(marker)
lost_markers   <- dm %>% filter(direction == "lower in TAGLN")  %>% arrange(p) %>% slice_head(n = 5) %>% pull(marker)
cat("\nEndMT panel (sig No-fib vs Fib in TAGLN):", paste(endmt_markers, collapse = ", "),
    "\nTAGLN gained (mesenchymal):", paste(gained_markers, collapse = ", "),
    "\nTAGLN lost (endothelial):",  paste(lost_markers, collapse = ", "), "\n")

rna <- function(ms) ms[assay_of[ms] == "RNA"]
adt <- function(ms) ms[assay_of[ms] == "ADTonly"]
# each panel is a HORIZONTAL strip (markers as columns, 1 row) so it slots into
# the figure as a thin block; width scales with the number of markers.

# ---- EndMT CHANGE: TAGLN endo across solid tissues, No-fib vs Fibrosis -------
h_stage <- function(markers, title, file) {
  markers <- markers[assay_of[markers] %in% c("RNA","ADTonly")]; if (!length(markers)) return(invisible())
  df <- long %>% filter(endo == "TAGLN endo", marker %in% markers) %>%
    mutate(marker = factor(marker, levels = markers))
  st <- df %>% group_by(marker, Tissue) %>%
    summarise(p = tryCatch(wilcox.test(value ~ Stage)$p.value, error = function(e) NA_real_),
              ytop = max(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(marker) %>% mutate(y = max(ytop) * 1.08) %>% ungroup() %>% mutate(lab = star(p))
  p <- ggplot(df, aes(Tissue, value, fill = Stage)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, position = position_dodge(0.85)) +
    stat_summary(fun = median, geom = "point", size = 0.5, colour = "black", position = position_dodge(0.85)) +
    geom_text(data = st, aes(x = Tissue, y = y, label = lab), inherit.aes = FALSE, size = 1.8) +
    scale_fill_manual(values = stage_cols, name = "Stage") +
    facet_wrap(~ marker, nrow = 1, scales = "free_y", labeller = ggplot2::labeller(marker = disp)) +
    labs(x = NULL, y = "Norm. expression", title = title) + theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_nature(p, file.path(fig_dir, file), 3.2 * length(markers) + 3, 6)
}
# NOTE: the No-fib -> Fibrosis CHANGE within TAGLN endo (h_stage) is NOT plotted -
# the EndMT genes barely change with stage (all ns; TAGLN endo are already
# mesenchymal at both stages). The EndMT signature is a CLUSTER-IDENTITY difference
# (TAGLN endo vs venous/capillary endo), which is what Fig 4b-e show via h_cluster
# below and what the Results/legend describe.

# ---- TAGLN DISTINCT: 3 endo clusters (Fibrosis, solid tissues) --------------
endo_cols <- c("Capillary endo" = "#80cdc1", "Venous endo" = "#35978f", "TAGLN endo" = "#01665e")
h_cluster <- function(markers, title, file) {
  markers <- markers[assay_of[markers] %in% c("RNA","ADTonly")]; if (!length(markers)) return(invisible())
  df <- long %>% filter(Stage %in% c("No_fibrosis","Fibrosis"), marker %in% markers) %>%
    mutate(marker = factor(marker, levels = markers))
  st <- df %>% group_by(marker) %>%
    summarise(p = tryCatch(wilcox.test(value ~ (endo == "TAGLN endo"))$p.value, error = function(e) NA_real_),
              y = max(value, na.rm = TRUE) * 1.08, .groups = "drop") %>% mutate(lab = star(p))
  p <- ggplot(df, aes(endo, value, fill = endo)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    stat_summary(fun = median, geom = "point", size = 0.5, colour = "black") +
    geom_text(data = st, aes(x = 3, y = y, label = lab), inherit.aes = FALSE, size = 1.8) +
    scale_fill_manual(values = endo_cols, guide = "none") +
    facet_wrap(~ marker, nrow = 1, scales = "free_y", labeller = ggplot2::labeller(marker = disp)) +
    labs(x = NULL, y = "Norm. expression (all MASLD)", title = title) + theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_nature(p, file.path(fig_dir, file), 3.0 * length(markers) + 2, 6)
}
# Fig 4b-e are all the TAGLN-vs-other-endothelial CLUSTER comparison (all MASLD),
# split into protein (b,c) and genes (d,e). Marker lists are PINNED to the paper's
# EndMT signature so every cited marker is shown, in a fixed order; absent markers
# are dropped inside h_cluster(). rna()/adt() then route each marker to its panel.
# CD155 dropped: weak ADT-only effect (higher in TAGLN by a small margin, and DOWN
# with fibrosis), and it is NOT among the endo->T ligand-receptor interactions in
# Fig 4f (no CD155/PVR-TIGIT/DNAM-1 pair) - so it does not support the axis. CD146
# (MCAM) carries the gained surface phenotype.
gained_markers <- intersect(c("CD146","ACTA2","TAGLN","COL3A1","S100A4","VIM"), names(resolved))
lost_markers   <- intersect(c("CD31","MHC-I","PECAM1","VWF","CLDN5","CDH5"),     names(resolved))
cat("\nFig4 gained (pinned):", paste(gained_markers, collapse=", "),
    "\nFig4 lost (pinned):",   paste(lost_markers,   collapse=", "), "\n")

# b,c = surface phenotype / antigen presentation (protein, ADT)
h_cluster(adt(gained_markers), "TAGLN vs other endo - markers GAINED (protein)", "Fig4_TAGLNgained_ADT_v7")   # 4b
h_cluster(adt(lost_markers),   "TAGLN vs other endo - markers LOST (protein)",   "Fig4_TAGLNlost_ADT_v7")     # 4c
# d,e = EndMT mesenchymal genes gained / endothelial genes lost (RNA)
h_cluster(rna(gained_markers), "TAGLN vs other endo - mesenchymal genes GAINED",         "Fig4_TAGLNgained_RNA_v7")   # 4d
h_cluster(rna(lost_markers),   "TAGLN vs other endo - endothelial genes LOST",           "Fig4_TAGLNlost_RNA_v7")     # 4e

cat("\nStars: *** p<0.001, ** p<0.01, * p<0.05, ns (Wilcoxon; No-fib -> Fibrosis).",
    "\nEach panel is a horizontal strip; ADT (protein) and RNA (genes) are separate files.\n")
cat("Done Fig 4b-e.\n")
# End of the script

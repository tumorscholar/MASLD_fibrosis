###############################################################################
# 31 - Figure 5 (FNA) - PBMC object: SoupX -> doublet -> QC -> cluster,
#      with baseline Fig-1 cell types transferred on (split UMAP by visit).
#
# One PBMC object holding BOTH visits (4 baseline + 4 follow-up patients):
#   baseline 10738 / 11040 / 11183   PBMC = demux from their SVP capture (reuse SoupX)
#   baseline 10291-1                 PBMC = separate capture GC-WL-10291-1-PBMC (fresh SoupX)
#   follow-up 11596 / 11303 / 11816 / 11937   standalone PBMC (fresh SoupX)
#
# Cell types LABEL-TRANSFERRED from SeuObjV7 (PBMC cells) - same vocabulary as
# Fig 1a; split UMAP (baseline vs follow-up) shows which cell types change.
#
# Output: SeuObjFNA_PBMCV7.rds  + Fig5 PBMC split-UMAP / composition
# NOTE: the FINAL clinical-annotated Fig 5b PCA is 32_Figure5b_PCA_clinical.R (liver);
#       this script's own PCA (Fig5b_pbmc_PCA_by_visit) is the quick per-sample view.
###############################################################################

library(Seurat); library(SoupX); library(DropletUtils)
library(scDblFinder); library(SingleCellExperiment)
library(harmony); library(dplyr); library(ggplot2); library(ggrepel)
source("00_nature_style.R")

# ---- Fig-1a cell-type colours (name -> colour), rebuilt to match Fig 1 ------
fig1_celltype_cols <- function() {
 ann  <- read.csv("/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure1/cluster_annotation_V7.csv",
                  stringsAsFactors = FALSE)                       # cluster (C#), cell_type
 cl2t <- setNames(ann$cell_type, sub("^C", "", ann$cluster))
 lineage_order_ids <- c(0,2,1,6,20, 4, 8,19,29, 3,11,23, 7,9,26, 14,25, 13,18,
                        12,27,24,28, 5,10,15, 17,21, 16,22)
 my_colors <- c("#7f0000","#b30000","#d7301f","#ef6548","#fc8d59","#6a51a3",
                "#08519c","#3182bd","#6baed6","#00441b","#238b45","#66c2a4",
                "#41ab5d","#74c476","#a1d99b","#addd8e","#78c679","#fec44f",
                "#fe9929","#8c510a","#bf812d","#dfc27d","#f6e8c3","#01665e",
                "#35978f","#80cdc1","#d94801","#fd8d3c","#737373","#bdbdbd")
 nm <- unname(cl2t[as.character(lineage_order_ids)]); keep <- !duplicated(nm)
 setNames(my_colors[keep], nm[keep])                             # first colour per cell-type name
}
ct_cols <- fig1_celltype_cols()

# ---- CONFIG -----------------------------------------------------------------
data_dir   <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_R2_all50_cite"
soupx_main <- "/data/Blizard-AlazawiLab/rk/soupx_v7"
soupx_fna  <- "/data/Blizard-AlazawiLab/rk/soupxFNA_v7"
out_dir    <- "/data/Blizard-AlazawiLab/rk/seurat"
ref_rds    <- file.path(out_dir, "SeuObjV7.rds")
fig_dir    <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
out_rds    <- file.path(out_dir, "SeuObjFNA_PBMCV7.rds")

hto_tags <- c(SAT = "C0251", VAT = "C0252", PBMC = "C0253")
fu_map   <- c("10738" = "11596", "10291-1" = "11303", "11040" = "11816", "11183" = "11937")  # 10738 f/u capture = 11596 (matches the cap table below and scripts 34/35/40)

set.seed(1234)
for (d in c(soupx_fna, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# base = patient; visit; mode = "demux" (SVP -> keep PBMC) or "standalone"; cap = folder.
# All 4 follow-up PBMC were sequenced standalone (blood only - no adipose at follow-up).
# 11303's follow-up capture is named -SVP but carries only PBMC, so it is standalone;
# any residual hashtag features are stripped from ADTonly in the standalone branch.
pbmc_tbl <- data.frame(
 base  = c("10738","11040","11183","10291-1",   "10738","10291-1","11040","11183"),
 visit = c(rep("Baseline", 4),                   rep("Followup", 4)),
 mode  = c("demux","demux","demux","standalone", "standalone","standalone","standalone","standalone"),
 cap   = c("GC-WL-10738-SVP","GC-WL-11040-SVP","GC-WL-11183-SVP","GC-WL-10291-1-PBMC",
           "GC-WL-11596-PBMC","GC-WL-11303-SVP","GC-WL-11816-PBMC","GC-WL-11937-PBMC"),
 stringsAsFactors = FALSE)
all_caps <- basename(dirname(dirname(
 list.files(data_dir, "filtered_feature_bc_matrix.h5", recursive = TRUE, full.names = TRUE))))
miss <- setdiff(pbmc_tbl$cap, all_caps)
if (length(miss)) stop("captures not found in ", data_dir, ": ", paste(miss, collapse = ", "))
cat("PBMC captures:\n"); print(pbmc_tbl, row.names = FALSE)

# ---- helpers ----------------------------------------------------------------
read10x_gex <- function(outs_dir, which) {
 m <- Read10X_h5(file.path(outs_dir, paste0(which, "_feature_bc_matrix.h5")))
 if (is.list(m)) m[["Gene Expression"]] else m
}
run_soupx_one <- function(outs_dir) {
 tod <- read10x_gex(outs_dir, "raw"); toc <- read10x_gex(outs_dir, "filtered")
 tod <- tod[rownames(toc), , drop = FALSE]
 o <- CreateSeuratObject(toc)
 o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, verbose = FALSE)
 o <- ScaleData(o, verbose = FALSE);     o <- RunPCA(o, npcs = 20, verbose = FALSE)
 o <- FindNeighbors(o, dims = 1:20, verbose = FALSE); o <- FindClusters(o, resolution = 0.5, verbose = FALSE)
 sc <- SoupChannel(tod, toc); sc <- setClusters(sc, setNames(as.character(Idents(o)), colnames(o)))
 sc <- autoEstCont(sc); adjustCounts(sc, roundToInt = TRUE)
}
run_dbl <- function(obj) {
 sce  <- scDblFinder(SingleCellExperiment(list(counts = GetAssayData(obj, "RNA", layer = "counts"))),
                     BPPARAM = BiocParallel::SerialParam(RNGseed = 1234))
 meta <- data.frame(scDblFinder.class = as.character(sce$scDblFinder.class),
                    scDblFinder.score = sce$scDblFinder.score,
                    row.names = colnames(sce))[colnames(obj), , drop = FALSE]
 AddMetaData(obj, metadata = meta)
}
soupx_rna_dir <- function(cap) {
 reuse <- file.path(soupx_main, paste0(cap, "_strainedCounts")); if (dir.exists(reuse)) return(reuse)
 fresh <- file.path(soupx_fna,  paste0(cap, "_strainedCounts"))
 if (!dir.exists(fresh)) DropletUtils:::write10xCounts(fresh, run_soupx_one(file.path(data_dir, cap, "outs")), overwrite = TRUE)
 fresh
}
h5_path <- function(cap) file.path(data_dir, cap, "outs", "filtered_feature_bc_matrix.h5")

# ---- build one PBMC object per capture (demux SVP, or standalone) -----------
build_one <- function(cap, base, visit, mode) {
 message("PBMC: ", cap, " (", visit, ", ", mode, ")")
 obj <- CreateSeuratObject(counts = Seurat::Read10X(soupx_rna_dir(cap)), project = cap, assay = "RNA")
 obj$Sample <- cap
 adt  <- Read10X_h5(h5_path(cap))[["Antibody Capture"]]
 core <- function(b) sub("-\\d+$", "", sub("^.*_", "", b))
 idx  <- match(core(colnames(obj)), core(colnames(adt))); keep <- !is.na(idx)
 cat("  RNA:", ncol(obj), "| ADT:", ncol(adt), "| overlap:", sum(keep), "\n")
 obj <- subset(obj, cells = colnames(obj)[keep]); adt <- adt[, idx[keep], drop = FALSE]
 obj <- RenameCells(obj, new.names = colnames(adt))
 obj[["ADT"]] <- CreateAssayObject(counts = adt)
 obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
 obj$percent.hb <- PercentageFeatureSet(obj, features = c("HBA1","HBA2","HBB"))
 
 if (mode == "demux") {                                        # SVP capture -> keep PBMC
  ac <- GetAssayData(obj, assay = "ADT", layer = "counts")
  obj[["HTO"]]     <- CreateAssayObject(counts = ac[hto_tags, , drop = FALSE])
  obj[["ADTonly"]] <- CreateAssayObject(counts = ac[setdiff(rownames(ac), hto_tags), , drop = FALSE])
  # drop cells with no hashtag signal (they break HTODemux' k-means as an all-zero cluster)
  obj <- subset(obj, cells = colnames(obj)[colSums(as.matrix(ac[hto_tags, , drop = FALSE])) > 0])
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2)
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.95)
  obj <- subset(obj, subset = HTO_classification.global == "Singlet")
  obj <- subset(obj, cells = colnames(obj)[obj$HTO_classification == hto_tags[["PBMC"]]])
 } else {                                                      # standalone PBMC (no demux)
  ac <- GetAssayData(obj, assay = "ADT", layer = "counts")    # strip hashtag features if present
  obj[["ADTonly"]] <- CreateAssayObject(counts = ac[setdiff(rownames(ac), hto_tags), , drop = FALSE])
 }
 obj$Tissue <- "PBMC"; obj$Patient_ID <- base; obj$Visit <- visit
 run_dbl(obj)
}
objs <- Map(build_one, pbmc_tbl$cap, pbmc_tbl$base, pbmc_tbl$visit, pbmc_tbl$mode)

# ---- merge + QC (Fig-1 PBMC thresholds) -------------------------------------
SeuObj <- merge(objs[[1]], y = objs[-1], add.cell.ids = paste(pbmc_tbl$cap, pbmc_tbl$visit, sep = "_"))
SeuObj <- subset(SeuObj, subset =
                  scDblFinder.class == "singlet" & nFeature_RNA > 200 & percent.mt < 20)
SeuObj$Visit <- factor(SeuObj$Visit, levels = c("Baseline","Followup"))
cat("\nCells per (Patient_ID,Visit):\n"); print(table(SeuObj$Patient_ID, SeuObj$Visit))

# ---- normalise, Harmony, cluster, UMAP --------------------------------------
DefaultAssay(SeuObj) <- "ADTonly"; SeuObj <- NormalizeData(SeuObj, normalization.method = "CLR", margin = 2)
DefaultAssay(SeuObj) <- "RNA"
SeuObj <- NormalizeData(SeuObj) |> FindVariableFeatures() |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)

# NOTE: the quick per-sample PCA figure (Fig5b_pbmc_PCA_by_visit) is NOT plotted
# here - the FINAL Fig 5b PCA is 32_Figure5b_PCA_clinical.R, which reads the "pca"
# reduction created by RunPCA above. This script only builds the object.

SeuObj <- RunHarmony(SeuObj, group.by.vars = "Sample")
SeuObj <- FindNeighbors(SeuObj, reduction = "harmony", dims = 1:30)
SeuObj <- FindClusters(SeuObj, resolution = 0.5, random.seed = 1234)
SeuObj <- RunUMAP(SeuObj, reduction = "harmony", dims = 1:30, seed.use = 1234)

# ---- transfer baseline Fig-1 cell types (PBMC reference) --------------------
ref <- readRDS(ref_rds)
if (is.null(ref$cell_type)) stop("SeuObjV7 has no $cell_type")
ref <- subset(ref, subset = Tissue == "PBMC"); DefaultAssay(ref) <- "RNA"
ref <- NormalizeData(ref) |> FindVariableFeatures() |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)
anchors <- FindTransferAnchors(reference = ref, query = SeuObj, dims = 1:30, reference.reduction = "pca")
SeuObj$cell_type <- TransferData(anchorset = anchors, refdata = ref$cell_type, dims = 1:30)$predicted.id
cat("\nTransferred cell types:\n"); print(table(SeuObj$cell_type))

saveRDS(SeuObj, out_rds)

# ---- Fig 5: split UMAP (baseline vs follow-up) + composition ----------------
p_umap <- DimPlot(SeuObj, group.by = "cell_type", split.by = "Visit", reduction = "umap",
                  raster = FALSE, label = FALSE) +
  scale_color_manual(values = ct_cols, na.value = "grey85") +
  ggtitle("PBMC FNA - baseline vs follow-up (Fig-1 cell types)") + theme_nature(7) +
  theme(legend.position = "none")
lab <- p_umap$data |> dplyr::group_by(cell_type) |>
  dplyr::summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")
p_umap <- p_umap + geom_text_repel(
  data = lab, aes(umap_1, umap_2, label = cell_type), inherit.aes = FALSE,
  size = 2.0, colour = "black", bg.color = "white", bg.r = 0.12,
  box.padding = 0.4, point.padding = 0.3, segment.color = "grey40", segment.size = 0.25,
  segment.alpha = 0.6, max.overlaps = Inf, min.segment.length = 0.3, force = 3,
  force_pull = 0.5, seed = 42)
save_nature(p_umap, file.path(fig_dir, "Fig5_pbmc_splitUMAP_celltype_v7"), 22, 11)

# composition CSV only (full table) - the plotted composition panel is 33b
# (all-immune, cluster-labelled). The raw all-cell bar and the within-T bar are
# no longer plotted here; 33/33b read this CSV.
comp <- SeuObj@meta.data |> dplyr::count(Visit, cell_type) |>
 dplyr::group_by(Visit) |> dplyr::mutate(pct = 100 * n / sum(n)) |> dplyr::ungroup()
write.csv(comp, file.path(fig_dir, "Fig5_pbmc_composition_v7.csv"), row.names = FALSE)  # full table (input to 33b)

cat("\nDone Fig 5 PBMC. Saved:", out_rds, "\n")
# End of the script

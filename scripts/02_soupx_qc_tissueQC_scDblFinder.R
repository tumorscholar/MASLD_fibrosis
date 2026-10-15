###############################################################################
# Quality control - SoupX ambient-RNA removal, doublet detection, filtering
#
# SoupX ambient-RNA removal  ->  scDblFinder doublets  ->  simple QC.
#   * Object 1 (all cells, doublet-scored) saved as SeuObjV7_preQC.rds
#   * Object 2 (QC-filtered main object)   saved as SeuObjV7.rds
#
# QC: keep singlets, nFeature_RNA > 200, and percent.mt < 30 (liver) / < 20
#     (PBMC/SAT/VAT).
###############################################################################

library(Seurat)
library(SoupX)
library(DropletUtils)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)

# ---- CONFIG -----------------------------------------------------------------
data_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_R2_all50_cite"
soupx_dir    <- "/data/Blizard-AlazawiLab/rk/soupx_v7"
qc_plot_dir  <- "/data/home/hdx044/seurat/QCplots_v7"
out_dir      <- "/data/Blizard-AlazawiLab/rk/seurat"

liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-1","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 21 liver (incl new 10291-1)
svp_pids   <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 20 SAT-VAT-PBMC

soupx_skip <- c("GC-WL-9991-LIVER")                # no raw matrix -> use raw counts
hto_tags   <- c(SAT = "C0251", VAT = "C0252", PBMC = "C0253")
run_scdblfinder <- TRUE

set.seed(1234)
dir.create(soupx_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)

# ---- Select the 41 main captures (21 liver + 20 SVP) ------------------------
all_h5    <- list.files(data_dir, pattern = "filtered_feature_bc_matrix.h5",
                        recursive = TRUE, full.names = TRUE)
all_names <- basename(dirname(dirname(all_h5)))
tissue_of <- ifelse(grepl("-LIVER$", all_names), "LIVER",
                    ifelse(grepl("-(SAT-VAT-PBMC|SVP)$", all_names), "SVP", "OTHER"))
pid_of    <- sub("^GC-WL-(.*)-(LIVER|SAT-VAT-PBMC|SVP)$", "\\1", all_names)
keep      <- (tissue_of == "LIVER" & pid_of %in% liver_pids) |
  (tissue_of == "SVP"   & pid_of %in% svp_pids)

files     <- all_h5[keep]
samples   <- all_names[keep]
soup_path <- gsub("/filtered_feature_bc_matrix.h5", "", files)
names(files) <- samples; names(soup_path) <- samples
cat("Selected", length(samples), "captures\n"); print(samples)

# ---- 1. SoupX ---------------------------------------------------------------
# Build the SoupChannel by hand (load10X is unreliable on CITE-seq data).
# These /outs have no analysis/ folder, so clusters are computed internally.
read10x_gex <- function(outs_dir, which) {
  m <- Read10X_h5(file.path(outs_dir, paste0(which, "_feature_bc_matrix.h5")))
  if (is.list(m)) m[["Gene Expression"]] else m
}
run_soupx_one <- function(outs_dir) {
  tod <- read10x_gex(outs_dir, "raw")
  toc <- read10x_gex(outs_dir, "filtered")
  tod <- tod[rownames(toc), , drop = FALSE]
  o <- CreateSeuratObject(toc)
  o <- NormalizeData(o, verbose = FALSE); o <- FindVariableFeatures(o, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE);     o <- RunPCA(o, npcs = 20, verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:20, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.5, verbose = FALSE)
  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc, setNames(as.character(Idents(o)), colnames(o)))
  sc <- autoEstCont(sc)
  adjustCounts(sc, roundToInt = TRUE)
}
setwd(soupx_dir)
for (x in samples) {
  if (x %in% soupx_skip) next
  out <- run_soupx_one(soup_path[x])
  DropletUtils:::write10xCounts(file.path(soupx_dir, paste0(x, "_strainedCounts")),
                                out, overwrite = TRUE)
}

# ---- Load corrected RNA (SoupX) + raw ADT -----------------------------------
soupx_files   <- list.files(soupx_dir, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
soupx_samples <- gsub(paste0(soupx_dir, "/?"), "", soupx_files)
soupx_samples <- gsub("_strainedCounts/matrix.mtx", "", soupx_samples)
soupx_files   <- gsub("/matrix.mtx", "", soupx_files)
names(soupx_files) <- soupx_samples

rawDat <- lapply(files, Read10X_h5); names(rawDat) <- samples

# ---- 2. Build objects, assign Tissue, run scDblFinder -----------------------
run_dbl <- function(obj) {
  if (!run_scdblfinder) {
    return(AddMetaData(obj, data.frame(scDblFinder.class = "singlet",
                                       scDblFinder.score = NA_real_,
                                       row.names = colnames(obj))))
  }
  # reproducible doublet calls: seed via a serial BiocParallel backend
  # (set.seed alone does NOT control scDblFinder's parallel RNG streams).
  sce  <- scDblFinder(SingleCellExperiment(list(counts = GetAssayData(obj, "RNA", layer = "counts"))),
                      BPPARAM = BiocParallel::SerialParam(RNGseed = 1234))
  meta <- data.frame(scDblFinder.class = as.character(sce$scDblFinder.class),
                     scDblFinder.score = sce$scDblFinder.score,
                     row.names = colnames(sce))[colnames(obj), , drop = FALSE]
  AddMetaData(obj, metadata = meta)
}

datObjs <- list()
for (x in samples) {
  message("Processing: ", x)
  rna <- if (x %in% soupx_skip) rawDat[[x]][["Gene Expression"]] else Seurat::Read10X(soupx_files[x])
  obj <- CreateSeuratObject(counts = rna, project = x, assay = "RNA")
  obj$Sample <- x
  
  # align RNA <-> ADT barcodes by the 16-bp barcode CORE, tolerating any
  # suffix/prefix difference (e.g. "-1" present on one but not the other, as in
  # some SoupX outputs). Match on core, then rename ADT columns to the RNA cells.
  adt  <- rawDat[[x]][["Antibody Capture"]]
  core <- function(b) sub("-\\d+$", "", sub("^.*_", "", b))
  idx  <- match(core(colnames(obj)), core(colnames(adt)))   # ADT column per RNA cell
  keep <- !is.na(idx)
  cat("  RNA:", ncol(obj), "| ADT:", ncol(adt), "| overlap:", sum(keep), "\n")
  if (sum(keep) == 0) stop("No overlapping RNA/ADT barcodes for ", x)
  obj <- subset(obj, cells = colnames(obj)[keep])
  adt <- adt[, idx[keep], drop = FALSE]
  obj <- RenameCells(obj, new.names = colnames(adt))        # standardise to raw-h5 barcodes ("-1")
  obj[["ADT"]] <- CreateAssayObject(counts = adt)           # so TCR matching stays consistent later
  
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.hb <- PercentageFeatureSet(obj, features = c("HBA1","HBA2","HBB"))
  
  if (grepl("(SAT-VAT-PBMC|SVP)$", x)) {                 # multiplexed: demux to tissue
    ac <- GetAssayData(obj, assay = "ADT", layer = "counts")
    obj[["HTO"]]     <- CreateAssayObject(counts = ac[hto_tags, , drop = FALSE])
    obj[["ADTonly"]] <- CreateAssayObject(counts = ac[setdiff(rownames(ac), hto_tags), , drop = FALSE])
    obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2)
    obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.95)
    obj <- subset(obj, subset = HTO_classification.global == "Singlet")
    tissue <- rep(NA_character_, ncol(obj))
    for (tt in names(hto_tags)) tissue[obj$HTO_classification == hto_tags[[tt]]] <- tt
    obj$Tissue     <- tissue
    obj$Patient_ID <- gsub("GC-WL-|-(SAT-VAT-PBMC|SVP)$", "", x)
  } else {                                                # liver: 130 proteins, no tags
    obj[["ADTonly"]] <- CreateAssayObject(counts = GetAssayData(obj, assay = "ADT", layer = "counts"))
    obj$Tissue     <- "LIVER"
    obj$Patient_ID <- gsub("GC-WL-|-LIVER$", "", x)
  }
  
  obj <- run_dbl(obj)
  datObjs[[x]] <- obj
}
cat("Objects created:", length(datObjs), "\n")

# ---- 3. Merge + pre-QC object + diagnostics ---------------------------------
# add.cell.ids prefixes each cell with its capture name -> "<Sample>_<barcode>",
# so TCR barcodes ("<capture>_<bc>") line up 1:1 downstream (Fig 2a/b, builder).
SeuObj <- merge(datObjs[[1]], y = datObjs[-1], add.cell.ids = names(datObjs))
SeuObj$Tissue <- factor(SeuObj$Tissue, levels = c("LIVER","PBMC","SAT","VAT"))

setwd(qc_plot_dir)
ggsave("preQC_byTissue.tiff",
       VlnPlot(SeuObj, group.by = "Tissue",
               features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hb"),
               pt.size = 0, ncol = 4),
       width = 14, height = 5, dpi = 300, compression = "lzw")
ggsave("preQC_mt_histograms.tiff",
       ggplot(SeuObj@meta.data, aes(percent.mt)) + geom_histogram(bins = 60) +
         facet_wrap(~Tissue, scales = "free_y") + theme_bw() +
         labs(title = "percent.mt by tissue"),
       width = 9, height = 6, dpi = 300, compression = "lzw")

cat("\nCells per tissue BEFORE QC:\n"); print(table(SeuObj$Tissue))
cat("\nscDblFinder calls:\n");          print(table(SeuObj$scDblFinder.class))
saveRDS(SeuObj, file.path(out_dir, "SeuObjV7_preQC.rds"))

# ---- 4. QC -> main object ---------------------------------------------------
SeuObj <- subset(
  SeuObj,
  subset =
    scDblFinder.class == "singlet" &
    nFeature_RNA > 200 &
    ((Tissue == "LIVER" & percent.mt < 30) |
       (Tissue %in% c("PBMC","SAT","VAT") & percent.mt < 20))
)

cat("\nCells per tissue AFTER QC:\n"); print(table(SeuObj$Tissue))
cat("\nCells per sample AFTER QC:\n"); print(table(SeuObj$Sample))
saveRDS(SeuObj, file.path(out_dir, "SeuObjV7.rds"))

# End of the script
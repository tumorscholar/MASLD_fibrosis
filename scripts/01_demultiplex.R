###############################################################################
# Demultiplexing of multiplexed captures (pre-QC object for TCR analysis)
#
# OBJECT 1 : pre-QC object for TCR / clonotype analysis.
#
# Purpose:
#  - Demultiplex the multiplexed SAT-VAT-PBMC captures (HTODemux) and assign a
#    Tissue label (SAT / VAT / PBMC) to every cell; liver captures are labelled
#    directly.
#  - NO cell-level QC filtering is applied here, so that clonotype-bearing cells
#    are NOT lost before TCR analysis.
#  - Export per-sample barcode->Tissue tables used to demultiplex the TCR contigs.
#
# Data source: the CITE-seq captures; the main liver and SAT-VAT-PBMC captures
#              are selected for analysis.
#
# Output (objects do NOT overwrite previous ones):
#  - SeuObjForDemuxV7.rds        (all cells, tissue-labelled)
#  - <barcode_dir>/<sample>_metadata.csv   (cell_barcode, sample, Tissue)
###############################################################################

library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)

# =============================================================================
# CONFIG  -- edit these for the re-run
# =============================================================================
data_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_R2_all50_cite"
out_dir      <- "/data/Blizard-AlazawiLab/rk/seurat"
barcode_dir  <- "/data/Blizard-AlazawiLab/rk/barcodes_v7"
preqc_rds    <- "SeuObjForDemuxV7.rds"

liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-1","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 21 liver (incl new 10291-1)
svp_pids   <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 20 SAT-VAT-PBMC

hto_tags   <- c(SAT = "C0251", VAT = "C0252", PBMC = "C0253")

dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(barcode_dir, showWarnings = FALSE, recursive = TRUE)
# =============================================================================

# ---- Locate & select the 41 main captures ----------------------------------
all_h5    <- list.files(data_dir, pattern = "filtered_feature_bc_matrix.h5",
                        recursive = TRUE, full.names = TRUE)
all_names <- basename(dirname(dirname(all_h5)))
tissue_of <- ifelse(grepl("-LIVER$", all_names), "LIVER",
              ifelse(grepl("-(SAT-VAT-PBMC|SVP)$", all_names), "SVP", "OTHER"))
pid_of    <- sub("^GC-WL-(.*)-(LIVER|SAT-VAT-PBMC|SVP)$", "\\1", all_names)
keep      <- (tissue_of == "LIVER" & pid_of %in% liver_pids) |
             (tissue_of == "SVP"   & pid_of %in% svp_pids)

files   <- all_h5[keep]
samples <- all_names[keep]
names(files) <- samples

miss_liver <- setdiff(liver_pids, pid_of[keep][tissue_of[keep]=="LIVER"])
miss_svp   <- setdiff(svp_pids,   pid_of[keep][tissue_of[keep]=="SVP"])
cat("Selected", length(samples), "captures\n")
if (length(miss_liver) | length(miss_svp))
  stop("Missing -> liver: ", paste(miss_liver, collapse=","),
       " | SVP: ", paste(miss_svp, collapse=","))

# ---- Load all data ---------------------------------------------------------
rawDat <- lapply(files, Read10X_h5)
names(rawDat) <- samples

# ---- Build objects, demultiplex, assign Tissue (NO QC) ---------------------
datObjs <- list()
for (x in samples) {
  obj <- CreateSeuratObject(counts = rawDat[[x]][["Gene Expression"]], project = x)
  obj[["ADT"]] <- CreateAssayObject(counts = rawDat[[x]][["Antibody Capture"]])
  obj$sample   <- x

  if (grepl("-(SAT-VAT-PBMC|SVP)$", x)) {
    adt <- obj@assays$ADT$counts
    obj[["HTO"]]     <- CreateAssayObject(counts = adt[hto_tags, ])
    obj[["ADTonly"]] <- CreateAssayObject(counts = adt[setdiff(rownames(adt), hto_tags), ])
    obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2)
    obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.95)
    obj <- subset(obj, subset = HTO_classification.global == "Singlet")
    tissue <- rep(NA_character_, ncol(obj))
    for (tt in names(hto_tags)) tissue[obj$HTO_classification == hto_tags[tt]] <- tt
    obj$Tissue     <- tissue
    obj$Patient_ID <- gsub("GC-WL-|-SAT-VAT-PBMC|-SVP", "", x)
  } else {
    obj[["ADTonly"]] <- CreateAssayObject(counts = obj@assays$ADT$counts)  # liver: 130 proteins, no tags
    obj$Tissue     <- "LIVER"
    obj$Patient_ID <- gsub("GC-WL-|-LIVER", "", x)
  }
  datObjs[[x]] <- obj
}

# ---- Merge -> pre-QC object ------------------------------------------------
SeuObjForDemux <- merge(x = datObjs[[1]], y = datObjs[-1])
saveRDS(SeuObjForDemux, file.path(out_dir, preqc_rds))   # -> SeuObjForDemuxV7.rds

# ---- Export per-sample barcode -> Tissue tables (for TCR demux) ------------

SeuObjForDemux <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjForDemuxV7.rds')

meta <- SeuObjForDemux@meta.data
meta$cell_barcode <- sub("_[0-9]+$", "", rownames(meta))  # strip merge suffix _<n>
for (s in unique(meta$sample)) {
  sub_df <- meta[meta$sample == s, c("cell_barcode", "sample", "Tissue")]
  write.csv(sub_df, file.path(barcode_dir, paste0(s, "_metadata.csv")), row.names = FALSE)
}
cat("Wrote barcode tables for", length(unique(meta$sample)), "captures to", barcode_dir, "\n")

# End of the script

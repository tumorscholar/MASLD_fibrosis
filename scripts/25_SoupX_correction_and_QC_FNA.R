###############################################################################
# Script name: SoupX_correction_and_QC_FNA.R
#
# Purpose:
#   - Perform SoupX ambient RNA contamination correction on FNA scRNA-seq data
#   - Reconstruct SoupX-corrected 10X matrices for each sample
#   - Load corrected RNA data together with raw ADT data
#   - Build Seurat objects and compute standard QC metrics
#   - Generate QC plots before and after filtering
#   - Save final Seurat objects for downstream analysis
#
# Input data:
#   - CellRanger outputs:
#       * filtered_feature_bc_matrix.h5
#       * Expected under:
#           /data/Blizard-AlazawiLab/rk/cellrangerFNA/
#           /data/Blizard-AlazawiLab/rk/cellranger/
#
# Main processing steps:
#   1. Identify CellRanger output directories and samples
#   2. Run SoupX (autoEstCont + adjustCounts) per sample
#   3. Write SoupX-corrected counts in 10X format
#   4. Load SoupX-corrected RNA and raw H5 (RNA + ADT)
#   5. Create Seurat objects with matched RNA and ADT barcodes
#   6. Calculate QC metrics (nFeature_RNA, nCount_RNA, mt%, hb%)
#   7. Generate QC plots before and after filtering
#   8. Filter cells based on QC thresholds
#
# Output files:
#   - SoupX-corrected 10X directories:
#       /data/Blizard-AlazawiLab/rk/soupxFNA/*_strainedCounts/
#   - QC plots (TIFF):
#       *_qcPlot.tiff
#       *_lmPlot.tiff
#       *_afterQCplot.tiff
#   - R object:
#       SeuObjFNA.rds
#
###############################################################################

# Load packages
library(Seurat)
library(hdf5r)
library(ggplot2)
library(SoupX)
library(DropletUtils)

# Locate files 
dataDir <- '/data/Blizard-AlazawiLab/rk/cellrangerFNA'
files <- list.files(path = dataDir,
                    pattern = 'filtered_feature_bc_matrix.h5',
                    recursive = T,
                    full.names = T)

# Extract sample names from file names
samples <- gsub(
  "/outs/filtered_feature_bc_matrix.h5$",
  "",
  gsub(paste0("^", dataDir, "/"), "", files)
)

# Create output directory if it does not exist
outDir  <- "/data/Blizard-AlazawiLab/rk/soupxFNA"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#Paths passed to SoupX (directory level)
soup_path <- gsub("/filtered_feature_bc_matrix.h5$", "", files)

#Run SoupX per sample

for (x in seq_along(soup_path)) {
  
  message("Processing sample: ", samples[x])
  
  # Load 10X data
  sc <- load10X(
    soup_path[x],
    includeFeatures = "Gene Expression"
  )
  
  # Estimate ambient RNA contamination
  sc <- autoEstCont(sc)
  
  # Adjust counts
  out <- adjustCounts(sc, roundToInt = TRUE)
  
  # Write SoupX-corrected counts (10x format)
  write10xCounts(
    path = file.path(outDir, paste0(samples[x], "_strainedCounts")),
    x = out
  )
}

# Locate SoupX output 
#Directories
dataDirSoupx <- "/data/Blizard-AlazawiLab/rk/soupxFNA"
dataDirFNA   <- "/data/Blizard-AlazawiLab/rk/cellrangerFNA"
dataDirCR    <- "/data/Blizard-AlazawiLab/rk/cellranger"

setwd("~/seurat/QCplots/FNA")
dir.create(".", showWarnings = FALSE, recursive = TRUE)

#Locate SoupX output (matrix.mtx)
matrix_files <- list.files(
  path       = dataDirSoupx,
  pattern    = "matrix.mtx",
  recursive  = TRUE,
  full.names = TRUE
)

# Folders containing SoupX-corrected data
soupx_dirs <- dirname(matrix_files)

# Sample names (derived from folder names)
sample_names_soupx <- basename(soupx_dirs)
sample_names_soupx <- sub("_strainedCounts$", "", sample_names_soupx)

cat("SoupX samples detected:\n")
print(sample_names_soupx)

# Build matching raw H5 paths
h5_files <- sapply(sample_names_soupx, function(s) {
  
  path_FNA <- file.path(
    dataDirFNA, s, "outs", "filtered_feature_bc_matrix.h5"
  )
  
  path_other <- file.path(
    dataDirCR, s, "outs", "filtered_feature_bc_matrix.h5"
  )
  
  if (file.exists(path_FNA))   return(path_FNA)
  if (file.exists(path_other)) return(path_other)
  NA_character_
})

# Verify existence
check_df <- data.frame(
  sample = sample_names_soupx,
  h5_path = h5_files,
  found = !is.na(h5_files),
  stringsAsFactors = FALSE
)

cat("\nH5 file check:\n")
print(check_df)

if (any(is.na(h5_files))) {
  stop(
    "Missing H5 files for: ",
    paste(sample_names_soupx[is.na(h5_files)], collapse = ", ")
  )
}

# Load data
cat("\nLoading SoupX-corrected RNA matrices...\n")
rawDat_after_soupx <- lapply(soupx_dirs, Seurat::Read10X)
names(rawDat_after_soupx) <- sample_names_soupx

cat("Loading raw H5 data (RNA + ADT)...\n")
rawDat <- lapply(h5_files, Read10X_h5)
names(rawDat) <- sample_names_soupx

# Build Seurat objects + QC
datObjs <- list()

for (x in sample_names_soupx) {
  
  cat("\n====================================\n")
  cat("Processing sample:", x, "\n")
  
  ## Create Seurat object (RNA from SoupX)
  datObjs[[x]] <- CreateSeuratObject(
    counts  = rawDat_after_soupx[[x]],
    project = x,
    assay   = "RNA"
  )
  
  # Match RNA and ADT barcodes
  rna_barcodes <- colnames(datObjs[[x]])
  adt_counts   <- rawDat[[x]][["Antibody Capture"]]
  
  shared_barcodes <- intersect(
    rna_barcodes,
    colnames(adt_counts)
  )
  
  cat(
    "RNA:", length(rna_barcodes),
    "| ADT:", ncol(adt_counts),
    "| Shared:", length(shared_barcodes), "\n"
  )
  
  # Subset to shared barcodes
  datObjs[[x]] <- subset(datObjs[[x]], cells = shared_barcodes)
  datObjs[[x]][["ADT"]] <- CreateAssayObject(
    counts = adt_counts[, shared_barcodes]
  )
  
  # Metadata
  datObjs[[x]]$Sample <- x
  
  # QC metrics
  datObjs[[x]][["percent.mt"]] <- PercentageFeatureSet(
    datObjs[[x]], pattern = "^MT-"
  )
  
  datObjs[[x]][["percent.hb"]] <- PercentageFeatureSet(
    datObjs[[x]],
    features = c("HBA1", "HBA2", "HBB")
  )
  
  # QC plots (before filtering)
  p1 <- VlnPlot(
    datObjs[[x]],
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.hb"
    ),
    ncol = 4
  )
  
  ggsave(
    plot    = p1,
    filename = paste0(x, "_qcPlot.tiff"),
    width   = 8,
    height  = 8,
    dpi     = 300,
    units   = "in"
  )
  
  p2 <- FeatureScatter(
    datObjs[[x]],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ) + geom_smooth(method = "lm")
  
  ggsave(
    plot    = p2,
    filename = paste0(x, "_lmPlot.tiff"),
    width   = 8,
    height  = 8,
    dpi     = 300,
    units   = "in"
  )
  
  ## Filtering
  datObjs[[x]] <- subset(
    datObjs[[x]],
    subset =
      nFeature_RNA > 200 &
      nFeature_RNA < 5000 &
      percent.mt < 10
  )
  
  # QC plots (after filtering)
  p3 <- VlnPlot(
    datObjs[[x]],
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.hb"
    ),
    ncol = 4
  )
  
  ggsave(
    plot    = p3,
    filename = paste0(x, "_afterQCplot.tiff"),
    width   = 8,
    height  = 8,
    dpi     = 300,
    units   = "in"
  )
  
  cat("Completed:", x, "\n")
}


setwd("~/seurat")
saveRDS(datObjs, 'SeuObjFNA.rds')

# End of the script

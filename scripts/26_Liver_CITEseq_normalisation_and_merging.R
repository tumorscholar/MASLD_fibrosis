###############################################################################
# Script name: Liver_CITEseq_normalisation_and_merging.R
#
# Purpose:
#   - Prepare liver CITEseq data for downstream analysis
#   - Add liver-specific metadata (tissue and patient ID)
#   - Perform RNA and ADT normalisation, variable feature selection, and scaling
#   - Integrate liver samples with matched control tissues processed
#     under the same experimental workflow
#
# Input:
#   - Seurat object list generated after QC and SoupX correction:
#       ~/seurat/SeuObjFNA.rds
#
# Main processing steps:
#   1. Subset liver samples from the Seurat object list
#   2. Annotate liver cells with tissue type and patient identifiers
#   3. Perform RNA normalisation, variable feature selection, and scaling
#   4. Normalise ADT signals (CLR) where available
#   5. Harmonise liver samples with other tissue datasets processed in parallel
#   6. Merge all processed objects into a single Seurat object
#
# Output:
#   - Combined Seurat object containing liver data (with matched reference tissues):
#       SeuObjFNAx.rds
#
#
###############################################################################


library(Seurat)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("~/seurat")

# Load data list
mylist <- readRDS("~/seurat/SeuObjFNA.rds")

# Split list into LIVER and SVP objects
names(mylist)

# LIVER
LIVER_idx <- grep("-LIVER$", names(mylist))
LIVERlist <- mylist[LIVER_idx]

# PBMC (pure PBMC only; exclude SAT-VAT)
PBMC_idx <- grep("-PBMC$", names(mylist))               # ends with -PBMC
PBMC_idx <- PBMC_idx[!grepl("SAT-VAT", names(mylist)[PBMC_idx], ignore.case = TRUE)]
PBMClist <- mylist[PBMC_idx]

# SVP (any with SAT-VAT in the name)
SVP_idx <- grep("SAT-VAT", names(mylist), ignore.case = TRUE)
SVPlist <- mylist[SVP_idx]

# --- Optional: sanity checks to ensure no overlap ---
stopifnot(length(intersect(names(LIVERlist), names(PBMClist))) == 0)
stopifnot(length(intersect(names(LIVERlist), names(SVPlist))) == 0)
stopifnot(length(intersect(names(PBMClist), names(SVPlist))) == 0)

# Show final selections
names(LIVERlist)
names(PBMClist)
names(SVPlist)

#### Add metadata, normalise, FindVariableFeatures and scale LIVER data ####

for (x in names(LIVERlist)) {
 
 # --- Metadata ---
 LIVERlist[[x]]$Tissue <- "LIVER"
 
 id <- gsub("-LIVER", "", x)
 id <- gsub("GC-WL-", "", id)
 LIVERlist[[x]]$Patient_ID <- id
 
 # --- ADT normalisation ---
 if ("ADT" %in% names(LIVERlist[[x]]@assays)) {
  DefaultAssay(LIVERlist[[x]]) <- "ADT"
  
  # ADTonly contains only the antibody features (1:130)
  LIVERlist[[x]][["ADTonly"]] <- CreateAssayObject(
   counts = LIVERlist[[x]]@assays$ADT$counts[1:130, ]
  )
  
  LIVERlist[[x]] <- NormalizeData(
   LIVERlist[[x]],
   assay = "ADTonly",
   normalization.method = "CLR",
   margin = 2
  )
 }
 
 # RNA normalisation
 DefaultAssay(LIVERlist[[x]]) <- "RNA"
 LIVERlist[[x]] <- NormalizeData(LIVERlist[[x]])
 LIVERlist[[x]] <- FindVariableFeatures(LIVERlist[[x]])
 LIVERlist[[x]] <- ScaleData(LIVERlist[[x]])
}

#### Add metadata, normalise, FindVariableFeatures and scale PBMC data ####

for (x in names(PBMClist)) {
 
 # Metadata
 PBMClist[[x]]$Tissue <- "PBMC"
 
 id <- gsub("-PBMC", "", x)
 id <- gsub("GC-WL-", "", id)
 PBMClist[[x]]$Patient_ID <- id
 
 # ADT if available (optional)
 if ("ADT" %in% names(PBMClist[[x]]@assays)) {
  DefaultAssay(PBMClist[[x]]) <- "ADT"
  
  PBMClist[[x]][["ADTonly"]] <- CreateAssayObject(
   counts = PBMClist[[x]]@assays$ADT$counts[1:130, ]
  )
  
  PBMClist[[x]] <- NormalizeData(
   PBMClist[[x]], assay = "ADTonly",
   normalization.method = "CLR", margin = 2
  )
 }
 
 # RNA normalisation
 DefaultAssay(PBMClist[[x]]) <- "RNA"
 PBMClist[[x]] <- NormalizeData(PBMClist[[x]])
 PBMClist[[x]] <- FindVariableFeatures(PBMClist[[x]])
 PBMClist[[x]] <- ScaleData(PBMClist[[x]])
}

#### Add metadata, normalise, FindVariableFeatures, scale and demultiplex SVP data ####

for (x in names(SVPlist)) {
 
 # Create HTO and ADTonly assays
 DefaultAssay(SVPlist[[x]]) <- "ADT"
 
 SVPlist[[x]][["HTO"]] <- CreateAssayObject(
  counts = SVPlist[[x]]@assays$ADT$counts[c("C0251", "C0252", "C0253"), ]
 )
 
 SVPlist[[x]][["ADTonly"]] <- CreateAssayObject(
  counts = SVPlist[[x]]@assays$ADT$counts[1:130, ]
 )
 
 # CLR normalization
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "HTO",
                               normalization.method = "CLR", margin = 2)
 
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "ADTonly",
                               normalization.method = "CLR", margin = 2)
 
 # HTO demultiplexing
 SVPlist[[x]] <- HTODemux(
  SVPlist[[x]],
  assay = "HTO",
  positive.quantile = 0.95
 )
 
 # Keep only singlets
 SVPlist[[x]] <- subset(
  SVPlist[[x]],
  subset = HTO_classification.global == "Singlet"
 )
 
 # Tissue metadata from HTO
 tissue_vector <- rep(NA, ncol(SVPlist[[x]]))
 
 tissue_vector[SVPlist[[x]]$HTO_classification == "C0251"] <- "SAT"
 tissue_vector[SVPlist[[x]]$HTO_classification == "C0252"] <- "VAT"
 tissue_vector[SVPlist[[x]]$HTO_classification == "C0253"] <- "PBMC"
 
 SVPlist[[x]]$Tissue <- tissue_vector
 
 # Patient ID
 id <- gsub("-SAT-VAT-PBMC", "", x)
 id <- gsub("GC-WL-", "", id)
 SVPlist[[x]]$Patient_ID <- id
 
 # Normalize RNA
 DefaultAssay(SVPlist[[x]]) <- "RNA"
 SVPlist[[x]] <- NormalizeData(SVPlist[[x]])
 SVPlist[[x]] <- FindVariableFeatures(SVPlist[[x]])
 SVPlist[[x]] <- ScaleData(SVPlist[[x]])
}

# Merge the liver and SAT, VAT and PBMC list
all.tissue.list <- c(LIVERlist, PBMClist, SVPlist)

SeuObj <- merge(x = all.tissue.list[[1]], 
                y = all.tissue.list[2:16])

saveRDS(SeuObj, 'SeuObjFNAx.rds')

# End of the script

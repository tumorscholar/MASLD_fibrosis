###############################################################################
# Script name: 02_soupx_qc_and_normalisation.R
#
# Description:
# This script performs ambient RNA removal, quality control (QC), and initial
# preprocessing of multi‑tissue single‑cell RNA‑seq data generated using
# 10x Genomics technology.
#
# The primary purpose of this script is to remove ambient RNA contamination
# using SoupX, apply cell‑level QC filtering, and generate cleaned Seurat
# objects for downstream analysis. This script operates on Cell Ranger output
# files (filtered_feature_bc_matrix.h5) and integrates both gene expression
# and antibody‑derived tag (ADT) data.
#
# Ambient RNA contamination is estimated on a per‑sample basis using SoupX,
# followed by correction of gene expression counts. Corrected count matrices
# are then used to construct Seurat objects, while ADT data are loaded from
# the original Cell Ranger outputs.
#
# Key steps performed in this script include:
#  - SoupX‑based ambient RNA estimation and correction
#  - Construction of Seurat objects from corrected RNA counts
#  - Integration of antibody capture (ADT) data
#  - Calculation of QC metrics (nFeature_RNA, nCount_RNA, mitochondrial and
#    hemoglobin gene percentages)
#  - Visualization of QC metrics using violin and scatter plots
#  - Filtering of low‑quality cells based on gene complexity and mitochondrial
#    content
#  - Assignment of tissue and patient metadata for liver and multiplexed
#    SAT‑VAT‑PBMC samples
#  - Normalization, variable feature selection, and scaling of RNA and ADT data
#
# Hashtag oligonucleotide (HTO) data are used to demultiplex multiplexed
# SAT‑VAT‑PBMC samples, retaining only singlet cells and assigning tissue
# identity at the cell level.
#
# Outputs:
#  - QC‑filtered and SoupX‑corrected Seurat objects
#  - Diagnostic QC plots for each sample
#
# Purpose:
# To generate clean, QC‑filtered, and ambient RNA‑corrected Seurat objects
# suitable for downstream clustering, cell type annotation, and integrative
# single‑cell analyses.
#

###############################################################################

# Load packages
library(Seurat)
library(hdf5r)
library(ggplot2)
library(SoupX)
library(DropletUtils)

# Set working directory
setwd("/data/Blizard-AlazawiLab/rk/soupx")

# Locate files 
dataDir <- '/data/Blizard-AlazawiLab/rk/cellranger'
files <- list.files(path = dataDir,
                    pattern = 'filtered_feature_bc_matrix.h5',
                    recursive = T,
                    full.names = T)

# Extract sample names from file names
pattern1 <- '/data/Blizard-AlazawiLab/rk/cellranger/'
samples <- gsub(pattern1,
                '',
                files)

pattern2 <- '/outs/filtered_feature_bc_matrix.h5'
samples <- gsub(pattern2,
                '',
                samples)

pattern3 <- '/filtered_feature_bc_matrix.h5'

# SoupX input file list
soup_path <- gsub(pattern3, '', files)

# Load all data to SoupX except 9991-LIVER
sample_indexes = c(1:length(soup_path))
remove_sample <- grep("GC-WL-9991-LIVER", soup_path) 
sample_indexes <- sample_indexes[-remove_sample]

for (x in sample_indexes) {
 sc = load10X(soup_path[x], includeFeatures = c("Gene Expression"))
 
 # Estimate rho
 sc = autoEstCont(sc)
 
 # Clean the data
 out = adjustCounts(sc, roundToInt = TRUE)
 
 # Save file
 DropletUtils:::write10xCounts(paste0("./", samples[x], "_strainedCounts"), out)
} 

# Locate files of SoupX output
dataDirSoupx <- '/data/Blizard-AlazawiLab/rk/soupx'
files_after_soupx <- list.files(path = dataDirSoupx,
                    pattern = 'matrix.mtx',
                    recursive = T,
                    full.names = T)

# Extract sample names from soupx output
pattern4 <- '/data/Blizard-AlazawiLab/rk/soupx/'
samples_after_soupx <- gsub(pattern4,
                            '',
                            files_after_soupx)

pattern5 <- '_strainedCounts/matrix.mtx'
samples_after_soupx <- gsub(pattern5,
                            '',
                samples_after_soupx)

# Remove /matrix.mtx from end of files_after_soupx
files_after_soupx <- gsub('/matrix.mtx', '', files_after_soupx)

# Load all data from soupx output for RNA
rawDat_after_soupx <- lapply(files_after_soupx, function(x){
 Seurat::Read10X(x)
})

# Specify list element names 
names(rawDat_after_soupx) <- samples_after_soupx

#Load all data without soupx for ADT
rawDat <- lapply(files, function(x){
 Read10X_h5(x)
}) 

# Specify list element names 
names(rawDat) <- samples

# Set working directory
setwd("~/seurat/QCplots")

# Make seurat objects for RNA in for loop
datObjs <- list() 
for (x in samples_after_soupx) {
 
 # Create seurat object for RNA
 datObjs[[x]] <- Seurat::CreateSeuratObject(counts = rawDat_after_soupx[[x]], project = x, assay = 'RNA')
 
 # Add antibody to new assay object
 datObjs[[x]][['ADT']] <- CreateAssayObject(counts = rawDat[[x]][['Antibody Capture']])
 
 # Add data of non soupx corrected sample
 datObjs[["GC-WL-9991-LIVER"]] <- CreateSeuratObject(rawDat[["GC-WL-9991-LIVER"]][['Gene Expression']], project = x, assay = 'RNA')
 datObjs[["GC-WL-9991-LIVER"]][['ADT']] <- CreateAssayObject(counts = rawDat[["GC-WL-9991-LIVER"]][['Antibody Capture']])
 
 # Add sample name in metadata
 datObjs[[x]]@meta.data$Sample <- x
 datObjs[["GC-WL-9991-LIVER"]]@meta.data$Sample <- "GC-WL-9991-LIVER"

 # Plot mitochondria and hemoglobin genes to remove dead and doublet cells
 datObjs[[x]][["percent.mt"]] <- PercentageFeatureSet(datObjs[[x]], pattern = "^MT-")
 features <- c('HBA1','HBA2','HBB')
 datObjs[[x]][["percent.hb"]] <- PercentageFeatureSet(datObjs[[x]], features = features)
 plt <- VlnPlot(datObjs[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
 ggsave(
  plot = plt,
  filename = paste0(x, '_qcPlot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
 plt <- FeatureScatter(datObjs[[x]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
 ggsave(
  plot = plt,
  filename = paste0(x, '_lmPlot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
 # Filtering
 datObjs[[x]] <- subset(datObjs[[x]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
 plt <-VlnPlot(datObjs[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 4)
 ggsave(
  plot = plt,
  filename = paste0(x, '_afterQCplot.tiff'),
  height = 8,
  width = 8,
  dpi = 300,
  units = 'in'
 )
}

setwd("~/seurat")
saveRDS(datObjs, 'SeuObj.rds')

# Set working directory
setwd("~/seurat")

# Load data list
mylist <- readRDS("~/seurat/SeuObj.rds")

# Split list into LIVER and SVP objects
names(mylist)

# Remove longitudinal samples
mylist[["GC-WL-11303-LIVER"]] <- NULL
mylist[["GC-WL-11303-PBMC"]] <- NULL


# Subset LIVER
grep('LIVER', names(mylist))
LIVERlist <- mylist[grep('LIVER', names(mylist))]

# Subset SVP
grep('SAT-VAT-PBMC', names(mylist))
SVPlist <- mylist[grep('SAT-VAT-PBMC', names(mylist))]

#### Add metadata, normalise, FindVariableFeatures and scale LIVER data ####
# Name LIVER list
names(LIVERlist)
LIVER <- names(LIVERlist)[1]
LIVERlist[[LIVER]]

for (x in names(LIVERlist)) {
  
  # Create metadata columns for Tissue type and Patient_ID 
  LIVERlist[[x]]$Tissue <- 'LIVER'
  id <- gsub('-LIVER', '', x)
  id <- gsub('GC-WL-', '', id)
  LIVERlist[[x]]$Patient_ID <- id
  
  # Normalise ADT data
  DefaultAssay(LIVERlist[[x]]) <- 'ADT'
  LIVERlist[[x]][['ADTonly']] <- CreateAssayObject(counts = LIVERlist[[x]]@assays$ADT$counts[c(1:130), ])
  LIVERlist[[x]] <- NormalizeData(LIVERlist[[x]], assay = "ADTonly", normalization.method = "CLR", margin = 2)
  
  # Normalise and scale RNA data
  DefaultAssay(LIVERlist[[x]]) <- "RNA"
  LIVERlist[[x]] <- NormalizeData(LIVERlist[[x]])
  LIVERlist[[x]] <- FindVariableFeatures(LIVERlist[[x]])
  LIVERlist[[x]] <- ScaleData(LIVERlist[[x]])
}

#### Add metadata, normalise, FindVariableFeatures, scale and demultiplex SVP data ####
# Name SVP list
names(SVPlist)
SVP <- names(SVPlist)[1]
SVPlist[[SVP]]

# Check ADT list
rownames(SVPlist$`GC-WL-10113-1-SAT-VAT-PBMC`@assays$ADT$counts)

for (x in names(SVPlist)) {
  # Split ADT and HTO data into separate assays
  SVPlist[[x]][['HTO']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c('C0251', 'C0252', 'C0253'), ])
  SVPlist[[x]][['ADTonly']] <- CreateAssayObject(counts = SVPlist[[x]]@assays$ADT$counts[c(1:130), ])
  
  # Normalise ADT and HTO data
  SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "HTO", normalization.method = "CLR", margin = 2)
  SVPlist[[x]] <- NormalizeData(SVPlist[[x]], assay = "ADTonly", normalization.method = "CLR", margin = 2)
  
  # Demultiplex cells based on HTO enrichment
  SVPlist[[x]] <- HTODemux(SVPlist[[x]], assay = "HTO", positive.quantile = 0.95)
  
  # Keep only singlets
  SVPlist[[x]] <- subset(SVPlist[[x]], subset = HTO_classification.global == 'Singlet')
  
  # Demux results are saved in HTO classification column under meta data
  # Create tissue vector using HTO classification
  tissue_vector <- c(rep(NA, length(SVPlist[[x]]$HTO_classification)))
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0251'] <- 'SAT'
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0252'] <- 'VAT'
  tissue_vector[SVPlist[[x]]$HTO_classification == 'C0253'] <- 'PBMC'
  
  # Create metadata columns for Tissue type and Patient_ID and add tissue vector in tissue type
  SVPlist[[x]]$Tissue <- tissue_vector
  id <- gsub('-SAT-VAT-PBMC', '', x)
  id <- gsub('GC-WL-', '', id)
  SVPlist[[x]]$Patient_ID <- id
  
  # Normalise RNA data
  DefaultAssay(SVPlist[[x]]) <- "RNA"
  SVPlist[[x]] <- NormalizeData(SVPlist[[x]])
  SVPlist[[x]] <- FindVariableFeatures(SVPlist[[x]])
  SVPlist[[x]] <- ScaleData(SVPlist[[x]])
}

# Merge the liver and SAT, VAT and PBMC list
all.tissue.list <- c(LIVERlist, SVPlist)

SeuObj <- merge(x = all.tissue.list[[1]], 
                y = all.tissue.list[2:40])

saveRDS(SeuObj, 'SeuObjx.rds')

# End of the script
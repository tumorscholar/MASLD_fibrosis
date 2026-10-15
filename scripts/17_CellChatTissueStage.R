###############################################################################
# 17 - run CellChat per tissue x stage on shared T cells + non-T clusters
#
# One script, parameterised by tissue via the CC_TISSUE env var (LIVER/SAT/VAT),
# so it runs as three jobs (see 15_cellchat_run.slurm). Groups are the cell-type
# NAMES from 13; rare groups are dropped by a min-cell filter (>= min_cells in
# the tissue) rather than hardcoded C-numbers, so the group set is consistent.
#
# Input : SharedTcells_NonTcell_V7.rds
# Output: cellchat_<tissue>_<stage>_V7.rds  (per stage)
###############################################################################
library(Seurat); library(CellChat); library(dplyr); library(future)
options(stringsAsFactors = FALSE)
plan(sequential)                                # sequential = memory-safe (liver)
options(future.globals.maxSize = 20 * 1024^3)

tissue    <- toupper(Sys.getenv("CC_TISSUE", "LIVER"))    # LIVER / SAT / VAT
min_cells <- 10                                           # drop groups smaller than this
in_rds    <- "/data/Blizard-AlazawiLab/rk/seurat/SharedTcells_NonTcell_V7.rds"
out_dir   <- file.path("/data/Blizard-AlazawiLab/rk/cellchat_v7", tolower(tissue), "sharedTcells")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)
cat("CellChat tissue:", tissue, " | out:", out_dir, "\n")

SeuObj <- readRDS(in_rds)
if ("Sample" %in% colnames(SeuObj@meta.data))            # CellChat wants 'samples'
  colnames(SeuObj@meta.data)[colnames(SeuObj@meta.data) == "Sample"] <- "samples"
DefaultAssay(SeuObj) <- "RNA"

# subset the tissue, set idents by cell-type name
obj_t <- subset(SeuObj, subset = Tissue == tissue)
Idents(obj_t) <- "cluster_cellchat"
obj_t$cluster_cellchat <- droplevels(factor(obj_t$cluster_cellchat))

# COMMON cell-type set across ALL CellChat tissues: keep a group only if it has
# >= min_cells in BOTH disease stages (No_fibrosis & Fibrosis) in EVERY tissue
# (LIVER, SAT, VAT), so all three share an identical group set for direct
# comparison. Healthy is not used for the cutoff (one patient, not in the
# differential). Technical clusters dropped; "T cell" (focus) always kept.
drop_technical <- c("Cycling cells", "Erythrocytes")
cc_tissues     <- c("LIVER", "SAT", "VAT")
sub_all <- subset(SeuObj, subset = Tissue %in% cc_tissues)
sub_all$cluster_cellchat <- droplevels(factor(sub_all$cluster_cellchat))
tab3 <- table(sub_all$cluster_cellchat, sub_all$Tissue, sub_all$Stage)
dis  <- grep("fib", dimnames(tab3)[[3]], ignore.case = TRUE, value = TRUE)   # No_fibrosis & Fibrosis
common <- rownames(tab3)[vapply(rownames(tab3), function(g)
  all(vapply(cc_tissues, function(tt) all(tab3[g, tt, dis] >= min_cells), logical(1))), logical(1))]
common <- union("T cell", setdiff(common, drop_technical))
cat("COMMON groups across", paste(cc_tissues, collapse = "/"), "(>=", min_cells,
    "in", paste(dis, collapse = " & "), ", all tissues):\n"); print(sort(common))

keep <- intersect(common, levels(droplevels(factor(obj_t$cluster_cellchat))))
obj_t <- subset(obj_t, idents = keep)
lev  <- c(if ("T cell" %in% keep) "T cell", sort(setdiff(keep, "T cell")))
obj_t$cluster_cellchat <- droplevels(factor(obj_t$cluster_cellchat, levels = lev))
Idents(obj_t) <- "cluster_cellchat"

objList    <- SplitObject(obj_t, split.by = "Stage")
CellChatDB <- CellChatDB.human

run_cellchat <- function(seu, outfile) {
  Idents(seu) <- "cluster_cellchat"
  cc <- createCellChat(seu, group.by = "ident", assay = "RNA")
  cc@DB <- CellChatDB
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  saveRDS(cc, file = outfile)
  message("saved: ", outfile)
}

for (st in c("Healthy", "No_fibrosis", "Fibrosis")) {
  if (st %in% names(objList) && ncol(objList[[st]]) > 0) {
    cat("\n== ", tissue, "/", st, " (", ncol(objList[[st]]), " cells ) ==\n", sep = "")
    run_cellchat(objList[[st]], paste0("cellchat_", tolower(tissue), "_", st, "_V7.rds"))
  } else message("skip ", st, " (absent)")
}
cat("\nDone:", tissue, "\n")
# End of the script
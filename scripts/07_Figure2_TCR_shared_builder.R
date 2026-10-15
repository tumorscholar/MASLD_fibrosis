###############################################################################
# Figure 2 (d-i) - shared TCR repertoire inputs
#
# Builds the objects the Figure 2 d-i panels read, so they share one repertoire
# and one cell object:
#   combinedTCR_V7.rds   - combined TCR repertoire (all captures)
#   TcellObjV7.rds       - T cells with clonal-expansion annotation (2e-h)
#   SeuObjV7_clonal.rds  - full object with clone-size annotation (2i)
#   unique CDR3aa table  - unique clonotypes for external annotation (2d)
#
# TCR is per capture (liver = -LIVER-TCR, SVP tissues share one -SVP-TCR). We
# combineTCR with the CAPTURE names, so TCR barcodes become "<capture>_<bc>" and
# line up 1:1 with the Seurat cell names (which are also "<capture>_<bc>"); tissue
# is then inherited from the (HTO-demuxed) Seurat metadata. Clonality (Expanded =
# clonotype seen in >1 cell) is defined on the repertoire and mapped to cells.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(scRepertoire)
library(Seurat)
library(tidyverse)

# ---- CONFIG (v7) ------------------------------------------------------------
tcr_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
seurat_rds  <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"        # annotated main object (Fig 1)
out_rds     <- "/data/Blizard-AlazawiLab/rk/seurat"
out_files   <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"
dir.create(out_files, showWarnings = FALSE, recursive = TRUE)

liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-1","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 21
svp_pids   <- setdiff(liver_pids, "10291-1")                        # 20

stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

tcr_file_for <- function(capture) {
  cands <- c(file.path(tcr_dir, paste0(capture, "-TCR"), "outs", "filtered_contig_annotations.csv"),
             file.path(tcr_dir, paste0(capture, "-TCR"),          "filtered_contig_annotations.csv"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else stop("No v7 TCR file for capture: ", capture)
}

# ---- 1. combined.TCR (barcodes prefixed with capture = match Seurat names) --
captures    <- c(paste0("GC-WL-", liver_pids, "-LIVER"),
                 paste0("GC-WL-", svp_pids,   "-SVP"))
contig_list <- lapply(captures, function(cap) read.csv(tcr_file_for(cap), stringsAsFactors = FALSE))
combined.TCR <- combineTCR(contig_list, samples = captures,
                           removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
cap_pid <- sub("^GC-WL-(.*)-(LIVER|SVP)$", "\\1", captures)
combined.TCR <- addVariable(combined.TCR, variable.name = "Type", variables = unname(stage_map[cap_pid]))
saveRDS(combined.TCR, file.path(out_rds, "combinedTCR_V7.rds"))

# ---- 2. unique clonotypes for Fig 2d annotation -----------------------------
tcr_all <- bind_rows(combined.TCR)
uniq_ctaa <- tcr_all %>% filter(!is.na(CTaa), CTaa != "") %>% distinct(CTaa)
write.csv(uniq_ctaa, file.path(out_files, "unique_TRA_TRB_CTaa_V7.csv"), row.names = FALSE)
cat("Total TCR cells:", nrow(tcr_all), "| unique CDR3aa clonotypes:", nrow(uniq_ctaa), "\n")

# ---- where do the TCR cells come from? per-capture counts (v7) --------------
# TCR barcode is "<capture>_<10x-bc>"; strip the last field to get the capture.
tcr_all$capture <- sub("_[^_]*$", "", tcr_all$barcode)
tcr_all$pid     <- sub("^GC-WL-(.*)-(LIVER|SVP)$", "\\1", tcr_all$capture)
tcr_all$libtype <- sub("^GC-WL-.*-(LIVER|SVP)$", "\\1", tcr_all$capture)
per_cap <- tcr_all %>% count(capture, libtype, pid, name = "tcr_cells") %>% arrange(desc(tcr_cells))
cat("\n== TCR cells per capture (v7) ==\n"); print(as.data.frame(per_cap), row.names = FALSE)
cat("\n== TCR cells per library type (LIVER captures vs SVP captures) ==\n")
print(tcr_all %>% count(libtype, name = "tcr_cells"))
write.csv(per_cap, file.path(out_files, "tcr_cells_per_capture_V7.csv"), row.names = FALSE)

# ---- 3. integrate with Seurat; cloneSize (Fig 2i) + Clonality --------------
# Clone size pooled across the cohort, as in the original analysis (one added
# sample doesn't change the approach). Bins match the published Fig 2i legend
# (Small 1<X<=5, Medium 5<X<=10, Large 10<X<=20, Hyperexpanded 20<X<=max); the
# top bound is the observed max so no clone is silently dropped (the old '21'
# would have binned off everything larger).
gene_freq  <- tcr_all %>% filter(!is.na(CTgene), CTgene != "") %>% count(CTgene, name = "n")
max_clone  <- max(gene_freq$n)
clone_bins <- c(Single = 1, Small = 5, Medium = 10, Large = 20, Hyperexpanded = max_clone)
cat("Largest clone (pooled, by gene):", max_clone, "-> Hyperexpanded 20 < X <=", max_clone, "\n")

obj <- readRDS(seurat_rds)
obj <- combineExpression(combined.TCR, obj, cloneCall = "gene", group.by = "Type",
                         proportion = FALSE, cloneSize = clone_bins)

# Clonality (Expanded = CDR3aa seen in >1 cell across the cohort) - original method
clone_freq    <- tcr_all %>% filter(!is.na(CTaa), CTaa != "") %>% count(CTaa, name = "n")
expanded_ctaa <- clone_freq$CTaa[clone_freq$n > 1]

# ---- how many expanded / hyperexpanded clones (pooled cohort, by CDR3aa) ----
cs <- clone_freq$n
clone_summary <- data.frame(
  category   = c("Singletons (=1)", "Expanded (>1)", "  Small (1<X<=5)",
                 "  Medium (5<X<=10)", "  Large (10<X<=20)", "  Hyperexpanded (>20)"),
  clonotypes = c(sum(cs == 1), sum(cs > 1), sum(cs > 1 & cs <= 5),
                 sum(cs > 5 & cs <= 10), sum(cs > 10 & cs <= 20), sum(cs > 20)),
  cells      = c(sum(cs[cs == 1]), sum(cs[cs > 1]), sum(cs[cs > 1 & cs <= 5]),
                 sum(cs[cs > 5 & cs <= 10]), sum(cs[cs > 10 & cs <= 20]), sum(cs[cs > 20])))
cat("\n== Clone-size summary (pooled cohort, CDR3aa) ==\n"); print(clone_summary, row.names = FALSE)
write.csv(clone_summary, file.path(out_files, "clone_size_summary_V7.csv"), row.names = FALSE)
ctaa_by_bc    <- tcr_all %>% distinct(barcode, .keep_all = TRUE) %>% select(barcode, CTaa)
obj$CTaa_tcr  <- ctaa_by_bc$CTaa[match(colnames(obj), ctaa_by_bc$barcode)]
obj$Clonality <- factor(ifelse(is.na(obj$CTaa_tcr), NA_character_,
                               ifelse(obj$CTaa_tcr %in% expanded_ctaa, "Expanded", "Non-expanded")),
                        levels = c("Non-expanded", "Expanded"))

saveRDS(obj, file.path(out_rds, "SeuObjV7_clonal.rds"))          # full object -> Fig 2i

# ---- 4. T cells only (have a TCR) -> Fig 2e/f/g/h ---------------------------
TcellObj <- subset(obj, cells = colnames(obj)[!is.na(obj$Clonality)])
saveRDS(TcellObj, file.path(out_rds, "TcellObjV7.rds"))

cat("\nT cells:", ncol(TcellObj), "\n")
print(table(Tissue = TcellObj$Tissue, Clonality = TcellObj$Clonality))

# End of the script
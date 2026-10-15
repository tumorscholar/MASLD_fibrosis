###############################################################################
# 4 - Build baseline_expanded_clones_tracking.csv  (canonical, reproducible)
#
# Regenerates the WHOLE clone-tracking table read by 36 (donut), 37 (contraction)
# and 40 (gene usage), for the four paired FNA donors in BOTH tissues.
#
# SOURCE OF CELLS: the SAME QC'd FNA objects used by the heatmap (39) and the
# per-clone signature (38) - SeuObjFNA_LIVERV7.rds / SeuObjFNA_PBMCV7.rds - and
# the SAME CTaa attachment (Sample + 16bp-core match, with the 11303 alias so
# 10291-1's follow-up PBMC is included). This keeps the tracking table, the
# heatmaps and the per-clone analysis on one identical cell population.
#
# For each donor x tissue: clonotypes with >=2 cells at baseline are baseline-
# expanded and tracked into the matched follow-up sample.
#
#   Presence         = Followup_Count > 0 ? "Present in Both" : "Lost in Followup"
#   Expansion_Status = Followup_Count == 0 ? "Lost in Followup"
#                      : Percent_Change > 0 ? "Expanded" : "Contracted"
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
suppressMessages({library(Seurat); library(scRepertoire); library(dplyr)})

tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
ann_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"     # dir 36/37/40 read
dir.create(ann_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(ann_dir, "baseline_expanded_clones_tracking.csv")

tissues <- list(
  LIVER = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVERV7.rds",
  PBMC  = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_PBMCV7.rds")

# TCR library named differently from the GEX capture (obj$Sample): map so CTaa
# attaches. 10291-1 follow-up PBMC: Sample = GC-WL-11303-SVP, TCR = GC-WL-11303-PBMC-TCR.
sample_alias <- c("GC-WL-11303-PBMC" = "GC-WL-11303-SVP")
core16 <- function(x) sub(".*?([ACGT]{16}).*", "\\1", x)

# ---- per-cell (Donor, Visit, CTaa) for one FNA object -----------------------
cells_of <- function(TISSUE, rds) {
  cat("\n=====================", TISSUE, "=====================\n")
  if (!file.exists(rds)) { cat("  object not found:", rds, "- skipping.\n"); return(NULL) }
  seu <- readRDS(rds)
  if (!"CTaa" %in% colnames(seu@meta.data)) {
    caps     <- basename(dirname(dirname(
      list.files(tcr_dir, "filtered_contig_annotations.csv", recursive = TRUE, full.names = TRUE))))
    cap_base <- sub("-TCR$", "", caps)
    cap_samp <- ifelse(cap_base %in% names(sample_alias), sample_alias[cap_base], cap_base)
    take     <- which(cap_samp %in% unique(seu$Sample) & !duplicated(cap_samp))
    use_caps <- caps[take]; use_samp <- cap_samp[take]
    cat(TISSUE, "TCR capture -> Sample:\n"); print(setNames(use_samp, use_caps))
    if (!length(use_caps)) { cat("  no matching TCR captures - skipping.\n"); return(NULL) }
    cl <- lapply(use_caps, function(cap)
      read.csv(file.path(tcr_dir, cap, "outs", "filtered_contig_annotations.csv"), stringsAsFactors = FALSE))
    combined <- combineTCR(cl, samples = use_samp,
                           removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
    tcr_key <- do.call(rbind, lapply(names(combined), function(s)
      data.frame(key = paste(s, core16(combined[[s]]$barcode)),
                 CTaa = combined[[s]]$CTaa, stringsAsFactors = FALSE)))
    seu$CTaa <- tcr_key$CTaa[match(paste(seu$Sample, core16(colnames(seu))), tcr_key$key)]
  }
  if (!"Visit" %in% colnames(seu@meta.data))
    seu$Visit <- ifelse(grepl("FNA", seu$Sample), "Followup", "Baseline")
  if (!"Patient_ID" %in% colnames(seu@meta.data)) seu$Patient_ID <- "unknown"
  df <- data.frame(Tissue = TISSUE,
                   Donor  = as.character(seu$Patient_ID),
                   Visit  = ifelse(as.character(seu$Visit) == "Followup", "Followup", "Baseline"),
                   CTaa   = seu$CTaa, stringsAsFactors = FALSE)
  df <- df[!is.na(df$CTaa) & df$CTaa != "", , drop = FALSE]
  cat("  clonotype+ cells per Donor x Visit:\n"); print(table(df$Donor, df$Visit))
  df
}

cells <- bind_rows(Map(cells_of, names(tissues), tissues))
if (is.null(cells) || !nrow(cells)) stop("no clonotype+ cells - check object paths / TCR match.")

# ---- tracking rows for one donor x tissue -----------------------------------
track_one <- function(donor, tissue) {
  d  <- cells[cells$Donor == donor & cells$Tissue == tissue, ]
  bl <- d$CTaa[d$Visit == "Baseline"]; fu <- d$CTaa[d$Visit == "Followup"]
  bl_tot <- length(bl); fu_tot <- length(fu)
  if (!bl_tot) return(NULL)
  bcnt <- as.data.frame(table(CTaa = bl), stringsAsFactors = FALSE); names(bcnt)[2] <- "Baseline_Count"
  fcnt <- as.data.frame(table(CTaa = fu), stringsAsFactors = FALSE); names(fcnt)[2] <- "Followup_Count"
  bcnt %>%
    filter(Baseline_Count >= 2) %>%                          # baseline-expanded only
    left_join(fcnt, by = "CTaa") %>%
    mutate(
      Followup_Count        = ifelse(is.na(Followup_Count), 0L, Followup_Count),
      Donor                 = donor,
      Tissue                = tissue,
      Baseline_Total_Tcells = bl_tot,
      Followup_Total_Tcells = fu_tot,
      Baseline_Percent      = round(100 * Baseline_Count / bl_tot, 3),
      Followup_Percent      = round(100 * Followup_Count / fu_tot, 3),
      Change_Absolute       = Followup_Count - Baseline_Count,
      Percent_Change        = round(Followup_Percent - Baseline_Percent, 3),
      Presence              = ifelse(Followup_Count > 0, "Present in Both", "Lost in Followup"),
      Expansion_Status      = ifelse(Followup_Count == 0, "Lost in Followup",
                                     ifelse(Percent_Change > 0, "Expanded", "Contracted"))) %>%
    select(Donor, Tissue, CTaa, Baseline_Count, Baseline_Total_Tcells, Baseline_Percent,
           Followup_Count, Followup_Total_Tcells, Followup_Percent, Change_Absolute,
           Percent_Change, Presence, Expansion_Status)
}

donors   <- sort(unique(cells$Donor))
tracking <- bind_rows(lapply(c("LIVER", "PBMC"), function(ti)
  bind_rows(lapply(donors, function(d) track_one(d, ti)))))
write.csv(tracking, out_csv, row.names = FALSE)

cat("\nDonor x Tissue clone counts (tracked / re-detected at follow-up):\n")
tracking %>% group_by(Tissue, Donor) %>%
  summarise(tracked = n(), redetected = sum(Followup_Count > 0), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows). Re-run the annotation build, then 36, 37 and 40.\n",
            out_csv, nrow(tracking)))
# End of the script
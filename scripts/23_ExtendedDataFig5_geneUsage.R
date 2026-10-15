###############################################################################
# 23 - Extended Data Figure 5 - TCR repertoire gene usage & CDR3 diversity per tissue
#
# Rebuilt from 31_ExtendedDataFigure4_PBMC.R and generalised to all four tissues.
#   (a) percentGenes  - TRB V-gene and J-gene usage, by disease stage
#   (b) positionalEntropy - TRB positional amino-acid entropy, by stage
#
# Split into 4 (LIVER / PBMC / SAT / VAT). LIVER has 21 patients (incl. 10291-1,
# liver-only); SAT/VAT/PBMC have 20 - the per-tissue file glob picks up whatever
# contig files exist for each tissue.
#
# Input : per-tissue TCR contig CSVs (demux_contig/TCR)
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(scRepertoire); library(ggplot2); library(dplyr); library(grid)
source("00_nature_style.R")

tcr_dir <- "/data/home/hdx044/files/screpertoire/demux_contig/TCR"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

stage_lev <- c("Healthy", "No_fibrosis", "Fibrosis")     # matches stage_cols names
# 21-patient baseline map (F0 = No_fibrosis; F1/F2/F3 all = Fibrosis; 10742 Healthy)
stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

leg_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1),
                   axis.line = element_line(colour = "black"),
                   legend.position = "right", legend.direction = "vertical",
                   legend.title = element_text(size = 8), legend.text = element_text(size = 8),
                   legend.key.size = unit(1.2, "lines"))

# ---- one tissue -> V usage, J usage, positional entropy ---------------------
run_tissue <- function(tissue) {
  if (tissue == "LIVER") {
    # resolve ONE liver contig per cohort patient: demux (LIVER_LIVER) first,
    # else the non-demux (LIVER_) file. Iterating the 21-patient map guarantees
    # 21 files, brings in 10291-1 (which only has the non-demux file), and can
    # NEVER pick up 11303 (a follow-up sample, absent from the map).
    resolve <- function(p) {
      f1 <- file.path(tcr_dir, paste0("GC-WL-", p, "-LIVER_LIVER_TCR_contig.csv"))
      f2 <- file.path(tcr_dir, paste0("GC-WL-", p, "-LIVER_TCR_contig.csv"))
      if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NA_character_
    }
    pid   <- names(stage_map)
    files <- vapply(pid, resolve, "")
    keep  <- nzchar(files) & !is.na(files)
    files <- unname(files[keep]); pid <- pid[keep]
  } else {
    files <- list.files(tcr_dir, pattern = paste0("-SAT-VAT-PBMC_", tissue, "_TCR_contig\\.csv$"),
                        full.names = TRUE)
    pid   <- sub("^GC-WL-(.*)-SAT-VAT-PBMC_.*$", "\\1", basename(files))
    keep  <- pid %in% names(stage_map)
    files <- files[keep]; pid <- pid[keep]
  }
  cat(tissue, "- contig files:", length(files), "(", paste(pid, collapse = ","), ")\n")
  if (!length(files)) { cat("  no contig files for", tissue, "- skipping.\n"); return(invisible(NULL)) }
  
  combined.TCR <- combineTCR(lapply(files, read.csv, stringsAsFactors = FALSE),
                             samples = pid, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
  combined.TCR <- addVariable(combined.TCR, variable.name = "Type", variables = unname(stage_map[pid]))
  combined.TCR <- lapply(combined.TCR, function(df) {
    df$Type <- factor(df$Type, levels = stage_lev); df })
  ord <- intersect(stage_lev, unique(unlist(lapply(combined.TCR, function(d) as.character(d$Type)))))
  
  # (a) TRB V-gene usage  - 3 x 12 in. coord_flip() puts GENES on the (tall) y-axis
  # and STAGE on the x-axis, matching the original 31_ExtendedDataFigure4 layout.
  # After coord_flip the gene aesthetic (x) is titled with xlab(), the stage
  # aesthetic (y) with ylab() - hence xlab("Genes") / ylab("Stage").
  pv <- percentGenes(combined.TCR, chain = "TRB", gene = "Vgene", group.by = "Type", order.by = ord) +
    coord_flip() + ggtitle(tissue) + xlab("Genes") + ylab("Stage") + theme_nature(7) + leg_theme
  save_nature(pv, file.path(fig_dir, paste0("EDFig4_", tolower(tissue), "_percentVgenes_v7")), 3*2.54, 12*2.54)
  
  # (a) TRB J-gene usage  - 3 x 5 in (GENES on y-axis, STAGE on x-axis)
  pj <- percentGenes(combined.TCR, chain = "TRB", gene = "Jgene", group.by = "Type", order.by = ord) +
    coord_flip() + ggtitle(tissue) + xlab("Genes") + ylab("Stage") + theme_nature(7) + leg_theme
  save_nature(pj, file.path(fig_dir, paste0("EDFig4_", tolower(tissue), "_percentJgenes_v7")), 3*2.54, 5*2.54)
  
  # (b) TRB positional amino-acid entropy
  pe <- positionalEntropy(combined.TCR, chain = "TRB", group.by = "Type", order.by = ord, aa.length = 20) +
    scale_colour_manual(name = "Type", values = stage_cols) +
    ggtitle(tissue) + theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(colour = "black"),
          legend.position = "right", legend.title = element_text(size = 8),
          legend.text = element_text(size = 8), legend.key.size = unit(1.2, "lines"))
  save_nature(pe, file.path(fig_dir, paste0("EDFig4_", tolower(tissue), "_positionalEntropy_v7")), 5*2.54, 5*2.54)
  
  cat("  done", tissue, "\n")
}

sel <- toupper(Sys.getenv("CC_TISSUE", ""))
tissues <- if (nzchar(sel)) sel else c("LIVER", "PBMC", "SAT", "VAT")
for (tt in tissues) run_tissue(tt)
cat("\nDone Extended Data Fig 4 (", paste(tissues, collapse = ", "), ").\n")
# End of the script
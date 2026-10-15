###############################################################################
# 38 - Figure 5 - how did the SAME baseline-expanded clones change at follow-up?
#
# Paired, same-clone comparison of the expanded-clone CYTOTOXIC signature (the
# marker array used to separate expanded vs non-expanded T cells in Fig 2d/e and
# ED10), between baseline and follow-up.
#
# CLONE SET (important): restricted to the SAME clones as the contraction figure
# (37_Figure5_contraction.R) - i.e. the BASELINE-EXPANDED clones that were
# re-detected at follow-up ("Present in Both" in baseline_expanded_clones_tracking.csv).
# This is why n here matches the contraction n per tissue, instead of counting
# every clonotype seen at both visits (expanded or not).
#
# Runs for LIVER and PBMC. For each tissue:
#   Cell level : cytotoxic module score + per-gene dot plot, baseline vs follow-up.
#   Clone level: mean score per clonotype per visit -> paired plot + Wilcoxon,
#                plus contraction (frequency change) linked to the score change.
#
# Inputs : SeuObjFNA_<TISSUE>V7.rds, v7_tcr_TR captures, tracking CSV.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(Seurat); library(scRepertoire); library(dplyr); library(tidyr); library(ggplot2)
library(readr)
suppressWarnings(suppressMessages(library(SeuratExtend)))
source("00_nature_style.R")

tcr_dir  <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
ann_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
fig_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# per-tissue config: GEX object + the Tissue label used in the tracking CSV
tissues <- list(
  LIVER = list(rds = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVERV7.rds", trk_tissue = "LIVER"),
  PBMC  = list(rds = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_PBMCV7.rds",  trk_tissue = "PBMC")
)

# EXACT array from V7 Fig 2e (08_Figure2ef_V7.R) - expanded-clone signature
expanded_genes <- c("FGFBP2","GZMH","GNLY","ADGRG1","FCRL6","GZMB","CX3CR1","S1PR5",
                    "NKG7","KLRD1","PRF1","TBX21","PLEK","FCGR3A","ZEB2","CST7",
                    "CCL5","EFHD2","CTSW","SPON2","CD8A","CD4")
# for the cytotoxicity MODULE SCORE drop the lineage markers (CD4/CD8A) so the
# score reflects effector/cytotoxic programme, not CD4-vs-CD8 identity.
score_genes <- setdiff(expanded_genes, c("CD4", "CD8A"))

# baseline-expanded, re-detected clones per tissue (SAME set as the contraction fig)
trk <- readr::read_csv(file.path(ann_dir, "baseline_expanded_clones_tracking.csv"),
                       show_col_types = FALSE)

cols <- c("Contracted" = "#2166ac", "Expanded" = "#b2182b")

###############################################################################
run_tissue <- function(TISSUE, cfg) {
  cat("\n=====================", TISSUE, "=====================\n")
  if (!file.exists(cfg$rds)) { cat("  object not found:", cfg$rds, "- skipping.\n"); return(invisible(NULL)) }
  seu <- readRDS(cfg$rds)
  
  # ---- attach CTaa per cell (map TCR barcodes onto GEX cells) ----------------
  # Some captures have the TCR library named differently from the GEX capture
  # (obj$Sample). e.g. patient 10291-1 follow-up PBMC: GEX Sample = GC-WL-11303-SVP
  # but its TCR library is GC-WL-11303-PBMC-TCR. Map TCR-base -> GEX Sample so the
  # follow-up cells still receive their CTaa (otherwise they are silently dropped).
  sample_alias <- c("GC-WL-11303-PBMC" = "GC-WL-11303-SVP")
  if (!"CTaa" %in% colnames(seu@meta.data)) {
    caps     <- basename(dirname(dirname(
      list.files(tcr_dir, "filtered_contig_annotations.csv", recursive = TRUE, full.names = TRUE))))
    cap_base <- sub("-TCR$", "", caps)
    cap_samp <- ifelse(cap_base %in% names(sample_alias), sample_alias[cap_base], cap_base)
    take     <- which(cap_samp %in% unique(seu$Sample) & !duplicated(cap_samp))
    use_caps <- caps[take]; use_samp <- cap_samp[take]
    cat(TISSUE, "TCR capture -> GEX Sample:\n"); print(setNames(use_samp, use_caps))
    if (!length(use_caps)) { cat("  no matching TCR captures - skipping.\n"); return(invisible(NULL)) }
    cl <- lapply(use_caps, function(cap)
      read.csv(file.path(tcr_dir, cap, "outs", "filtered_contig_annotations.csv"),
               stringsAsFactors = FALSE))
    combined <- combineTCR(cl, samples = use_samp,
                           removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
    # Robust match on (Sample, 16-bp core barcode). GEX cell names differ by
    # tissue - liver is "<Sample>_<16bp>-1", but PBMC is "<cap>_<Visit>_<16bp>"
    # (no -1), so an exact-string match gives 0%. Pulling the 16-bp core from
    # both sides and pairing it with the sample is format-agnostic and cannot
    # collide across samples.
    core16  <- function(x) sub(".*?([ACGT]{16}).*", "\\1", x)
    tcr_key <- do.call(rbind, lapply(names(combined), function(s)
      data.frame(key  = paste(s, core16(combined[[s]]$barcode)),
                 CTaa = combined[[s]]$CTaa, stringsAsFactors = FALSE)))
    g_key <- paste(seu$Sample, core16(colnames(seu)))
    seu$CTaa <- tcr_key$CTaa[match(g_key, tcr_key$key)]
    cat("  cells with a TCR (CTaa):", sum(!is.na(seu$CTaa)), "of", ncol(seu),
        sprintf("(%.1f%% - only T cells carry a TCR)\n", 100 * mean(!is.na(seu$CTaa))))
  }
  
  # ---- Visit (infer from FNA capture name if the column is missing) ----------
  if (!"Visit" %in% colnames(seu@meta.data))
    seu$Visit <- ifelse(grepl("FNA", seu$Sample), "Followup", "Baseline")
  seu$Visit <- factor(seu$Visit, levels = c("Baseline", "Followup"))
  
  # ---- join RNA layers on the FULL object BEFORE any subsetting --------------
  # (a joined assay stays 2D for any later subset; joining after subsetting to a
  #  few clones can leave a 1-cell layer that collapses to a vector and errors.)
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
  if (!"data" %in% Layers(seu[["RNA"]])) seu <- NormalizeData(seu, verbose = FALSE)
  
  # ---- restrict to the SAME baseline-expanded clones as the contraction fig --
  exp_clones <- trk %>% filter(Tissue == cfg$trk_tissue, Presence == "Present in Both") %>%
    pull(CTaa) %>% unique()
  cat("  baseline-expanded clones re-detected at follow-up (tracking CSV):", length(exp_clones), "\n")
  
  keep_cells <- colnames(seu)[!is.na(seu$CTaa) & seu$CTaa %in% exp_clones]
  if (!length(keep_cells)) {
    cat("  0 cells map to these clones (check the TCR<->GEX match rate above) -",
        "skipping", TISSUE, "\n"); return(invisible(NULL)) }
  tc <- subset(seu, cells = keep_cells)
  # keep only clones that actually have GEX cells at BOTH visits (needed for a
  # paired score); the per-clone table below drops single-visit clones anyway.
  both <- tc@meta.data %>% distinct(CTaa, Visit) %>% count(CTaa) %>% filter(n == 2) %>% pull(CTaa)
  tc <- subset(tc, subset = CTaa %in% both)
  cat("  of those, with GEX cells at both visits:", length(both), "\n")
  cat("  cells in these clones per visit:\n"); print(table(tc$Visit))
  if (length(both) < 5) { cat("  <5 paired clones - too few to plot for", TISSUE, "\n"); return(invisible(NULL)) }
  
  # ---- cytotoxic module score ------------------------------------------------
  DefaultAssay(tc) <- "RNA"   # layers already joined + normalized above
  eg <- intersect(expanded_genes, rownames(tc))
  sg <- intersect(score_genes,    rownames(tc))
  tc <- AddModuleScore(tc, features = list(sg), name = "Cytotoxic", seed = 1234)
  tc$Cytotoxic <- tc$Cytotoxic1
  
  Idents(tc) <- "Visit"
  # Panel A (per-gene dot plot) not plotted - not used in the final figure.
  
  # ---- per clone: cytotoxic score (from GEX) + contraction status taken from
  #      the SAME frequency table as the contraction figure (37), so 'Contracted'
  #      means the identical thing in both panels (repertoire %, not GEX counts).
  frq <- trk %>% filter(Tissue == cfg$trk_tissue, Presence == "Present in Both") %>%
    transmute(CTaa,
              Baseline_Percent = as.numeric(Baseline_Percent),
              Followup_Percent = as.numeric(Followup_Percent)) %>% distinct()
  clone <- tc@meta.data %>% group_by(CTaa, Visit) %>%
    summarise(score = mean(Cytotoxic), n = dplyr::n(), .groups = "drop") %>%
    pivot_wider(names_from = Visit, values_from = c(score, n)) %>%
    filter(!is.na(score_Baseline), !is.na(score_Followup)) %>%
    left_join(frq, by = "CTaa") %>%
    filter(!is.na(Baseline_Percent), !is.na(Followup_Percent)) %>%
    mutate(d_score = score_Followup - score_Baseline,
           d_freq  = Followup_Percent - Baseline_Percent,   # repertoire %, same metric as fig 37
           status  = factor(ifelse(Followup_Percent < Baseline_Percent, "Contracted", "Expanded"),
                            levels = c("Contracted", "Expanded")))
  write.csv(clone, file.path(tab_dir, sprintf("Fig5_%s_sharedClone_cyto_contraction_v7.csv", TISSUE)),
            row.names = FALSE)
  pw <- wilcox.test(clone$score_Baseline, clone$score_Followup, paired = TRUE)
  dir_word <- if (mean(clone$d_score) < 0) "LOWER (less cytotoxic)" else "HIGHER (more cytotoxic)"
  # explicit reconciliation with the contraction figure's n and % contracted
  cat(sprintf("  RECONCILE: %d baseline-expanded clones re-detected (contraction fig 37); %d of them have\n",
              length(exp_clones), nrow(clone)))
  cat(sprintf("             a transcriptome at BOTH visits (rest lost a cell in GEX QC). Of those %d, %d (%.0f%%)\n",
              nrow(clone), sum(clone$status == "Contracted"), 100 * mean(clone$status == "Contracted")))
  cat("             are contracted - SAME frequency definition as fig 37 (which reports over all re-detected clones).\n")
  cat(sprintf("  Paired cytotoxic score (n=%d): p=%.3g, mean d=%+.3f -> follow-up is %s | Contracted=%d Expanded=%d\n",
              nrow(clone), pw$p.value, mean(clone$d_score), dir_word,
              sum(clone$status == "Contracted"), sum(clone$status == "Expanded")))
  
  # ---- Panel B: paired per-clone score, coloured by contraction status -------
  longs <- clone %>% select(CTaa, status, score_Baseline, score_Followup) %>%
    pivot_longer(c(score_Baseline, score_Followup), names_to = "Visit", values_to = "score") %>%
    mutate(Visit = factor(sub("score_", "", Visit), levels = c("Baseline", "Followup")))
  p_pair <- ggplot(longs, aes(Visit, score, group = CTaa, colour = status)) +
    geom_line(alpha = 0.5, linewidth = 0.3) + geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = cols, name = NULL) +
    annotate("text", x = 1.5, y = max(longs$score), vjust = 1,
             label = sprintf("n=%d clones, paired p=%.3g", nrow(clone), pw$p.value), size = 2.3) +
    labs(x = NULL, y = "Cytotoxic score (per clone)",
         title = sprintf("%s: same clones at follow-up", TISSUE)) +
    theme_nature(7) + theme(legend.position = "bottom")
  # per-clone paired cytoscore saved for BOTH tissues (ED12a: liver + PBMC, matching)
  save_nature(p_pair, file.path(fig_dir, sprintf("Fig5_%s_sharedClone_cytoscore_paired_v7", TISSUE)), 7, 7)
  
  # ---- contraction-vs-cytotoxicity stats (console only; panels C/D not plotted)
  rho <- suppressWarnings(cor.test(clone$d_freq, clone$d_score, method = "spearman"))
  pw2 <- tryCatch(wilcox.test(d_score ~ status, data = clone)$p.value, error = function(e) NA_real_)
  cat(sprintf("  contraction vs cytotoxicity: Spearman rho=%.2f p=%.3g | d_score Contracted vs Expanded p=%.3g\n",
              rho$estimate, rho$p.value, pw2))
  
  invisible(data.frame(Tissue = TISSUE, n_clones = nrow(clone), paired_p = pw$p.value,
                       mean_d_score = mean(clone$d_score),
                       n_contracted = sum(clone$status == "Contracted"),
                       n_expanded   = sum(clone$status == "Expanded")))
}

# ---- run both tissues + a combined summary ---------------------------------
res <- do.call(rbind, lapply(names(tissues), function(t) run_tissue(t, tissues[[t]])))
cat("\n================ Fig 5 clone-expression summary ================\n")
if (!is.null(res)) { print(res, row.names = FALSE)
  write.csv(res, file.path(tab_dir, "Fig5_sharedClone_cyto_summary_v7.csv"), row.names = FALSE) }
cat("\nDone Fig 5 clone expression (LIVER + PBMC, baseline-expanded clones).\n")
# End of the script
###############################################################################
# 39 - Figure 5 - whole-repertoire expanded-clone signature, baseline vs FU
#
# Companion to 38 (the PAIRED same-clone analysis). This is the REPERTOIRE-LEVEL
# view, in the spirit of Extended Data Fig 10 / Fig 2e: it averages the
# expanded-clone signature (the SAME 22-gene array used to separate expanded vs
# non-expanded T cells) over the repertoire and compares baseline vs follow-up
# in LIVER and PBMC.
#
# TWO populations, side by side:
#   (1) All clonotype+ T cells  - every T cell with a recovered TCR (any size)
#   (2) Expanded clones (size>1) - clone observed in >1 cell (Fig 2e definition)
#
# OUTPUTS
#   (A) POOLED heatmap  : genes x {LIVER BL/FU, PBMC BL/FU}, faceted by population
#   (B) PER-PATIENT heatmaps (SUPPLEMENTARY): one figure per patient x tissue,
#       genes x {BL, FU}, faceted total-T vs expanded-T. Colour is the within-gene
#       z (computed across all patients+visits in that population, so patient
#       panels are directly comparable); number = mean log-norm expression;
#       star on FU tile = per-patient Wilcoxon baseline vs follow-up (cell level).
#
# NOTE on interpretation (kept honest): this repertoire-level average can move
# DOWN even when individual persisting clones move UP (script 38) - because the
# large cytotoxic clones contract in FREQUENCY, so fewer cytotoxic cells remain.
#
# Inputs : SeuObjFNA_<TISSUE>V7.rds + v7_tcr_TR captures.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(Seurat); library(scRepertoire); library(dplyr); library(tidyr); library(ggplot2)
source("00_nature_style.R")

tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

tissues <- list(
  LIVER = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVERV7.rds",
  PBMC  = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_PBMCV7.rds"
)

# EXACT array from V7 Fig 2e (08_Figure2ef_V7.R)
expanded_genes <- c("FGFBP2","GZMH","GNLY","ADGRG1","FCRL6","GZMB","CX3CR1","S1PR5",
                    "NKG7","KLRD1","PRF1","TBX21","PLEK","FCGR3A","ZEB2","CST7",
                    "CCL5","EFHD2","CTSW","SPON2","CD8A","CD4")
score_genes <- setdiff(expanded_genes, c("CD4", "CD8A"))   # cytotoxic module (no lineage)

core16 <- function(x) sub(".*?([ACGT]{16}).*", "\\1", x)

###############################################################################
# returns cell-level data.frame: Tissue, Patient_ID, Visit, expanded, Cytotoxic + gene cols
per_tissue <- function(TISSUE, rds) {
  cat("\n=====================", TISSUE, "=====================\n")
  if (!file.exists(rds)) { cat("  object not found:", rds, "- skipping.\n"); return(NULL) }
  seu <- readRDS(rds)
  
  # ---- attach CTaa per cell (robust Sample + 16bp-core match) ----------------
  # Some captures have the TCR library named differently from the GEX capture
  # (obj$Sample). Map TCR-capture base name -> GEX Sample so CTaa still attaches.
  # e.g. patient 10291-1 follow-up PBMC: GEX Sample = GC-WL-11303-SVP but its TCR
  # library is GC-WL-11303-PBMC-TCR (the cellranger --id is GC-WL-11303-PBMC).
  # Add any further GEX/TCR naming mismatches here as "TCR-base" = "GEX-Sample".
  sample_alias <- c("GC-WL-11303-PBMC" = "GC-WL-11303-SVP")
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
  seu$Visit <- factor(seu$Visit, levels = c("Baseline", "Followup"))
  if (!"Patient_ID" %in% colnames(seu@meta.data)) seu$Patient_ID <- "unknown"
  # diagnostic: clonotype recovery per Patient x Visit (spot a dropped follow-up,
  # e.g. 10291-1 PBMC follow-up before the TCR-name alias was added)
  cat("  CTaa recovery (TRUE = cell has a clonotype), Patient_ID x Visit:\n")
  print(table(Patient = seu$Patient_ID, Visit = seu$Visit, hasCTaa = !is.na(seu$CTaa)))
  
  # ---- join + normalize, keep ALL clonotype+ T cells -------------------------
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu)
  if (!"data" %in% Layers(seu[["RNA"]])) seu <- NormalizeData(seu, verbose = FALSE)
  keep <- colnames(seu)[!is.na(seu$CTaa)]
  if (!length(keep)) { cat("  0 clonotype+ cells (check TCR match) - skipping.\n"); return(NULL) }
  tc <- subset(seu, cells = keep)
  
  # clone size within this tissue object -> expanded = observed in >1 cell (Fig2e)
  sz <- table(tc$CTaa); tc$expanded <- as.logical(sz[as.character(tc$CTaa)] > 1)
  cat("  clonotype+ T cells per visit:\n"); print(table(tc$Visit))
  cat("  expanded (size>1) per visit:\n");  print(table(tc$Visit, tc$expanded))
  
  genes <- intersect(expanded_genes, rownames(tc))
  expr  <- FetchData(tc, vars = genes, layer = "data")
  # pad any gene missing in this tissue so tissues bind cleanly
  for (g in setdiff(expanded_genes, genes)) expr[[g]] <- NA_real_
  expr <- expr[, expanded_genes, drop = FALSE]
  tc <- AddModuleScore(tc, features = list(intersect(score_genes, rownames(tc))),
                       name = "Cytotoxic", seed = 1234)
  data.frame(Tissue = TISSUE, Patient_ID = as.character(tc$Patient_ID),
             Visit = tc$Visit, expanded = tc$expanded,
             Cytotoxic = tc$Cytotoxic1, expr, check.names = FALSE)
}

cells <- bind_rows(lapply(names(tissues), function(t) per_tissue(t, tissues[[t]])))
if (is.null(cells) || !nrow(cells)) stop("No tissue produced data - check paths / TCR match.")

star <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", ""))))

###############################################################################
# (A) POOLED heatmap  (genes x tissue/visit, faceted by population)
###############################################################################
agg <- function(df, pop_label) {
  long <- df %>% select(Tissue, Visit, all_of(expanded_genes)) %>%
    pivot_longer(-c(Tissue, Visit), names_to = "gene", values_to = "e") %>% filter(!is.na(e))
  m <- long %>% group_by(Tissue, gene, Visit) %>% summarise(mean_expr = mean(e), .groups = "drop")
  p <- long %>% group_by(Tissue, gene) %>%
    summarise(p = tryCatch(wilcox.test(e ~ Visit)$p.value, error = function(x) NA_real_), .groups = "drop")
  for (tt in unique(df$Tissue)) { s <- df[df$Tissue == tt, ]
  mu <- tapply(s$Cytotoxic, s$Visit, mean)
  cat(sprintf("  [%s | %s] cytotoxic module score  BL=%.3f  FU=%.3f  d=%+.3f -> %s\n",
              pop_label, tt, mu["Baseline"], mu["Followup"], mu["Followup"] - mu["Baseline"],
              ifelse(mu["Followup"] < mu["Baseline"], "LESS cytotoxic at FU", "MORE cytotoxic at FU"))) }
  left_join(m, p, by = c("Tissue", "gene")) %>% mutate(pop = pop_label)
}

cat("\n--- module-score summary by population (pooled) ---\n")
dat <- bind_rows(
  agg(cells,                    "All clonotype+ T cells"),
  agg(cells[cells$expanded, ],  "Expanded clones (size>1)")
) %>%
  mutate(col = factor(paste(Tissue, ifelse(Visit == "Baseline", "BL", "FU")),
                      levels = c("LIVER BL","LIVER FU","PBMC BL","PBMC FU")),
         gene = factor(gene, levels = rev(expanded_genes)),
         pop  = factor(pop, levels = c("All clonotype+ T cells","Expanded clones (size>1)"))) %>%
  filter(!is.na(col)) %>%
  group_by(pop, gene) %>%
  mutate(z = { s <- sd(mean_expr); if (is.na(s) || s == 0) 0 else (mean_expr - mean(mean_expr)) / s }) %>%
  ungroup() %>%
  mutate(lab = ifelse(grepl("FU$", col), sprintf("%.2f%s", mean_expr, star(p)), sprintf("%.2f", mean_expr)))
write.csv(dat, file.path(tab_dir, "Fig5_repertoire_signature_heatmap_v7.csv"), row.names = FALSE)

hm <- ggplot(dat, aes(col, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 1.8) +
  facet_wrap(~ pop) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       name = "z-score\n(within gene)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Repertoire-wide expanded-clone signature, baseline vs follow-up",
       subtitle = "number = mean log-norm expression; colour = within-gene z; * baseline vs follow-up (Wilcoxon)") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5.5),
        panel.grid = element_blank(), plot.subtitle = element_text(size = 5.5),
        strip.background = element_rect(fill = "grey92", colour = NA))
save_nature(hm, file.path(fig_dir, "Fig5_repertoire_signature_heatmap_v7"), 18, 14)   # 2-tissue (supp)

# ---- MAIN Fig 5d: LIVER-only bulk heatmap (all clonotype+ | expanded) --------
# whole-repertoire signature falls at follow-up (contraction of the cytotoxic pool);
# the per-clone persistence is the supplementary nuance (script 38 / ED12).
datL <- dplyr::filter(dat, col %in% c("LIVER BL","LIVER FU"))
hmL <- ggplot(datL, aes(col, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 1.9) +
  facet_wrap(~ pop) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       name = "z-score\n(within gene)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Liver clonotype+ T cells: cytotoxic/effector signature, baseline vs follow-up",
       subtitle = "number = mean log-norm expression; colour = within-gene z; * baseline vs follow-up (Wilcoxon)") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
        panel.grid = element_blank(), plot.subtitle = element_text(size = 5.2),
        strip.background = element_rect(fill = "grey92", colour = NA))
save_nature(hmL, file.path(fig_dir, "Fig5d_repertoire_signature_heatmap_liver_v7"), 11, 14)

###############################################################################
# (B) PER-PATIENT heatmaps (supplementary) - total-T and expanded-T
###############################################################################
agg_patient <- function(df, pop_label) {
  long <- df %>% select(Tissue, Patient_ID, Visit, all_of(expanded_genes)) %>%
    pivot_longer(-c(Tissue, Patient_ID, Visit), names_to = "gene", values_to = "e") %>% filter(!is.na(e))
  m <- long %>% group_by(Tissue, Patient_ID, gene, Visit) %>%
    summarise(mean_expr = mean(e), n_cell = n(), .groups = "drop")
  p <- long %>% group_by(Tissue, Patient_ID, gene) %>%
    summarise(p = tryCatch(wilcox.test(e ~ Visit)$p.value, error = function(x) NA_real_), .groups = "drop")
  left_join(m, p, by = c("Tissue", "Patient_ID", "gene")) %>% mutate(pop = pop_label)
}

datP <- bind_rows(
  agg_patient(cells,                   "All clonotype+ T cells"),
  agg_patient(cells[cells$expanded, ], "Expanded clones (size>1)")
) %>%
  mutate(pop  = factor(pop, levels = c("All clonotype+ T cells","Expanded clones (size>1)")),
         gene = factor(gene, levels = rev(expanded_genes)),
         Visit = factor(ifelse(Visit == "Baseline", "BL", "FU"), levels = c("BL","FU"))) %>%
  # z within gene across ALL patients+visits of that population -> panels comparable
  group_by(pop, gene) %>%
  mutate(z = { s <- sd(mean_expr); if (is.na(s) || s == 0) 0 else (mean_expr - mean(mean_expr)) / s }) %>%
  ungroup() %>%
  mutate(lab = ifelse(Visit == "FU", sprintf("%.2f%s", mean_expr, star(p)), sprintf("%.2f", mean_expr)))
write.csv(datP, file.path(tab_dir, "Fig5_repertoire_signature_heatmap_perPatient_v7.csv"), row.names = FALSE)

zlim <- max(abs(datP$z), na.rm = TRUE)
plot_patient <- function(TIS, pid) {
  d <- filter(datP, Tissue == TIS, Patient_ID == pid)
  if (!nrow(d)) return(invisible(NULL))
  p <- ggplot(d, aes(Visit, gene, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 1.7) +
    facet_wrap(~ pop) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                         limits = c(-zlim, zlim), name = "z-score\n(within gene)") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         title = sprintf("%s FNA - patient %s: signature baseline vs follow-up", TIS, pid),
         subtitle = "number = mean log-norm expression; * baseline vs follow-up (Wilcoxon, cell level)") +
    theme_nature(7) +
    theme(axis.text.x = element_text(size = 6), panel.grid = element_blank(),
          plot.subtitle = element_text(size = 5.2),
          strip.background = element_rect(fill = "grey92", colour = NA))
  save_nature(p, file.path(fig_dir, sprintf("Fig5_repertoire_signature_heatmap_Patient_%s_%s_v7", pid, TIS)),
              9, 12)
  cat("  wrote per-patient heatmap:", TIS, pid, "\n")
}
# per-patient heatmaps for BOTH tissues (4 liver + 4 PBMC) -> supplementary.
# The pooled LIVER bulk heatmap is main Fig 5d; these per-patient panels show the
# signature is consistent across individuals in each tissue.
for (TIS in unique(datP$Tissue))
  for (pid in sort(unique(datP$Patient_ID[datP$Tissue == TIS]))) plot_patient(TIS, pid)

cat("\nDone Fig 5 repertoire signature heatmaps (pooled + per-patient; all clonotype+ and expanded).\n")
# End of the script
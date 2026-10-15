###############################################################################
# Figure 3a-c - cross-tissue sharing of expanded T-cell clonotypes
#
# Fig 3a : 4-way tissue Venn (LIVER/SAT/VAT/PBMC) of EXPANDED clonotypes,
#          one panel per fibrosis stage (Healthy / No Fibrosis / Fibrosis).
# Fig 3b : % of hyperexpanded T cells (clones > 20 cells) among all TCR-defined
#          T cells, per tissue.
# Fig 3c : 4-way tissue Venn of HYPEREXPANDED clonotypes (No Fibrosis / Fibrosis).
#
# It also classifies every expanded clonotype by its tissue footprint
# (PBMC-only / one solid tissue / shared between solid tissues) and writes the
# shared-class-annotated T-cell object that Fig 3 d-i and the SCENIC panels read:
#       expandedTcellFinalV7.rds
#
# TCR handling matches Fig 2 (06): TCR is per capture (liver = -LIVER, the three
# SVP tissues share one -SVP library); the SVP library is split to SAT/VAT/PBMC
# via the HTO-demuxed GEX barcode->tissue map. Clone size is pooled across the
# cohort, as in the original analysis.
###############################################################################

.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))   # scRepertoire library location
library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ggVennDiagram)
library(patchwork)
source("00_nature_style.R")            # theme_nature(), stage_cols, save_nature()

# ---- v7 data locations ------------------------------------------------------
tcr_dir     <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"        # per capture: GC-WL-<cap>-TCR/outs/
barcode_dir <- "/data/Blizard-AlazawiLab/rk/barcodes_v7"             # 01_demultiplex bc->Tissue tables (NO QC)
tcell_rds<- "/data/Blizard-AlazawiLab/rk/seurat/TcellObjV7.rds"       # T cells w/ CTaa_tcr + Clonality
obj_out  <- "/data/Blizard-AlazawiLab/rk/seurat/expandedTcellFinalV7.rds"
fig_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure3"
tab_dir  <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure3"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---- patients per tissue + stage (same as 06/07) ----------------------------
liver_pids <- c("9680","9932","9961","9991","9999","10113-1","10113-2","10202",
                "10203","10205","10291-1","10291-2","10380","10634","10738","10742",
                "11040","11051","11183","11327","11471")            # 21
svp_pids   <- setdiff(liver_pids, "10291-1")                        # 20

stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")
stage_lev <- c("Healthy","No_fibrosis","Fibrosis")
disp      <- c("Healthy"="Healthy","No_fibrosis"="No Fibrosis","Fibrosis"="Fibrosis")

# ---- barcode -> tissue map (from 01_demultiplex: NO cell-level QC) ----------
# Use the pre-QC demux tables, NOT the QC'd object (02): tissue only needs the
# HTO call, so no clonotype-bearing cell should be lost to RNA/ADT QC. (Sourcing
# from SeuObjV7_preQC dropped ~33% of TCR cells - up to 99% in captures with an
# ADT/RNA barcode-match failure, e.g. 10634-LIVER - which is not real signal.)
bfiles <- list.files(barcode_dir, pattern = "_metadata\\.csv$", full.names = TRUE)
if (!length(bfiles)) stop("No barcode tables in ", barcode_dir, " - run 01_demultiplex.R first.")
gex_map <- do.call(rbind, lapply(bfiles, read.csv, stringsAsFactors = FALSE))
gex_map <- data.frame(Sample = gex_map$sample, Tissue = gex_map$Tissue,
                      bc = gex_map$cell_barcode, stringsAsFactors = FALSE)
gex_map <- gex_map[!is.na(gex_map$Tissue), , drop = FALSE]        # keep HTO-assigned cells
cat("demux barcode->tissue rows:", nrow(gex_map), "\n")

capture_of <- function(tissue, pid)
  if (tissue == "LIVER") paste0("GC-WL-", pid, "-LIVER") else paste0("GC-WL-", pid, "-SVP")

tcr_file_for <- function(capture) {
  cands <- c(file.path(tcr_dir, paste0(capture, "-TCR"), "outs", "filtered_contig_annotations.csv"),
             file.path(tcr_dir, paste0(capture, "-TCR"),          "filtered_contig_annotations.csv"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else stop("No v7 TCR file for capture: ", capture)
}

# contigs for one (tissue, patient):
#  LIVER -> unambiguous tissue, keep ALL productive-TCR cells (no GEX needed);
#           the manuscript reads the liver contig file whole. GEX-filtering liver
#           wrongly collapsed captures whose GEX/TCR barcodes disagree (e.g. 10634).
#  SVP   -> pooled library, demux SAT/VAT/PBMC via the HTO barcode->tissue map.
read_contig <- function(tissue, pid) {
  cap <- capture_of(tissue, pid)
  d   <- read.csv(tcr_file_for(cap), stringsAsFactors = FALSE)
  if (tissue == "LIVER") return(d)
  keep <- gex_map$bc[gex_map$Sample == cap & gex_map$Tissue == tissue]
  d[d$barcode %in% keep, , drop = FALSE]
}

# one long TCR table per tissue: combineTCR gives CTaa, then tag Tissue + patient
build_tcr_df <- function(tissue, pids) {
  cl <- lapply(pids, function(p) read_contig(tissue, p))
  ct <- combineTCR(cl, samples = pids, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
  ct <- addVariable(ct, variable.name = "Type", variables = unname(stage_map[pids]))
  df <- bind_rows(ct)
  df$Tissue     <- tissue
  df$Patient_ID <- df$sample                 # combineTCR sets 'sample' = pid
  df
}

tissues    <- list(LIVER = liver_pids, PBMC = svp_pids, SAT = svp_pids, VAT = svp_pids)
all_tcells <- bind_rows(Map(build_tcr_df, names(tissues), tissues)) %>%
  filter(!is.na(CTaa), CTaa != "")
all_tcells$Stage <- factor(stage_map[all_tcells$Patient_ID], levels = stage_lev)

# ---- pooled clone sizes (cohort-wide, by CDR3aa) ----------------------------
clone_counts <- all_tcells %>% count(CTaa, name = "n_cells")
expanded_clones      <- clone_counts %>% filter(n_cells > 1)
hyperexpanded_clones <- clone_counts %>% filter(n_cells > 20)

clone_summary <- data.frame(
  category   = c("Total TCR cells", "Unique clonotypes", "Singletons (=1)",
                 "Expanded (>1)", "Hyperexpanded (>20)"),
  clonotypes = c(NA, nrow(clone_counts), sum(clone_counts$n_cells == 1),
                 nrow(expanded_clones), nrow(hyperexpanded_clones)),
  cells      = c(nrow(all_tcells), NA, sum(clone_counts$n_cells == 1),
                 sum(expanded_clones$n_cells), sum(hyperexpanded_clones$n_cells)))
cat("\n== Fig 3 clone summary (pooled cohort, CDR3aa) ==\n"); print(clone_summary, row.names = FALSE)
write.csv(clone_summary, file.path(tab_dir, "Fig3_clone_summary_v7.csv"), row.names = FALSE)

# long CTaa x Stage x tissue table (one row per clone per tissue per stage)
clone_long <- all_tcells %>%
  distinct(CTaa, Stage, Tissue) %>%
  mutate(Stage = as.character(Stage))

# ============================ Fig 3a / 3c: Venns =============================
tissue_order <- c("LIVER", "SAT", "VAT", "PBMC")
venn_stage <- function(ctaa_keep, stage_key, title_word) {
  d <- clone_long %>% filter(Stage == stage_key, CTaa %in% ctaa_keep)
  if (!nrow(d)) { message("No clones for ", title_word, " / ", stage_key); return(NULL) }
  sets <- d %>% group_by(Tissue) %>% summarise(cl = list(unique(CTaa)), .groups = "drop") %>% deframe()
  sets <- sets[intersect(tissue_order, names(sets))]         # fixed panel order
  n_st <- n_distinct(d$CTaa)
  ggVennDiagram(sets, label_alpha = 0, edge_size = 0.6, set_size = 2.6, label_size = 2.2) +
    scale_fill_gradient(low = "#FFFFFF", high = unname(stage_cols[stage_key])) +
    scale_colour_manual(values = rep("grey30", length(sets))) +
    # just enough drawing area so the LIVER / PBMC set labels are not clipped and
    # the region-count text top/bottom has headroom - tightened to remove the
    # large empty padding around and between the Venns.
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    theme_void(base_family = "Arial") +                    # remove the X / Y axes entirely
    theme(legend.position = "none",
          plot.title  = element_text(hjust = 0.5, size = 7, family = "Arial"),
          plot.margin = margin(1, 2, 1, 2)) +              # small margin -> Venns sit closer
    labs(title = paste0(title_word, " - ", disp[stage_key], " (n = ", n_st, ")"))
}

# Fig 3a - expanded clonotypes, one Venn per stage
a_plots <- lapply(stage_lev, function(s) venn_stage(expanded_clones$CTaa, s, "Expanded clonotypes"))
a_plots <- a_plots[!sapply(a_plots, is.null)]
fig3a <- wrap_plots(a_plots, nrow = 1)
# Venns are fixed-aspect (~square); a 30-cm-wide canvas for 3 panels left a lot
# of empty space on either side of each Venn. Narrow the canvas so each panel is
# about as wide as the Venn is tall (~6-7 cm), removing the outside-Venn space.
save_nature(fig3a, file.path(fig_dir, "Fig3a_expanded_venn_v7"), w_cm = 21, h_cm = 6)

# Fig 3c - hyperexpanded clonotypes (No Fibrosis / Fibrosis; none in Healthy)
c_plots <- lapply(c("No_fibrosis", "Fibrosis"),
                  function(s) venn_stage(hyperexpanded_clones$CTaa, s, "Hyperexpanded clonotypes"))
c_plots <- c_plots[!sapply(c_plots, is.null)]
fig3c <- wrap_plots(c_plots, nrow = 1)
# 2 fixed-aspect Venns - narrow the canvas to match (removes outside-Venn space).
save_nature(fig3c, file.path(fig_dir, "Fig3c_hyperexpanded_venn_v7"), w_cm = 14, h_cm = 6)

# ---- per-stage 4-way hyperexpanded sharing counts (for Results text [85]) ----
# For each stage: total hyperexpanded clonotypes, how many are shared across ALL
# FOUR tissues (LIVER/SAT/VAT/PBMC), and how many are absent from the liver.
hyper_venn_counts <- do.call(rbind, lapply(stage_lev, function(s) {
  d <- clone_long %>% filter(Stage == s, CTaa %in% hyperexpanded_clones$CTaa)
  if (!nrow(d)) return(NULL)
  memb  <- d %>% group_by(CTaa) %>% summarise(tis = list(sort(unique(Tissue))), .groups = "drop")
  total <- nrow(memb)
  all4  <- sum(vapply(memb$tis, function(t) all(c("LIVER","SAT","VAT","PBMC") %in% t), logical(1)))
  noliv <- sum(vapply(memb$tis, function(t) !("LIVER" %in% t), logical(1)))
  data.frame(Stage = s, n_hyperexpanded = total,
             shared_all4 = all4, pct_shared_all4 = round(100 * all4 / total, 1),
             absent_from_liver = noliv, pct_absent_from_liver = round(100 * noliv / total, 1))
}))
cat("\n== Fig 3c hyperexpanded 4-way sharing per stage (Results [85]) ==\n")
print(hyper_venn_counts, row.names = FALSE)
write.csv(hyper_venn_counts, file.path(tab_dir, "Fig3c_hyperexpanded_venn_counts_v7.csv"), row.names = FALSE)


# ============================ Fig 3b: % hyperexpanded (stacked) ==============
# Bar per tissue (x = LIVER/PBMC/SAT/VAT); height = % of that tissue's TCR-defined
# T cells that are hyperexpanded (>20-cell clone); each bar STACKED by stage
# (Healthy / No fibrosis / Fibrosis). Denominator is the tissue total, so the
# stage segments add up to the tissue's overall hyperexpanded %. Because
# hyperexpanded clones appear only in MASLD, the Healthy segment is ~0.
hyper_ctaa <- hyperexpanded_clones$CTaa
tissue_totals <- all_tcells %>% count(Tissue, name = "n_tissue")
b_stack <- all_tcells %>%
  filter(CTaa %in% hyper_ctaa) %>%
  count(Tissue, Stage, name = "n_hyper") %>%
  tidyr::complete(Tissue = c("LIVER","PBMC","SAT","VAT"), Stage = stage_lev,
                  fill = list(n_hyper = 0)) %>%
  left_join(tissue_totals, by = "Tissue") %>%
  mutate(pct_hyper = 100 * n_hyper / n_tissue,
         Tissue = factor(Tissue, levels = c("LIVER","PBMC","SAT","VAT")),
         Stage  = factor(Stage,  levels = stage_lev))
cat("\n== Fig 3b: % hyperexpanded cells per tissue x stage (stacked) ==\n")
print(as.data.frame(b_stack), row.names = FALSE)
write.csv(b_stack, file.path(tab_dir, "Fig3b_pct_hyperexpanded_stacked_v7.csv"), row.names = FALSE)

stage_fill_b <- scale_fill_manual(name = "Stage", values = stage_cols,
                                  labels = c("Healthy", "No fibrosis", "Fibrosis"))
fig3b <- ggplot(b_stack, aes(x = Tissue, y = pct_hyper, fill = Stage)) +
  geom_col(width = 0.7) +                                   # stacked (default position)
  stage_fill_b +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hyperexpanded T cells\n(% of TCR-defined T cells)") +
  theme_nature(7) +
  theme(
    legend.position = "right", legend.key.size = unit(0.7, "lines"))
save_nature(fig3b, file.path(fig_dir, "Fig3b_pct_hyperexpanded_stacked_v7"), w_cm = 7, h_cm = 5)

# ============ tissue-footprint classification of EXPANDED clones =============
# per clonotype: which tissues it is seen in (cohort-wide), then classify.
classify <- function(tv) {
  solid <- intersect(tv, c("LIVER","SAT","VAT"))
  if (all(tv == "PBMC"))                            "a_PBMC"
  else if (length(solid) == 1 && !("PBMC" %in% tv)) "b_One_solid_tissue"
  else if (length(solid) >= 2 && !("PBMC" %in% tv)) "c_Shared_between_solid_tissues"
  else                                              "d_Other"
}
clone_class <- all_tcells %>%
  filter(CTaa %in% expanded_clones$CTaa) %>%
  distinct(CTaa, Tissue) %>%
  group_by(CTaa) %>%
  summarise(class = classify(unique(Tissue)), .groups = "drop")

class_tab <- clone_class %>%
  left_join(clone_counts, by = "CTaa") %>%
  group_by(class) %>%
  summarise(n_clonotypes = n(), n_cells = sum(n_cells), .groups = "drop")
cat("\n== Expanded-clone tissue footprint (cohort-wide) ==\n"); print(as.data.frame(class_tab), row.names = FALSE)
write.csv(class_tab, file.path(tab_dir, "Fig3_shared_class_summary_v7.csv"), row.names = FALSE)

# ---- annotate the V7 T-cell object and save the shared subset ---------------
obj <- readRDS(tcell_rds)
obj$Stage <- factor(unname(stage_map[as.character(obj$Patient_ID)]), levels = stage_lev)

obj$expanded_shared_class <- NA_character_
ex <- which(obj$Clonality == "Expanded")
obj$expanded_shared_class[ex] <- clone_class$class[match(obj$CTaa_tcr[ex], clone_class$CTaa)]
cat("\n== expanded_shared_class on V7 T cells ==\n")
print(table(obj$expanded_shared_class, useNA = "ifany"))

keep_classes <- c("a_PBMC", "b_One_solid_tissue", "c_Shared_between_solid_tissues")
obj_sub <- subset(obj, cells = colnames(obj)[obj$expanded_shared_class %in% keep_classes])
# saved with the ORIGINAL class labels (a_/b_/c_) - Fig 3 d-i recode to friendly
# names in-memory, and the SCENIC scripts (14/17) key off these original labels.
obj_sub$expanded_shared_class <- factor(obj_sub$expanded_shared_class, levels = keep_classes)
cat("\n== shared-class T cells kept (a/b/c) ==\n"); print(table(obj_sub$expanded_shared_class))

saveRDS(obj_sub, obj_out)
cat("\nSaved:", obj_out, "(", ncol(obj_sub), "cells )\n")
# End of the script
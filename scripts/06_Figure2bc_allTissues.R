###############################################################################
# Figure 2b,c - TCR repertoire diversity across tissues and fibrosis stages
#
# Fig 2b : percentage of unique clonotypes (unique CDR3aa / T cells) per fibrosis
#          stage, faceted by tissue.
# Fig 2c : TCR diversity (Shannon and inverse Simpson) per tissue.
#
# Input : demultiplexed gene-expression object and per-capture TCR contigs
# Output: repertoire tables and Figures 2b, 2c (PDF + SVG + TIFF)
###############################################################################

.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))   # scRepertoire library location
library(scRepertoire)
library(patchwork)
library(tidyverse)
library(Seurat)                        # for the GEX barcode -> tissue map
source("00_nature_style.R")            # theme_nature(), stage_cols, ...

# ---- v7 data locations ------------------------------------------------------
tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"           # one dir per capture: GC-WL-<cap>-TCR/outs/
gex_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"   # demuxed GEX -> barcode->tissue map
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"   # figure panels
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"    # data tables
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---- theme + saver (kept from your scripts) --------------------------------
my_theme <- theme_nature(7) +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        strip.text = element_text(size = 8, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey92", color = NA),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7), legend.key.size = unit(0.7, "lines"))
save_figure <- function(plot, filename, width_cm = 18, height_cm = 6, path = fig_dir) {
  grDevices::cairo_pdf(file.path(path, paste0(filename, ".pdf")), width = width_cm/2.54, height = height_cm/2.54)
  print(plot); dev.off()
  svglite::svglite(file.path(path, paste0(filename, ".svg")), width = width_cm/2.54, height = height_cm/2.54)
  print(plot); dev.off()
  tiff(file.path(path, paste0(filename, ".tiff")), width = width_cm, height = height_cm, units = "cm", res = 600, compression = "lzw")
  print(plot); dev.off()
  message("Saved: ", filename)
}

# ---- CONFIG: patients per tissue + stage --------------------------------------
liver_pids <- c("10742","9961","9999","10113-2","9680","10203","10113-1","10380",
                "10291-2","10202","10634","10738","10205","9991","9932","11040",
                "11051","11183","11471","11327","10291-1")            # 21
svp_pids   <- setdiff(liver_pids, "10291-1")                          # 20

# Stage per patient (F0 -> No fibrosis, F1-F3 -> Fibrosis, HC -> Healthy)
stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")
disp <- c("Healthy"="Healthy","No_fibrosis"="No Fibrosis","Fibrosis"="Fibrosis")
lev  <- c("Healthy","No Fibrosis","Fibrosis")

# ---- GEX barcode -> tissue map (from the HTO-demuxed object) ----------------
obj <- readRDS(gex_rds)
gex_map <- data.frame(Sample = obj$Sample, Tissue = obj$Tissue,
                      bc = sub("^.*_", "", rownames(obj@meta.data)),
                      stringsAsFactors = FALSE)
rm(obj); gc()

# capture name for a (tissue, patient): liver = -LIVER, others share one -SVP VDJ library
capture_of <- function(tissue, pid)
  if (tissue == "LIVER") paste0("GC-WL-", pid, "-LIVER") else paste0("GC-WL-", pid, "-SVP")

tcr_file_for <- function(capture) {
  cands <- c(file.path(tcr_dir, paste0(capture, "-TCR"), "outs", "filtered_contig_annotations.csv"),
             file.path(tcr_dir, paste0(capture, "-TCR"),          "filtered_contig_annotations.csv"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else stop("No v7 TCR file for capture: ", capture)
}

# read a capture's contigs and keep only TCR cells present in the pre-QC demuxed
# GEX object for this tissue - applied to ALL tissues (liver -> 12,454), so the
# inclusion rule is identical everywhere and matches Fig 2a.
read_contig <- function(tissue, pid) {
  cap  <- capture_of(tissue, pid)
  d    <- read.csv(tcr_file_for(cap), stringsAsFactors = FALSE)
  keep <- gex_map$bc[gex_map$Sample == cap & gex_map$Tissue == tissue]
  d[d$barcode %in% keep, , drop = FALSE]
}

# build combined.TCR for one tissue (v7, demuxed per patient)
build_tcr <- function(tissue, pids) {
  cl <- lapply(pids, function(p) read_contig(tissue, p))
  ct <- combineTCR(cl, samples = pids, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
  ct <- addVariable(ct, variable.name = "Type", variables = unname(stage_map[pids]))
  lapply(ct, function(df) { df$Type <- factor(unname(disp[as.character(df$Type)]), levels = lev); df })
}

tissues <- list(LIVER = liver_pids, PBMC = svp_pids, SAT = svp_pids, VAT = svp_pids)
cts     <- Map(build_tcr, names(tissues), tissues)

stage_fill <- scale_fill_manual(name = "Stage",
                                values = c("Healthy" = unname(stage_cols["Healthy"]),
                                           "No Fibrosis" = unname(stage_cols["No_fibrosis"]),
                                           "Fibrosis" = unname(stage_cols["Fibrosis"])),
                                breaks = lev)

# ============================ Fig 2b =========================================
# % unique clonotypes per fibrosis stage via scRepertoire::clonalQuant
# (cloneCall = "aa", scale = TRUE) - the validated default. group.by = "Type"
# pools all of a stage's cells, so clonal expansion shows (expanded samples pull
# it down).
#
# Layout matches Fig 2a EXACTLY: a single facet_wrap(~Tissue) plot, so the four
# tissue columns are EQUAL width and each carries a grey strip box with the
# tissue name (same style as 2a / 2g / 2h). We reuse scRepertoire's own computed
# values (each tissue's clonalQuant plot data) and pool them into one data frame,
# then facet - no re-implementation of the clonalQuant maths.
# same labels + strip style as Fig 2a (non-bold "Liver T cells ...", grey box)
tissue_lab_b <- as_labeller(c(LIVER = "Liver T cells", PBMC = "PBMC T cells",
                              SAT = "SAT T cells", VAT = "VAT T cells"))

# Pull the per-stage unique-clonotype % straight from scRepertoire (exportTable),
# one clean row per stage per tissue, then draw a single faceted bar plot. This
# avoids reusing clonalQuant's own ggplot (whose stat_summary layer is tied to
# its original data and cannot simply be re-faceted).
quant_one <- function(ct, tis) {
  d <- as.data.frame(clonalQuant(ct, cloneCall = "aa", chain = "both",
                                 group.by = "Type", scale = TRUE, exportTable = TRUE))
  num <- vapply(d, is.numeric, logical(1))
  grp <- if ("Type"   %in% names(d)) "Type"   else if ("values" %in% names(d)) "values" else names(d)[!num][1]
  val <- if ("scaled" %in% names(d)) "scaled" else names(d)[num][sum(num)]     # scaled % (last numeric)
  data.frame(Type = as.character(d[[grp]]), value = as.numeric(d[[val]]),
             Tissue = tis, stringsAsFactors = FALSE)
}
quant_df <- bind_rows(Map(quant_one, cts, names(cts)))
quant_df$Type   <- factor(quant_df$Type,   levels = lev)
quant_df$Tissue <- factor(quant_df$Tissue, levels = names(cts))

fig2b <- ggplot(quant_df, aes(x = Type, y = value, fill = Type)) +
  geom_col(width = 0.65) +
  facet_wrap(~ Tissue, nrow = 1, labeller = tissue_lab_b) +
  stage_fill +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = NULL, x = NULL, y = "Unique clonotypes (% of T cells)") +
  theme_nature(7) +                                   # identical strip style to Fig 2a
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
save_figure(fig2b, "Fig2b_uniqueClones_V7", width_cm = 12, height_cm = 5)   # same footprint as Fig 2a

# ============================ Fig 2c =========================================
# Shannon AND inverse Simpson diversity per tissue, INTERLEAVED side by side:
#   Liver Shannon | Liver Inv.Simpson | PBMC Shannon | PBMC Inv.Simpson | ...
# Each panel is scRepertoire's default clonalDiversity plot for ONE metric (the
# same call the previous version used, just rearranged). Every panel keeps its
# own y-scale (Shannon and inv.Simpson are on different ranges); the tissue and
# metric are named in the panel title, stage shown by fill + one shared legend.
mk_div <- function(ct, tis, metric, mlab) {
  # title = tissue only; the metric (Shannon / inv.Simpson) is already shown in
  # the clonalDiversity panel itself, so naming it here would be redundant.
  clonalDiversity(ct, cloneCall = "aa", group.by = "Type", metric = metric) +
    stage_fill +
    labs(title = tis, x = NULL, y = NULL) +
    my_theme +
    theme(plot.title = element_text(size = 7, face = "plain", hjust = 0.5),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

panels_c <- list()
for (i in seq_along(cts)) {
  tis <- names(cts)[i]; ct <- cts[[i]]
  panels_c <- c(panels_c, list(
    mk_div(ct, tis, "shannon",     "Shannon"),
    mk_div(ct, tis, "inv.simpson", "Inv. Simpson")))
}
fig2c <- wrap_plots(panels_c, nrow = 1, guides = "collect") & theme(legend.position = "right")
save_figure(fig2c, "Fig2c_shannon_invSimpson_perTissue_V7", width_cm = 26, height_cm = 6)

# End of the script
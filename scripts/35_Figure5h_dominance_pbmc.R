###############################################################################
# 35 - Figure 5h - PBMC FNA TCR top-100 clonal dominance (baseline vs follow-up)
#
# Same analysis as Fig 5d (liver, script 34), for the paired PBMC follow-up FNAs.
#
# ---------------------------------------------------------------------------
# BASELINE PBMC SOURCE (unified with scripts 38/39/4_Build) -- IMPORTANT:
# Earlier this script read the raw v7_tcr_TR contigs and, when a donor had no
# standalone -PBMC capture at baseline, fell back to the MULTIPLEXED SVP capture
# (SAT-VAT-PBMC pooled). That mixed adipose T-cell clones into the "PBMC"
# repertoire for 10738/11040/11183, so their baseline top-100 dominance was
# computed on a contaminated pool, while 10291-1 (deep standalone -PBMC library)
# was not - i.e. the four donors were not on the same PBMC population.
#
# This version derives BOTH baseline and follow-up PBMC clonotypes from the
# QC'd FNA PBMC object (SeuObjFNA_PBMCV7.rds) - the SAME demultiplexed,
# tissue-QC'd PBMC cells used by scripts 38, 39 and 4_Build. CTaa is attached
# per cell with the identical robust Sample + 16 bp-core match (incl. the
# GC-WL-11303-PBMC -> GC-WL-11303-SVP alias that recovers 10291-1's follow-up).
# Top-100 dominance is then the top-100 CTaa clones' share of the clonotype+
# T cells in each Patient x Visit. Clone identity is CTaa (matching 38/39),
# not CTstrict as before.
# ---------------------------------------------------------------------------
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(Seurat); library(scRepertoire); library(dplyr); library(tidyr); library(ggplot2)
source("00_nature_style.R")

tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
rds_pbmc <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_PBMCV7.rds"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

core16 <- function(x) sub(".*?([ACGT]{16}).*", "\\1", x)

# ---- load PBMC FNA object + attach CTaa (verbatim pattern from script 39) ----
if (!file.exists(rds_pbmc)) stop("PBMC FNA object not found: ", rds_pbmc)
seu <- readRDS(rds_pbmc)

# Some captures name the TCR library differently from the GEX capture (obj$Sample).
# Map TCR-capture base name -> GEX Sample so CTaa still attaches. e.g. patient
# 10291-1 follow-up PBMC: GEX Sample = GC-WL-11303-SVP but its TCR library is
# GC-WL-11303-PBMC-TCR (cellranger --id GC-WL-11303-PBMC).
sample_alias <- c("GC-WL-11303-PBMC" = "GC-WL-11303-SVP")
if (!"CTaa" %in% colnames(seu@meta.data)) {
  caps     <- basename(dirname(dirname(
    list.files(tcr_dir, "filtered_contig_annotations.csv", recursive = TRUE, full.names = TRUE))))
  cap_base <- sub("-TCR$", "", caps)
  cap_samp <- ifelse(cap_base %in% names(sample_alias), sample_alias[cap_base], cap_base)
  take     <- which(cap_samp %in% unique(seu$Sample) & !duplicated(cap_samp))
  use_caps <- caps[take]; use_samp <- cap_samp[take]
  cat("PBMC TCR capture -> Sample:\n"); print(setNames(use_samp, use_caps))
  if (!length(use_caps)) stop("no matching TCR captures for PBMC FNA object")
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
if (!"Patient_ID" %in% colnames(seu@meta.data)) stop("PBMC FNA object has no Patient_ID")

# diagnostic: clonotype recovery per Patient x Visit (spot a dropped follow-up)
cat("  CTaa recovery (TRUE = cell has a clonotype), Patient_ID x Visit:\n")
print(table(Patient = seu$Patient_ID, Visit = seu$Visit, hasCTaa = !is.na(seu$CTaa)))

# ---- per Patient x Visit top-100 dominance over clonotype+ T cells -----------
md <- seu@meta.data[!is.na(seu$CTaa), c("Patient_ID", "Visit", "CTaa")]
md$Patient_ID <- as.character(md$Patient_ID)

top100 <- do.call(rbind, lapply(split(md, list(md$Patient_ID, md$Visit), drop = TRUE), function(s) {
  f <- sort(table(s$CTaa), decreasing = TRUE)
  data.frame(Donor = s$Patient_ID[1],
             Timepoint = as.character(s$Visit[1]),
             Tissue = "PBMC",
             n_cells = sum(f),
             Top100_Percent = round(100 * sum(f[1:min(100, length(f))]) / sum(f), 2),
             stringsAsFactors = FALSE)
}))
rownames(top100) <- NULL
top100$Timepoint <- factor(top100$Timepoint, levels = c("Baseline", "Followup"))
cat("\n== Top-100 dominance (PBMC, FNA object) ==\n"); print(top100, row.names = FALSE)
write.csv(top100, file.path(tab_dir, "Fig5h_top100_dominance_pbmc_v7.csv"), row.names = FALSE)

wide <- top100 %>% select(Donor, Timepoint, Top100_Percent) %>%
  pivot_wider(names_from = Timepoint, values_from = Top100_Percent) %>%
  filter(!is.na(Baseline) & !is.na(Followup))
if (nrow(wide) > 1) cat(sprintf("\nPaired Wilcoxon (PBMC): p = %.4f | mean change = %.2f%%\n",
    wilcox.test(wide$Baseline, wide$Followup, paired = TRUE)$p.value,
    mean(wide$Followup - wide$Baseline)))

fig5h <- ggplot(top100, aes(Timepoint, Top100_Percent, group = Donor, colour = Donor)) +
  geom_line(linewidth = 0.8, alpha = 0.85) + geom_point(size = 2.5) + facet_wrap(~ Tissue) +
  scale_colour_manual(values = npg_pal) + scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60)) +  # shared scale liver+PBMC
  labs(title = "Top-100 TCR clonal dominance - PBMC FNA", subtitle = "Paired baseline vs follow-up",
       x = NULL, y = "Top-100 clones (% of repertoire)", colour = "Donor") +
  theme_nature(7) + theme(legend.position = "bottom",
                          strip.background = element_rect(fill = "grey92", colour = NA))
save_nature(fig5h, file.path(fig_dir, "Fig5h_FNA_top100_dominance_pbmc_v7"), 6, 6)
cat("\nDone Fig 5h (PBMC, unified on FNA object).\n")
# End of the script

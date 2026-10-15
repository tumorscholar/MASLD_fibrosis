###############################################################################
# Figure 2a - T-cell fraction across tissues and fibrosis stages
#
# Quantifies the proportion of sequenced cells that are T cells (carry a
# productive TCR), per patient and tissue, and plots the mean per fibrosis
# stage across the four tissues.
#
# Each cell inherits its tissue label from HTO demultiplexing; a cell is scored
# as a T cell when its barcode has a productive TRA/TRB contig.
#   numerator   = cells with a productive TCR   (per patient and tissue)
#   denominator = qc passed cells               (per patient and tissue)
#   value       = 100 * numerator / denominator
#
# Input : demultiplexed gene-expression object and TCR contig annotations
# Output: T-cell fraction table and Figure 2a (PDF + TIFF)
###############################################################################

library(Seurat)
library(tidyverse)
source("00_nature_style.R")            # theme_nature(), stage_cols, save_nature()


# ============================ CONFIG =========================================
gex_rds <- "/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"
tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"   # <<< your v7 TCR folder
out_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"   # figure panels
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"    # data tables
productive_only <- TRUE                # count cells with a productive TRA/TRB

# Given a GEX capture name (obj$Sample) return its filtered_contig_annotations.csv.
# EDIT the patterns to match your v7_tcr directory naming, then check the
# "Captures with TCR" report the script prints.
tcr_file_for <- function(sample) {
  cands <- c(
    file.path(tcr_dir, paste0(sample, "-TCR"), "outs", "filtered_contig_annotations.csv"),
    file.path(tcr_dir, paste0(sample, "-TCR"),          "filtered_contig_annotations.csv"),
    file.path(tcr_dir, sample, "outs",                  "filtered_contig_annotations.csv"),
    file.path(tcr_dir, sample,                          "filtered_contig_annotations.csv")
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

# Stage per patient 
stage_map <- c(
  "9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
  "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
  "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
  "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
  "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
# =============================================================================

# ---- GEX cells: barcode, capture, patient, tissue ---------------------------
obj <- readRDS(gex_rds)
# Cell names are "<Sample>_<barcode>" (e.g. GC-WL-...-LIVER_AAACCTGAGAGGGCTT-1);
# take the 10x barcode after the last underscore so it matches the TCR barcode.
gex <- obj@meta.data %>%
  mutate(bc = sub("^.*_", "", rownames(.))) %>%
  transmute(Sample, Patient_ID, Tissue, bc)

# ---- TCR: set of T-cell barcodes per capture --------------------------------
captures <- unique(gex$Sample)
tcr_bc <- lapply(captures, function(s) {
  f <- tcr_file_for(s)
  if (is.na(f)) { message("No TCR file for capture: ", s); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (productive_only && "productive" %in% names(d))
    d <- d[tolower(as.character(d$productive)) %in% c("true", "1"), , drop = FALSE]
  if ("chain" %in% names(d)) d <- d[d$chain %in% c("TRA", "TRB"), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  data.frame(Sample = s, bc = unique(d$barcode), stringsAsFactors = FALSE)
}) %>% bind_rows()
cat("Captures with a TCR file:", length(unique(tcr_bc$Sample)), "of", length(captures), "\n")

# ---- Flag GEX cells that have a TCR (match within capture) ------------------
gex$hasTCR <- paste(gex$Sample, gex$bc) %in% paste(tcr_bc$Sample, tcr_bc$bc)

# ---- Per patient x tissue: T-cell fraction ----------------------------------
frac <- gex %>%
  group_by(Patient_ID, Tissue) %>%
  summarise(nTotal = n(), nT = sum(hasTCR), .groups = "drop") %>%
  mutate(pct    = 100 * nT / nTotal,
         Stage  = factor(stage_map[Patient_ID], levels = c("Healthy","No_fibrosis","Fibrosis")),
         Tissue = factor(Tissue, levels = c("LIVER","PBMC","SAT","VAT")))
write.csv(frac, file.path(tab_dir, "Fig2a_Tcell_fraction_perPatientTissue_V7.csv"), row.names = FALSE)

# ---- Per-patient x tissue: T cells vs total sequenced (printed BEFORE plotting) ----
cat("\n=== Per patient x tissue (T cells / total sequenced / %) ===\n")
frac %>%
  arrange(Tissue, Patient_ID) %>%
  select(Patient_ID, Tissue, Stage, nT, nTotal, pct) %>%
  print(n = Inf)

# ---- Totals per tissue: T cells vs total sequenced (printed + saved BEFORE plotting) ----
tissue_totals <- frac %>%
  group_by(Tissue) %>%
  summarise(total_Tcells    = sum(nT),
            total_sequenced = sum(nTotal),
            pct_overall     = round(100 * sum(nT) / sum(nTotal), 1),
            n_patients      = n(), .groups = "drop")
cat("\n=== Totals per tissue (T cells / total sequenced cells) ===\n")
print(tissue_totals)
cat("\nGrand total: T cells =", sum(frac$nT), "| sequenced cells =", sum(frac$nTotal), "\n")
write.csv(tissue_totals, file.path(tab_dir, "Fig2a_tissue_totals_V7.csv"), row.names = FALSE)

# ---- Plot: mean bar (coloured by Stage), per fibrosis stage ----
if (!exists("stage_pts"))
  stage_pts <- c("Healthy" = "#9E82C4", "No_fibrosis" = "#6E2E8A", "Fibrosis" = "#8E1044")
tissue_lab <- as_labeller(c(LIVER = "Liver T cells", PBMC = "PBMC T cells",
                            SAT = "SAT T cells", VAT = "VAT T cells"))

p <- ggplot(frac, aes(x = Stage, y = pct)) +
  stat_summary(aes(fill = Stage), fun = mean, geom = "col", width = 0.65) +   # no bar outline
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_discrete(labels = c("Healthy", "No fibrosis", "Fibrosis")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ Tissue, nrow = 1, labeller = tissue_lab) +
  labs(x = NULL, y = "T cells % compared to\ntotal cell") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_nature(p, file.path(out_dir, "Fig2a_Tcell_fraction_V7"), w_cm = 12, h_cm = 5)

# ---- per-sample dot version (mean bar + SEM error bar + each patient as a point) --
# Saved as a SEPARATE file; the no-dot version above is kept unchanged.
p_dots <- p +
  stat_summary(data = frac, aes(x = Stage, y = pct), fun.data = mean_se,
               geom = "errorbar", width = 0.3, linewidth = 0.4, colour = "black",
               inherit.aes = FALSE) +
  geom_jitter(data = frac, aes(x = Stage, y = pct),
              width = 0.15, height = 0, size = 0.5, colour = "black", alpha = 0.75,
              inherit.aes = FALSE)
save_nature(p_dots, file.path(out_dir, "Fig2a_Tcell_fraction_dots_V7"), w_cm = 12, h_cm = 5)

# End of the script
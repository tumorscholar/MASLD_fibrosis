###############################################################################
# Figure 2g,h - CD8 T-cell phenotypes within expanded clones
#
# Fig 2g : CD8+ CD69+ CXCR6+ T cells
# Fig 2h : CD8+ CD69+ CD103- T cells        (CD103 = ITGAE, negative)
#
#   - marker positive = log-normalised RNA expression > 0.25 (CD8 = CD8A)
#   - value           = percentage of EXPANDED T cells (per patient and tissue)
#                       matching each phenotype
#   - bar             = mean across patients, per fibrosis stage
#   - patients with fewer than 10 expanded cells in a tissue are excluded, as a
#     percentage from very few cells is unstable
#
# Input : T-cell Seurat object with clonal-expansion annotation
# Output: per-patient and per-tissue tables and Figures 2g, 2h (PDF + TIFF)
###############################################################################

library(Seurat)
library(tidyverse)
source("00_nature_style.R")            # theme_nature(), stage_cols, stage_pts, save_nature()

# ---- CONFIG (v7) ------------------------------------------------------------
tcell_rds <- "/data/Blizard-AlazawiLab/rk/seurat/TcellObjV7.rds"
out_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure2"   # figure panels
tab_dir   <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure2"    # data tables
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

stage_map <- c("9680"="Fibrosis","9932"="Fibrosis","9961"="No_fibrosis","9991"="Fibrosis",
               "9999"="No_fibrosis","10113-1"="Fibrosis","10113-2"="No_fibrosis","10202"="Fibrosis",
               "10203"="Fibrosis","10205"="Fibrosis","10291-1"="Fibrosis","10291-2"="Fibrosis",
               "10380"="Fibrosis","10634"="Fibrosis","10738"="Fibrosis","10742"="Healthy",
               "11040"="Fibrosis","11051"="Fibrosis","11183"="Fibrosis","11327"="Fibrosis","11471"="Fibrosis")

# ---- load T cells; gate markers on RNA (log-normalised > 0.25) --------------
obj <- readRDS(tcell_rds)
DefaultAssay(obj) <- "RNA"

need <- c("CD8A", "CD69", "CXCR6", "ITGAE")            # CD8A = CD8; ITGAE = CD103
missing <- need[!need %in% rownames(obj[["RNA"]])]
if (length(missing)) stop("Missing genes in RNA assay: ", paste(missing, collapse = ", "))

# Use the LOG-NORMALISED data layer so the 0.25 threshold is a real cutoff.
# (On raw counts, >0.25 == >0, i.e. any single detected molecule counts as
# positive - that made almost every expanded cell "positive" and inflated the %.)
expr <- FetchData(obj, vars = need, layer = "data")    # cells x genes, log-normalised
expr <- expr[colnames(obj), , drop = FALSE]            # keep cell order aligned

cd8   <- expr[, "CD8A"]  > 0.25
cd69  <- expr[, "CD69"]  > 0.25
cxcr6 <- expr[, "CXCR6"] > 0.25
itgae <- expr[, "ITGAE"] > 0.25

obj$is_autoaggr <- cd8 & cd69 &  cxcr6                  # Fig 2g: CD8+ CD69+ CXCR6+
obj$is_trm      <- cd8 & cd69 & !itgae                  # Fig 2h: CD8+ CD69+ CD103- (ITGAE-)

# ---- sanity read-out (double-check the numbers) -----------------------------
cat("\n== marker positivity (RNA data > 0.25), across ALL T cells ==\n")
print(round(100 * colMeans(expr > 0.25, na.rm = TRUE), 2))
cat("CD8+CD69+CXCR6+ :", round(100 * mean(obj$is_autoaggr, na.rm = TRUE), 2),
    "%  |  CD8+CD69+CD103- :", round(100 * mean(obj$is_trm, na.rm = TRUE), 2),
    "% of all T cells\n")

# ---- Stage (from patient) + keep EXPANDED T cells only ----------------------
md <- obj@meta.data %>%
  mutate(Patient_ID = as.character(Patient_ID),
         Stage  = factor(stage_map[Patient_ID], levels = c("Healthy","No_fibrosis","Fibrosis")),
         Tissue = factor(Tissue, levels = c("LIVER","PBMC","SAT","VAT"))) %>%
  filter(Clonality == "Expanded")

cat("Expanded T cells used:", nrow(md), "\n")

# ---- per patient x tissue: % of expanded that are autoaggressive / Trm ------
per_pt <- md %>%
  group_by(Patient_ID, Tissue) %>%
  summarise(nExpanded = n(),
            nAuto     = sum(is_autoaggr),
            nTrm      = sum(is_trm),
            .groups   = "drop") %>%
  mutate(pct_auto = 100 * nAuto / nExpanded,
         pct_trm  = 100 * nTrm  / nExpanded,
         Stage    = factor(stage_map[Patient_ID], levels = c("Healthy","No_fibrosis","Fibrosis")))
write.csv(per_pt, file.path(tab_dir, "Fig2gh_perPatientTissue_V7.csv"), row.names = FALSE)

# ---- print numerators / denominators BEFORE plotting ------------------------
cat("\n=== Per patient x tissue (expanded T cells | autoaggr | Trm | %) ===\n")
per_pt %>% arrange(Tissue, Patient_ID) %>%
  select(Patient_ID, Tissue, Stage, nExpanded, nAuto, pct_auto, nTrm, pct_trm) %>%
  print(n = Inf)

# ---- drop patients with too few expanded cells (statistical, not result-driven);
#      a % from <10 cells is noise (1 expanded cell reads 0% or 100%). ----------
MIN_EXP_CELLS <- 10
dropped <- per_pt %>% filter(nExpanded < MIN_EXP_CELLS)
cat("\nPatients dropped (<", MIN_EXP_CELLS, "expanded cells in a tissue):", nrow(dropped), "\n")
if (nrow(dropped)) print(as.data.frame(dropped[, c("Patient_ID","Tissue","nExpanded")]), row.names = FALSE)
per_pt_plot <- per_pt %>% filter(nExpanded >= MIN_EXP_CELLS)

# ---- Totals per tissue x stage = MEAN of per-patient % (+/- SEM) = your Prism -
tissue_stage_tot <- per_pt_plot %>%
  group_by(Tissue, Stage) %>%
  summarise(n_patients = n(),
            CD8_CD69_CXCR6_mean = round(mean(pct_auto), 2), CD8_CD69_CXCR6_sem = round(sd(pct_auto)/sqrt(n()), 2),
            CD8_CD69_CD103neg_mean = round(mean(pct_trm), 2), CD8_CD69_CD103neg_sem = round(sd(pct_trm)/sqrt(n()), 2),
            .groups = "drop")
cat("\n=== Totals per tissue x stage (mean of patients +/- SEM) ===\n")
print(as.data.frame(tissue_stage_tot), row.names = FALSE)
write.csv(tissue_stage_tot, file.path(tab_dir, "Fig2gh_tissueStage_totals_V7.csv"), row.names = FALSE)

# ---- plot helper (same look as Fig 2a) --------------------------------------
tissue_lab_g <- as_labeller(c(LIVER="Liver", PBMC="PBMC", SAT="SAT", VAT="VAT"))
mk_panel <- function(df, yvar, ylab) {
  ggplot(df, aes(x = Stage, y = .data[[yvar]])) +
    stat_summary(aes(fill = Stage), fun = mean, geom = "col", width = 0.65) +   # no outline
    scale_fill_manual(values = stage_cols, guide = "none") +
    scale_x_discrete(labels = c("Healthy", "No fibrosis", "Fibrosis")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~ Tissue, nrow = 1, labeller = tissue_lab_g) +
    labs(x = NULL, y = ylab) +
    theme_nature(7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_g <- mk_panel(per_pt_plot, "pct_auto", "CD8⁺CD69⁺CXCR6⁺ T cells\n(% of expanded T cells)")
p_h <- mk_panel(per_pt_plot, "pct_trm",  "CD8⁺CD69⁺CD103⁻ T cells\n(% of expanded T cells)")

save_nature(p_g, file.path(out_dir, "Fig2g_CD8_CD69_CXCR6_expanded_V7"),    w_cm = 12, h_cm = 5)
save_nature(p_h, file.path(out_dir, "Fig2h_CD8_CD69_CD103neg_expanded_V7"), w_cm = 12, h_cm = 5)

# ---- per-sample dot versions (each patient = one point on top of the mean bar) --
# Saved as SEPARATE files; the no-dot versions above are kept unchanged. Points
# use per_pt_plot (the same >=10-expanded-cell patients that build the bars).
add_dots <- function(base, df, yvar)
  base +
  stat_summary(data = df, aes(x = Stage, y = .data[[yvar]]), fun.data = mean_se,
               geom = "errorbar", width = 0.3, linewidth = 0.4, colour = "black",
               inherit.aes = FALSE) +
  geom_jitter(data = df, aes(x = Stage, y = .data[[yvar]]),
              width = 0.15, height = 0, size = 0.5, colour = "black", alpha = 0.75,
              inherit.aes = FALSE)
save_nature(add_dots(p_g, per_pt_plot, "pct_auto"),
            file.path(out_dir, "Fig2g_CD8_CD69_CXCR6_expanded_dots_V7"),    w_cm = 12, h_cm = 5)
save_nature(add_dots(p_h, per_pt_plot, "pct_trm"),
            file.path(out_dir, "Fig2h_CD8_CD69_CD103neg_expanded_dots_V7"), w_cm = 12, h_cm = 5)

# End of the script
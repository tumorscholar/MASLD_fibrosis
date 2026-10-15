###############################################################################
# 32 - Figure 5b - FNA PCA, baseline -> follow-up, clinically annotated
#      LIVER and PBMC (one clean panel each + a combined two-panel figure).
#
# What moves the samples after surgery? Per-sample (patient x visit) centroid in
# each tissue's own PC space (same construction as the quick PCA in
# 30_Figure5_liverFNA.R / 31_Figure5_pbmcFNA.R). PC1 vs PC2, coloured by a
# clinical variable, with a paired baseline -> follow-up arrow per patient.
#
# Kept honest (n = 8 per tissue): the PC1 <-> clinical correlations are computed
# per tissue and shown in each subtitle, so the reader sees what PC1 tracks. Any
# sample that is a clear PC2 outlier (robust |z| > 2.5) is flagged automatically
# rather than interpreted - the axis story is not hardcoded to one tissue.
#
# Inputs : SeuObjFNA_LIVERV7.rds , SeuObjFNA_PBMCV7.rds  (Patient_ID / Visit +
#          the "pca" reduction; written by scripts 30 / 31)
# Output : Fig5b_PCA_clinical_LIVER_v7 , Fig5b_PCA_clinical_PBMC_v7 ,
#          Fig5b_PCA_clinical_bothTissues_v7  (pdf+tiff) + centroid CSVs
###############################################################################
suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(ggrepel); library(patchwork) })
source("00_nature_style.R")

COLOUR_BY <- "BMI"     # <- switch to "CAP" (liver-specific steatosis), "Weight" or "Stiffness";
#    the colourbar re-labels automatically.
lab_of <- c(BMI = "BMI (kg/m2)", CAP = "CAP (dB/m)", Weight = "Weight (kg)",
            Stiffness = "Liver stiffness (kPa)")

tissues <- c(LIVER = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVERV7.rds",
             PBMC  = "/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_PBMCV7.rds")
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- clinical metadata per patient x visit (Table 1 / Table 2) --------------
clin <- read.csv(text = "
Patient_ID,Visit,BaselineStage,Diabetes,CAP,Stiffness,BMI,Weight
10291-1,Baseline,F1,Yes,400,21.5,61.4,199
10291-1,Followup,F1,No,376,6.6,55.6,180
10738,Baseline,F1,No,388,6.6,37.4,128
10738,Followup,F1,No,247,3.3,28.3,97
11040,Baseline,F3,Yes,384,19.3,46.7,135
11040,Followup,F3,No,228,9.3,30.1,87
11183,Baseline,F3,No,400,46.5,57.5,170
11183,Followup,F3,No,230,6.5,32.4,96
", stringsAsFactors = FALSE, strip.white = TRUE)

cor_vars <- c("BMI","Weight","CAP","Stiffness")

# ---- build one tissue's clinically-annotated centroid PCA panel -------------
make_panel <- function(TISSUE, rds) {
  cat("\n=====================", TISSUE, "=====================\n")
  if (!file.exists(rds)) { cat("  object not found:", rds, "- skipping.\n"); return(NULL) }
  obj <- readRDS(rds)
  
  pc <- as.data.frame(Embeddings(obj, "pca")[, 1:2]); colnames(pc) <- c("PC1","PC2")
  pc$Patient_ID <- obj$Patient_ID; pc$Visit <- as.character(obj$Visit)
  cent <- pc %>% group_by(Patient_ID, Visit) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  
  df <- cent %>% left_join(clin, by = c("Patient_ID","Visit")) %>%
    mutate(Visit = factor(Visit, levels = c("Baseline","Followup")), Tissue = TISSUE)
  write.csv(df, file.path(tab_dir, sprintf("Fig5b_PCA_centroids_clinical_%s_v7.csv", TISSUE)),
            row.names = FALSE)
  
  # PC1 <-> clinical correlations (descriptive)
  cors <- sapply(cor_vars, function(v) suppressWarnings(cor(df$PC1, df[[v]])))
  cat("  PC1 vs clinical (Pearson r, n =", nrow(df), "):\n"); print(round(cors, 2))
  cat(sprintf("  PC1 baseline mean %.2f  ->  follow-up mean %.2f\n",
              mean(df$PC1[df$Visit == "Baseline"]), mean(df$PC1[df$Visit == "Followup"])))
  sub <- paste0("PC1 vs clinical (n=", nrow(df), ", descriptive): ",
                paste(sprintf("%s r=%.2f", cor_vars, cors), collapse = "  |  "))
  
  # automatic outlier flag (robust z on EITHER PC1 or PC2, so a sample that
  # stretches either axis - e.g. 10738 on PC2 in liver, 11040 on PC1 in PBMC -
  # is flagged and not over-interpreted)
  rz <- function(v) { m <- median(v); s <- mad(v); if (s == 0) s <- sd(v)
  if (is.na(s) || s == 0) rep(0, length(v)) else abs(v - m) / s }
  df$is_out <- rz(df$PC1) > 2.5 | rz(df$PC2) > 2.5
  outs <- df[df$is_out, , drop = FALSE]
  
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_path(aes(group = Patient_ID), colour = "grey65", linewidth = 0.5,   # baseline -> follow-up
              arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    geom_point(aes(colour = .data[[COLOUR_BY]], shape = Visit), size = 4, stroke = 0.4) +
    scale_shape_manual(values = c(Baseline = 16, Followup = 17)) +
    scale_colour_viridis_c(name = lab_of[[COLOUR_BY]]) +
    scale_x_continuous(expand = expansion(mult = 0.18)) +   # headroom so edge labels/arrows aren't clipped
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    coord_cartesian(clip = "off") +
    geom_text_repel(aes(label = Patient_ID), size = 2.3, max.overlaps = Inf,
                    box.padding = 0.4, segment.size = 0.2, segment.colour = "grey60") +
    labs(title = sprintf("%s FNA - PCA, baseline -> follow-up", TISSUE),
         subtitle = sub, x = "PC1", y = "PC2") +
    theme_nature(7) +
    theme(plot.subtitle = element_text(size = 5.0, colour = "grey35"),
          legend.key.size = unit(0.35, "cm"))
  if (nrow(outs))
    p <- p + geom_text(data = outs, aes(PC1, PC2, label = "outlier"),
                       size = 1.9, colour = "grey45", hjust = -0.2, vjust = 1.4)
  
  save_nature(p, file.path(fig_dir, sprintf("Fig5b_PCA_clinical_%s_v7", TISSUE)), 9, 8)
  list(panel = p, df = df)
}

res <- lapply(names(tissues), function(t) make_panel(t, tissues[[t]]))
names(res) <- names(tissues); res <- res[!vapply(res, is.null, logical(1))]

# liver (main) + PBMC (supp) are saved separately inside make_panel(); the
# combined two-panel figure is not plotted (kept the tidy centroid table only).
write.csv(bind_rows(lapply(res, `[[`, "df")),
          file.path(tab_dir, "Fig5b_PCA_centroids_clinical_v7.csv"), row.names = FALSE)
cat("\nsaved clean Fig 5b PCA for", paste(names(res), collapse = " + "),
    "(coloured by ", COLOUR_BY, ")\n", sep = "")
# End of the script
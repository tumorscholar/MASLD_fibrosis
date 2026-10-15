###############################################################################
# 37 - Figure 5 - clonal contraction of baseline-expanded clones at follow-up
#
# Uses FREQUENCY (% of repertoire), not raw cell count, because sequencing depth
# differs between visits. For the clones re-detected at follow-up (present in
# both visits) it draws a paired baseline->follow-up frequency plot per tissue
# and a paired Wilcoxon, plus an overlap/contraction summary table.
#
# Input : baseline_expanded_clones_tracking.csv
###############################################################################
library(dplyr); library(tidyr); library(ggplot2); library(readr)
if (file.exists("00_nature_style.R")) source("00_nature_style.R")
if (!exists("theme_nature")) theme_nature <- function(...) theme_bw(base_size = 8)
if (!exists("save_nature")) save_nature <- function(p, f, w, h) {
  ggsave(paste0(f, ".pdf"), p, width = w/2.54, height = h/2.54)
  ggsave(paste0(f, ".tiff"), p, width = w/2.54, height = h/2.54, dpi = 600, compression = "lzw", bg = "white") }

ann_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

trk <- readr::read_csv(file.path(ann_dir, "baseline_expanded_clones_tracking.csv"), show_col_types = FALSE) %>%
  mutate(Baseline_Percent = as.numeric(Baseline_Percent),
         Followup_Percent = as.numeric(Followup_Percent),
         Shared = Presence == "Present in Both")

# ---- overlap + contraction summary per tissue -------------------------------
summ <- trk %>% group_by(Tissue) %>%
  summarise(n_tracked = n(),
            n_shared  = sum(Shared),
            n_lost    = sum(!Shared),
            pct_shared = round(100 * n_shared / n_tracked, 1),
            .groups = "drop")
sh <- trk %>% filter(Shared) %>% group_by(Tissue) %>%
  summarise(n_shared = n(),
            pct_contracted = round(100 * mean(Followup_Percent < Baseline_Percent), 1),  # of re-detected
            pct_expanded   = round(100 * mean(Followup_Percent > Baseline_Percent), 1),
            median_baseline_pct = round(median(Baseline_Percent), 3),
            median_followup_pct = round(median(Followup_Percent), 3),
            wilcoxon_p = tryCatch(signif(wilcox.test(Baseline_Percent, Followup_Percent, paired = TRUE)$p.value, 3),
                                  error = function(e) NA),
            .groups = "drop")
summary_tbl <- left_join(summ, sh %>% select(-n_shared), by = "Tissue")
cat("=== Contraction / overlap summary ===\n"); print(as.data.frame(summary_tbl), row.names = FALSE)
write.csv(summary_tbl, file.path(tab_dir, "Fig5_contraction_summary_V7.csv"), row.names = FALSE)

# ---- paired frequency plot (shared clones only) -----------------------------
shared <- trk %>% filter(Shared) %>%
  mutate(dir = ifelse(Followup_Percent < Baseline_Percent, "Contracted", "Expanded"))
long <- shared %>% select(Tissue, CTaa, dir, Baseline_Percent, Followup_Percent) %>%
  pivot_longer(c(Baseline_Percent, Followup_Percent), names_to = "Visit", values_to = "Freq") %>%
  mutate(Visit = factor(ifelse(grepl("Baseline", Visit), "Baseline", "Followup"),
                        levels = c("Baseline", "Followup")))

# per-facet annotation (n, % contracted, p)
facet_lab <- sh %>% mutate(
  lab = sprintf("n=%d | %.0f%% contracted | p=%s", n_shared, pct_contracted,
                ifelse(is.na(wilcoxon_p), "NA", format(wilcoxon_p, scientific = TRUE))))

p <- ggplot(long, aes(Visit, Freq, group = CTaa, colour = dir)) +
  geom_line(alpha = 0.4, linewidth = 0.3) + geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~ Tissue) + scale_y_log10() +
  scale_colour_manual(values = c("Contracted" = "#2166ac", "Expanded" = "#b2182b"), name = NULL) +
  geom_text(data = facet_lab, aes(x = 1.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 2.2) +
  labs(x = NULL, y = "Clone frequency (% of repertoire, log10)",
       title = "Baseline-expanded clones re-detected at follow-up",
       subtitle = "paired frequency; blue = contracted, red = expanded") +
  theme_nature(7) + theme(legend.position = "bottom",
                          strip.background = element_rect(fill = "grey92", colour = NA))
save_nature(p, file.path(fig_dir, "Fig5_contraction_paired_frequency_V7"), 12, 8)

# ---- one-line numbers for the manuscript sentence ---------------------------
cat("\n=== sentence numbers ===\n")
for (i in seq_len(nrow(summary_tbl))) with(summary_tbl[i, ],
  cat(sprintf("%s: %d baseline-expanded clones tracked; %d (%.1f%%) re-detected at follow-up, of which %.0f%% contracted and %.0f%% expanded (paired Wilcoxon p=%s).\n",
      Tissue, n_tracked, n_shared, pct_shared, pct_contracted, pct_expanded,
      ifelse(is.na(wilcoxon_p), "NA", format(wilcoxon_p, scientific = TRUE)))))
cat("\nDone Fig 5 contraction.\n")
# End of the script

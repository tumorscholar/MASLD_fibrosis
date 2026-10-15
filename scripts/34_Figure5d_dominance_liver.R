###############################################################################
# 34 - Figure 5d - liver FNA TCR top-100 clonal dominance (baseline vs follow-up)
#
# Paired baseline/follow-up liver FNA. Top-100 clonal occupancy per sample
# (fraction of the repertoire held by the 100 most abundant clonotypes), with a
# paired Wilcoxon test. Line plot per donor.
#
# TCR source: v7 --chain=TR (v7_tcr_TR). Capture folders are auto-detected by
# patient id + LIVER/LIVERFNA, so it adapts to the real folder names.
###############################################################################
.libPaths(c("~/Rlibs/scRepertoire-2.2.1", .libPaths()))
library(scRepertoire); library(dplyr); library(tidyr); library(ggplot2)
source("00_nature_style.R")

tcr_dir <- "/data/Blizard-AlazawiLab/rk/cellranger/v7_tcr_TR"
fig_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/Figure5"
tab_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/Figure5"
for (d in c(fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# baseline patient -> follow-up capture id (edit here if an id differs)
fu_map <- c("10738" = "11596", "10291-1" = "11303", "11040" = "11816", "11183" = "11937")

caps_avail <- list.dirs(tcr_dir, recursive = FALSE, full.names = FALSE)
cat("TCR captures available (", length(caps_avail), "):\n"); print(caps_avail)
find_cap <- function(id, tis_rx) {
  hit <- grep(paste0("^GC-WL-", id, "-.*(", tis_rx, ")"), caps_avail, value = TRUE)
  if (length(hit) == 0) stop("no TCR capture for id=", id, " tissue~", tis_rx,
                             " (check the list printed above)"); hit[1]
}
read_contig <- function(cap) {
  f <- file.path(tcr_dir, cap, "outs", "filtered_contig_annotations.csv")
  if (!file.exists(f)) stop("missing contig file: ", f); read.csv(f, stringsAsFactors = FALSE)
}

# build paired baseline / follow-up liver TCR
contig_list <- list(); sample_names <- c(); type_vec <- c()
for (base in names(fu_map)) {
  b_cap <- find_cap(base,        "LIVER")     # matches -LIVER or -LIVERFNA
  f_cap <- find_cap(fu_map[base],"LIVER")
  for (k in 1:2) {
    vis <- c("Baseline", "Followup")[k]
    contig_list[[length(contig_list) + 1]] <- read_contig(c(b_cap, f_cap)[k])
    sample_names <- c(sample_names, paste0(vis, "-", base, "-LIVER"))
    type_vec     <- c(type_vec, paste0(vis, "-LIVER"))
  }
}
combined.TCR <- combineTCR(contig_list, samples = sample_names,
                           removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
for (i in seq_along(combined.TCR))
  combined.TCR[[i]]$Type <- factor(type_vec[i], levels = c("Baseline-LIVER", "Followup-LIVER"))

# ---- top-100 clonal occupancy per sample -----------------------------------
top100 <- bind_rows(lapply(seq_along(combined.TCR), function(i) {
  d  <- combined.TCR[[i]]; nm <- names(combined.TCR)[i]
  cf <- sort(table(d$CTstrict), decreasing = TRUE)
  data.frame(Sample = nm, Timepoint = ifelse(grepl("Baseline", nm), "Baseline", "Followup"),
             Tissue = "LIVER", Donor = gsub("(Baseline-|Followup-)(.+)-LIVER", "\\2", nm),
             Top100_Percent = round(100 * sum(cf[1:min(100, length(cf))]) / sum(cf), 2))
}))
top100$Timepoint <- factor(top100$Timepoint, levels = c("Baseline", "Followup"))
cat("\n== Top-100 dominance (liver) ==\n"); print(top100, row.names = FALSE)
write.csv(top100, file.path(tab_dir, "Fig5d_top100_dominance_liver_v7.csv"), row.names = FALSE)

wide <- top100 %>% select(Donor, Timepoint, Top100_Percent) %>%
  pivot_wider(names_from = Timepoint, values_from = Top100_Percent) %>%
  filter(!is.na(Baseline) & !is.na(Followup))
if (nrow(wide) > 1) cat(sprintf("\nPaired Wilcoxon (liver): p = %.4f | mean change = %.2f%%\n",
    wilcox.test(wide$Baseline, wide$Followup, paired = TRUE)$p.value,
    mean(wide$Followup - wide$Baseline)))

# ---- Fig 5d: paired line plot ----------------------------------------------
fig5d <- ggplot(top100, aes(Timepoint, Top100_Percent, group = Donor, colour = Donor)) +
  geom_line(linewidth = 0.8, alpha = 0.85) + geom_point(size = 2.5) + facet_wrap(~ Tissue) +
  scale_colour_manual(values = npg_pal) + scale_y_continuous(breaks = seq(0, 60, 10), limits = c(0, 60)) +  # shared scale liver+PBMC
  labs(title = "Top-100 TCR clonal dominance - liver FNA", subtitle = "Paired baseline vs follow-up",
       x = NULL, y = "Top-100 clones (% of repertoire)", colour = "Donor") +
  theme_nature(7) + theme(legend.position = "bottom",
                          strip.background = element_rect(fill = "grey92", colour = NA))
save_nature(fig5d, file.path(fig_dir, "Fig5d_FNA_top100_dominance_liver_v7"), 6, 6)
cat("\nDone Fig 5d (liver).\n")
# End of the script

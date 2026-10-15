################################################################################
# 00_nature_style.R
#
# Shared "Nature main-journal" figure style for every figure script in this
# project. Source this file at the top of a figure script:
#
#     source("00_nature_style.R")
#
# It provides:
#   npg_pal      - the official Nature Publishing Group (npg) categorical palette
#   stage_cols   - fixed colours for Healthy / No_fibrosis / Fibrosis
#   tissue_cols  - fixed colours for Liver / PBMC / SAT / VAT
#   theme_nature(base) - a clean Arial theme at Nature's 5-7 pt spec
#   save_nature(p, file, w_cm, h_cm) - exports PDF (vector) + 600-dpi TIFF
#
# Font: Nature requires a sans-serif face (Helvetica/Arial) at 5-7 pt final
# size. Register Arial once per R session before plotting (see setup_arial.R /
# fonts.R). On a machine without Arial, Liberation Sans is a metric-identical
# substitute.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

## ---- Font bootstrap: make Arial available on EVERY device -------------------
## No terminal and no system font install needed. We register Liberation Sans
## (metric-identical to Arial) from a .ttf found in the working dir / ~ / ~/.fonts,
## under BOTH the "Arial" and the default "sans" family, and turn on showtext so
## that ggplot AND grid-based packages (pheatmap, ComplexHeatmap) AND base
## graphics all render in Arial. Put LiberationSans-Regular.ttf (and -Bold.ttf)
## in your scripts folder if they are not already installed on the machine.
.nat_find_font <- function(fname) {
  cand <- c(fname,
            file.path(getwd(), fname),
            file.path(path.expand("~"), fname),
            file.path(path.expand("~/.fonts"), fname),
            file.path("/usr/share/fonts/truetype/liberation", fname),
            file.path("/usr/share/fonts/truetype/liberation2", fname))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
.nat_reg  <- .nat_find_font("LiberationSans-Regular.ttf")
.nat_regb <- .nat_find_font("LiberationSans-Bold.ttf")

# systemfonts registration (used by svglite / ragg / ggplot text metrics)
if (!is.na(.nat_reg) && requireNamespace("systemfonts", quietly = TRUE)) {
  try(systemfonts::register_font("Arial", plain = .nat_reg,
                                 bold = if (!is.na(.nat_regb)) .nat_regb else .nat_reg), silent = TRUE)
}

# showtext: renders text on cairo_pdf / tiff / png / grid / base as Arial.
# "sans" is overridden too so pheatmap / ComplexHeatmap / base default text
# (which use the sans family) also come out Arial - no per-package edits needed.
if (!is.na(.nat_reg) &&
    requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  try({
    sysfonts::font_add("Arial", regular = .nat_reg,
                       bold = if (!is.na(.nat_regb)) .nat_regb else .nat_reg)
    sysfonts::font_add("sans",  regular = .nat_reg,
                       bold = if (!is.na(.nat_regb)) .nat_regb else .nat_reg)
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 600)   # MUST match the 600-dpi TIFF exports below
    options(nature.arial = TRUE)
    message("nature_style: Arial (Liberation Sans) active on all devices via showtext (dpi = 600).")
  }, silent = TRUE)
} else {
  message("nature_style: showtext / Liberation Sans not found - falling back to the device default font.")
}
# NOTE on dpi: showtext sizes RASTER text by dpi, so every tiff()/png() you open
# directly should use res = 600 (as the figure scripts do). Vector PDFs are not
# affected by this dpi. If you ever export a raster at another dpi, wrap it with
# showtext::showtext_opts(dpi = <that dpi>) ... showtext::showtext_opts(dpi = 600).

## ---- Palettes ---------------------------------------------------------------
# Nature Publishing Group (npg) categorical palette
npg_pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
             "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# -----------------------------------------------------------------------------
# PAPER-WIDE COLOUR SCHEME. Stage and tissue use DIFFERENT colour families so
# they never clash: stage = purple/lavender/teal ; tissue = red/blue/green/navy.
# Use stage_cols for stage-grouped panels, tissue_cols for tissue-grouped panels.
# -----------------------------------------------------------------------------
# Disease-stage palette: lilac -> violet -> magenta. Three DISTINCT hues (easy to
# tell No fibrosis from Fibrosis) that still read light->saturated as severity
# rises, and stay clear of the tissue red/blue/green/navy.
stage_cols <- c("Healthy" = "#C3A9E0", "No_fibrosis" = "#8E44AD", "Fibrosis" = "#C2185B")
# darker shades for the jittered points that sit on the stage bars
stage_pts  <- c("Healthy" = "#9E82C4", "No_fibrosis" = "#6E2E8A", "Fibrosis" = "#8E1044")

# Tissue palette (npg red / blue / green / navy) - distinct hue family from stage
tissue_cols <- c("Liver" = "#E64B35", "LIVER" = "#E64B35",
                 "PBMC" = "#4DBBD5", "SAT" = "#00A087", "VAT" = "#3C5488")

## ---- Theme ------------------------------------------------------------------
theme_nature <- function(base = 7) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      axis.line       = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.3, colour = "black"),
      axis.text       = element_text(colour = "black", size = base - 1),
      axis.title      = element_text(colour = "black", size = base),
      plot.title      = element_text(hjust = 0.5, size = base, face = "plain"),
      plot.subtitle   = element_text(hjust = 0.5, size = base - 1),
      strip.text      = element_text(size = base, face = "plain"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.text     = element_text(size = base - 1),
      legend.title    = element_text(size = base),
      legend.key.size = unit(3, "mm")
    )
}

## ---- Export helper: PDF (vector) + high-res TIFF ----------------------------
# w_cm, h_cm = final panel size in centimetres (no scaling, so fonts stay true).
save_nature <- function(p, file, w_cm, h_cm, dpi = 600) {
  grDevices::cairo_pdf(paste0(file, ".pdf"), width = w_cm / 2.54, height = h_cm / 2.54)
  print(p); grDevices::dev.off()
  ggplot2::ggsave(paste0(file, ".tiff"), plot = p,
                  width = w_cm, height = h_cm, units = "cm",
                  dpi = dpi, device = "tiff", compression = "lzw", bg = "white")
}
################################################################################
# =============================================================================
# 44 - hLICA pipeline, STEP 3 of 3: the three composition figures
#
#   From the STEP 2 transfer table, keeps confident lymphoid calls (score>=0.5)
#   and plots each of the 8 lymphoid subtypes as a % of that donor's lymphocytes.
#   Our 10742 (lymphoid only, myeloid dropped) is the red diamond.
#   Three panels: all donors / BMI-verified lean / documented-obese excluded.
#
#   INPUT : /gpfs/scratch/hdx044/hlica/transfer_hlica_cells.csv   (from STEP 2)
#   OUTPUT: FIGDIR/ExtFig_detailed_{allSamples,leanOnly,excludeObese}.{pdf,png}
#   Run   : Rscript 41_hlica_make_healthy_composition_figs.R
# =============================================================================
suppressMessages({library(data.table); library(ggplot2)})

## ---- paths: transfer table in extended TABLES, plots in extended FIGURES ----
CSV    <- "/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData/transfer_hlica_cells.csv"  # from STEP 2
FIGDIR <- "/data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData"
dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)

dt <- fread(CSV)            # was: dt <- fread(CSV)[score >= 0.5]
types <- c("Cytotoxic CD8 T cells","Activated CD8 T cells","Activated T cells","Memory T cells",
           "Naive/TCM T cells","Cytotoxic NK cells","Naive B cells","Plasma cells")
short <- c("Cyto CD8","Act CD8","Act T","Mem T","Naive/TCM","Cyto NK","Naive B","Plasma"); names(short)<-types
dt <- dt[pred %in% types]; dt[, pred := factor(pred, levels=types)]

## 10742 reference (LIVER Healthy cell counts, your annotation)
ocnt <- c(261,485,211,221,217,332,105,901); names(ocnt) <- types
ours <- data.table(pred=factor(types,levels=types), prop=100*ocnt/sum(ocnt))

bmi <- dt[, .(BMI=BMI[1]), by=donor]
npg <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2")

makefig <- function(keep, title, subtitle, fname){
  comp <- dt[donor %in% keep, .N, by=.(donor,study,pred)]; comp[, prop:=100*N/sum(N), by=donor]
  nD <- uniqueN(comp$donor); nS <- uniqueN(comp$study)
  p <- ggplot(comp, aes(pred, prop)) +
    geom_boxplot(outlier.shape=NA, width=0.62, fill="#EEEEEE", colour="#999999", linewidth=0.3) +
    geom_jitter(aes(colour=study), width=0.12, size=1.1, alpha=0.8) +
    geom_point(data=ours, aes(pred, prop), shape=23, size=3.4, fill="#E64B35", colour="black", stroke=0.6) +
    scale_colour_manual(values=npg, name="Study") + scale_x_discrete(labels=short) +
    labs(x=NULL, y="% of lymphocytes (per donor)", title=title, subtitle=sprintf(subtitle, nD, nS)) +
    theme_classic(base_size=11) +
    theme(axis.text.x=element_text(angle=25,hjust=1), plot.title=element_text(size=12),
          plot.subtitle=element_text(size=8.5, colour="grey30"))
  ggsave(file.path(FIGDIR,paste0(fname,".pdf")), p, width=9.2, height=4.9)
  ggsave(file.path(FIGDIR,paste0(fname,".png")), p, width=9.2, height=4.9, dpi=300)
  cat(fname,":",nD,"donors,",nS,"studies\n")
}

all      <- unique(dt$donor)
lean     <- bmi[BMI < 25, donor]
notobese <- bmi[is.na(BMI) | BMI < 25, donor]

makefig(all,      "Healthy liver lymphoid composition - all donors",
        "Red diamond = 10742 vs %d healthy donors (%d studies); label-transfer annotation",
        "ExtFig_detailed_allSamples")
makefig(lean,     "Healthy liver lymphoid composition - BMI-verified lean donors",
        "Red diamond = 10742 vs %d lean donors, BMI<25 (%d study); label-transfer annotation",
        "ExtFig_detailed_leanOnly")
makefig(notobese, "Healthy liver lymphoid composition - documented obese excluded",
        "Red diamond = 10742 vs %d donors, BMI>=25 excluded (%d studies); label-transfer annotation",
        "ExtFig_detailed_excludeObese")

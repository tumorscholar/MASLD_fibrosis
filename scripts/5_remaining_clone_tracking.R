###############################################################################
# 5 - Find clone-tracking clones that are NOT yet in the antigen-annotation file
#     and export their CDR3-beta sequences, ready to push through TCRMatch/VDJdb.
#
# The donut (36) joins the tracking table to v7_all_clones_annotation.csv; any
# clone absent from that file defaults to "Unannotated". After a tracking
# rebuild (new donors / cells), some clones - e.g. patient 10291-1's PBMC - were
# never run through the databases, which deflates the annotated %.
#
# This script writes the unique missing CDR3-beta list so you can re-run only the
# gap (not the whole cohort). Feed the output to your usual TCRMatch / VDJdb /
# TCRdb step, then rbind the new rows into v7_all_clones_annotation.csv and re-run 36.
###############################################################################
ann_dir <- "/data/Blizard-AlazawiLab/rk/Tcell/CloneAnnotation"
trk_csv <- file.path(ann_dir, "baseline_expanded_clones_tracking.csv")
ann_csv <- file.path(ann_dir, "v7_all_clones_annotation.csv")
out_txt <- file.path(ann_dir, "v7_missing_CDR3b_to_annotate.txt")   # one CDR3b per line
out_map <- file.path(ann_dir, "v7_missing_clones_to_annotate.csv")  # Tissue, CTaa, CDR3b

trk <- read.csv(trk_csv, stringsAsFactors = FALSE, check.names = FALSE)
ann <- read.csv(ann_csv, stringsAsFactors = FALSE, check.names = FALSE)

# CDR3-beta = the chain(s) after the LAST underscore in CTaa (scRepertoire packs
# it as "CDR3a_CDR3b"; either side may carry ";"-separated dual chains).
beta_of <- function(ctaa) {
 b <- sub("^.*_", "", ctaa)                 # keep the beta group
 unlist(strsplit(b, ";", fixed = TRUE))     # split dual-beta
}

# already-processed beta sequences (from the annotation file's CTaaB column)
have_beta <- unique(unlist(lapply(ann$CTaaB, function(x)
 if (is.na(x) || x == "") character(0) else unlist(strsplit(x, ";", fixed = TRUE)))))
# ... and, as a fallback, betas parsed straight from the annotation CTaa
have_beta <- unique(c(have_beta, unlist(lapply(ann$CTaa, beta_of))))

# tracking clones whose beta is not yet annotated
rows <- do.call(rbind, lapply(seq_len(nrow(trk)), function(i) {
 bs <- beta_of(trk$CTaa[i]); miss <- setdiff(bs, have_beta)
 if (!length(miss)) return(NULL)
 data.frame(Tissue = trk$Tissue[i], Donor = trk$Donor[i],
            CTaa = trk$CTaa[i], CDR3b = miss, stringsAsFactors = FALSE)
}))

if (is.null(rows) || !nrow(rows)) {
 cat("All tracking clones are already covered by the annotation file - nothing to do.\n")
} else {
 cdr3b <- sort(unique(rows$CDR3b))
 writeLines(cdr3b, out_txt)
 write.csv(unique(rows), out_map, row.names = FALSE)
 cat(sprintf("Missing clones: %d rows across %d unique CDR3-beta.\n", nrow(rows), length(cdr3b)))
 cat("By Tissue x Donor:\n")
 print(as.data.frame(table(rows$Tissue, rows$Donor)))
 cat(sprintf("\nWrote:\n  %s   (CDR3b list for TCRMatch/VDJdb)\n  %s   (Tissue/CTaa/CDR3b map)\n",
             out_txt, out_map))
 cat("Next: run your TCRMatch/VDJdb/TCRdb step on the CDR3b list, add the new\n",
     "rows to v7_all_clones_annotation.csv (Tissue, CTaa, CTaaB, ... , Annotated),\n",
     "then re-run 36_Figure5_TCRdonut.R.\n")
}
# End of the script
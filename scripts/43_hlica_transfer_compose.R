# =============================================================================
# 43 - hLICA pipeline, STEP 2 of 3: label transfer (our SeuObjV7 -> hLICA query)
#
#   Seurat FindTransferAnchors + TransferData, using our FULL cell_type
#   annotation as reference. Writes one row per query cell:
#       transfer_<name>_cells.csv  (barcode, pred, score, study, donor, BMI,
#                                   author_cell_type)
#   Also prints a quick coarse-lineage (T/NK/B/Plasma) sanity check vs our 10742.
#
#   INPUT : a query folder from STEP 1 (matrix.mtx.gz, barcodes.tsv, genes.tsv,
#           metadata.csv) + /data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds
#   USAGE : Rscript 40_hlica_transfer_compose.R <query_dir> <name>
#             Rscript 40_hlica_transfer_compose.R /gpfs/scratch/hdx044/hlica/query_hlica hlica
#             Rscript 40_hlica_transfer_compose.R /gpfs/scratch/hdx044/hlica/query_guilliams_lean guilliams_lean
#   NEXT  : 41_hlica_make_healthy_composition_figs.R
# =============================================================================
suppressMessages({library(Seurat);library(Matrix);library(data.table);library(ggplot2)})
args=commandArgs(trailingOnly=TRUE); QDIR=args[1]; QNAME=args[2]
OUT="/data/Blizard-AlazawiLab/rk/Tcell/tables/ExtendedData"   # transfer_*_cells.csv -> extended TABLES
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)           # (QDIR = query matrices, in scratch)
ref=readRDS("/data/Blizard-AlazawiLab/rk/seurat/SeuObjV7.rds"); DefaultAssay(ref)="RNA"
ctcol=intersect(c("cell_type","celltype"),colnames(ref@meta.data))[1]
ref=NormalizeData(ref,verbose=F); ref=FindVariableFeatures(ref,verbose=F); ref=ScaleData(ref,verbose=F); ref=RunPCA(ref,verbose=F)
m=ReadMtx(mtx=file.path(QDIR,"matrix.mtx.gz"),cells=file.path(QDIR,"barcodes.tsv"),features=file.path(QDIR,"genes.tsv"),feature.column=1)
q=CreateSeuratObject(m); q=NormalizeData(q,verbose=F)
cat("shared genes ref/query:",length(intersect(rownames(ref),rownames(q))),"\n")
anc=FindTransferAnchors(reference=ref,query=q,dims=1:30,reference.reduction="pca",verbose=F)
pred=TransferData(anchorset=anc,refdata=ref@meta.data[[ctcol]],dims=1:30,verbose=F); q=AddMetaData(q,pred)
md=fread(file.path(QDIR,"metadata.csv"))
dt=data.table(barcode=colnames(q),pred=q$predicted.id,score=q$prediction.score.max); dt=merge(dt,md,by="barcode")
fwrite(dt,file.path(OUT,paste0("transfer_",QNAME,"_cells.csv")))
cat("median pred score:",median(dt$score)," | cells>=0.5:",round(mean(dt$score>=0.5)*100,1),"%\n")
dt=dt[score>=0.5]
lin=function(x) fifelse(x %in% c("Cytotoxic CD8 T cells","Activated CD8 T cells","Activated T cells","Memory T cells","Naive/TCM T cells"),"T",
              fifelse(x=="Cytotoxic NK cells","NK",fifelse(x=="Naive B cells","B",fifelse(x=="Plasma cells","Plasma",NA_character_))))
dt[,l:=lin(pred)]; d=dt[!is.na(l)]; comp=d[,.N,by=.(donor,study,l)]; comp[,prop:=100*N/sum(N),by=donor]
rr=data.table(ct=ref@meta.data[[ctcol]],Tissue=as.character(ref$Tissue),Stage=gsub("No_Fibrosis","No_fibrosis",as.character(ref$Stage)))
rh=rr[Tissue=="LIVER"&Stage=="Healthy"]; rh[,l:=lin(ct)]; rh=rh[!is.na(l)]; oc=rh[,.N,by=l]; oc[,prop:=100*N/sum(N)]
ord=c("T","NK","B","Plasma")
cat("\n10742 vs transferred",QNAME,"(YOUR labels):\n")
for(k in ord){s=comp[l==k]$prop; oo=oc[l==k]$prop
 cat(sprintf(" %-7s ours=%.1f min/med/max=%.1f/%.1f/%.1f within=%s\n",k,oo,min(s),median(s),max(s),oo>=min(s)&oo<=max(s)))}
cat("donors:",uniqueN(comp$donor)," studies:",uniqueN(comp$study),"\n")

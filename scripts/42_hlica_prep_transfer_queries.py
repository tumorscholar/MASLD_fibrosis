# =============================================================================
# 42 - hLICA pipeline, STEP 1 of 3: prepare label-transfer query matrices
#
#   Reads the hLICA lymphocyte atlas, keeps fresh cells (suspension_type=="cell",
#   drops nuclei), parses BMI, and dumps two 10x-style query folders for Seurat:
#       query_hlica           - all studies except Narmada_2024
#       query_guilliams_lean  - Guilliams_2022 donors with BMI < 25
#   Each folder: matrix.mtx.gz, barcodes.tsv, genes.tsv (feature_name symbols),
#                metadata.csv (barcode, study, donor, BMI, author_cell_type)
#
#   INPUT : /gpfs/scratch/hdx044/hlica/hlica_lymphocyte.h5ad
#   NEXT  : 40_hlica_transfer_compose.R  <query_dir>  <name>
#   Run   : python 39_hlica_prep_transfer_queries.py
# =============================================================================
import anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, os
H5AD="/gpfs/scratch/hdx044/hlica/hlica_lymphocyte.h5ad"; OUT="/gpfs/scratch/hdx044/hlica"   # data (atlas + query dirs) in SCRATCH
a=ad.read_h5ad(H5AD); o=a.obs.copy()
o["BMI"]=pd.to_numeric(o["donor_BMI_at_collection"].astype(str).replace("unknown",np.nan),errors="coerce")
cell=o["suspension_type"].astype(str)=="cell"
if "counts" in a.layers: X=a.layers["counts"]; genes=a.var.get("feature_name",pd.Series(a.var_names,index=a.var_names))
elif a.raw is not None:  X=a.raw.X;            genes=a.raw.var.get("feature_name",pd.Series(a.raw.var_names,index=a.raw.var_names))
else:                    X=a.X;               genes=a.var.get("feature_name",pd.Series(a.var_names,index=a.var_names))
X=sp.csr_matrix(X); genes=[str(g) for g in list(genes)]
print("counts:",X.shape,"| genes:",len(genes),"| eg:",genes[:5])
dbmi=o.groupby("donor_id")["BMI"].transform("first")
def dump(mask,name):
    idx=np.where(mask.values)[0]; sub=X[idx,:]; d=os.path.join(OUT,name); os.makedirs(d,exist_ok=True)
    sio.mmwrite(os.path.join(d,"matrix.mtx"),sub.T.tocoo()); os.system("gzip -f "+os.path.join(d,"matrix.mtx"))
    bc=[str(x) for x in a.obs_names[idx]]
    open(os.path.join(d,"barcodes.tsv"),"w").write("\n".join(bc)+"\n")
    open(os.path.join(d,"genes.tsv"),"w").write("\n".join(genes)+"\n")
    md=o.iloc[idx][["STUDY","donor_id","BMI","author_cell_type"]].copy(); md.insert(0,"barcode",bc)
    md.columns=["barcode","study","donor","BMI","author_cell_type"]; md.to_csv(os.path.join(d,"metadata.csv"),index=False)
    print(name,":",sub.shape[0],"cells,",md.donor.nunique(),"donors")
dump(cell & (o["STUDY"].astype(str)!="Narmada_2024"), "query_hlica")
dump(cell & (o["STUDY"].astype(str)=="Guilliams_2022") & (dbmi<25),  "query_guilliams_lean")
print("PREP DONE")

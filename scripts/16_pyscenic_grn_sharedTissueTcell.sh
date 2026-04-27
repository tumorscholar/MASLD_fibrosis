#!/bin/bash

############################################################
# Slurm job configuration (QMUL Apocrita / Andrena)
############################################################
#SBATCH --job-name=SCENIC_GRN_Tcell
#SBATCH --output=SCENIC_GRN_Tcell.out
#SBATCH --error=SCENIC_GRN_Tcell.err
#SBATCH --mail-user=r.kumar@qmul.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=240:00:00

############################################################
# Load required software
############################################################
module load miniforge/24.7.1

############################################################
# Activate PySCENIC environment
############################################################
conda activate scenic_py38

############################################################
# Move to working directory
############################################################
cd /data/Blizard-AlazawiLab/rk/scenicTcell || exit 1

############################################################
# Run PySCENIC GRN inference
############################################################
pyscenic grn \
  /data/Blizard-AlazawiLab/rk/seurat/TissueExpandedTcells.loom \
  /data/Blizard-AlazawiLab/rk/scenicTcell/config/tfs.txt \
  -o /data/Blizard-AlazawiLab/rk/scenicTcell/results/adj.csv \
  --num_workers 4

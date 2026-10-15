#!/bin/bash
###############################################################################
# 20_cellchat_run.slurm - run 15_CellChatTissueStage.R for LIVER / SAT / VAT
# Submit:  sbatch 15_CellChatTissueStage.slurm
###############################################################################
#SBATCH --job-name=cellchat_v7
#SBATCH --output=/gpfs/scratch/hdx044/cellchat_v7_%A_%a.log
#SBATCH --mail-user=r.kumar@qmul.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-2

module load R/4.5.1

TISSUES=(LIVER SAT VAT)
export CC_TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}
echo "== CellChat for $CC_TISSUE =="

cd /data/home/hdx044/MASLD_fibrosis/scripts || exit 1
export R_PROFILE_USER="$(pwd)/.Rprofile"

Rscript /data/home/hdx044/MASLD_fibrosis/scripts/15_CellChatTissueStage.R
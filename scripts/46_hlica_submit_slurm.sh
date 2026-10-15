#!/usr/bin/env bash
#SBATCH --job-name=hlica_pipe
#SBATCH --partition=compute            # <-- SET to your partition (run: sinfo -s)
# #SBATCH --account=YOUR_ACCOUNT       # <-- uncomment/set if your cluster requires it
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G                      # atlas load + Seurat anchors are memory-heavy
#SBATCH --time=06:00:00
#SBATCH --output=hlica_pipe_%j.log     # stdout+stderr
# =============================================================================
# 46 - submit the hLICA pipeline (45) to a COMPUTE node. NEVER run on login-01.
#   sbatch 46_hlica_submit_slurm.sh
#   squeue -u $USER ; tail -f hlica_pipe_*.log
# =============================================================================
set -euo pipefail

# make python + R available (edit to match your cluster) -----------------------
source /data/Blizard-AlazawiLab/rk/miniconda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate base 2>/dev/null || true
# module load R/4.3.1        # <-- if R is an environment module, not in your conda env

bash "$HOME/MASLD_fibrosis/scripts/45_hlica_run_pipeline.sh"

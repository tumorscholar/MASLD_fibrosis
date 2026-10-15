#!/usr/bin/env bash
# =============================================================================
# 45 - hLICA composition pipeline: RUN ORDER (this is the start script)
#        42 prep queries (Python) -> 43 label transfer (R) -> 44 figures (R)
#   Data stays in scratch; transfer table -> Tcell/tables/ExtendedData;
#   figures -> Tcell/figures/ExtendedData.
#   Run on a COMPUTE node (see 46_hlica_submit_slurm.sh), not login.
#     bash 45_hlica_run_pipeline.sh
# =============================================================================
set -euo pipefail

HLICA=/gpfs/scratch/hdx044/hlica            # DATA (atlas + query matrices) in scratch
SDIR="$HOME/MASLD_fibrosis/scripts"         # scripts live here
cd "$HLICA"

echo "== 1/3  prepare query matrices =="
python  "$SDIR/42_hlica_prep_transfer_queries.py"

echo "== 2/3  label transfer (hLICA) =="
Rscript "$SDIR/43_hlica_transfer_compose.R" "$HLICA/query_hlica" hlica
# optional lean-only query:
# Rscript "$SDIR/43_hlica_transfer_compose.R" "$HLICA/query_guilliams_lean" guilliams_lean

echo "== 3/3  composition figures =="
Rscript "$SDIR/44_hlica_make_healthy_composition_figs.R"

echo "DONE -> /data/Blizard-AlazawiLab/rk/Tcell/figures/ExtendedData/ExtFig_detailed_*"

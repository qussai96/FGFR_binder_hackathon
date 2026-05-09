#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_SCRIPT="$SCRIPT_DIR/slurm_round0_pipeline.sh"

NUM_DESIGNS=5
ALL_STRATEGIES=0
RELAX_FILTERS=0
PREP_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --num-designs)
      NUM_DESIGNS="$2"
      shift 2
      ;;
    --all-strategies)
      ALL_STRATEGIES=1
      shift
      ;;
    --relaxed-filters)
      RELAX_FILTERS=1
      shift
      ;;
    --prep-only)
      PREP_ONLY=1
      shift
      ;;
    -h|--help)
      cat <<'USAGE'
Usage:
  bash submit_round0_slurm.sh [--num-designs N] [--all-strategies] [--relaxed-filters] [--prep-only]

Behavior:
  - Always submits to Slurm GPU-A40.
  - No local prep or local design execution.
  - --prep-only is accepted for compatibility and performs a no-op submission dry intent.
USAGE
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "ERROR: Missing Slurm pipeline script: $PIPELINE_SCRIPT" >&2
  exit 3
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch not found. This workflow requires Slurm submission." >&2
  exit 4
fi

echo "=== FGFR2 Round0 Slurm Submission ==="
echo "Campaign dir : $CAMPAIGN_DIR"
echo "Designs/strat: $NUM_DESIGNS"
echo "All strategies: $ALL_STRATEGIES"
echo "Relaxed filt.: $RELAX_FILTERS"

if [[ "$PREP_ONLY" -eq 1 ]]; then
  echo "INFO: --prep-only requested in submit-only mode; no local prep is run."
  echo "INFO: Skipping submission by request."
  exit 0
fi

JOB_ID=$(sbatch \
  --job-name "fgfr2_round0" \
  --export "ALL,CAMPAIGN_DIR=$CAMPAIGN_DIR,NUM_DESIGNS=$NUM_DESIGNS,ALL_STRATEGIES=$ALL_STRATEGIES,RELAX_FILTERS=$RELAX_FILTERS" \
  "$PIPELINE_SCRIPT" | awk '{print $NF}')

echo "Submitted Round0 pipeline to GPU-A40: job $JOB_ID"
echo "Monitor with: squeue -u $USER"
echo "Logs from Slurm: fgfr2_round0_${JOB_ID}.out / .err"

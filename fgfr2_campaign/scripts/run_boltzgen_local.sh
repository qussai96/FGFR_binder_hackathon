#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_boltzgen_local.sh --config <path> --strategy <name> --num-designs <N> --outdir <path> [--dry-run]

Cluster-only wrapper for BoltzGen on Slurm GPU-A40.
No local execution path.
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$CAMPAIGN_DIR/.." && pwd)"
LOG_DIR="$CAMPAIGN_DIR/out/logs"
mkdir -p "$LOG_DIR"
SBATCH_SCRIPT="$SCRIPT_DIR/slurm_boltzgen.sh"

CONFIG=""
STRATEGY=""
NUM_DESIGNS=""
OUTDIR=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --strategy) STRATEGY="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$CONFIG" || -z "$STRATEGY" || -z "$NUM_DESIGNS" || -z "$OUTDIR" ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage
  exit 2
fi

if [[ "$NUM_DESIGNS" -gt 10 ]]; then
  echo "WARNING: Requested $NUM_DESIGNS designs; recommended 5-10 for initial A40 batches." >&2
fi

if [[ ! -f "$SBATCH_SCRIPT" ]]; then
  echo "ERROR: Missing Slurm runner script: $SBATCH_SCRIPT" >&2
  exit 3
fi

mkdir -p "$OUTDIR"
LOG_FILE="$LOG_DIR/boltzgen_${STRATEGY}_$(date +%Y%m%d_%H%M%S).log"

NATIVE_SPEC="$CONFIG"

SBATCH_CMD=(
  sbatch
  --job-name "fgfr2_${STRATEGY}"
  --export "ALL,STRATEGY=$STRATEGY,NUM_DESIGNS=$NUM_DESIGNS,CAMPAIGN_DIR=$CAMPAIGN_DIR,CONFIG_PATH=$NATIVE_SPEC,OUTDIR=$OUTDIR"
  "$SBATCH_SCRIPT"
  --strategy "$STRATEGY"
  --num-designs "$NUM_DESIGNS"
)

echo "Strategy: $STRATEGY"
echo "Config: $CONFIG"
echo "Outdir: $OUTDIR"
echo "Submission log: $LOG_FILE"
echo "Command shape: ${SBATCH_CMD[*]}"
echo "Required env hints: Slurm access and GPU-A40 partition availability"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY RUN: not submitting Slurm job."
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch not found. Slurm submission is required in cluster-only mode." >&2
  exit 4
fi

{
  echo "[INFO] Submitting BoltzGen Slurm job at $(date --iso-8601=seconds)"
  "${SBATCH_CMD[@]}"
  echo "[INFO] Submission complete at $(date --iso-8601=seconds)"
} 2>&1 | tee "$LOG_FILE"

#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_rfdiffusion_local.sh --config <path> --strategy <name> --num-designs <N> --outdir <path> [--num-batches M] [--dry-run]

Local-only runner for RFdiffusion3 on NVIDIA A40.
No Lyceum, no S3, no Modal.

Supports three execution modes (auto-detected):
1) `rfdiffusion3` CLI (legacy wrapper mode)
2) `rfd3 design` (Rosetta Foundry RFdiffusion3)
3) `apptainer exec docker://rosettacommons/foundry rfd3 design`
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$CAMPAIGN_DIR/.." && pwd)"
LOG_DIR="$CAMPAIGN_DIR/out/logs"
mkdir -p "$LOG_DIR"

CONFIG=""
STRATEGY=""
NUM_DESIGNS=""
OUTDIR=""
NUM_BATCHES="${NUM_BATCHES:-1}"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --strategy) STRATEGY="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --num-batches) NUM_BATCHES="$2"; shift 2 ;;
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

if [[ -f "$ROOT_DIR/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$ROOT_DIR/.venv/bin/activate"
elif [[ -f "$CAMPAIGN_DIR/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$CAMPAIGN_DIR/.venv/bin/activate"
else
  echo "INFO: No local venv found; using current environment."
fi

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "ERROR: nvidia-smi not found. Local GPU environment is required." >&2
  exit 3
fi

nvidia-smi -L || true

mkdir -p "$OUTDIR"
LOG_FILE="$LOG_DIR/rfdiffusion3_${STRATEGY}_$(date +%Y%m%d_%H%M%S).log"

RFDIFFUSION_CMD="${RFDIFFUSION_CMD:-rfdiffusion3}"
CMD=("$RFDIFFUSION_CMD" "--config" "$CONFIG" "--num_designs" "$NUM_DESIGNS" "--output_dir" "$OUTDIR")

echo "Strategy: $STRATEGY"
echo "Config: $CONFIG"
echo "Outdir: $OUTDIR"
echo "Batches: $NUM_BATCHES"
echo "Log: $LOG_FILE"
echo "Required env hints: CUDA_VISIBLE_DEVICES (optional), RFDIFFUSION_CMD (optional override), APPTAINER_BINDPATH (optional)"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY RUN: not executing RFdiffusion3."
  exit 0
fi

run_legacy_wrapper() {
  local cmd=("$RFDIFFUSION_CMD" "--config" "$CONFIG" "--num_designs" "$NUM_DESIGNS" "--output_dir" "$OUTDIR")
  echo "Mode: legacy wrapper"
  echo "Command shape: ${cmd[*]}"
  "${cmd[@]}"
}

run_rfd3_cli() {
  local cmd=(
    rfd3 design
    "out_dir=$OUTDIR"
    "inputs=$CONFIG"
    "diffusion_batch_size=$NUM_DESIGNS"
    "n_batches=$NUM_BATCHES"
  )
  echo "Mode: rfd3 CLI"
  echo "Command shape: ${cmd[*]}"
  "${cmd[@]}"
}

run_apptainer_foundry() {
  local cfg_abs out_abs
  cfg_abs="$(python - <<PY
import os
print(os.path.abspath('$CONFIG'))
PY
)"
  out_abs="$(python - <<PY
import os
print(os.path.abspath('$OUTDIR'))
PY
)"

  mkdir -p "$out_abs"
  local cmd=(
    apptainer exec --nv
    "docker://rosettacommons/foundry"
    rfd3 design
    "out_dir=$out_abs"
    "inputs=$cfg_abs"
    "diffusion_batch_size=$NUM_DESIGNS"
    "n_batches=$NUM_BATCHES"
  )
  echo "Mode: apptainer+foundry"
  echo "Command shape: ${cmd[*]}"
  "${cmd[@]}"
}

{
  echo "[INFO] Starting RFdiffusion3 local run at $(date --iso-8601=seconds)"
  if command -v "$RFDIFFUSION_CMD" >/dev/null 2>&1; then
    run_legacy_wrapper
  elif command -v rfd3 >/dev/null 2>&1; then
    run_rfd3_cli
  elif command -v apptainer >/dev/null 2>&1; then
    run_apptainer_foundry
  else
    echo "ERROR: No RFdiffusion executable found (rfdiffusion3/rfd3/apptainer)." >&2
    exit 4
  fi
  echo "[INFO] RFdiffusion3 local run complete at $(date --iso-8601=seconds)"
} 2>&1 | tee "$LOG_FILE"

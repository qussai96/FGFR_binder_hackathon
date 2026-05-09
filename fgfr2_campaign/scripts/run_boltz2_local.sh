#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_boltz2_local.sh --config <path> --strategy <name> --num-designs <N> --outdir <path> [--dry-run]

Local-only runner for Boltz-2/validation jobs on NVIDIA A40.
No Lyceum, no S3, no Modal.
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
LOG_FILE="$LOG_DIR/boltz2_${STRATEGY}_$(date +%Y%m%d_%H%M%S).log"

BOLTZ2_CMD="${BOLTZ2_CMD:-boltz2}"
CMD=("$BOLTZ2_CMD" "--config" "$CONFIG" "--num_designs" "$NUM_DESIGNS" "--out" "$OUTDIR")

echo "Strategy: $STRATEGY"
echo "Config: $CONFIG"
echo "Outdir: $OUTDIR"
echo "Log: $LOG_FILE"
echo "Command shape: ${CMD[*]}"
echo "Required env hints: CUDA_VISIBLE_DEVICES (optional), BOLTZ2_CMD (optional override)"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY RUN: not executing Boltz-2 job."
  exit 0
fi

if ! command -v "$BOLTZ2_CMD" >/dev/null 2>&1; then
  echo "ERROR: Boltz-2 command '$BOLTZ2_CMD' not found." >&2
  echo "TODO: install/configure Boltz-2 locally, then rerun this script." >&2
  echo "Expected command: ${CMD[*]}" >&2
  exit 4
fi

{
  echo "[INFO] Starting Boltz-2 local run at $(date --iso-8601=seconds)"
  "${CMD[@]}"
} 2>&1 | tee "$LOG_FILE"

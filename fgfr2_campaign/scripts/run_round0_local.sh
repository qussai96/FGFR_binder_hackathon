#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_SCRIPT="$SCRIPT_DIR/submit_round0_slurm.sh"

DRY_RUN=0
PASS_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      cat <<'USAGE'
Usage: run_round0_local.sh [--dry-run] [submit_round0_slurm.sh options]

This wrapper enforces cluster-only execution and forwards to Slurm submission.
Examples:
  bash fgfr2_campaign/scripts/run_round0_local.sh
  bash fgfr2_campaign/scripts/run_round0_local.sh --num-designs 5 --all-strategies
  bash fgfr2_campaign/scripts/run_round0_local.sh --dry-run
USAGE
      exit 0
      ;;
    *)
      PASS_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
  echo "ERROR: Missing Slurm submit helper: $SUBMIT_SCRIPT" >&2
  exit 2
fi

echo "INFO: Local execution is disabled. Submitting via Slurm GPU-A40 only."

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "DRY RUN: bash $SUBMIT_SCRIPT --prep-only ${PASS_ARGS[*]}"
  exit 0
fi

exec bash "$SUBMIT_SCRIPT" "${PASS_ARGS[@]}"

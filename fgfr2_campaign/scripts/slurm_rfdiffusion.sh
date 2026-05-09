#!/usr/bin/env bash
#SBATCH --job-name=fgfr2_rfd3
#SBATCH --partition=GPU-A40
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

if [[ -n "${CAMPAIGN_DIR:-}" ]]; then
  CAMPAIGN_DIR="$(cd "$CAMPAIGN_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/fgfr2_campaign" ]]; then
  CAMPAIGN_DIR="${SLURM_SUBMIT_DIR}/fgfr2_campaign"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

SCRIPT_DIR="$CAMPAIGN_DIR/scripts"

STRATEGY="${STRATEGY:-strategy1_ligand_site_blocker}"
NUM_DESIGNS="${NUM_DESIGNS:-5}"
NUM_BATCHES="${NUM_BATCHES:-1}"
CONFIG_PATH="${CONFIG_PATH:-$CAMPAIGN_DIR/configs/rfd3/${STRATEGY}_design_spec.json}"
OUTDIR="${OUTDIR:-$CAMPAIGN_DIR/out/rfdiffusion3/${STRATEGY}_round0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strategy) STRATEGY="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --num-batches) NUM_BATCHES="$2"; shift 2 ;;
    *) echo "WARNING: Unknown argument: $1 (ignored)" >&2; shift ;;
  esac
done

CONDA_ENV=/home/students/q.abbas/.conda/envs/rfd3_py312
export PATH="$CONDA_ENV/bin:$PATH"

echo "=== FGFR2 RFdiffusion3 Slurm Job ==="
echo "Job ID   : ${SLURM_JOB_ID:-local}"
echo "Node     : $(hostname)"
echo "Strategy : $STRATEGY"
echo "Designs  : $NUM_DESIGNS"
echo "Batches  : $NUM_BATCHES"
echo "Config   : $CONFIG_PATH"
echo "Outdir   : $OUTDIR"
echo "Python   : $(python --version 2>&1)"
echo "Start    : $(date --iso-8601=seconds)"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "ERROR: RFdiffusion config not found: $CONFIG_PATH" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

nvidia-smi -L || true

bash "$SCRIPT_DIR/run_rfdiffusion_local.sh" \
  --config "$CONFIG_PATH" \
  --strategy "$STRATEGY" \
  --num-designs "$NUM_DESIGNS" \
  --num-batches "$NUM_BATCHES" \
  --outdir "$OUTDIR"

echo "=== RFdiffusion job complete: $(date --iso-8601=seconds) ==="

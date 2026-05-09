#!/usr/bin/env bash
#SBATCH --job-name=fgfr2_boltzgen
#SBATCH --partition=GPU-A40
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ---------------------------------------------------------------------------
# Slurm GPU-A40 runner for BoltzGen — FGFR2 binder campaign
# Submit with:
#   sbatch slurm_boltzgen.sh --strategy strategy1 --num-designs 5
# Optional overrides via environment variables before sbatch:
#   STRATEGY, NUM_DESIGNS, CONFIG_PATH, OUTDIR
# ---------------------------------------------------------------------------
set -euo pipefail

# In Slurm, scripts can run from a spooled copy under /var/spool/slurm.
# Prefer explicit CAMPAIGN_DIR or SLURM_SUBMIT_DIR to keep paths in workspace.
if [[ -n "${CAMPAIGN_DIR:-}" ]]; then
  CAMPAIGN_DIR="$(cd "$CAMPAIGN_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/fgfr2_campaign" ]]; then
  CAMPAIGN_DIR="${SLURM_SUBMIT_DIR}/fgfr2_campaign"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

LOG_DIR="$CAMPAIGN_DIR/out/logs"
mkdir -p "$LOG_DIR"
  CONVERT_SCRIPT="$CAMPAIGN_DIR/scripts/04_convert_strategy_to_boltzgen_spec.py"

# ---- parse args (passed after sbatch slurm_boltzgen.sh ...) ----------------
STRATEGY="${STRATEGY:-strategy1_ligand_site_blocker}"
NUM_DESIGNS="${NUM_DESIGNS:-5}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strategy)    STRATEGY="$2";     shift 2 ;;
    --num-designs) NUM_DESIGNS="$2";  shift 2 ;;
    *) echo "WARNING: Unknown argument: $1 (ignored)" >&2; shift ;;
  esac
done

CONFIG_PATH="${CONFIG_PATH:-$CAMPAIGN_DIR/configs/${STRATEGY}.yaml}"
OUTDIR="${OUTDIR:-$CAMPAIGN_DIR/out/boltzgen/${STRATEGY}_round0}"
RELAX_FILTERS="${RELAX_FILTERS:-0}"
mkdir -p "$OUTDIR"

CONFIG_BASENAME="$(basename "$CONFIG_PATH")"
if [[ "$CONFIG_BASENAME" == *_design_spec.yaml || "$CONFIG_BASENAME" == *_design_spec.yml ]]; then
  NATIVE_CONFIG_PATH="$CONFIG_PATH"
else
  NATIVE_CONFIG_PATH="$CAMPAIGN_DIR/configs/native/${STRATEGY}_design_spec.yaml"
fi

# ---- activate conda env -----------------------------------------------------
CONDA_ENV=/home/students/q.abbas/anaconda3/envs/FGFR_hack
export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH=""

echo "=== FGFR2 BoltzGen Slurm Job ==="
echo "Job ID   : ${SLURM_JOB_ID:-local}"
echo "Node     : $(hostname)"
echo "Strategy : $STRATEGY"
echo "Designs  : $NUM_DESIGNS"
echo "Config   : $CONFIG_PATH"
echo "Native   : $NATIVE_CONFIG_PATH"
echo "Outdir   : $OUTDIR"
echo "Relaxed filters: $RELAX_FILTERS"
echo "Python   : $(python --version 2>&1)"
echo "Start    : $(date --iso-8601=seconds)"

# ---- GPU diagnostics --------------------------------------------------------
echo ""
echo "--- GPU Info ---"
nvidia-smi -L
python -c "
import torch
print(f'PyTorch: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA device: {torch.cuda.get_device_name(0)}')
    print(f'VRAM: {torch.cuda.get_device_properties(0).total_memory // 1024**3} GB')
"
echo ""

# ---- validate config exists -------------------------------------------------
if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "ERROR: Config not found: $CONFIG_PATH" >&2
  echo "Run the prep steps first: bash submit_round0_slurm.sh --prep-only" >&2
  exit 1
fi

if [[ "$NATIVE_CONFIG_PATH" != "$CONFIG_PATH" ]]; then
  if [[ ! -f "$CONVERT_SCRIPT" ]]; then
    echo "ERROR: Missing converter script: $CONVERT_SCRIPT" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$NATIVE_CONFIG_PATH")"
  echo "Converting strategy YAML to native BoltzGen design-spec..."
  python "$CONVERT_SCRIPT" \
    --strategy-config "$CONFIG_PATH" \
    --output "$NATIVE_CONFIG_PATH"
fi

# ---- run boltzgen -----------------------------------------------------------
LOG_FILE="$LOG_DIR/boltzgen_${STRATEGY}_${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}.log"
echo "Log: $LOG_FILE"

{
  echo "[INFO] BoltzGen start: $(date --iso-8601=seconds)"
  if [[ "$RELAX_FILTERS" == "1" ]]; then
    boltzgen run \
      "$NATIVE_CONFIG_PATH" \
      --num_designs "$NUM_DESIGNS" \
      --output "$OUTDIR" \
      --config filtering filter_designfolding=false filter_biased=false
  else
    boltzgen run \
      "$NATIVE_CONFIG_PATH" \
      --num_designs "$NUM_DESIGNS" \
      --output "$OUTDIR"
  fi
  echo "[INFO] BoltzGen done: $(date --iso-8601=seconds)"
} 2>&1 | tee "$LOG_FILE"

echo ""
echo "=== Job complete: $(date --iso-8601=seconds) ==="
echo "Outputs in: $OUTDIR"

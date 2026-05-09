#!/usr/bin/env bash
#SBATCH --job-name=fgfr2_round0
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
ROOT_DIR="$(cd "$CAMPAIGN_DIR/.." && pwd)"
LOG_DIR="$CAMPAIGN_DIR/out/logs"
mkdir -p "$LOG_DIR"

CONDA_ENV=/home/students/q.abbas/anaconda3/envs/FGFR_hack
export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH=""

NUM_DESIGNS="${NUM_DESIGNS:-5}"
ALL_STRATEGIES="${ALL_STRATEGIES:-0}"
RELAX_FILTERS="${RELAX_FILTERS:-0}"

RECEPTOR_CHAIN="A"
LIGAND_CHAIN="B"
OFFTARGET_CHAIN="C"

echo "=== FGFR2 Round0 Slurm Pipeline ==="
echo "Job ID     : ${SLURM_JOB_ID:-n/a}"
echo "Campaign   : $CAMPAIGN_DIR"
echo "Designs    : $NUM_DESIGNS"
echo "All strats : $ALL_STRATEGIES"
echo "Relaxed    : $RELAX_FILTERS"
echo "Start      : $(date --iso-8601=seconds)"

echo "--- GPU Info ---"
nvidia-smi -L || true

python - <<'PY'
import torch
print(f"PyTorch: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"CUDA device: {torch.cuda.get_device_name(0)}")
PY

echo "[1/7] Ensuring structures..."
bash "$SCRIPT_DIR/00_download_target.sh"

echo "[2/7] Building FGFR2-FGF1 interface hotspots..."
python "$SCRIPT_DIR/01_find_interface_residues.py" \
  --structure "$CAMPAIGN_DIR/structures/1DJS.pdb" \
  --receptor-chain "$RECEPTOR_CHAIN" \
  --ligand-chain "$LIGAND_CHAIN" \
  --cutoff 5.0 \
  --output "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json"

echo "[3/7] Building specificity filter..."
if [[ -f "$CAMPAIGN_DIR/structures/1CVS_offtarget_fgfr1.pdb" ]]; then
  python "$SCRIPT_DIR/03_build_specificity_filter.py" \
    --hotspots-json "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json" \
    --target-structure "$CAMPAIGN_DIR/structures/1DJS.pdb" \
    --target-chain "$RECEPTOR_CHAIN" \
    --offtarget-structure "$CAMPAIGN_DIR/structures/1CVS_offtarget_fgfr1.pdb" \
    --offtarget-chain "$OFFTARGET_CHAIN" \
    --output "$CAMPAIGN_DIR/configs/fgfr2_vs_fgfr1_specificity.json"
else
  echo "WARNING: Missing 1CVS_offtarget_fgfr1.pdb; skipping specificity filter."
fi

echo "[4/7] Generating strategy configs..."
python "$SCRIPT_DIR/02_generate_strategy_configs.py" \
  --campaign-dir "$CAMPAIGN_DIR" \
  --hotspots-json "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json" \
  --specificity-json "$CAMPAIGN_DIR/configs/fgfr2_vs_fgfr1_specificity.json"

STRATEGIES=(strategy1_ligand_site_blocker)
if [[ "$ALL_STRATEGIES" == "1" ]]; then
  STRATEGIES=(
    strategy1_ligand_site_blocker
    strategy2_dimer_or_activation_blocker
    strategy3_surface_explorer
  )
fi

echo "[5/7] Running BoltzGen for selected strategy set..."
for STRAT in "${STRATEGIES[@]}"; do
  STRAT_CFG="$CAMPAIGN_DIR/configs/${STRAT}.yaml"
  STRAT_OUT="$CAMPAIGN_DIR/out/boltzgen/${STRAT}_round0"
  mkdir -p "$STRAT_OUT"

  NATIVE_CFG="$CAMPAIGN_DIR/configs/native/${STRAT}_design_spec.yaml"
  python "$SCRIPT_DIR/04_convert_strategy_to_boltzgen_spec.py" \
    --strategy-config "$STRAT_CFG" \
    --output "$NATIVE_CFG"

  if [[ "$RELAX_FILTERS" == "1" ]]; then
    boltzgen run \
      "$NATIVE_CFG" \
      --num_designs "$NUM_DESIGNS" \
      --output "$STRAT_OUT" \
      --config filtering filter_designfolding=false filter_biased=false
  else
    boltzgen run \
      "$NATIVE_CFG" \
      --num_designs "$NUM_DESIGNS" \
      --output "$STRAT_OUT"
  fi
done

echo "[6/7] Syncing campaign artifacts..."
python "$CAMPAIGN_DIR/sync_designs_local.py" --campaign-dir "$CAMPAIGN_DIR"

echo "[7/7] Regenerating dashboard..."
python "$CAMPAIGN_DIR/generate_pages_local.py" \
  --campaign-dir "$CAMPAIGN_DIR" \
  --docs-root "$ROOT_DIR/docs"

echo "Done at $(date --iso-8601=seconds)"
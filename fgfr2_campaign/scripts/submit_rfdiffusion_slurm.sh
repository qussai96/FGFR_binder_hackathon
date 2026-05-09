#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# FGFR2 RFdiffusion Round0 - Slurm submission helper
#
# Usage:
#   bash submit_rfdiffusion_slurm.sh [--num-designs N] [--num-batches M] [--all-strategies] [--prep-only]
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CONDA_ENV=/home/students/q.abbas/anaconda3/envs/FGFR_hack
PYTHON="$CONDA_ENV/bin/python"

NUM_DESIGNS=5
NUM_BATCHES=1
ALL_STRATEGIES=0
PREP_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --num-batches) NUM_BATCHES="$2"; shift 2 ;;
    --all-strategies) ALL_STRATEGIES=1; shift ;;
    --prep-only) PREP_ONLY=1; shift ;;
    -h|--help)
      echo "Usage: submit_rfdiffusion_slurm.sh [--num-designs N] [--num-batches M] [--all-strategies] [--prep-only]"
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1" >&2; exit 2 ;;
  esac
done

RECEPTOR_CHAIN="A"
LIGAND_CHAIN="B"
OFFTARGET_CHAIN="C"

echo "=== FGFR2 RFdiffusion Round0 - Slurm Submission ==="
echo "Campaign dir : $CAMPAIGN_DIR"
echo "Python env   : $PYTHON"
echo "Designs/strat: $NUM_DESIGNS"
echo "Batches      : $NUM_BATCHES"
echo ""

echo "[1/4] Checking structures..."
bash "$SCRIPT_DIR/00_download_target.sh"

echo "[2/4] Computing interface hotspots (FGFR2:$RECEPTOR_CHAIN / FGF1:$LIGAND_CHAIN)..."
"$PYTHON" "$SCRIPT_DIR/01_find_interface_residues.py" \
  --structure "$CAMPAIGN_DIR/structures/1DJS.pdb" \
  --receptor-chain "$RECEPTOR_CHAIN" \
  --ligand-chain "$LIGAND_CHAIN" \
  --cutoff 5.0 \
  --output "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json"

echo "[3/4] Building FGFR1 specificity filter (off-target chain: $OFFTARGET_CHAIN)..."
if [[ -f "$CAMPAIGN_DIR/structures/1CVS_offtarget_fgfr1.pdb" ]]; then
  "$PYTHON" "$SCRIPT_DIR/03_build_specificity_filter.py" \
    --hotspots-json "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json" \
    --target-structure "$CAMPAIGN_DIR/structures/1DJS.pdb" \
    --target-chain "$RECEPTOR_CHAIN" \
    --offtarget-structure "$CAMPAIGN_DIR/structures/1CVS_offtarget_fgfr1.pdb" \
    --offtarget-chain "$OFFTARGET_CHAIN" \
    --output "$CAMPAIGN_DIR/configs/fgfr2_vs_fgfr1_specificity.json"
else
  echo "WARNING: 1CVS_offtarget_fgfr1.pdb not found; skipping specificity filter."
fi

echo "[4/4] Generating strategy YAML configs..."
"$PYTHON" "$SCRIPT_DIR/02_generate_strategy_configs.py" \
  --campaign-dir "$CAMPAIGN_DIR" \
  --hotspots-json "$CAMPAIGN_DIR/configs/fgfr2_fgf1_interface_hotspots.json" \
  --specificity-json "$CAMPAIGN_DIR/configs/fgfr2_vs_fgfr1_specificity.json"

echo ""
echo "Prep steps complete."

if [[ "$PREP_ONLY" -eq 1 ]]; then
  echo "--prep-only: skipping Slurm submission."
  exit 0
fi

STRATEGIES=(strategy1_ligand_site_blocker)
if [[ "$ALL_STRATEGIES" -eq 1 ]]; then
  STRATEGIES=(
    strategy1_ligand_site_blocker
    strategy2_dimer_or_activation_blocker
    strategy3_surface_explorer
  )
fi

mkdir -p "$CAMPAIGN_DIR/configs/rfd3"

JOB_IDS=()
for STRAT in "${STRATEGIES[@]}"; do
  CFG="$CAMPAIGN_DIR/configs/${STRAT}.yaml"
  if [[ ! -f "$CFG" ]]; then
    echo "WARNING: Config not found, skipping: $CFG" >&2
    continue
  fi

  RFD3_CFG="$CAMPAIGN_DIR/configs/rfd3/${STRAT}_design_spec.json"
  "$PYTHON" "$SCRIPT_DIR/06_convert_strategy_to_rfd3_spec.py" \
    --strategy-config "$CFG" \
    --output "$RFD3_CFG"

  JOB_ID=$(sbatch \
    --job-name="fgfr2_rfd3_${STRAT}" \
    --export=ALL,STRATEGY="$STRAT",NUM_DESIGNS="$NUM_DESIGNS",NUM_BATCHES="$NUM_BATCHES",CAMPAIGN_DIR="$CAMPAIGN_DIR",CONFIG_PATH="$RFD3_CFG" \
    "$SCRIPT_DIR/slurm_rfdiffusion.sh" \
    --strategy "$STRAT" --num-designs "$NUM_DESIGNS" --num-batches "$NUM_BATCHES" \
    | awk '{print $NF}')
  JOB_IDS+=("$JOB_ID")
  echo "Submitted $STRAT -> job $JOB_ID"
done

echo ""
echo "=== Submitted ${#JOB_IDS[@]} RFdiffusion job(s) to GPU-A40 ==="
echo "Monitor with:  squeue -u $USER"
echo "Logs in:       $CAMPAIGN_DIR/out/logs/"
echo "Outputs in:    $CAMPAIGN_DIR/out/rfdiffusion3/"

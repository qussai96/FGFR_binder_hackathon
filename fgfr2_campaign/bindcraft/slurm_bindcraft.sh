#!/usr/bin/env bash
#SBATCH --job-name=bc_fgfr2_ipsae
#SBATCH --partition=GPU-A40
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# ============================================================================
# Slurm wrapper for run_bindcraft_fgfr2.py
# ipSAE-maximizing BindCraft pipeline for FGFR2 selective binder design
# ============================================================================

set -euo pipefail

# --- conda env (BindCraft + PyRosetta + ColabDesign + AlphaFold2 params) ---
# Adjust the conda root path if yours differs.
CONDA_ROOT="${CONDA_ROOT:-/home/students/q.abbas/anaconda3}"
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate "${BINDCRAFT_ENV:-FGFR_hack}"

# --- working directory ---
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# --- inputs (the only two things the user changes) ---
export FGFR2_INPUT_PDB="${FGFR2_INPUT_PDB:-fgfr2_campaign/structures/FGFR2_1DJS_chainA_clean.pdb}"
export FGFR2_HOTSPOTS="${FGFR2_HOTSPOTS:-A283,A251,A281,A346,A317,A173}"

# --- off-target for selectivity gap ---
export FGFR1_OFFTARGET_PDB="${FGFR1_OFFTARGET_PDB:-fgfr2_campaign/structures/FGFR1_1CVS_chainC_clean.pdb}"

# --- output path (timestamped to avoid overwriting prior runs) ---
TS="$(date +%Y%m%d_%H%M%S)"
export FGFR2_DESIGN_PATH="${FGFR2_DESIGN_PATH:-fgfr2_campaign/out/bindcraft/fgfr2_ipsae_${TS}/}"

# --- BindCraft installation path ---
export BINDCRAFT_FOLDER="${BINDCRAFT_FOLDER:-$HOME/bindcraft}"

# --- log header ---
echo "=== BindCraft FGFR2 ipSAE run ==="
echo "Job ID:           ${SLURM_JOB_ID:-local}"
echo "Node:             $(hostname)"
echo "Conda env:        $(conda info --envs | grep '*' || true)"
echo "Input PDB:        $FGFR2_INPUT_PDB"
echo "Hotspots:         $FGFR2_HOTSPOTS"
echo "Off-target PDB:   $FGFR1_OFFTARGET_PDB"
echo "Design path:      $FGFR2_DESIGN_PATH"
echo "BindCraft folder: $BINDCRAFT_FOLDER"
echo "Start:            $(date --iso-8601=seconds)"

nvidia-smi -L || true

# --- run the pipeline ---
python fgfr2_campaign/bindcraft/run_bindcraft_fgfr2.py

echo "End:              $(date --iso-8601=seconds)"

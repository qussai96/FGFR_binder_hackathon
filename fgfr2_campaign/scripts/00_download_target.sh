#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STRUCT_DIR="$CAMPAIGN_DIR/structures"
ROOT_DIR="$(cd "$CAMPAIGN_DIR/.." && pwd)"
TARGET_URL="https://files.rcsb.org/download/1DJS.pdb"
OFFTARGET_URL="https://files.rcsb.org/download/1CVS.pdb"
RAW_PDB="$STRUCT_DIR/1DJS.pdb"
CLEAN_PDB="$STRUCT_DIR/fgfr2_target_clean.pdb"
OFFTARGET_PDB="$STRUCT_DIR/1CVS_offtarget_fgfr1.pdb"
LOCAL_TARGET_PDB="$ROOT_DIR/1DJS.pdb"
LOCAL_OFFTARGET_PDB="$ROOT_DIR/1CVS.pdb"
INSPECT_SCRIPT="$SCRIPT_DIR/inspect_cif.py"

mkdir -p "$STRUCT_DIR"

if [[ -f "$LOCAL_TARGET_PDB" ]]; then
  echo "Using local FGFR2 target PDB: $LOCAL_TARGET_PDB"
  cp "$LOCAL_TARGET_PDB" "$RAW_PDB"
else
  if command -v curl >/dev/null 2>&1; then
    echo "Downloading 1DJS from RCSB with curl..."
    curl -fsSL "$TARGET_URL" -o "$RAW_PDB"
  elif command -v wget >/dev/null 2>&1; then
    echo "Downloading 1DJS from RCSB with wget..."
    wget -q "$TARGET_URL" -O "$RAW_PDB"
  else
    echo "ERROR: Neither curl nor wget is available, and local 1DJS.pdb is missing." >&2
    exit 2
  fi
fi

if [[ -f "$LOCAL_OFFTARGET_PDB" ]]; then
  echo "Using local FGFR1 off-target PDB: $LOCAL_OFFTARGET_PDB"
  cp "$LOCAL_OFFTARGET_PDB" "$OFFTARGET_PDB"
elif command -v curl >/dev/null 2>&1; then
  echo "Downloading 1CVS from RCSB with curl..."
  curl -fsSL "$OFFTARGET_URL" -o "$OFFTARGET_PDB"
elif command -v wget >/dev/null 2>&1; then
  echo "Downloading 1CVS from RCSB with wget..."
  wget -q "$OFFTARGET_URL" -O "$OFFTARGET_PDB"
else
  echo "WARNING: Could not stage FGFR1 off-target structure (1CVS). Specificity filtering will be limited." >&2
fi

echo "Saved FGFR2 target structure to: $RAW_PDB"
if [[ -f "$OFFTARGET_PDB" ]]; then
  echo "Saved FGFR1 off-target structure to: $OFFTARGET_PDB"
fi

# Cleaning can be tool-dependent; use a safe fallback copy so downstream paths are stable.
cp "$RAW_PDB" "$CLEAN_PDB"
echo "Prepared cleaned target placeholder at: $CLEAN_PDB"

if command -v python3 >/dev/null 2>&1; then
  echo "Inspecting FGFR2 chain/entity information..."
  python3 "$INSPECT_SCRIPT" "$RAW_PDB" --json "$STRUCT_DIR/1DJS_inspection.json" || {
    echo "WARNING: Detailed inspection failed. Install BioPython for best results: pip install biopython" >&2
  }
  if [[ -f "$OFFTARGET_PDB" ]]; then
    echo "Inspecting FGFR1 off-target chain/entity information..."
    python3 "$INSPECT_SCRIPT" "$OFFTARGET_PDB" --json "$STRUCT_DIR/1CVS_inspection.json" || true
  fi
else
  echo "WARNING: python3 not found; skipping detailed structure inspection." >&2
  exit 2
fi

echo "Done. Do not assume receptor/ligand chain IDs until inspection output is reviewed."

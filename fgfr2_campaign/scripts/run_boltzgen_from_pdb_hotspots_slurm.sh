#!/usr/bin/env bash
#SBATCH --job-name=hotspots
#SBATCH --partition=GPU-A40
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash fgfr2_campaign/scripts/run_boltzgen_from_pdb_hotspots_slurm.sh \
    --pdb <cleaned_target.pdb> \
    --receptor-chain <CHAIN_ID> \
    --hotspots <A:160,A:163,... | hotspots.json | hotspots.txt> \
    [--num-designs 50] \
    [--min-length 30] \
    [--max-length 55] \
    [--outdir fgfr2_campaign/out/boltzgen/custom_run] \
    [--job-name fgfr2_hotspots] \
    [--partition GPU-A40] \
    [--dry-run]

What it does:
- Submits a Slurm GPU-A40 job (no local BoltzGen execution).
- Builds a native BoltzGen design-spec from your cleaned PDB + hotspots.
- Runs BoltzGen and writes final candidate metrics.
- Produces filtered candidate outputs for sequences with length < max-length:
  - filtered_candidates_len_lt_<max>.csv
  - filtered_candidates_len_lt_<max>.json
USAGE
}

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
# In Slurm, scripts can run from a spooled copy under /var/spool/slurm.
# Prefer explicit CAMPAIGN_DIR or SLURM_SUBMIT_DIR to keep paths in workspace.
if [[ -n "${CAMPAIGN_DIR:-}" ]]; then
    CAMPAIGN_DIR="$(cd "$CAMPAIGN_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/fgfr2_campaign" ]]; then
    CAMPAIGN_DIR="${SLURM_SUBMIT_DIR}/fgfr2_campaign"
else
    CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

PDB_PATH=""
RECEPTOR_CHAIN="A"
HOTSPOTS_INPUT=""
NUM_DESIGNS="50"
MIN_LENGTH="30"
MAX_LENGTH="55"
OUTDIR=""
JOB_NAME="fgfr2_hotspots"
PARTITION="GPU-A40"
DRY_RUN=0
RUN_INTERNAL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pdb) PDB_PATH="$2"; shift 2 ;;
    --receptor-chain) RECEPTOR_CHAIN="$2"; shift 2 ;;
    --hotspots) HOTSPOTS_INPUT="$2"; shift 2 ;;
    --num-designs) NUM_DESIGNS="$2"; shift 2 ;;
    --min-length) MIN_LENGTH="$2"; shift 2 ;;
    --max-length) MAX_LENGTH="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --run-internal) RUN_INTERNAL=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$PDB_PATH" || -z "$HOTSPOTS_INPUT" ]]; then
  echo "ERROR: --pdb and --hotspots are required." >&2
  usage
  exit 2
fi

if ! [[ "$NUM_DESIGNS" =~ ^[0-9]+$ ]] || [[ "$NUM_DESIGNS" -lt 1 ]]; then
  echo "ERROR: --num-designs must be a positive integer." >&2
  exit 2
fi
if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MIN_LENGTH" -lt 1 ]]; then
  echo "ERROR: --min-length must be a positive integer." >&2
  exit 2
fi
if ! [[ "$MAX_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MAX_LENGTH" -lt 2 ]]; then
  echo "ERROR: --max-length must be an integer >= 2." >&2
  exit 2
fi
if [[ "$MIN_LENGTH" -ge "$MAX_LENGTH" ]]; then
  echo "ERROR: --min-length must be < --max-length." >&2
  exit 2
fi

PDB_ABS="$(readlink -f "$PDB_PATH")"
if [[ ! -f "$PDB_ABS" ]]; then
  echo "ERROR: PDB file not found: $PDB_ABS" >&2
  exit 2
fi

if [[ "$HOTSPOTS_INPUT" == *":"* || "$HOTSPOTS_INPUT" == *,* ]]; then
  HOTSPOTS_ABS="$HOTSPOTS_INPUT"
else
  HOTSPOTS_ABS="$(readlink -f "$HOTSPOTS_INPUT")"
  if [[ ! -f "$HOTSPOTS_ABS" ]]; then
    echo "ERROR: Hotspots file not found: $HOTSPOTS_ABS" >&2
    exit 2
  fi
fi

if [[ -z "$OUTDIR" ]]; then
  TS="$(date +%Y%m%d_%H%M%S)"
  OUTDIR="$CAMPAIGN_DIR/out/boltzgen/custom_hotspots_${TS}"
fi
OUTDIR_ABS="$(readlink -m "$OUTDIR")"

if [[ "$RUN_INTERNAL" -eq 0 ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found. This script is Slurm-only." >&2
    exit 3
  fi

  SBATCH_CMD=(
    sbatch
    --job-name "$JOB_NAME"
    --partition "$PARTITION"
    --gres gpu:1
    --cpus-per-task 8
    --mem 48G
    --time 12:00:00
    --output "%x_%j.out"
    --error "%x_%j.err"
    "$SCRIPT_PATH"
    --run-internal
    --pdb "$PDB_ABS"
    --receptor-chain "$RECEPTOR_CHAIN"
    --hotspots "$HOTSPOTS_ABS"
    --num-designs "$NUM_DESIGNS"
    --min-length "$MIN_LENGTH"
    --max-length "$MAX_LENGTH"
    --outdir "$OUTDIR_ABS"
    --job-name "$JOB_NAME"
    --partition "$PARTITION"
  )

  echo "Submitting Slurm job (GPU-A40 workflow)..."
  echo "PDB: $PDB_ABS"
  echo "Hotspots: $HOTSPOTS_INPUT"
  echo "Num designs: $NUM_DESIGNS"
  echo "Length window: $MIN_LENGTH..$MAX_LENGTH (final report filters length < $MAX_LENGTH)"
  echo "Outdir: $OUTDIR_ABS"
  echo "Command: ${SBATCH_CMD[*]}"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY RUN: not submitting."
    exit 0
  fi

  "${SBATCH_CMD[@]}"
  exit 0
fi

# -------- Slurm internal execution path --------
CONDA_ENV="/home/students/q.abbas/anaconda3/envs/FGFR_hack"
export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH=""

mkdir -p "$OUTDIR_ABS"
mkdir -p "$CAMPAIGN_DIR/configs/native"

RUN_TAG="custom_hotspots_${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}"
SPEC_PATH="$CAMPAIGN_DIR/configs/native/${RUN_TAG}_design_spec.yaml"

printf "=== BoltzGen from PDB + hotspots ===\n"
printf "Job ID: %s\n" "${SLURM_JOB_ID:-local}"
printf "Node: %s\n" "$(hostname)"
printf "PDB: %s\n" "$PDB_ABS"
printf "Receptor chain: %s\n" "$RECEPTOR_CHAIN"
printf "Hotspots: %s\n" "$HOTSPOTS_INPUT"
printf "Designs: %s\n" "$NUM_DESIGNS"
printf "Length range: %s..%s\n" "$MIN_LENGTH" "$MAX_LENGTH"
printf "Output: %s\n" "$OUTDIR_ABS"
printf "Spec: %s\n" "$SPEC_PATH"
printf "Start: %s\n" "$(date --iso-8601=seconds)"

nvidia-smi -L || true

python - <<'PY' "$PDB_ABS" "$RECEPTOR_CHAIN" "$HOTSPOTS_ABS" "$MIN_LENGTH" "$MAX_LENGTH" "$SPEC_PATH"
import json
import os
import random
import re
import sys
from typing import Dict, List

import yaml
from Bio.PDB import PDBParser

pdb_path = sys.argv[1]
receptor_chain = sys.argv[2]
hotspots_input = sys.argv[3]
min_len = int(sys.argv[4])
max_len = int(sys.argv[5])
spec_path = sys.argv[6]


def is_aa_residue(residue) -> bool:
    return residue.id[0] == " "


def map_auth_to_pos(structure_path: str, chain_id: str) -> Dict[int, int]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", structure_path)
    model = next(structure.get_models())
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in {structure_path}")
    chain = model[chain_id]
    mapping: Dict[int, int] = {}
    pos = 0
    for residue in chain:
        if not is_aa_residue(residue):
            continue
        pos += 1
        auth_seq = int(residue.id[1])
        if auth_seq not in mapping:
            mapping[auth_seq] = pos
    return mapping


def parse_hotspots_token(token: str, default_chain: str) -> List[str]:
    vals = []
    aa_one_letter = set("ACDEFGHIKLMNPQRSTVWY")
    for raw in token.split(","):
        item = raw.strip()
        if not item:
            continue
        if ":" in item:
            vals.append(item)
        elif re.fullmatch(r"([A-Za-z])(\d+)", item):
            m = re.fullmatch(r"([A-Za-z])(\d+)", item)
            assert m is not None
            prefix = m.group(1).upper()
            # Support both chain-style (A283) and residue-style (V317) tokens.
            # If prefix is an amino-acid code and not the receptor chain, treat
            # it as residue annotation and map to the default receptor chain.
            if prefix in aa_one_letter and prefix != str(default_chain).upper():
                vals.append(f"{default_chain}:{m.group(2)}")
            else:
                vals.append(f"{prefix}:{m.group(2)}")
        elif re.fullmatch(r"[0-9]+", item):
            vals.append(f"{default_chain}:{item}")
    return vals


def parse_hotspots(value: str, default_chain: str) -> List[str]:
    if os.path.isfile(value):
        if value.lower().endswith(".json"):
            with open(value, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            out: List[str] = []
            if isinstance(data, dict) and "interface_residues" in data:
                for item in data.get("interface_residues", []):
                    chain = item.get("auth_asym_id") or item.get("label_asym_id") or default_chain
                    seq = item.get("auth_seq_id")
                    if seq is None:
                        seq = item.get("label_seq_id")
                    if chain and seq is not None:
                        out.append(f"{chain}:{int(seq)}")
                return sorted(set(out))
            if isinstance(data, dict) and "hotspots" in data and isinstance(data["hotspots"], list):
                for item in data["hotspots"]:
                    if isinstance(item, str):
                        out.extend(parse_hotspots_token(item, default_chain))
                return sorted(set(out))
            if isinstance(data, list):
                for item in data:
                    if isinstance(item, str):
                        out.extend(parse_hotspots_token(item, default_chain))
                    elif isinstance(item, dict):
                        c = item.get("chain", default_chain)
                        s = item.get("seq")
                        if s is not None:
                            out.append(f"{c}:{int(s)}")
                return sorted(set(out))
            raise ValueError("Unsupported hotspots JSON format")
        # Plain text file with one hotspot per line or comma-separated.
        with open(value, "r", encoding="utf-8") as handle:
            text = handle.read().strip()
        tokens = re.split(r"[\n\s]+", text.replace(";", ","))
        out: List[str] = []
        for token in tokens:
            if not token:
                continue
            out.extend(parse_hotspots_token(token, default_chain))
        return sorted(set(out))

    return sorted(set(parse_hotspots_token(value, default_chain)))


all_hotspots = parse_hotspots(hotspots_input, receptor_chain)
selected_auth = []
for h in all_hotspots:
    if ":" not in h:
        continue
    c, s = h.split(":", 1)
    if c.strip() != receptor_chain:
        continue
    try:
        selected_auth.append(int(s))
    except ValueError:
        continue
selected_auth = sorted(set(selected_auth))

if len(selected_auth) > 5:
    selected_auth = sorted(random.sample(selected_auth, 5))

seq_map = map_auth_to_pos(pdb_path, receptor_chain)
mapped_positions = [seq_map[x] for x in selected_auth if x in seq_map]

if not mapped_positions:
    raise ValueError("No hotspots mapped to receptor chain sequence positions")

spec = {
    "entities": [
        {
            "file": {
                "path": os.path.abspath(pdb_path),
                "include": [{"chain": {"id": receptor_chain}}],
                "not_design": [{"chain": {"id": receptor_chain}}],
                "binding_types": [
                    {
                        "chain": {
                            "id": receptor_chain,
                            "binding": ",".join(str(x) for x in sorted(set(mapped_positions))),
                        }
                    }
                ],
            }
        },
        {
            "protein": {
                "id": "X",
                "sequence": f"{min_len}..{max_len}",
                "msa": "empty",
            }
        },
    ]
}

os.makedirs(os.path.dirname(spec_path), exist_ok=True)
with open(spec_path, "w", encoding="utf-8") as handle:
    yaml.safe_dump(spec, handle, sort_keys=False)

print(f"Wrote native BoltzGen spec: {spec_path}")
print(f"Hotspots parsed: {len(all_hotspots)}")
print(f"Hotspots used on chain {receptor_chain}: {len(selected_auth)}")
print(f"Selected auth hotspots: {selected_auth}")
print(f"Hotspots mapped to positions: {len(mapped_positions)}")
PY

boltzgen run "$SPEC_PATH" --num_designs "$NUM_DESIGNS" --output "$OUTDIR_ABS"

python - <<'PY' "$OUTDIR_ABS" "$MAX_LENGTH"
import csv
import glob
import json
import os
import sys
from typing import Dict, List

outdir = sys.argv[1]
max_len = int(sys.argv[2])


def _to_float(value: str) -> float | None:
    text = str(value).strip()
    if text in {"", "None", "null", "nan", "NaN"}:
        return None
    try:
        return float(text)
    except Exception:
        return None

candidates = sorted(glob.glob(os.path.join(outdir, "final_ranked_designs", "final_designs_metrics*.csv")))
if not candidates:
    raise FileNotFoundError(f"No final metrics CSV found under {outdir}/final_ranked_designs")

metrics_csv = candidates[0]
rows_out: List[Dict[str, str]] = []

with open(metrics_csv, "r", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        seq = (row.get("designed_sequence") or row.get("sequence") or "").strip()
        if not seq:
            continue
        seq_len = len(seq)
        # User requested below 55 aa style filtering: strict less-than max_len.
        if seq_len >= max_len:
            continue
        rows_out.append(
            {
                "id": row.get("id", ""),
                "final_rank": row.get("final_rank", ""),
                "sequence": seq,
                "sequence_length": str(seq_len),
                "quality_score": row.get("quality_score", ""),
                "design_to_target_iptm": row.get("design_to_target_iptm", ""),
                "min_design_to_target_pae": row.get("min_design_to_target_pae", ""),
                "design_ptm": row.get("design_ptm", ""),
                "design_to_target_ipsae": row.get("design_to_target_ipsae", ""),
                "target_to_design_ipsae": row.get("target_to_design_ipsae", ""),
                "pass_filters": row.get("pass_filters", ""),
            }
        )

def _rank_key(item: Dict[str, str]):
    ipsae = _to_float(item.get("design_to_target_ipsae", ""))
    iptm = _to_float(item.get("design_to_target_iptm", ""))
    ptm = _to_float(item.get("design_ptm", ""))
    quality = _to_float(item.get("quality_score", ""))
    final_rank = _to_float(item.get("final_rank", ""))
    return (
        -(ipsae if ipsae is not None else -1e9),
        -(iptm if iptm is not None else -1e9),
        -(ptm if ptm is not None else -1e9),
        -(quality if quality is not None else -1e9),
        (final_rank if final_rank is not None else 1e9),
    )


rows_out.sort(key=_rank_key)

for idx, row in enumerate(rows_out, start=1):
    row["rank_ipsae_first"] = str(idx)

csv_path = os.path.join(outdir, f"filtered_candidates_len_lt_{max_len}.csv")
json_path = os.path.join(outdir, f"filtered_candidates_len_lt_{max_len}.json")

fieldnames = [
    "rank_ipsae_first",
    "id",
    "final_rank",
    "sequence",
    "sequence_length",
    "quality_score",
    "design_to_target_iptm",
    "min_design_to_target_pae",
    "design_ptm",
    "design_to_target_ipsae",
    "target_to_design_ipsae",
    "pass_filters",
]

with open(csv_path, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows_out)

with open(json_path, "w", encoding="utf-8") as handle:
    json.dump(
        {
            "source_metrics_csv": metrics_csv,
            "ranking_policy": [
                "higher design_to_target_ipsae",
                "higher design_to_target_iptm",
                "higher design_ptm",
                "higher quality_score",
                "lower final_rank",
            ],
            "length_filter": f"sequence_length < {max_len}",
            "num_candidates": len(rows_out),
            "candidates": rows_out,
        },
        handle,
        indent=2,
    )

print(f"Source metrics: {metrics_csv}")
print(f"Filtered candidates: {len(rows_out)}")
print(f"Wrote: {csv_path}")
print(f"Wrote: {json_path}")
PY

printf "Done: %s\n" "$(date --iso-8601=seconds)"

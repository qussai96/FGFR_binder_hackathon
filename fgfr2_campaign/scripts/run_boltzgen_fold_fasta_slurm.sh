#!/usr/bin/env bash
#SBATCH --job-name=fold_fasta
#SBATCH --partition=GPU-A40
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash fgfr2_campaign/scripts/run_boltzgen_fold_fasta_slurm.sh \
    --fasta <input.fasta> \
        [--target-pdb <clean_target.pdb>] \
        [--target-chain A] \
    [--outdir fgfr2_campaign/out/boltzgen/fasta_fold_run] \
    [--partition GPU-A40] \
    [--job-name fgfr2_fold_fasta] \
    [--max-seqs 0] \
    [--min-length 1] \
    [--max-length 2000] \
    [--allow-x] \
    [--dry-run]

What it does:
- Slurm-only submission and execution on GPU-A40.
- Reads a FASTA list of proteins (e.g. RFdiffusion3 + ProteinMPNN outputs).
- Optionally folds each sequence in complex with a target PDB chain.
- Ranks all fold results and writes:
  - fold_ranking.csv
  - fold_ranking.json
USAGE
}

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
if [[ -n "${CAMPAIGN_DIR:-}" ]]; then
    CAMPAIGN_DIR="$(cd "$CAMPAIGN_DIR" && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/fgfr2_campaign" ]]; then
    CAMPAIGN_DIR="${SLURM_SUBMIT_DIR}/fgfr2_campaign"
else
    CAMPAIGN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

FASTA_PATH=""
TARGET_PDB=""
TARGET_CHAIN="A"
OUTDIR=""
PARTITION="GPU-A40"
JOB_NAME="fgfr2_fold_fasta"
MAX_SEQS="0"
MIN_LENGTH="1"
MAX_LENGTH="2000"
ALLOW_X=0
DRY_RUN=0
RUN_INTERNAL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA_PATH="$2"; shift 2 ;;
        --target-pdb) TARGET_PDB="$2"; shift 2 ;;
        --target-chain) TARGET_CHAIN="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --max-seqs) MAX_SEQS="$2"; shift 2 ;;
    --min-length) MIN_LENGTH="$2"; shift 2 ;;
    --max-length) MAX_LENGTH="$2"; shift 2 ;;
    --allow-x) ALLOW_X=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --run-internal) RUN_INTERNAL=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$FASTA_PATH" ]]; then
  echo "ERROR: --fasta is required." >&2
  usage
  exit 2
fi

if ! [[ "$MAX_SEQS" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --max-seqs must be an integer >= 0." >&2
  exit 2
fi
if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MIN_LENGTH" -lt 1 ]]; then
  echo "ERROR: --min-length must be an integer >= 1." >&2
  exit 2
fi
if ! [[ "$MAX_LENGTH" =~ ^[0-9]+$ ]] || [[ "$MAX_LENGTH" -lt "$MIN_LENGTH" ]]; then
  echo "ERROR: --max-length must be an integer >= --min-length." >&2
  exit 2
fi

FASTA_ABS="$(readlink -f "$FASTA_PATH")"
if [[ ! -f "$FASTA_ABS" ]]; then
  echo "ERROR: FASTA file not found: $FASTA_ABS" >&2
  exit 2
fi

TARGET_PDB_ABS=""
if [[ -n "$TARGET_PDB" ]]; then
    TARGET_PDB_ABS="$(readlink -f "$TARGET_PDB")"
    if [[ ! -f "$TARGET_PDB_ABS" ]]; then
        echo "ERROR: Target PDB file not found: $TARGET_PDB_ABS" >&2
        exit 2
    fi
    if [[ -z "$TARGET_CHAIN" ]]; then
        echo "ERROR: --target-chain cannot be empty when --target-pdb is set." >&2
        exit 2
    fi
fi

if [[ -z "$OUTDIR" ]]; then
  TS="$(date +%Y%m%d_%H%M%S)"
  OUTDIR="$CAMPAIGN_DIR/out/boltzgen/fasta_fold_${TS}"
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
    --time 24:00:00
    --output "%x_%j.out"
    --error "%x_%j.err"
    "$SCRIPT_PATH"
    --run-internal
    --fasta "$FASTA_ABS"
    --outdir "$OUTDIR_ABS"
    --partition "$PARTITION"
    --job-name "$JOB_NAME"
    --max-seqs "$MAX_SEQS"
    --min-length "$MIN_LENGTH"
    --max-length "$MAX_LENGTH"
  )

  if [[ "$ALLOW_X" -eq 1 ]]; then
    SBATCH_CMD+=(--allow-x)
  fi
    if [[ -n "$TARGET_PDB_ABS" ]]; then
        SBATCH_CMD+=(--target-pdb "$TARGET_PDB_ABS" --target-chain "$TARGET_CHAIN")
    fi

  echo "Submitting FASTA folding job to Slurm..."
  echo "FASTA: $FASTA_ABS"
  echo "Outdir: $OUTDIR_ABS"
  echo "max-seqs: $MAX_SEQS (0 means all)"
  echo "length filter for accepted sequences: $MIN_LENGTH..$MAX_LENGTH"
  echo "allow X residues: $ALLOW_X"
    if [[ -n "$TARGET_PDB_ABS" ]]; then
        echo "Target PDB: $TARGET_PDB_ABS"
        echo "Target chain: $TARGET_CHAIN"
    fi
  echo "Command: ${SBATCH_CMD[*]}"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY RUN: not submitting."
    exit 0
  fi

  "${SBATCH_CMD[@]}"
  exit 0
fi

CONDA_ENV="/home/students/q.abbas/anaconda3/envs/FGFR_hack"
export PATH="$CONDA_ENV/bin:$PATH"
export PYTHONPATH=""

mkdir -p "$OUTDIR_ABS"
mkdir -p "$OUTDIR_ABS/specs"
mkdir -p "$OUTDIR_ABS/per_sequence"

MANIFEST_JSON="$OUTDIR_ABS/parsed_sequences_manifest.json"
RANKING_CSV="$OUTDIR_ABS/fold_ranking.csv"
RANKING_JSON="$OUTDIR_ABS/fold_ranking.json"

printf "=== BoltzGen FASTA folding ===\n"
printf "Job ID: %s\n" "${SLURM_JOB_ID:-local}"
printf "Node: %s\n" "$(hostname)"
printf "FASTA: %s\n" "$FASTA_ABS"
printf "Outdir: %s\n" "$OUTDIR_ABS"
if [[ -n "$TARGET_PDB_ABS" ]]; then
    printf "Target PDB: %s\n" "$TARGET_PDB_ABS"
    printf "Target chain: %s\n" "$TARGET_CHAIN"
fi
printf "Start: %s\n" "$(date --iso-8601=seconds)"

nvidia-smi -L || true

python - <<'PY' "$FASTA_ABS" "$MANIFEST_JSON" "$MAX_SEQS" "$MIN_LENGTH" "$MAX_LENGTH" "$ALLOW_X"
import json
import re
import random
import sys
from typing import Dict, List

fasta_path = sys.argv[1]
manifest_path = sys.argv[2]
max_seqs = int(sys.argv[3])
min_len = int(sys.argv[4])
max_len = int(sys.argv[5])
allow_x = bool(int(sys.argv[6]))

valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
if allow_x:
    valid_chars.add("X")

entries: List[Dict[str, str]] = []
current_header = None
seq_parts: List[str] = []

with open(fasta_path, "r", encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                entries.append({"header": current_header, "sequence": "".join(seq_parts)})
            current_header = line[1:].strip() or f"seq_{len(entries)+1}"
            seq_parts = []
        else:
            seq_parts.append(line.replace(" ", "").upper())

if current_header is not None:
    entries.append({"header": current_header, "sequence": "".join(seq_parts)})

parsed = []
seen_ids = {}
for i, item in enumerate(entries, start=1):
    header = item["header"]
    seq = item["sequence"].upper()
    seq_len = len(seq)
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", header).strip("_")
    # BoltzGen reserves '_native' in YAML names; avoid it in generated IDs.
    safe = re.sub(r"_native", "_nat", safe, flags=re.IGNORECASE)
    if not safe:
        safe = f"seq_{i}"
    safe = safe[:120]
    seen_ids[safe] = seen_ids.get(safe, 0) + 1
    if seen_ids[safe] > 1:
        safe = f"{safe}_{seen_ids[safe]}"

    invalid = sorted({c for c in seq if c not in valid_chars})
    status = "accepted"
    reason = ""
    if seq_len < min_len or seq_len > max_len:
        status = "skipped"
        reason = f"length_out_of_range_{min_len}_{max_len}"
    elif invalid:
        status = "skipped"
        reason = f"invalid_residues:{''.join(invalid)}"

    parsed.append(
        {
            "index": i,
            "id": safe,
            "header": header,
            "sequence": seq,
            "length": seq_len,
            "status": status,
            "reason": reason,
        }
    )

accepted = [x for x in parsed if x["status"] == "accepted"]
if max_seqs > 0:
    accepted = accepted[:max_seqs]
    accepted_ids = {x["id"] for x in accepted}
    for row in parsed:
        if row["status"] == "accepted" and row["id"] not in accepted_ids:
            row["status"] = "skipped"
            row["reason"] = f"max_seqs_limit:{max_seqs}"

with open(manifest_path, "w", encoding="utf-8") as handle:
    json.dump(
        {
            "fasta_path": fasta_path,
            "max_seqs": max_seqs,
            "min_length": min_len,
            "max_length": max_len,
            "allow_x": allow_x,
            "total_entries": len(parsed),
            "accepted_entries": sum(1 for x in parsed if x["status"] == "accepted"),
            "entries": parsed,
        },
        handle,
        indent=2,
    )

print(f"Wrote manifest: {manifest_path}")
print(f"Total FASTA entries: {len(parsed)}")
print(f"Accepted for folding: {sum(1 for x in parsed if x['status'] == 'accepted')}")
PY

ACCEPTED_IDS=$(python - <<'PY' "$MANIFEST_JSON"
import json
import sys
with open(sys.argv[1], 'r', encoding='utf-8') as h:
    data = json.load(h)
ids = [e['id'] for e in data.get('entries', []) if e.get('status') == 'accepted']
print('\n'.join(ids))
PY
)

if [[ -z "$ACCEPTED_IDS" ]]; then
  echo "No accepted sequences after validation."
  echo "Manifest: $MANIFEST_JSON"
  exit 0
fi

while IFS= read -r SEQ_ID; do
  [[ -z "$SEQ_ID" ]] && continue
  SEQ_OUT="$OUTDIR_ABS/per_sequence/$SEQ_ID"
  SPEC_PATH="$OUTDIR_ABS/specs/${SEQ_ID}_design_spec.yaml"
  mkdir -p "$SEQ_OUT"

  python - <<'PY' "$MANIFEST_JSON" "$SEQ_ID" "$SPEC_PATH" "$TARGET_PDB_ABS" "$TARGET_CHAIN"
import json
import os
import sys
import yaml

manifest_path = sys.argv[1]
seq_id = sys.argv[2]
spec_path = sys.argv[3]
target_pdb = sys.argv[4]
target_chain = sys.argv[5]

with open(manifest_path, 'r', encoding='utf-8') as h:
    data = json.load(h)
entry = next((e for e in data.get('entries', []) if e.get('id') == seq_id), None)
if entry is None:
    raise ValueError(f"Sequence id not found in manifest: {seq_id}")

entities = []
if target_pdb:
    entities.append(
        {
            'file': {
                'path': os.path.abspath(target_pdb),
                'include': [{'chain': {'id': target_chain}}],
                'not_design': [{'chain': {'id': target_chain}}],
            }
        }
    )

entities.append(
    {
        'protein': {
            'id': 'X',
            'sequence': entry['sequence'],
            'msa': 'empty',
        }
    }
)

spec = {'entities': entities}

with open(spec_path, 'w', encoding='utf-8') as h:
    yaml.safe_dump(spec, h, sort_keys=False)
print(f"Wrote spec: {spec_path}")
PY

  echo "[Fold] $SEQ_ID"
    if ! boltzgen run "$SPEC_PATH" --num_designs 1 --output "$SEQ_OUT"; then
        echo "[WARN] BoltzGen failed for $SEQ_ID; marking as failed and continuing." >&2
    fi
done <<< "$ACCEPTED_IDS"

python - <<'PY' "$MANIFEST_JSON" "$OUTDIR_ABS" "$RANKING_CSV" "$RANKING_JSON"
import csv
import glob
import json
import os
import sys
from typing import Dict, List, Optional

manifest_path = sys.argv[1]
outdir = sys.argv[2]
ranking_csv = sys.argv[3]
ranking_json = sys.argv[4]


def _to_float(v) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip()
    if s in {'', 'None', 'null', 'nan', 'NaN'}:
        return None
    try:
        return float(s)
    except Exception:
        return None


def _to_bool(v) -> Optional[bool]:
    s = str(v).strip().lower()
    if s in {'true', '1', 'yes'}:
        return True
    if s in {'false', '0', 'no'}:
        return False
    return None


with open(manifest_path, 'r', encoding='utf-8') as h:
    manifest = json.load(h)

rows: List[Dict[str, object]] = []
for e in manifest.get('entries', []):
    if e.get('status') != 'accepted':
        rows.append(
            {
                'id': e.get('id', ''),
                'header': e.get('header', ''),
                'length': e.get('length', ''),
                'status': e.get('status', 'skipped'),
                'reason': e.get('reason', ''),
                'sequence': e.get('sequence', ''),
            }
        )
        continue

    seq_id = e['id']
    seq_out = os.path.join(outdir, 'per_sequence', seq_id)
    candidates = sorted(glob.glob(os.path.join(seq_out, 'final_ranked_designs', 'final_designs_metrics*.csv')))

    if not candidates:
        rows.append(
            {
                'id': seq_id,
                'header': e.get('header', ''),
                'length': e.get('length', ''),
                'status': 'failed',
                'reason': 'missing_final_metrics_csv',
                'sequence': e.get('sequence', ''),
            }
        )
        continue

    metrics_csv = candidates[0]
    picked = None
    with open(metrics_csv, 'r', encoding='utf-8') as h:
        reader = csv.DictReader(h)
        for row in reader:
            picked = row
            break

    if picked is None:
        rows.append(
            {
                'id': seq_id,
                'header': e.get('header', ''),
                'length': e.get('length', ''),
                'status': 'failed',
                'reason': 'empty_final_metrics_csv',
                'sequence': e.get('sequence', ''),
            }
        )
        continue

    rows.append(
        {
            'id': seq_id,
            'header': e.get('header', ''),
            'length': e.get('length', ''),
            'status': 'done',
            'reason': '',
            'sequence': e.get('sequence', ''),
            'source_metrics_csv': metrics_csv,
            'final_rank': picked.get('final_rank', ''),
            'quality_score': picked.get('quality_score', ''),
            'pass_filters': picked.get('pass_filters', ''),
            'design_ptm': picked.get('design_ptm', picked.get('ptm', '')),
            'design_iptm': picked.get('design_iptm', picked.get('iptm', '')),
            'complex_plddt': picked.get('complex_plddt', ''),
            'design_to_target_iptm': picked.get('design_to_target_iptm', ''),
            'min_design_to_target_pae': picked.get('min_design_to_target_pae', ''),
            'design_to_target_ipsae': picked.get('design_to_target_ipsae', ''),
        }
    )


def rank_key(item: Dict[str, object]):
    status = str(item.get('status', ''))
    if status != 'done':
        return (1, 1, 1e9, 1e9, 1e9, 1e9, 1e9)
    pass_filters = _to_bool(item.get('pass_filters', ''))
    pass_sort = 0 if pass_filters is True else 1
    quality = _to_float(item.get('quality_score'))
    ptm = _to_float(item.get('design_ptm'))
    iptm = _to_float(item.get('design_iptm'))
    plddt = _to_float(item.get('complex_plddt'))
    final_rank = _to_float(item.get('final_rank'))
    return (
        0,
        pass_sort,
        -(quality if quality is not None else -1e9),
        -(ptm if ptm is not None else -1e9),
        -(iptm if iptm is not None else -1e9),
        -(plddt if plddt is not None else -1e9),
        (final_rank if final_rank is not None else 1e9),
    )


rankable = [r for r in rows if str(r.get('status', '')) == 'done']
rankable.sort(key=rank_key)
for i, row in enumerate(rankable, start=1):
    row['rank_overall'] = i

failed = [r for r in rows if str(r.get('status', '')) != 'done']
rows_sorted = rankable + failed

fieldnames = [
    'rank_overall',
    'id',
    'header',
    'status',
    'reason',
    'length',
    'quality_score',
    'pass_filters',
    'design_ptm',
    'design_iptm',
    'complex_plddt',
    'final_rank',
    'design_to_target_iptm',
    'min_design_to_target_pae',
    'design_to_target_ipsae',
    'source_metrics_csv',
    'sequence',
]

os.makedirs(os.path.dirname(ranking_csv), exist_ok=True)
with open(ranking_csv, 'w', encoding='utf-8', newline='') as h:
    writer = csv.DictWriter(h, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows_sorted:
        writer.writerow({k: row.get(k, '') for k in fieldnames})

with open(ranking_json, 'w', encoding='utf-8') as h:
    json.dump(
        {
            'fasta_path': manifest.get('fasta_path'),
            'total_entries': manifest.get('total_entries'),
            'accepted_entries': manifest.get('accepted_entries'),
            'completed_entries': len(rankable),
            'failed_or_skipped_entries': len(failed),
            'ranking_policy': [
                'status=done first',
                'pass_filters=true first',
                'higher quality_score',
                'higher design_ptm',
                'higher design_iptm',
                'higher complex_plddt',
                'lower final_rank',
            ],
            'results': rows_sorted,
        },
        h,
        indent=2,
    )

print(f"Wrote ranking CSV: {ranking_csv}")
print(f"Wrote ranking JSON: {ranking_json}")
print(f"Completed entries: {len(rankable)}")
print(f"Failed/skipped entries: {len(failed)}")
PY

printf "Done: %s\n" "$(date --iso-8601=seconds)"
printf "Manifest: %s\n" "$MANIFEST_JSON"
printf "Ranking CSV: %s\n" "$RANKING_CSV"
printf "Ranking JSON: %s\n" "$RANKING_JSON"

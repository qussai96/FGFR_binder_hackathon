#!/usr/bin/env python3
"""Local design synchronization and ranking for FGFR2 campaign.

This script scans local output folders, standardizes design artifacts, and writes
fgfr2_campaign/designs/index.json without requiring S3, Lyceum, or Modal.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from glob import glob
from typing import Any, Dict, List, Optional, Tuple

MAX_BINDER_LENGTH_AA = 50

METRIC_ALIASES = {
    "iptm": ["iptm", "ipTM", "interface_ptm"],
    "target_ipTM": ["target_ipTM", "target_iptm", "fgfr2_ipTM", "fgfr2_iptm"],
    "ptm": ["ptm", "pTM"],
    "target_pTM": ["target_pTM", "target_ptm", "fgfr2_pTM", "fgfr2_ptm"],
    "plddt": ["plddt", "pLDDT", "mean_plddt", "avg_plddt"],
    "min_interaction_pae": ["min_interaction_pae", "interaction_pae_min", "min_pae"],
    "rmsd": ["rmsd", "sc_rmsd", "self_consistency_rmsd"],
    "clash_flag": ["clash", "has_clash", "clash_flag"],
    "binder_length": ["binder_length", "binder_len", "length"],
    "binder_sequence": ["binder_sequence", "binder_seq", "sequence"],
    "off_target_ipTM": ["off_target_ipTM", "offtarget_iptm", "fgfr1_iptm"],
    "ipsae": ["ipsae", "ipSAE", "interface_pSAE", "interface_psae"],
    "target_ipSAE": ["target_ipSAE", "target_ipsae", "fgfr2_ipSAE", "fgfr2_ipsae", "ipsae", "ipSAE"],
    "off_target_ipSAE": ["off_target_ipSAE", "offtarget_ipsae", "fgfr1_ipSAE", "fgfr1_ipsae"],
    "pyrosetta_total_score": ["pyrosetta_total_score", "rosetta_total_score", "total_score"],
    "pyrosetta_interface_dG": ["pyrosetta_interface_dG", "rosetta_interface_dg", "interface_dG", "interface_dg"],
    "pyrosetta_interface_delta_sasa": [
        "pyrosetta_interface_delta_sasa",
        "rosetta_interface_delta_sasa",
        "interface_delta_sasa",
        "dSASA_int",
    ],
    "pyrosetta_packstat": ["pyrosetta_packstat", "rosetta_packstat", "packstat"],
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text in {"", ".", "?", "None", "null", "nan", "NaN"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _safe_bool(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    return None


def _get_first_metric(metrics: Dict[str, Any], aliases: List[str]) -> Any:
    for key in aliases:
        if key in metrics:
            return metrics[key]
        for m_key, m_value in metrics.items():
            if m_key.lower() == key.lower():
                return m_value
    return None


def _normalize_key(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", key)


def _detect_strategy(path: str) -> str:
    lower = path.lower()
    if "strategy1" in lower:
        return "strategy1_ligand_site_blocker"
    if "strategy2" in lower:
        return "strategy2_dimer_or_activation_blocker"
    if "strategy3" in lower:
        return "strategy3_surface_explorer"
    return "unknown"


def _make_design_id(tool: str, relative_path: str) -> str:
    stem = os.path.splitext(os.path.basename(relative_path))[0]
    slug = _normalize_key(stem)[:40]
    h = hashlib.sha1(f"{tool}:{relative_path}".encode("utf-8")).hexdigest()[:10]
    return f"{tool}_{slug}_{h}"


def _read_json(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if isinstance(data, dict):
            return data
    except Exception:
        return None
    return None


def _read_csv_first_row(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                return dict(row)
    except Exception:
        return None
    return None


def _find_metric_candidates(structure_path: str) -> List[str]:
    base = os.path.dirname(structure_path)
    stem = os.path.splitext(os.path.basename(structure_path))[0]
    candidates = []
    for ext in ("json", "csv"):
        sidecar = os.path.join(base, f"{stem}.{ext}")
        if os.path.isfile(sidecar):
            candidates.append(sidecar)
    for ext in ("json", "csv"):
        candidates.extend(glob(os.path.join(base, f"*metrics*.{ext}")))
        candidates.extend(glob(os.path.join(base, f"*score*.{ext}")))
    deduped = []
    seen = set()
    for c in candidates:
        if c not in seen:
            seen.add(c)
            deduped.append(c)
    return deduped


def _load_metrics_for_structure(structure_path: str) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    for metric_path in _find_metric_candidates(structure_path):
        if metric_path.endswith(".json"):
            data = _read_json(metric_path)
        else:
            data = _read_csv_first_row(metric_path)
        if not data:
            continue
        merged.update(data)
    return merged


def _extract_chain_lengths_pdb(path: str) -> Dict[str, int]:
    residues = defaultdict(set)
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("ATOM"):
                    continue
                chain = line[21].strip() or "_"
                resseq = line[22:26].strip()
                icode = line[26].strip()
                key = (resseq, icode)
                residues[chain].add(key)
    except Exception:
        return {}
    return {k: len(v) for k, v in residues.items()}


def _extract_chain_lengths_cif(path: str) -> Dict[str, int]:
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except Exception:
        return {}

    try:
        d = MMCIF2Dict(path)
    except Exception:
        return {}

    def to_list(v):
        return v if isinstance(v, list) else [v]

    chain = to_list(d.get("_atom_site.auth_asym_id", []))
    seq = to_list(d.get("_atom_site.auth_seq_id", []))
    grp = to_list(d.get("_atom_site.group_PDB", []))

    n = min(len(chain), len(seq), len(grp))
    residues = defaultdict(set)
    for i in range(n):
        if str(grp[i]).upper() != "ATOM":
            continue
        c = str(chain[i]).strip()
        s = str(seq[i]).strip()
        if c and s and s not in {".", "?"}:
            residues[c].add(s)
    return {k: len(v) for k, v in residues.items()}


def _infer_binder_length(structure_path: str, metrics: Dict[str, Any]) -> Optional[int]:
    direct = _safe_float(_get_first_metric(metrics, METRIC_ALIASES["binder_length"]))
    if direct is not None:
        return int(direct)

    ext = os.path.splitext(structure_path)[1].lower()
    if ext == ".pdb":
        lens = _extract_chain_lengths_pdb(structure_path)
    elif ext == ".cif":
        lens = _extract_chain_lengths_cif(structure_path)
    else:
        lens = {}

    if not lens:
        return None

    sorted_lens = sorted(lens.values())
    return sorted_lens[0] if len(sorted_lens) > 1 else sorted_lens[0]


def _collect_structures(tool_dir: str) -> List[str]:
    paths = []
    for ext in ("*.pdb", "*.cif", "*.cif.gz"):
        paths.extend(glob(os.path.join(tool_dir, "**", ext), recursive=True))
    return sorted([p for p in paths if os.path.isfile(p)])


def _normalized(values: List[Optional[float]], inverse: bool = False) -> List[Optional[float]]:
    present = [v for v in values if v is not None]
    if not present:
        return [None for _ in values]
    lo = min(present)
    hi = max(present)
    if abs(hi - lo) < 1e-12:
        base = [1.0 if v is not None else None for v in values]
    else:
        base = [None if v is None else (v - lo) / (hi - lo) for v in values]
    if inverse:
        return [None if b is None else 1.0 - b for b in base]
    return base


def _binder_within_limit(length: Any) -> Optional[bool]:
    flen = _safe_float(length)
    if flen is None:
        return None
    return flen <= float(MAX_BINDER_LENGTH_AA)


def _first_present(*values: Any) -> Any:
    for value in values:
        if value is not None:
            return value
    return None


def _merge_evaluation_metrics(record: Dict[str, Any], evaluation: Dict[str, Any]) -> None:
    eval_metrics = evaluation.get("metrics") or {}
    if not isinstance(eval_metrics, dict):
        return

    merged_metrics = {
        "ipTM": _safe_float(
            _first_present(
                _get_first_metric(eval_metrics, METRIC_ALIASES["target_ipTM"]),
                _get_first_metric(eval_metrics, METRIC_ALIASES["iptm"]),
            )
        ),
        "pTM": _safe_float(
            _first_present(
                _get_first_metric(eval_metrics, METRIC_ALIASES["target_pTM"]),
                _get_first_metric(eval_metrics, METRIC_ALIASES["ptm"]),
            )
        ),
        "min_interaction_pae": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["min_interaction_pae"])),
        "off_target_ipTM": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["off_target_ipTM"])),
        "ipSAE": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["target_ipSAE"])),
        "off_target_ipSAE": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["off_target_ipSAE"])),
        "pyrosetta_total_score": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["pyrosetta_total_score"])),
        "pyrosetta_interface_dG": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["pyrosetta_interface_dG"])),
        "pyrosetta_interface_delta_sasa": _safe_float(
            _get_first_metric(eval_metrics, METRIC_ALIASES["pyrosetta_interface_delta_sasa"])
        ),
        "pyrosetta_packstat": _safe_float(_get_first_metric(eval_metrics, METRIC_ALIASES["pyrosetta_packstat"])),
    }

    for key, value in merged_metrics.items():
        if value is not None:
            record["metrics"][key] = value

    target_ipsae = _safe_float(record["metrics"].get("ipSAE"))
    off_target_ipsae = _safe_float(record["metrics"].get("off_target_ipSAE"))
    if target_ipsae is not None and off_target_ipsae is not None:
        record["metrics"]["ipSAE_specificity_gap"] = round(target_ipsae - off_target_ipsae, 6)


def _compute_scores(records: List[Dict[str, Any]]) -> None:
    iptm_vals = [_safe_float(r["metrics"].get("ipTM")) for r in records]
    ptm_vals = [_safe_float(r["metrics"].get("pTM")) for r in records]
    min_pae_vals = [_safe_float(r["metrics"].get("min_interaction_pae")) for r in records]
    rmsd_vals = [_safe_float(r["metrics"].get("rmsd")) for r in records]

    n_iptm = _normalized(iptm_vals, inverse=False)
    n_ptm = _normalized(ptm_vals, inverse=False)
    n_inv_pae = _normalized(min_pae_vals, inverse=True)
    n_inv_rmsd = _normalized(rmsd_vals, inverse=True)

    combo_counts = defaultdict(int)
    for r in records:
        blen = r["metrics"].get("binder_length")
        bin_id = None if blen is None else int(blen) // 10
        combo_counts[(r.get("strategy"), bin_id)] += 1

    diversity_raw = []
    for r in records:
        blen = r["metrics"].get("binder_length")
        bin_id = None if blen is None else int(blen) // 10
        count = combo_counts[(r.get("strategy"), bin_id)]
        diversity_raw.append(1.0 / count if count > 0 else None)
    diversity_norm = _normalized(diversity_raw, inverse=False)

    for i, r in enumerate(records):
        clash = r["metrics"].get("clash_flag")
        no_clash_bonus = None
        if clash is not None:
            no_clash_bonus = 0.0 if bool(clash) else 1.0

        length_ok = _binder_within_limit(r["metrics"].get("binder_length"))
        r["metrics"]["passes_binder_length_filter"] = length_ok

        components = {
            "normalized_ipTM": n_iptm[i],
            "normalized_pTM": n_ptm[i],
            "normalized_inverse_min_interaction_pae": n_inv_pae[i],
            "normalized_inverse_rmsd": n_inv_rmsd[i],
            "no_clash_bonus": no_clash_bonus,
            "diversity_bonus": diversity_norm[i],
        }
        r["score_components"] = components

        weights = {
            "normalized_ipTM": 0.45,
            "normalized_pTM": 0.20,
            "normalized_inverse_min_interaction_pae": 0.25,
            "normalized_inverse_rmsd": 0.05,
            "no_clash_bonus": 0.03,
            "diversity_bonus": 0.02,
        }
        required_keys = (
            "normalized_ipTM",
            "normalized_pTM",
            "normalized_inverse_min_interaction_pae",
        )

        if length_ok is False:
            r["composite_score"] = None
            r["specificity_adjusted_score"] = None
            r.setdefault("warnings", []).append(
                f"Binder length exceeds {MAX_BINDER_LENGTH_AA} aa filter; excluded from ranking."
            )
            continue

        missing_required = [key for key in required_keys if components.get(key) is None]
        if missing_required:
            r["composite_score"] = None
            r["specificity_adjusted_score"] = None
            r.setdefault("warnings", []).append(
                "Insufficient core metrics for ranking: " + ", ".join(missing_required)
            )
            continue

        weighted = 0.0
        weight_sum = 0.0
        missing = []
        for key, weight in weights.items():
            value = components.get(key)
            if value is None:
                missing.append(key)
                continue
            weighted += weight * float(value)
            weight_sum += weight

        base_score = round(weighted / weight_sum, 6) if weight_sum > 0 else None
        r["composite_score"] = base_score

        off_target_iptm = _safe_float(r["metrics"].get("off_target_ipTM"))
        target_iptm = _safe_float(r["metrics"].get("ipTM"))
        specificity_gap = None
        specificity_penalty = 0.0
        if target_iptm is not None and off_target_iptm is not None:
            specificity_gap = round(target_iptm - off_target_iptm, 6)
            if off_target_iptm >= 0.60:
                specificity_penalty = 0.15
            elif off_target_iptm >= 0.45:
                specificity_penalty = 0.08

        r["metrics"]["specificity_gap"] = specificity_gap
        if base_score is not None:
            r["specificity_adjusted_score"] = round(max(0.0, base_score - specificity_penalty), 6)
        else:
            r["specificity_adjusted_score"] = None

        if missing:
            r.setdefault("warnings", []).append(
                "Missing metrics for composite score components: " + ", ".join(missing)
            )
        if off_target_iptm is None:
            r.setdefault("warnings", []).append("Missing FGFR1 off-target metric (off_target_ipTM); specificity risk not quantified.")
        elif off_target_iptm >= 0.60:
            r.setdefault("warnings", []).append("High predicted FGFR1 off-target binding risk (off_target_ipTM >= 0.60).")


def sync_designs(campaign_dir: str) -> Dict[str, Any]:
    out_dir = os.path.join(campaign_dir, "out")
    designs_dir = os.path.join(campaign_dir, "designs")
    os.makedirs(designs_dir, exist_ok=True)

    tool_dirs = {
        "boltzgen": os.path.join(out_dir, "boltzgen"),
        "rfdiffusion3": os.path.join(out_dir, "rfdiffusion3"),
        "boltz2": os.path.join(out_dir, "boltz2"),
    }

    records: List[Dict[str, Any]] = []
    for tool, tool_dir in tool_dirs.items():
        if not os.path.isdir(tool_dir):
            continue
        for structure_path in _collect_structures(tool_dir):
            rel = os.path.relpath(structure_path, tool_dir)
            design_id = _make_design_id(tool, rel)
            strategy = _detect_strategy(structure_path)
            raw_metrics = _load_metrics_for_structure(structure_path)

            metrics = {
                "ipTM": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["iptm"])),
                "pTM": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["ptm"])),
                "pLDDT": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["plddt"])),
                "min_interaction_pae": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["min_interaction_pae"])),
                "rmsd": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["rmsd"])),
                "clash_flag": _safe_bool(_get_first_metric(raw_metrics, METRIC_ALIASES["clash_flag"])),
                "off_target_ipTM": _safe_float(_get_first_metric(raw_metrics, METRIC_ALIASES["off_target_ipTM"])),
                "binder_length": None,
                "binder_sequence": _get_first_metric(raw_metrics, METRIC_ALIASES["binder_sequence"]),
            }
            metrics["binder_length"] = _infer_binder_length(structure_path, raw_metrics)

            warnings = [
                "Computational binder candidate only. No experimental efficacy claim."
            ]
            if metrics["binder_length"] is not None and metrics["binder_length"] > MAX_BINDER_LENGTH_AA:
                warnings.append(
                    f"Binder length > {MAX_BINDER_LENGTH_AA} aa; redesign or filter before prioritization."
                )

            record = {
                "design_id": design_id,
                "tool": tool,
                "strategy": strategy,
                "source_path": os.path.abspath(structure_path),
                "metrics": metrics,
                "warnings": warnings,
                "notes": [],
            }
            records.append(record)

    eval_jsons = glob(os.path.join(out_dir, "boltz2", "evaluations", "**", "evaluation_metrics.json"), recursive=True)
    eval_by_design = {}
    for epath in eval_jsons:
        data = _read_json(epath)
        if not data:
            continue
        did = data.get("design_id")
        if did:
            eval_by_design[did] = data

    for record in records:
        did = record["design_id"]
        if did in eval_by_design:
            record["evaluation"] = eval_by_design[did]
            _merge_evaluation_metrics(record, eval_by_design[did])

    _compute_scores(records)

    for record in records:
        did = record["design_id"]
        ddir = os.path.join(designs_dir, did)
        os.makedirs(ddir, exist_ok=True)

        src = record["source_path"]
        ext = os.path.splitext(src)[1].lower()
        dst_name = "designed.pdb" if ext == ".pdb" else "designed.cif"
        dst = os.path.join(ddir, dst_name)
        shutil.copy2(src, dst)

        metrics_out = {
            "design_id": did,
            "tool": record["tool"],
            "strategy": record["strategy"],
            "source_path": record["source_path"],
            "metrics": record["metrics"],
            "score_components": record.get("score_components"),
            "composite_score": record.get("composite_score"),
            "specificity_adjusted_score": record.get("specificity_adjusted_score"),
            "warnings": record.get("warnings", []),
            "evaluation": record.get("evaluation"),
            "updated_at": now_utc(),
        }
        with open(os.path.join(ddir, "metrics.json"), "w", encoding="utf-8") as handle:
            json.dump(metrics_out, handle, indent=2)

        record["design_path"] = os.path.abspath(dst)

    def sort_key(item: Dict[str, Any]) -> Tuple[int, float, float, float]:
        metrics = item.get("metrics") or {}

        ipsae = _safe_float(metrics.get("ipSAE"))
        if ipsae is not None:
            ipsae_gap = _safe_float(metrics.get("ipSAE_specificity_gap"))
            if ipsae_gap is None:
                ipsae_gap = -1e9

            score = item.get("specificity_adjusted_score")
            if score is None:
                score = item.get("composite_score")
            score_value = -float(score) if score is not None else 0.0
            return (0, -float(ipsae), -float(ipsae_gap), score_value)

        score = item.get("specificity_adjusted_score")
        if score is None:
            score = item.get("composite_score")
        if score is None:
            return (1, 0.0, 0.0, 0.0)
        return (1, 0.0, 0.0, -float(score))

    records = sorted(records, key=sort_key)

    index = {
        "campaign": "fgfr2",
        "generated_at": now_utc(),
        "num_designs": len(records),
        "designs": records,
        "global_warnings": [
            "All outputs are computational binder candidates for future experimental validation.",
            f"Binder ranking limit set to <= {MAX_BINDER_LENGTH_AA} aa for future prioritization.",
            "Dashboard ranking currently prioritizes higher target ipSAE (then ipSAE specificity gap) when available.",
            "Composite scores require core metrics (ipTM, pTM, min_interaction_pae); incomplete designs are not ranked.",
            "Specificity against FGFR1 should be evaluated before prioritizing candidates.",
        ],
    }

    with open(os.path.join(designs_dir, "index.json"), "w", encoding="utf-8") as handle:
        json.dump(index, handle, indent=2)

    return index


def main() -> int:
    parser = argparse.ArgumentParser(description="Sync local FGFR2 design outputs into standardized design index")
    parser.add_argument(
        "--campaign-dir",
        default=os.path.abspath(os.path.join(os.path.dirname(__file__), ".")),
        help="Path to fgfr2_campaign",
    )
    args = parser.parse_args()

    index = sync_designs(args.campaign_dir)
    print(f"Wrote index with {index['num_designs']} designs: {os.path.join(args.campaign_dir, 'designs', 'index.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

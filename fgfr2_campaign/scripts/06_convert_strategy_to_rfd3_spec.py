#!/usr/bin/env python3
"""Convert FGFR strategy YAML to RFdiffusion3 JSON design spec.

This creates a Foundry-compatible JSON config used by:
  rfd3 design out_dir=... inputs=<json>
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Dict, List, Tuple

MAX_BINDER_LENGTH_AA = 50
DEFAULT_MIN_BINDER_LENGTH_AA = 30

import yaml
from Bio.PDB import MMCIFParser, PDBParser


def _load_yaml(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Strategy config must be a mapping")
    return data


def _abs(path: str) -> str:
    return os.path.abspath(path)


def _parse_hotspots(values: List[str], receptor_chain: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in values or []:
        if not isinstance(item, str) or ":" not in item:
            continue
        chain, seq = item.split(":", 1)
        chain = chain.strip()
        seq = seq.strip()
        if not chain or not seq.isdigit():
            continue
        if chain != receptor_chain:
            continue
        out[f"{chain}{int(seq)}"] = "CA"
    return out


def _build_contig_string(structure_path: str, receptor_chain: str, min_len: int, max_len: int) -> str:
    """Build RFdiffusion contig string that skips missing residues.
    
    Example: if chain A has residues 147-161, 163-185, 187-362 (skipping 162, 186),
    returns: "60-120,/0,A147-161,A163-185,A187-362"
    """
    parser = MMCIFParser(QUIET=True) if structure_path.lower().endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("target", structure_path)
    model = next(structure.get_models())
    if receptor_chain not in model:
        raise ValueError(f"Receptor chain '{receptor_chain}' not found in {structure_path}")

    # Collect all standard residue numbers
    residues = set()
    for residue in model[receptor_chain].get_residues():
        hetfield, resseq, _icode = residue.id
        if hetfield.strip() == "":
            residues.add(int(resseq))

    if not residues:
        raise ValueError(f"No standard residues found for chain '{receptor_chain}' in {structure_path}")

    # Build continuous ranges from sorted residue list
    sorted_res = sorted(residues)
    ranges = []
    start = sorted_res[0]
    end = sorted_res[0]

    for res in sorted_res[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append((start, end))
            start = res
            end = res
    ranges.append((start, end))

    # Format: "binder_len,/0,chain_res1-res2,chain_res3-res4,..."
    receptor_ranges = ",".join(f"{receptor_chain}{s}-{e}" for s, e in ranges)
    return f"{min_len}-{max_len},/0,{receptor_ranges}"


def build_rfd3_spec(strategy_cfg: Dict) -> Tuple[str, Dict]:
    strategy_name = str(strategy_cfg.get("strategy_name") or "fgfr2_design")

    target = strategy_cfg.get("target") or {}
    binder = strategy_cfg.get("binder") or {}

    structure_path = str(target.get("structure_path") or "").strip()
    receptor_chain = str(target.get("receptor_chain") or "A").strip()
    min_len = int(binder.get("min_length") or DEFAULT_MIN_BINDER_LENGTH_AA)
    max_len = int(binder.get("max_length") or MAX_BINDER_LENGTH_AA)
    max_len = min(max_len, MAX_BINDER_LENGTH_AA)
    min_len = min(min_len, max_len)
    hotspots = _parse_hotspots(strategy_cfg.get("hotspots") or [], receptor_chain)

    if not structure_path:
        raise ValueError("Missing target.structure_path in strategy config")

    structure_path = _abs(structure_path)
    contig = _build_contig_string(structure_path, receptor_chain, min_len, max_len)

    spec: Dict = {
        "dialect": 2,
        "input": structure_path,
        "contig": contig,
        "is_non_loopy": True,
    }

    if hotspots:
        spec["infer_ori_strategy"] = "hotspots"
        spec["select_hotspots"] = hotspots

    return strategy_name, {strategy_name: spec}


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert FGFR strategy YAML to RFdiffusion3 JSON")
    parser.add_argument("--strategy-config", required=True, help="Path to strategy YAML")
    parser.add_argument("--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    cfg = _load_yaml(args.strategy_config)
    strategy_name, payload = build_rfd3_spec(cfg)

    out_path = _abs(args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    print(f"Wrote RFdiffusion3 JSON spec: {out_path}")
    print(f"Strategy: {strategy_name}")
    print(f"Contig: {list(payload.values())[0]['contig']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

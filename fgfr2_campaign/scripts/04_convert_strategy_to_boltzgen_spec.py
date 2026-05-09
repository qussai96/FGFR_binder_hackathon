#!/usr/bin/env python3
"""Convert FGFR strategy YAML into a native BoltzGen design-spec YAML."""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Dict, List, Tuple

MAX_BINDER_LENGTH_AA = 50
DEFAULT_MIN_BINDER_LENGTH_AA = 30

import yaml
from Bio.PDB import PDBParser


def _load_yaml(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return data


def _is_aa_residue(residue) -> bool:
    hetfield = residue.id[0]
    return hetfield == " "


def _map_auth_seq_to_chain_pos(structure_path: str, chain_id: str) -> Dict[int, int]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", structure_path)
    model = next(structure.get_models())
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in {structure_path}")

    chain = model[chain_id]
    mapping: Dict[int, int] = {}
    pos = 0
    for residue in chain:
        if not _is_aa_residue(residue):
            continue
        pos += 1
        auth_seq = int(residue.id[1])
        # Keep first occurrence for insertion-code duplicates.
        if auth_seq not in mapping:
            mapping[auth_seq] = pos
    return mapping


def _extract_hotspots(strategy_cfg: Dict, receptor_chain: str) -> List[int]:
    values = []
    for item in strategy_cfg.get("hotspots", []) or []:
        if isinstance(item, str):
            token = item
        elif isinstance(item, dict) and len(item) == 1:
            key = next(iter(item.keys()))
            token = f"{key}:{item[key]}"
        else:
            continue

        if ":" not in token:
            continue
        chain, seq = token.split(":", 1)
        if chain.strip() != receptor_chain:
            continue
        try:
            values.append(int(seq))
        except ValueError:
            continue
    return sorted(set(values))


def _to_range_string(positions: List[int]) -> str:
    return ",".join(str(p) for p in sorted(set(positions)))


def convert(strategy_path: str, output_path: str) -> Tuple[str, int, int]:
    cfg = _load_yaml(strategy_path)

    target = cfg.get("target", {})
    receptor_chain = str(target.get("receptor_chain", "A"))
    ligand_chain = target.get("ligand_reference_chain")

    binder = cfg.get("binder", {})
    min_len = int(binder.get("min_length", DEFAULT_MIN_BINDER_LENGTH_AA))
    max_len = int(binder.get("max_length", MAX_BINDER_LENGTH_AA))
    max_len = min(max_len, MAX_BINDER_LENGTH_AA)
    min_len = min(min_len, max_len)
    if min_len > max_len:
        min_len, max_len = max_len, min_len

    structure_path = str(target.get("structure_path", "")).strip()
    if not structure_path:
        raise ValueError("Missing target.structure_path in strategy config")

    # Resolve relative paths from workspace root (current working directory expected at repo root).
    abs_structure_path = os.path.abspath(structure_path)
    if not os.path.isfile(abs_structure_path):
        raise FileNotFoundError(f"Structure file not found: {abs_structure_path}")

    auth_hotspots = _extract_hotspots(cfg, receptor_chain)
    seq_map = _map_auth_seq_to_chain_pos(abs_structure_path, receptor_chain)
    pos_hotspots = [seq_map[h] for h in auth_hotspots if h in seq_map]

    include_chains = [{"chain": {"id": receptor_chain}}]
    if ligand_chain:
        include_chains.append({"chain": {"id": str(ligand_chain)}})

    entities = [
        {
            "file": {
                "path": abs_structure_path,
                "include": include_chains,
                "not_design": [{"chain": {"id": receptor_chain}}],
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

    if pos_hotspots:
        entities[0]["file"]["binding_types"] = [
            {
                "chain": {
                    "id": receptor_chain,
                    "binding": _to_range_string(pos_hotspots),
                }
            }
        ]

    native = {"entities": entities}
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(native, handle, sort_keys=False)

    return output_path, len(auth_hotspots), len(pos_hotspots)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert strategy YAML to native BoltzGen design-spec")
    parser.add_argument("--strategy-config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    out, n_auth, n_mapped = convert(args.strategy_config, args.output)
    print(f"Wrote BoltzGen design-spec: {out}")
    print(f"Hotspots (auth numbering): {n_auth} | mapped to chain positions: {n_mapped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

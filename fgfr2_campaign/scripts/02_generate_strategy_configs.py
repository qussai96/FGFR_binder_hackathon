#!/usr/bin/env python3
"""Generate FGFR2 strategy config files from interface hotspot output."""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from typing import Any, Dict, List


def _load_hotspots(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _hotspot_strings(hotspots_json: Dict[str, Any]) -> List[str]:
    results: List[str] = []
    for item in hotspots_json.get("interface_residues", []):
        chain = item.get("auth_asym_id") or item.get("label_asym_id") or item.get("receptor_chain_parser")
        seq = item.get("auth_seq_id")
        if seq is None:
            seq = item.get("label_seq_id")
        if chain is None or seq is None:
            continue
        results.append(f"{chain}:{seq}")
    return sorted(set(results))


def _yaml_list(values: List[str], indent: int = 2) -> str:
    if not values:
        return "[]"
    pad = " " * indent
    return "\n".join(f"{pad}- {v}" for v in values)


def _load_specificity(path: str) -> Dict[str, Any]:
  if not path or not os.path.isfile(path):
    return {}
  with open(path, "r", encoding="utf-8") as handle:
    return json.load(handle)


def _preferred_specific_hotspots(spec_json: Dict[str, Any], chain_id: Optional[str]) -> List[str]:
  if not spec_json:
    return []
  out = []
  for item in spec_json.get("preferred_hotspots", []):
    resseq = item.get("target_resseq")
    if resseq is None or chain_id is None:
      continue
    out.append(f"{chain_id}:{resseq}")
  return sorted(set(out))


def _conserved_hotspots(spec_json: Dict[str, Any], chain_id: Optional[str]) -> List[str]:
  if not spec_json:
    return []
  out = []
  for item in spec_json.get("conserved_hotspots", []):
    resseq = item.get("target_resseq")
    if resseq is None or chain_id is None:
      continue
    out.append(f"{chain_id}:{resseq}")
  return sorted(set(out))


def generate_configs(campaign_dir: str, hotspots_path: str, specificity_path: Optional[str] = None) -> None:
    cfg_dir = os.path.join(campaign_dir, "configs")
    os.makedirs(cfg_dir, exist_ok=True)

    hs = _load_hotspots(hotspots_path)
    spec = _load_specificity(specificity_path) if specificity_path else {}
    hotspot_residues = _hotspot_strings(hs)
    receptor_chain = hs.get("receptor_chain_input")
    preferred_specific = _preferred_specific_hotspots(spec, receptor_chain)
    conserved_specific = _conserved_hotspots(spec, receptor_chain)
    generated_at = datetime.now(timezone.utc).isoformat()

    warning_lines = []
    if hs.get("warning"):
        warning_lines.append(hs["warning"])
    if not hotspot_residues:
        warning_lines.append("No interface hotspots available; strategy1 confidence reduced until chain IDs are corrected.")
    if spec.get("warning"):
      warning_lines.append(spec["warning"])

    strategy1 = f"""campaign: fgfr2
strategy_name: strategy1_ligand_site_blocker
purpose: Design computational binder candidates that target the FGFR2 surface used by FGF1.
target:
  pdb_id: 1DJS
  structure_path: fgfr2_campaign/structures/fgfr2_target_clean.pdb
  receptor_chain: {hs.get('receptor_chain_input')}
  ligand_reference_chain: {hs.get('ligand_chain_input')}
off_target:
  pdb_id: 1CVS
  structure_path: fgfr2_campaign/structures/1CVS_offtarget_fgfr1.pdb
  receptor_chain: {spec.get('offtarget_chain')}
binder:
  min_length: 30
  max_length: 50
hotspot_source_json: fgfr2_campaign/configs/fgfr2_fgf1_interface_hotspots.json
hotspots:
{_yaml_list(hotspot_residues, indent=2)}
specificity_preferred_hotspots:
{_yaml_list(preferred_specific, indent=2) if preferred_specific else '[]'}
specificity_conserved_hotspots_to_avoid:
{_yaml_list(conserved_specific, indent=2) if conserved_specific else '[]'}
scoring_focus:
  - interface_confidence
  - min_interaction_pae
  - ipTM
  - fgfr1_off_target_avoidance
warnings:
{_yaml_list(warning_lines, indent=2) if warning_lines else '  - None'}
labels:
  - computational binder candidates
  - local A40 workflow
notes:
  generated_at_utc: {generated_at}
"""

    strategy2 = f"""campaign: fgfr2
strategy_name: strategy2_dimer_or_activation_blocker
purpose: Explore a plausible receptor activation/dimerization-relevant surface; fallback is secondary surface near ligand-binding domain.
target:
  pdb_id: 1DJS
  structure_path: fgfr2_campaign/structures/fgfr2_target_clean.pdb
  receptor_chain: {hs.get('receptor_chain_input')}
off_target:
  pdb_id: 1CVS
  structure_path: fgfr2_campaign/structures/1CVS_offtarget_fgfr1.pdb
  receptor_chain: {spec.get('offtarget_chain')}
binder:
  min_length: 30
  max_length: 50
mode:
  confidence: lower
  selection: secondary surface near ligand-binding domain
hotspot_source_json: fgfr2_campaign/configs/fgfr2_fgf1_interface_hotspots.json
hotspots:
{_yaml_list(hotspot_residues[:20], indent=2) if hotspot_residues else '[]'}
specificity_conserved_hotspots_to_avoid:
{_yaml_list(conserved_specific[:20], indent=2) if conserved_specific else '[]'}
warnings:
  - Secondary-surface strategy may not map true activation interface from 1DJS alone.
  - Computational binder candidates only; not validated inhibitors.
notes:
  generated_at_utc: {generated_at}
"""

    strategy3 = f"""campaign: fgfr2
strategy_name: strategy3_surface_explorer
purpose: Model-free surface binder exploration on FGFR2 extracellular region.
target:
  pdb_id: 1DJS
  structure_path: fgfr2_campaign/structures/fgfr2_target_clean.pdb
  receptor_chain: {hs.get('receptor_chain_input')}
off_target:
  pdb_id: 1CVS
  structure_path: fgfr2_campaign/structures/1CVS_offtarget_fgfr1.pdb
  receptor_chain: {spec.get('offtarget_chain')}
binder:
  min_length: 30
  max_length: 50
mode:
  hotspot_constrained: false
  broad_surface_sampling: true
specificity:
  prefer_nonconserved_residues: true
  avoid_conserved_fgfr1_positions: true
warnings:
  - No explicit hotspot constraints by default.
  - Computational binder candidates only; not validated inhibitors.
notes:
  generated_at_utc: {generated_at}
"""

    with open(os.path.join(cfg_dir, "strategy1_ligand_site_blocker.yaml"), "w", encoding="utf-8") as h:
        h.write(strategy1)
    with open(os.path.join(cfg_dir, "strategy2_dimer_or_activation_blocker.yaml"), "w", encoding="utf-8") as h:
        h.write(strategy2)
    with open(os.path.join(cfg_dir, "strategy3_surface_explorer.yaml"), "w", encoding="utf-8") as h:
        h.write(strategy3)

    print("Generated strategy config files:")
    print(f"- {os.path.join(cfg_dir, 'strategy1_ligand_site_blocker.yaml')}")
    print(f"- {os.path.join(cfg_dir, 'strategy2_dimer_or_activation_blocker.yaml')}")
    print(f"- {os.path.join(cfg_dir, 'strategy3_surface_explorer.yaml')}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate FGFR2 strategy YAML configs from hotspot JSON")
    parser.add_argument(
        "--campaign-dir",
        default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        help="Path to fgfr2_campaign directory",
    )
    parser.add_argument(
        "--hotspots-json",
        default=None,
        help="Path to hotspot JSON. Defaults to configs/fgfr2_fgf1_interface_hotspots.json",
    )
    parser.add_argument(
      "--specificity-json",
      default=None,
      help="Optional specificity JSON from 03_build_specificity_filter.py",
    )
    args = parser.parse_args()

    hotspots_path = args.hotspots_json or os.path.join(args.campaign_dir, "configs", "fgfr2_fgf1_interface_hotspots.json")
    generate_configs(campaign_dir=args.campaign_dir, hotspots_path=hotspots_path, specificity_path=args.specificity_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Inspect a structure file (PDB or mmCIF) for FGFR2 campaign setup."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Set


def _safe_text(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", ".", "?"}:
        return None
    return text


def _safe_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", ".", "?"}:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _to_list(value):
    if isinstance(value, list):
        return value
    return [value]


def _parse_pdb_compnd_descriptions(path: str) -> List[Dict[str, object]]:
    compnd_lines = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("COMPND"):
                    compnd_lines.append(line[10:].strip())
                elif compnd_lines and not line.startswith("COMPND"):
                    break
    except Exception:
        return []

    if not compnd_lines:
        return []

    text = " ".join(compnd_lines)
    blocks = [b.strip() for b in text.split("MOL_ID:") if b.strip()]

    entities = []
    for block in blocks:
        block = block.replace("  ", " ")
        parts = [p.strip() for p in block.split(";") if p.strip()]
        entity = {"entity_id": None, "description": None, "chains": []}
        for part in parts:
            if part.isdigit() and entity["entity_id"] is None:
                entity["entity_id"] = part
            elif part.startswith("MOLECULE:"):
                entity["description"] = part.split(":", 1)[1].strip()
            elif part.startswith("CHAIN:"):
                chain_text = part.split(":", 1)[1].strip()
                chains = [c.strip() for c in chain_text.split(",") if c.strip()]
                entity["chains"] = chains
        entities.append(entity)
    return entities


def inspect_structure(path: str) -> Dict[str, object]:
    try:
        from Bio.PDB import MMCIFParser, PDBParser
    except Exception as exc:
        print("ERROR: BioPython is required. Install with: pip install biopython", file=sys.stderr)
        raise SystemExit(2) from exc

    abs_path = os.path.abspath(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".cif":
        parser = MMCIFParser(QUIET=True)
        file_format = "mmCIF"
    else:
        parser = PDBParser(QUIET=True)
        file_format = "PDB"

    structure = parser.get_structure("target", abs_path)
    model = next(structure.get_models())

    chain_residues_auth: Dict[str, Set[int]] = defaultdict(set)
    non_protein_ligands = set()

    for chain in model:
        for residue in chain.get_residues():
            hetfield, resseq, _icode = residue.id
            if hetfield.strip() == "":
                chain_residues_auth[chain.id].add(int(resseq))
            elif hetfield.strip() not in {"W"}:
                resname = residue.get_resname().strip()
                if resname not in {"HOH", "WAT", "DOD"}:
                    non_protein_ligands.add(resname)

    residue_ranges = {
        c: {"min": min(vals), "max": max(vals)}
        for c, vals in chain_residues_auth.items()
        if vals
    }

    entities = []
    if file_format == "mmCIF":
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict

        mmcif = MMCIF2Dict(abs_path)
        entity_ids = _to_list(mmcif.get("_entity.id", []))
        entity_desc = _to_list(mmcif.get("_entity.pdbx_description", []))
        asym_ids = _to_list(mmcif.get("_struct_asym.id", []))
        asym_entity_ids = _to_list(mmcif.get("_struct_asym.entity_id", []))

        entity_map = {}
        for i, eid in enumerate(entity_ids):
            entity_map[str(eid)] = {
                "entity_id": str(eid),
                "description": _safe_text(entity_desc[i]) if i < len(entity_desc) else None,
                "chains": [],
            }
        for i, asym in enumerate(asym_ids):
            if i < len(asym_entity_ids):
                eid = str(asym_entity_ids[i])
                if eid in entity_map:
                    entity_map[eid]["chains"].append(str(asym))
        entities = list(entity_map.values())
    else:
        entities = _parse_pdb_compnd_descriptions(abs_path)

    candidate_fgfr2 = []
    candidate_fgf1 = []
    candidate_fgfr1 = []
    for ent in entities:
        desc = (ent.get("description") or "").lower()
        chains = ent.get("chains") or []
        if "fgfr2" in desc or "fibroblast growth factor receptor 2" in desc:
            candidate_fgfr2.extend(chains)
        if "fgf1" in desc or "fibroblast growth factor 1" in desc:
            candidate_fgf1.extend(chains)
        if "fgfr1" in desc or "fibroblast growth factor receptor 1" in desc:
            candidate_fgfr1.extend(chains)

    return {
        "file": abs_path,
        "file_format": file_format,
        "chain_ids": sorted(chain_residues_auth.keys()),
        "auth_residue_ranges": residue_ranges,
        "entities": entities,
        "non_protein_ligands": sorted(non_protein_ligands),
        "candidate_fgfr2_receptor_chains": sorted(set(candidate_fgfr2)),
        "candidate_fgf1_ligand_chains": sorted(set(candidate_fgf1)),
        "candidate_fgfr1_receptor_chains": sorted(set(candidate_fgfr1)),
    }


def _print_summary(summary: Dict[str, object]) -> None:
    print(f"Structure file: {summary['file']}")
    print(f"File format: {summary['file_format']}")

    print("\nChain IDs")
    print("- " + (", ".join(summary.get("chain_ids", [])) or "None"))

    print("\nResidue ranges by chain (auth numbering)")
    for chain, bounds in (summary.get("auth_residue_ranges") or {}).items():
        print(f"- {chain}: {bounds['min']}..{bounds['max']}")

    print("\nEntity descriptions")
    for ent in summary.get("entities", []):
        print(
            "- entity {eid}: chains={chains}, description={desc}".format(
                eid=ent.get("entity_id"),
                chains=ent.get("chains"),
                desc=ent.get("description"),
            )
        )

    print("\nNon-protein ligands")
    ligs = summary.get("non_protein_ligands", [])
    if ligs:
        for lig in ligs:
            print(f"- {lig}")
    else:
        print("- None detected")

    print("\nCandidate FGFR2 receptor chains")
    fgfr2 = summary.get("candidate_fgfr2_receptor_chains", [])
    print("- " + (", ".join(fgfr2) if fgfr2 else "None detected from entity descriptions"))

    print("\nCandidate FGF1 ligand chains")
    fgf1 = summary.get("candidate_fgf1_ligand_chains", [])
    print("- " + (", ".join(fgf1) if fgf1 else "None detected from entity descriptions"))

    print("\nCandidate FGFR1 receptor chains (off-target) ")
    fgfr1 = summary.get("candidate_fgfr1_receptor_chains", [])
    print("- " + (", ".join(fgfr1) if fgfr1 else "None detected from entity descriptions"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect PDB/mmCIF content for FGFR2 campaign")
    parser.add_argument("structure", help="Path to PDB or mmCIF")
    parser.add_argument("--json", dest="json_out", help="Optional JSON output path")
    args = parser.parse_args()

    summary = inspect_structure(args.structure)
    _print_summary(summary)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)
        print(f"\nWrote summary JSON: {args.json_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Find receptor residues at a protein-protein interface from PDB/mmCIF."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Dict, Optional


def _distance(a, b) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _is_heavy_atom(atom) -> bool:
    element = (getattr(atom, "element", "") or "").strip().upper()
    name = atom.get_name().strip().upper()
    if element == "H":
        return False
    if not element and name.startswith("H"):
        return False
    return True


def _choose_parser(structure_path: str):
    from Bio.PDB import MMCIFParser, PDBParser

    if structure_path.lower().endswith(".cif"):
        return MMCIFParser(QUIET=True), "mmCIF"
    return PDBParser(QUIET=True), "PDB"


def find_interface_residues(structure_path: str, receptor_chain: str, ligand_chain: str, cutoff: float):
    try:
        parser, _format = _choose_parser(structure_path)
    except Exception as exc:
        print("ERROR: BioPython is required. Install with: pip install biopython", file=sys.stderr)
        raise SystemExit(2) from exc

    from Bio.PDB.Polypeptide import is_aa

    structure = parser.get_structure("target", structure_path)
    model = next(structure.get_models())

    if receptor_chain not in model:
        raise ValueError(f"Receptor chain '{receptor_chain}' not found in parsed structure.")
    if ligand_chain not in model:
        raise ValueError(f"Ligand chain '{ligand_chain}' not found in parsed structure.")

    receptor = model[receptor_chain]
    ligand = model[ligand_chain]

    ligand_atoms = [a for a in ligand.get_atoms() if _is_heavy_atom(a)]
    if not ligand_atoms:
        raise ValueError("No ligand heavy atoms found. Check ligand chain ID.")

    contacting = []
    for residue in receptor.get_residues():
        hetfield, resseq, icode = residue.id
        # Keep only standard amino-acid residues in the receptor interface map.
        if hetfield.strip() != "":
            continue
        if not is_aa(residue, standard=True):
            continue

        receptor_atoms = [a for a in residue.get_atoms() if _is_heavy_atom(a)]
        if not receptor_atoms:
            continue

        min_dist = float("inf")
        in_contact = False
        for ra in receptor_atoms:
            for la in ligand_atoms:
                d = _distance(ra.coord, la.coord)
                if d < min_dist:
                    min_dist = d
                if d <= cutoff:
                    in_contact = True
                    break
            if in_contact:
                break

        if in_contact:
            contacting.append(
                {
                    "receptor_chain_parser": receptor_chain,
                    "resname": residue.get_resname().strip(),
                    "insertion_code": (icode or " ").strip() or None,
                    "auth_asym_id": receptor_chain,
                    "auth_seq_id": int(resseq),
                    "label_asym_id": None,
                    "label_seq_id": None,
                    "min_heavy_atom_distance": round(float(min_dist), 3),
                }
            )

    contacting.sort(key=lambda x: (x.get("auth_seq_id") is None, x.get("auth_seq_id") or 10**9, x["resname"]))
    return contacting


def main() -> int:
    parser = argparse.ArgumentParser(description="Find FGFR2 receptor residues at FGFR2-FGF1 interface")
    parser.add_argument("--structure", default=None, help="Input structure path (PDB or mmCIF)")
    parser.add_argument("--cif", default=None, help="Backward-compatible alias for --structure")
    parser.add_argument("--receptor-chain", required=True, help="Receptor chain ID")
    parser.add_argument("--ligand-chain", required=True, help="Ligand chain ID")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Heavy-atom distance cutoff in Angstrom")
    parser.add_argument(
        "--output",
        default=None,
        help="Output JSON path (default: fgfr2_campaign/configs/fgfr2_fgf1_interface_hotspots.json)",
    )
    args = parser.parse_args()

    structure_path = os.path.abspath(args.structure or args.cif or "")
    if not structure_path:
        raise ValueError("Provide --structure (or --cif for compatibility).")

    out_path = args.output
    if out_path is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        campaign_dir = os.path.abspath(os.path.join(script_dir, ".."))
        out_path = os.path.join(campaign_dir, "configs", "fgfr2_fgf1_interface_hotspots.json")

    contacting = find_interface_residues(
        structure_path=structure_path,
        receptor_chain=args.receptor_chain,
        ligand_chain=args.ligand_chain,
        cutoff=args.cutoff,
    )

    payload = {
        "target": "FGFR2 extracellular ligand-binding domain",
        "source_structure": os.path.basename(structure_path),
        "source_path": structure_path,
        "receptor_chain_input": args.receptor_chain,
        "ligand_chain_input": args.ligand_chain,
        "cutoff_angstrom": args.cutoff,
        "num_interface_residues": len(contacting),
        "warning": None,
        "interface_residues": contacting,
    }

    if not contacting:
        payload["warning"] = (
            "No interface residues found with the given chain IDs and cutoff. "
            "Verify chain IDs from inspect_cif.py output."
        )

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    print(f"Wrote hotspot JSON: {out_path}")
    print(f"Receptor chain: {args.receptor_chain} | Ligand chain: {args.ligand_chain} | Cutoff: {args.cutoff:.2f} A")
    if contacting:
        print("Receptor residues contacting ligand (FGF1 candidate):")
        for r in contacting:
            print(
                "- chain={chain} res={res} auth_seq={auth} label_seq={label} icode={icode} min_d={dist}A".format(
                    chain=r.get("auth_asym_id") or r.get("receptor_chain_parser"),
                    res=r.get("resname"),
                    auth=r.get("auth_seq_id"),
                    label=r.get("label_seq_id"),
                    icode=r.get("insertion_code") or "-",
                    dist=r.get("min_heavy_atom_distance"),
                )
            )
    else:
        print("WARNING: no interface residues found.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

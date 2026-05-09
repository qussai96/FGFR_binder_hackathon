#!/usr/bin/env python3
"""Build FGFR2-vs-FGFR1 specificity filter from interface hotspots.

Produces a JSON with preferred (less FGFR1-conserved) and conserved hotspot sets.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Dict, List, Optional, Tuple

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _load_json(path: str) -> Dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _extract_chain_residues(structure_path: str, chain_id: str) -> List[Tuple[int, str]]:
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if structure_path.lower().endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("x", structure_path)
    model = next(structure.get_models())
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in {structure_path}")

    out = []
    for residue in model[chain_id].get_residues():
        hetfield, resseq, _icode = residue.id
        if hetfield.strip() != "":
            continue
        aa = AA3_TO_1.get(residue.get_resname().strip().upper(), "X")
        out.append((int(resseq), aa))
    return out


def _align_map(target_res: List[Tuple[int, str]], off_res: List[Tuple[int, str]]) -> Dict[int, Optional[Tuple[int, str]]]:
    from Bio import Align

    target_seq = "".join(a for _i, a in target_res)
    off_seq = "".join(a for _i, a in off_res)

    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5

    aln = aligner.align(target_seq, off_seq)[0]
    aligned_target, aligned_off = str(aln[0]), str(aln[1])

    mapping: Dict[int, Optional[Tuple[int, str]]] = {}
    ti = 0
    oi = 0
    for t_char, o_char in zip(aligned_target, aligned_off):
        if t_char != "-":
            t_resseq = target_res[ti][0]
            if o_char != "-":
                mapping[t_resseq] = (off_res[oi][0], off_res[oi][1])
                oi += 1
            else:
                mapping[t_resseq] = None
            ti += 1
        elif o_char != "-":
            oi += 1

    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(description="Build FGFR1 off-target exclusion hotspots from FGFR2 interface")
    parser.add_argument("--hotspots-json", required=True, help="FGFR2 interface hotspot JSON")
    parser.add_argument("--target-structure", required=True, help="FGFR2 structure path (e.g., 1DJS.pdb)")
    parser.add_argument("--target-chain", required=True, help="FGFR2 receptor chain ID")
    parser.add_argument("--offtarget-structure", required=True, help="FGFR1 off-target structure path (e.g., 1CVS.pdb)")
    parser.add_argument("--offtarget-chain", required=True, help="FGFR1 receptor chain ID")
    parser.add_argument("--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    hotspot_json = _load_json(args.hotspots_json)

    target_res = _extract_chain_residues(args.target_structure, args.target_chain)
    off_res = _extract_chain_residues(args.offtarget_structure, args.offtarget_chain)
    mapping = _align_map(target_res, off_res)

    target_res_aa = {r: aa for r, aa in target_res}

    preferred = []
    conserved = []
    unmatched = []

    for h in hotspot_json.get("interface_residues", []):
        t_resseq = h.get("auth_seq_id")
        if t_resseq is None:
            continue
        t_resseq = int(t_resseq)
        t_aa = target_res_aa.get(t_resseq)
        mapped = mapping.get(t_resseq)
        if mapped is None:
            unmatched.append({"target_resseq": t_resseq, "target_aa": t_aa})
            preferred.append({"target_resseq": t_resseq, "target_aa": t_aa, "reason": "no FGFR1 alignment"})
            continue

        off_resseq, off_aa = mapped
        if t_aa == off_aa:
            conserved.append(
                {
                    "target_resseq": t_resseq,
                    "target_aa": t_aa,
                    "offtarget_resseq": off_resseq,
                    "offtarget_aa": off_aa,
                    "reason": "conserved FGFR2/FGFR1 residue",
                }
            )
        else:
            preferred.append(
                {
                    "target_resseq": t_resseq,
                    "target_aa": t_aa,
                    "offtarget_resseq": off_resseq,
                    "offtarget_aa": off_aa,
                    "reason": "non-conserved FGFR2/FGFR1 residue",
                }
            )

    warning = None
    if not preferred:
        warning = "No preferred non-conserved hotspots found. FGFR1 off-target avoidance confidence is low."

    output = {
        "target": "FGFR2",
        "off_target": "FGFR1",
        "target_structure": os.path.abspath(args.target_structure),
        "offtarget_structure": os.path.abspath(args.offtarget_structure),
        "target_chain": args.target_chain,
        "offtarget_chain": args.offtarget_chain,
        "preferred_hotspots": preferred,
        "conserved_hotspots": conserved,
        "unmatched_hotspots": unmatched,
        "warning": warning,
        "note": "Use preferred hotspots for positive targeting and down-weight conserved hotspots to improve FGFR2 specificity.",
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)

    print(f"Wrote specificity filter JSON: {args.output}")
    print(f"Preferred hotspots: {len(preferred)} | Conserved hotspots: {len(conserved)} | Unmatched: {len(unmatched)}")
    if warning:
        print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

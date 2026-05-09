#!/usr/bin/env python3
"""Aggregate per-strategy BoltzGen final metrics into one comparison table."""

from __future__ import annotations

import argparse
import csv
import glob
import os
from typing import Dict, List


def _to_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _to_boolish(value: str) -> str:
    v = str(value).strip().lower()
    if v in {"true", "1", "yes"}:
        return "true"
    if v in {"false", "0", "no"}:
        return "false"
    return value


def _latest_metrics_file(strategy_dir: str) -> str | None:
    p = os.path.join(strategy_dir, "final_ranked_designs", "final_designs_metrics_30.csv")
    return p if os.path.isfile(p) else None


def aggregate(campaign_dir: str, out_csv: str) -> int:
    boltz_root = os.path.join(campaign_dir, "out", "boltzgen")
    strategy_dirs = sorted(glob.glob(os.path.join(boltz_root, "*_round0")))

    rows: List[Dict[str, str]] = []
    for sdir in strategy_dirs:
        metrics_csv = _latest_metrics_file(sdir)
        if not metrics_csv:
            continue

        strategy_name = os.path.basename(sdir).replace("_round0", "")
        with open(metrics_csv, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for r in reader:
                rows.append(
                    {
                        "strategy": strategy_name,
                        "run_dir": sdir,
                        "id": r.get("id", ""),
                        "final_rank": r.get("final_rank", ""),
                        "pass_filters": _to_boolish(r.get("pass_filters", "")),
                        "num_filters_passed": r.get("num_filters_passed", ""),
                        "design_to_target_iptm": r.get("design_to_target_iptm", ""),
                        "min_design_to_target_pae": r.get("min_design_to_target_pae", ""),
                        "design_ptm": r.get("design_ptm", ""),
                        "filter_rmsd": r.get("filter_rmsd", ""),
                        "designfolding_filter_rmsd": r.get("designfolding-filter_rmsd", ""),
                        "plip_hbonds_refolded": r.get("plip_hbonds_refolded", ""),
                        "delta_sasa_refolded": r.get("delta_sasa_refolded", ""),
                        "quality_score": r.get("quality_score", ""),
                        "sequence": r.get("designed_sequence", ""),
                    }
                )

    # Sort by strategy then rank if possible.
    def sort_key(x: Dict[str, str]):
        try:
            rk = int(float(x.get("final_rank", "999999")))
        except Exception:
            rk = 999999
        return (x.get("strategy", ""), rk)

    rows.sort(key=sort_key)

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    headers = [
        "strategy",
        "run_dir",
        "id",
        "final_rank",
        "pass_filters",
        "num_filters_passed",
        "design_to_target_iptm",
        "min_design_to_target_pae",
        "design_ptm",
        "filter_rmsd",
        "designfolding_filter_rmsd",
        "plip_hbonds_refolded",
        "delta_sasa_refolded",
        "quality_score",
        "sequence",
    ]

    with open(out_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote strategy comparison: {out_csv}")
    print(f"Rows: {len(rows)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate strategy final metrics")
    parser.add_argument(
        "--campaign-dir",
        default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        help="Path to fgfr2_campaign",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (default: fgfr2_campaign/out/strategy_comparison.csv)",
    )
    args = parser.parse_args()

    out_csv = args.output or os.path.join(args.campaign_dir, "out", "strategy_comparison.csv")
    return aggregate(args.campaign_dir, out_csv)


if __name__ == "__main__":
    raise SystemExit(main())

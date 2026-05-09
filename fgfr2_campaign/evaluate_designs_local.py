#!/usr/bin/env python3
"""Prepare and optionally run local FGFR2 binder validation jobs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


METRIC_ALIASES = {
    "iptm": ["iptm", "ipTM", "interface_ptm"],
    "ptm": ["ptm", "pTM"],
    "min_interaction_pae": ["min_interaction_pae", "interaction_pae_min", "min_pae"],
    "ipsae": ["ipsae", "ipSAE", "interface_pSAE", "interface_psae"],
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


def _load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_sequence(seq_inline: Optional[str], fasta_path: Optional[str]) -> Optional[str]:
    if seq_inline:
        return seq_inline.strip().replace(" ", "")
    if fasta_path:
        seq_parts = []
        with open(fasta_path, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                seq_parts.append(line)
        return "".join(seq_parts)
    return None


def _score_value(design: Dict[str, Any]) -> Optional[float]:
    return _safe_float(design.get("specificity_adjusted_score")) or _safe_float(design.get("composite_score"))


def _select_candidates(index: Dict[str, Any], top_n: int, force: bool) -> List[Dict[str, Any]]:
    designs = index.get("designs", [])
    pending = []
    for design in designs:
        if not force:
            ev = design.get("evaluation") or {}
            if ev.get("status") == "done":
                continue
        pending.append(design)

    pending.sort(key=lambda item: (_score_value(item) is None, -(_score_value(item) or 0.0)))
    return pending[:top_n]


def _write_fasta(path: str, receptor_name: str, receptor_seq: str, binder_seq: str, design_id: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{receptor_name}_{design_id}\n")
        handle.write(receptor_seq + "\n")
        handle.write(f">binder_{design_id}\n")
        handle.write(binder_seq + "\n")


def _flatten_metrics(data: Dict[str, Any], prefix: str = "") -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    for key, value in data.items():
        key_text = str(key)
        prefixed = f"{prefix}{key_text}" if prefix else key_text
        if isinstance(value, dict):
            flat.update(_flatten_metrics(value, prefix=f"{prefixed}."))
        else:
            flat[key_text] = value
            flat[prefixed] = value
    return flat


def _load_metric_file(path: str) -> Dict[str, Any]:
    try:
        if path.endswith(".json"):
            return _flatten_metrics(_load_json(path))
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                return dict(row)
    except Exception:
        return {}
    return {}


def _collect_metrics_from_dir(root: str) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    if not root or not os.path.isdir(root):
        return merged
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            lower = filename.lower()
            if not lower.endswith((".json", ".csv")):
                continue
            if filename == "evaluation_metrics.json":
                continue
            merged.update(_load_metric_file(os.path.join(dirpath, filename)))
    return merged


def _get_first_metric(metrics: Dict[str, Any], aliases: List[str]) -> Any:
    lowered = {str(key).lower(): value for key, value in metrics.items()}
    for alias in aliases:
        if alias in metrics:
            return metrics[alias]
        found = lowered.get(alias.lower())
        if found is not None:
            return found
    return None


def _quote_command(parts: List[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def _build_conda_run_cmd(conda_env: Optional[str], command: List[str]) -> List[str]:
    if not conda_env:
        return command
    return ["conda", "run", "-n", conda_env, *command]


def _build_boltz_eval_command(eval_cmd_parts: List[str], input_fasta: str, output_dir: str) -> List[str]:
    """Build evaluator command with CLI-specific flags.

    Supports both legacy boltz2_eval wrappers and native boltz predict CLI.
    """
    if not eval_cmd_parts:
        return []

    cmd = list(eval_cmd_parts)
    executable = os.path.basename(cmd[0]).lower()
    has_predict_subcmd = len(cmd) > 1 and cmd[1].lower() == "predict"

    if executable == "boltz":
        if not has_predict_subcmd:
            cmd.append("predict")
        return [*cmd, input_fasta, "--out_dir", output_dir, "--model", "boltz2"]

    return [*cmd, "--input_fasta", input_fasta, "--output_dir", output_dir, "--max_jobs", "1"]


def _command_available(command_name: str, conda_env: Optional[str]) -> bool:
    if not conda_env:
        return shutil.which(command_name) is not None
    if shutil.which("conda") is None:
        return False
    probe_cmd = [
        "conda",
        "run",
        "-n",
        conda_env,
        "python",
        "-c",
        f"import shutil, sys; sys.exit(0 if shutil.which({command_name!r}) else 1)",
    ]
    completed = subprocess.run(probe_cmd, check=False)
    return completed.returncode == 0


def _extract_boltz2_metrics(target_dir: Optional[str], off_target_dir: Optional[str]) -> Tuple[Dict[str, Any], List[str]]:
    warnings: List[str] = []
    target_raw = _collect_metrics_from_dir(target_dir or "")
    off_target_raw = _collect_metrics_from_dir(off_target_dir or "")

    metrics = {
        "ipTM": _safe_float(_get_first_metric(target_raw, METRIC_ALIASES["iptm"])),
        "pTM": _safe_float(_get_first_metric(target_raw, METRIC_ALIASES["ptm"])),
        "min_interaction_pae": _safe_float(_get_first_metric(target_raw, METRIC_ALIASES["min_interaction_pae"])),
        "off_target_ipTM": _safe_float(_get_first_metric(off_target_raw, METRIC_ALIASES["iptm"])),
        "off_target_pTM": _safe_float(_get_first_metric(off_target_raw, METRIC_ALIASES["ptm"])),
        "ipSAE": _safe_float(_get_first_metric(target_raw, METRIC_ALIASES["ipsae"])),
        "off_target_ipSAE": _safe_float(_get_first_metric(off_target_raw, METRIC_ALIASES["ipsae"])),
    }

    if target_dir and os.path.isdir(target_dir) and metrics["ipTM"] is None:
        warnings.append("Boltz-2 target evaluation finished but target ipTM was not found in output files.")
    if off_target_dir and os.path.isdir(off_target_dir) and metrics["off_target_ipTM"] is None:
        warnings.append("Boltz-2 off-target evaluation finished but FGFR1 off-target ipTM was not found in output files.")
    if target_dir and os.path.isdir(target_dir) and metrics["ipSAE"] is None:
        warnings.append("Boltz-2 target evaluation did not expose ipSAE in parsed output files.")
    if off_target_dir and os.path.isdir(off_target_dir) and metrics["off_target_ipSAE"] is None:
        warnings.append("Boltz-2 off-target evaluation did not expose ipSAE in parsed output files.")

    return metrics, warnings


def _compute_ipsae_summary(metrics: Dict[str, Any]) -> Dict[str, Any]:
    target_ipsae = _safe_float(metrics.get("ipSAE"))
    off_target_ipsae = _safe_float(metrics.get("off_target_ipSAE"))

    summary = {
        "status": "missing",
        "source": "boltz2_outputs",
        "metrics": {},
        "warnings": [],
    }

    if target_ipsae is not None:
        summary["metrics"]["ipSAE"] = target_ipsae
    if off_target_ipsae is not None:
        summary["metrics"]["off_target_ipSAE"] = off_target_ipsae

    if target_ipsae is not None and off_target_ipsae is not None:
        summary["metrics"]["ipSAE_specificity_gap"] = round(target_ipsae - off_target_ipsae, 6)
        summary["status"] = "done"
    elif target_ipsae is not None or off_target_ipsae is not None:
        summary["status"] = "partial"
        summary["warnings"].append("Only one ipSAE value was available; specificity gap could not be computed.")
    else:
        summary["warnings"].append("ipSAE was not found in the current Boltz-2 evaluation outputs.")

    return summary


def _compute_pyrosetta_metrics(structure_path: str) -> Tuple[Dict[str, Any], List[str]]:
    import pyrosetta  # type: ignore

    warnings: List[str] = []
    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init("-mute all")

    pose = pyrosetta.Pose()
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, structure_path)

    scorefxn = pyrosetta.get_fa_scorefxn()
    metrics: Dict[str, Any] = {
        "pyrosetta_total_score": float(scorefxn(pose)),
        "pyrosetta_chain_count": int(pose.num_chains()),
    }

    chain_info = []
    pdb_info = pose.pdb_info()
    for chain_index in range(1, pose.num_chains() + 1):
        begin = pose.conformation().chain_begin(chain_index)
        end = pose.conformation().chain_end(chain_index)
        chain_letter = pdb_info.chain(begin).strip() if pdb_info is not None else ""
        if not chain_letter:
            chain_letter = chr(64 + chain_index) if chain_index <= 26 else f"C{chain_index}"
        chain_info.append((chain_index, chain_letter, int(end - begin + 1)))

    if len(chain_info) < 2:
        warnings.append("PyRosetta interface analysis skipped because the structure has fewer than 2 chains.")
        return metrics, warnings

    if any(len(chain_letter) != 1 for _, chain_letter, _ in chain_info):
        warnings.append("PyRosetta interface analysis skipped because chain IDs are not single-character PDB IDs.")
        return metrics, warnings

    binder_chain = min(chain_info, key=lambda item: item[2])
    receptor_group = "".join(chain_letter for _, chain_letter, _ in chain_info if chain_letter != binder_chain[1])
    if not receptor_group:
        warnings.append("PyRosetta interface analysis skipped because a receptor chain group could not be determined.")
        return metrics, warnings

    interface = f"{receptor_group}_{binder_chain[1]}"
    metrics["pyrosetta_interface"] = interface
    metrics["pyrosetta_binder_chain"] = binder_chain[1]

    iam_cls = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover
    mover = None
    builders = [
        lambda: iam_cls(interface),
        lambda: iam_cls(interface, False),
        lambda: iam_cls(interface, False, scorefxn, False, False, False, False),
    ]
    for builder in builders:
        try:
            mover = builder()
            break
        except TypeError:
            continue

    if mover is None:
        warnings.append("PyRosetta interface analyzer could not be constructed with the available API signature.")
        return metrics, warnings

    if hasattr(mover, "set_scorefunction"):
        mover.set_scorefunction(scorefxn)
    if hasattr(mover, "set_compute_packstat"):
        mover.set_compute_packstat(True)

    mover.apply(pose)

    getter_map = {
        "get_interface_dG": "pyrosetta_interface_dG",
        "get_interface_delta_sasa": "pyrosetta_interface_delta_sasa",
        "get_interface_packstat": "pyrosetta_packstat",
        "get_packstat": "pyrosetta_packstat",
    }
    for getter_name, metric_name in getter_map.items():
        getter = getattr(mover, getter_name, None)
        if not callable(getter):
            continue
        try:
            value = getter()
        except Exception:
            continue
        numeric = _safe_float(value)
        if numeric is not None:
            metrics[metric_name] = numeric

    return metrics, warnings


def _run_pyrosetta_validation(
    structure_path: Optional[str],
    output_path: str,
    pyrosetta_conda_env: Optional[str],
) -> Dict[str, Any]:
    result = {
        "status": "todo",
        "tool": "pyrosetta",
        "metrics": {},
        "warnings": [],
        "output_json": os.path.abspath(output_path),
    }

    if not structure_path or not os.path.isfile(structure_path):
        result["status"] = "missing"
        result["warnings"].append("Structure path is missing; PyRosetta validation was skipped.")
    else:
        if pyrosetta_conda_env:
            if shutil.which("conda") is None:
                result["status"] = "unavailable"
                result["warnings"].append(
                    f"Conda was not found on PATH; cannot run PyRosetta env '{pyrosetta_conda_env}'."
                )
                with open(output_path, "w", encoding="utf-8") as handle:
                    json.dump(result, handle, indent=2)
                return result

            script_path = os.path.abspath(__file__)
            cmd = [
                "conda",
                "run",
                "-n",
                pyrosetta_conda_env,
                "python",
                script_path,
                "--pyrosetta-only",
                "--structure-path",
                os.path.abspath(structure_path),
                "--pyrosetta-output-path",
                os.path.abspath(output_path),
            ]
            completed = subprocess.run(cmd, check=False)
            if completed.returncode == 0 and os.path.isfile(output_path):
                payload = _load_json(output_path)
                result.update(
                    {
                        "status": payload.get("status", "done"),
                        "metrics": payload.get("metrics", {}),
                        "warnings": payload.get("warnings", []),
                    }
                )
            else:
                result["status"] = "failed"
                result["warnings"].append(
                    f"PyRosetta conda run failed with code {completed.returncode}: {_quote_command(cmd)}"
                )
            with open(output_path, "w", encoding="utf-8") as handle:
                json.dump(result, handle, indent=2)
            return result

        external_cmd = os.environ.get("PYROSETTA_VALIDATE_CMD")
        if external_cmd:
            cmd = shlex.split(external_cmd) + ["--structure", os.path.abspath(structure_path), "--output", os.path.abspath(output_path)]
            completed = subprocess.run(cmd, check=False)
            if completed.returncode == 0 and os.path.isfile(output_path):
                payload = _load_json(output_path)
                result.update({
                    "status": payload.get("status", "done"),
                    "metrics": payload.get("metrics", {}),
                    "warnings": payload.get("warnings", []),
                })
            else:
                result["status"] = "failed"
                result["warnings"].append(
                    f"External PyRosetta validator exited with code {completed.returncode}: {_quote_command(cmd)}"
                )
        else:
            try:
                metrics, warnings = _compute_pyrosetta_metrics(structure_path)
                result["status"] = "done"
                result["metrics"] = metrics
                result["warnings"] = warnings
            except Exception as exc:
                result["status"] = "unavailable"
                result["warnings"].append(
                    "PyRosetta validation unavailable in the current environment. "
                    f"Install/configure PyRosetta or set PYROSETTA_VALIDATE_CMD. Details: {exc}"
                )

    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    return result


def _merge_metrics(dst: Dict[str, Any], src: Dict[str, Any]) -> None:
    for key, value in src.items():
        if value is not None:
            dst[key] = value


def _combine_statuses(*statuses: str) -> str:
    active = [status for status in statuses if status and status not in {"missing", "skipped"}]
    if any(status == "failed" for status in active):
        return "failed"
    if any(status == "done" for status in active):
        if any(status in {"partial", "todo", "dry-run", "unavailable"} for status in active):
            return "partial"
        return "done"
    if any(status == "partial" for status in active):
        return "partial"
    if any(status == "unavailable" for status in active):
        return "unavailable"
    if any(status == "dry-run" for status in active):
        return "dry-run"
    return "todo"


def _run_pyrosetta_only(structure_path: str, output_path: str) -> int:
    result = {
        "status": "unavailable",
        "tool": "pyrosetta",
        "metrics": {},
        "warnings": [],
        "output_json": os.path.abspath(output_path),
    }
    try:
        metrics, warnings = _compute_pyrosetta_metrics(structure_path)
        result["status"] = "done"
        result["metrics"] = metrics
        result["warnings"] = warnings
    except Exception as exc:
        result["status"] = "unavailable"
        result["warnings"].append(f"PyRosetta execution failed in conda env: {exc}")

    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)

    return 0 if result["status"] == "done" else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate top FGFR2 computational binder candidates locally")
    parser.add_argument(
        "--campaign-dir",
        default=os.path.abspath(os.path.join(os.path.dirname(__file__), ".")),
        help="Path to fgfr2_campaign",
    )
    parser.add_argument("--top-n", type=int, default=10, help="Top unevaluated candidates to process (max recommended 10)")
    parser.add_argument("--receptor-seq", default=None, help="FGFR2 receptor sequence for complex evaluation")
    parser.add_argument("--receptor-fasta", default=None, help="FASTA file with FGFR2 receptor sequence")
    parser.add_argument("--offtarget-seq", default=None, help="FGFR1 off-target receptor sequence")
    parser.add_argument("--offtarget-fasta", default=None, help="FASTA file with FGFR1 receptor sequence")
    parser.add_argument("--dry-run", action="store_true", help="Prepare commands and files without running evaluation")
    parser.add_argument("--force", action="store_true", help="Re-run validation even if evaluation status is already done")
    parser.add_argument("--no-sync", action="store_true", help="Do not call sync_designs_local.py at end")
    parser.add_argument(
        "--boltz2-conda-env",
        default=os.environ.get("BOLTZ2_CONDA_ENV", "boltz2"),
        help="Conda environment used to run Boltz-2 evaluation commands (set empty string to disable).",
    )
    parser.add_argument(
        "--pyrosetta-conda-env",
        default=os.environ.get("PYROSETTA_CONDA_ENV", "pyrosetta"),
        help="Conda environment used to run PyRosetta validation (set empty string to disable).",
    )
    parser.add_argument("--pyrosetta-only", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--structure-path", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--pyrosetta-output-path", default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.pyrosetta_only:
        if not args.structure_path or not args.pyrosetta_output_path:
            raise ValueError("--pyrosetta-only requires --structure-path and --pyrosetta-output-path")
        return _run_pyrosetta_only(args.structure_path, args.pyrosetta_output_path)

    if args.top_n > 10:
        print("WARNING: top-n > 10 requested. For A40 hackathon workflows, keep first runs <= 10.")

    campaign_dir = os.path.abspath(args.campaign_dir)
    index_path = os.path.join(campaign_dir, "designs", "index.json")
    if not os.path.isfile(index_path):
        raise FileNotFoundError(f"Missing index file: {index_path}. Run sync_designs_local.py first.")

    index = _load_json(index_path)
    selected = _select_candidates(index, args.top_n, args.force)

    receptor_seq = _read_sequence(args.receptor_seq, args.receptor_fasta)
    off_target_seq = _read_sequence(args.offtarget_seq, args.offtarget_fasta)

    out_eval_root = os.path.join(campaign_dir, "out", "boltz2", "evaluations")
    out_input_root = os.path.join(campaign_dir, "out", "boltz2", "eval_inputs")
    os.makedirs(out_eval_root, exist_ok=True)
    os.makedirs(out_input_root, exist_ok=True)

    todo_path = os.path.join(campaign_dir, "out", "boltz2", "eval_todos.sh")
    todo_lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]

    boltz2_conda_env = (args.boltz2_conda_env or "").strip() or None
    pyrosetta_conda_env = (args.pyrosetta_conda_env or "").strip() or None

    eval_cmd = os.environ.get("BOLTZ2_EVAL_CMD", "boltz2_eval")
    eval_cmd_parts = shlex.split(eval_cmd)
    if not eval_cmd_parts:
        raise ValueError("BOLTZ2_EVAL_CMD resolved to an empty command.")

    has_eval_cmd = _command_available(eval_cmd_parts[0], boltz2_conda_env)
    has_boltz_cli = _command_available("boltz", boltz2_conda_env)

    if receptor_seq is None:
        print("WARNING: FGFR2 receptor sequence not provided. Target evaluation commands will be TODO only.")
    if off_target_seq is None:
        print("WARNING: FGFR1 off-target sequence not provided. Specificity evaluation commands will be TODO only.")

    for design in selected:
        did = design["design_id"]
        metrics = design.get("metrics") or {}
        binder_seq = (metrics.get("binder_sequence") or "").strip()
        structure_path = design.get("design_path") or design.get("source_path")

        eval_dir = os.path.join(out_eval_root, did)
        os.makedirs(eval_dir, exist_ok=True)

        target_fasta = os.path.join(out_input_root, f"{did}.target.fasta")
        off_target_fasta = os.path.join(out_input_root, f"{did}.offtarget.fasta")
        target_out = os.path.join(eval_dir, "target_fgfr2")
        off_target_out = os.path.join(eval_dir, "offtarget_fgfr1")
        pyrosetta_out = os.path.join(eval_dir, "pyrosetta_metrics.json")

        warnings = ["Computational binder candidate only. Not experimentally validated."]
        aggregated_metrics: Dict[str, Any] = {}

        boltz2_result = {
            "status": "todo",
            "tool": "boltz2",
            "metrics": {},
            "warnings": [],
            "target_input_fasta": None,
            "offtarget_input_fasta": None,
            "target_output_dir": os.path.abspath(target_out),
            "offtarget_output_dir": os.path.abspath(off_target_out),
        }

        command_specs: List[Tuple[str, List[str]]] = []
        if not binder_seq:
            boltz2_result["warnings"].append("Missing binder sequence in metrics; cannot create Boltz-2 evaluation FASTA files.")
        else:
            if receptor_seq:
                _write_fasta(target_fasta, "FGFR2_receptor", receptor_seq, binder_seq, did)
                boltz2_result["target_input_fasta"] = os.path.abspath(target_fasta)
                target_cmd = _build_boltz_eval_command(eval_cmd_parts, target_fasta, target_out)
                command_specs.append(("target_fgfr2", _build_conda_run_cmd(boltz2_conda_env, target_cmd)))
            else:
                boltz2_result["warnings"].append("Missing FGFR2 receptor sequence for target evaluation.")

            if off_target_seq:
                _write_fasta(off_target_fasta, "FGFR1_offtarget", off_target_seq, binder_seq, did)
                boltz2_result["offtarget_input_fasta"] = os.path.abspath(off_target_fasta)
                off_target_cmd = _build_boltz_eval_command(eval_cmd_parts, off_target_fasta, off_target_out)
                command_specs.append(("offtarget_fgfr1", _build_conda_run_cmd(boltz2_conda_env, off_target_cmd)))
            else:
                boltz2_result["warnings"].append("Missing FGFR1 off-target sequence; specificity evaluation skipped.")

        ran_boltz2 = False
        for _, cmd in command_specs:
            todo_lines.append(_quote_command(cmd))
            if args.dry_run:
                print("DRY RUN:", _quote_command(cmd))
                boltz2_result["status"] = "dry-run"
                continue
            if not has_eval_cmd:
                if boltz2_conda_env:
                    boltz2_result["warnings"].append(
                        f"Evaluator command '{eval_cmd_parts[0]}' is unavailable in conda env '{boltz2_conda_env}'."
                    )
                    if eval_cmd_parts[0] == "boltz2_eval" and has_boltz_cli:
                        boltz2_result["warnings"].append(
                            "Found 'boltz' CLI in boltz2 env. If needed, set BOLTZ2_EVAL_CMD to your Boltz evaluation wrapper."
                        )
                else:
                    boltz2_result["warnings"].append(f"Local evaluator command not found: {eval_cmd_parts[0]}")
                boltz2_result["status"] = "todo"
                continue
            print("Running:", _quote_command(cmd))
            completed = subprocess.run(cmd, check=False)
            ran_boltz2 = True
            if completed.returncode == 0:
                boltz2_result["status"] = "done"
            else:
                boltz2_result["status"] = "failed"
                boltz2_result["warnings"].append(
                    f"Evaluation command exited with code {completed.returncode}: {_quote_command(cmd)}"
                )

        boltz_metrics, boltz_warnings = _extract_boltz2_metrics(target_out, off_target_out)
        boltz2_result["warnings"].extend(boltz_warnings)
        boltz2_result["metrics"] = boltz_metrics
        if not ran_boltz2 and any(value is not None for value in boltz_metrics.values()):
            boltz2_result["status"] = "done"
        _merge_metrics(aggregated_metrics, boltz_metrics)

        ipsae_result = _compute_ipsae_summary(aggregated_metrics)
        _merge_metrics(aggregated_metrics, ipsae_result["metrics"])

        pyrosetta_result = _run_pyrosetta_validation(structure_path, pyrosetta_out, pyrosetta_conda_env)
        _merge_metrics(aggregated_metrics, pyrosetta_result.get("metrics", {}))

        warnings.extend(boltz2_result["warnings"])
        warnings.extend(ipsae_result["warnings"])
        warnings.extend(pyrosetta_result.get("warnings", []))

        result = {
            "design_id": did,
            "status": _combine_statuses(
                boltz2_result["status"],
                ipsae_result["status"],
                pyrosetta_result["status"],
            ),
            "tool": "local_validation",
            "structure_path": os.path.abspath(structure_path) if structure_path else None,
            "output_dir": os.path.abspath(eval_dir),
            "metrics": aggregated_metrics,
            "validations": {
                "boltz2": boltz2_result,
                "ipsae": ipsae_result,
                "pyrosetta": pyrosetta_result,
            },
            "warnings": warnings,
            "updated_at": now_utc(),
        }
        with open(os.path.join(eval_dir, "evaluation_metrics.json"), "w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)

    with open(todo_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(todo_lines) + "\n")
    os.chmod(todo_path, 0o755)

    print(f"Wrote evaluation TODO commands: {todo_path}")

    if not args.no_sync:
        sync_script = os.path.join(campaign_dir, "sync_designs_local.py")
        if os.path.isfile(sync_script):
            subprocess.run([sys.executable, sync_script, "--campaign-dir", campaign_dir], check=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

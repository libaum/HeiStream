#!/usr/bin/env python3
"""
Utility to launch BuffCut experiments described by JSON configs.

Example:
    python tools/run_experiments.py --config experiments/configs/buffer_sweep.json
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Any, Optional


def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a batch of BuffCut experiments defined in a config file."
    )
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to the JSON configuration that defines the experiment.",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("experiments") / "results",
        help="Folder where run artifacts should be written "
        "(default: experiments/results).",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=None,
        help="Override the executable specified in the config (e.g., build/buffcut).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the commands that would be executed without running them.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue executing remaining runs even if one fails.",
    )
    return parser.parse_args()


def load_config(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)


def merge_arg_maps(maps: Iterable[Optional[Dict[str, Any]]]) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    for mapping in maps:
        if not mapping:
            continue
        for key, value in mapping.items():
            merged[key] = value
    return merged


def flatten_args(arg_map: Dict[str, Any]) -> List[str]:
    args: List[str] = []
    for key, value in arg_map.items():
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                args.append(str(key))
            continue
        if isinstance(value, list):
            for entry in value:
                args.append(str(key))
                args.append(str(entry))
            continue
        args.append(str(key))
        args.append(str(value))
    return args


def ensure_executable(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(f"Executable '{path}' does not exist. Build the project first.")
    if not os.access(path, os.X_OK):
        raise SystemExit(f"Executable '{path}' is not marked as executable.")
    return path


def snapshot_config(config: Dict[str, Any], experiment_root: Path) -> None:
    experiment_root.mkdir(parents=True, exist_ok=True)
    snapshot_path = experiment_root / "config_snapshot.json"
    snapshot_path.write_text(json.dumps(config, indent=2) + "\n")


def write_meta(meta: Dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(meta, indent=2) + "\n")


def run_experiments(cfg: Dict[str, Any], args: argparse.Namespace) -> None:
    experiment_name = cfg.get("experiment_name")
    if not experiment_name:
        raise SystemExit("Config must contain 'experiment_name'.")

    binary_path = args.binary or cfg.get("base_command")
    if not binary_path:
        raise SystemExit("Config must specify 'base_command' or pass --binary.")
    binary_path = ensure_executable(Path(binary_path)).resolve()

    results_root = args.results_root / experiment_name
    snapshot_config(cfg, results_root)

    graphs = cfg.get("graphs", [])
    if not graphs:
        raise SystemExit("Config must list at least one graph in 'graphs'.")

    variants = cfg.get("variants", [])
    if not variants:
        raise SystemExit("Config must list at least one entry under 'variants'.")

    k_values = cfg.get("k_values")
    if k_values is not None:
        if not isinstance(k_values, list) or not k_values:
            raise SystemExit("'k_values' must be a non-empty list when provided.")
        k_values = [int(k) for k in k_values]
    else:
        k_values = [None]

    for graph in graphs:
        graph_name = graph.get("name")
        graph_path = graph.get("path")
        if not graph_name or not graph_path:
            raise SystemExit("Each graph needs 'name' and 'path'.")
        graph_abs = Path(graph_path).resolve()
        if not graph_abs.exists():
            raise SystemExit(f"Graph file '{graph_abs}' does not exist.")
        graph_args = graph.get("args", {})

        for k_value in k_values:
            for variant in variants:
                variant_name = variant.get("name")
                if not variant_name:
                    raise SystemExit("Each variant needs a 'name'.")
                variant_args = variant.get("args", {})

                run_dir = results_root / graph_name
                if k_value is not None:
                    run_dir = run_dir / f"k_{k_value}"
                run_dir = run_dir / variant_name
                run_dir.mkdir(parents=True, exist_ok=True)
                combined_args = merge_arg_maps(
                    [cfg.get("common_args"), graph_args, variant_args]
                )
                if k_value is not None:
                    combined_args["--k"] = k_value
                if "--k" not in combined_args:
                    raise SystemExit("No value for --k supplied (set common_args['--k'] or provide k_values).")

                cmd = [str(binary_path)] + flatten_args(combined_args) + [str(graph_abs)]

                timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
                log_path = run_dir / "run.log"
                meta = {
                    "experiment": experiment_name,
                    "graph": graph_name,
                    "graph_path": str(graph_abs),
                    "variant": variant_name,
                    "k": combined_args.get("--k"),
                    "command": cmd,
                    "combined_args": combined_args,
                    "run_directory": str(run_dir.resolve()),
                    "timestamp_utc": timestamp,
                }

                if args.dry_run:
                    print("[DRY-RUN]", " ".join(cmd))
                    continue

                log_path.write_text(
                    f"# Command: {' '.join(cmd)}\n# Started: {timestamp} UTC\n\n"
                )
                start = time.time()
                with log_path.open("a") as log_file:
                    result = subprocess.run(
                        cmd,
                        cwd=run_dir,
                        stdout=log_file,
                        stderr=subprocess.STDOUT,
                        check=False,
                    )
                duration = time.time() - start
                bin_files = sorted(p.name for p in run_dir.glob("*.bin"))

                meta.update(
                    {
                        "duration_seconds": duration,
                        "returncode": result.returncode,
                        "status": "success" if result.returncode == 0 else "failed",
                        "bin_files": bin_files,
                    }
                )
                write_meta(meta, run_dir / "run_meta.json")

                if result.returncode != 0:
                    msg = f"Run {experiment_name}/{graph_name}/k_{combined_args.get('--k')}/{variant_name} failed."
                    if args.keep_going:
                        print(msg, file=sys.stderr)
                        continue
                    raise SystemExit(msg)


def main() -> None:
    cli_args = parse_cli()
    config = load_config(cli_args.config)
    run_experiments(config, cli_args)


if __name__ == "__main__":
    main()

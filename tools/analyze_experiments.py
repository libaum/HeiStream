#!/usr/bin/env python3
"""
Summarize results from experiments executed via tools/run_experiments.py.

Usage:
    python tools/analyze_experiments.py experiments/results/buffer_sweep
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from typing import Dict, Iterable, List, Tuple

from tools import parse_partition_log as ppl


def parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate metrics from run directories containing FlatBuffers."
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="Run directories, experiment folders, or specific run_meta.json files.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Optional path to write the tabular output as CSV.",
    )
    return parser.parse_args()


def discover_run_dirs(paths: Iterable[Path]) -> List[Path]:
    run_dirs: List[Path] = []
    for path in paths:
        if path.is_file():
            if path.name == "run_meta.json":
                run_dirs.append(path.parent.resolve())
            elif path.suffix == ".bin":
                run_dirs.append(path.parent.resolve())
        elif path.is_dir():
            for meta_file in path.rglob("run_meta.json"):
                run_dirs.append(meta_file.parent.resolve())
    return run_dirs


def load_meta(run_dir: Path) -> Dict:
    meta_path = run_dir / "run_meta.json"
    if not meta_path.exists():
        raise SystemExit(f"No run_meta.json found in {run_dir}")
    return json.loads(meta_path.read_text())


def pick_log_file(run_dir: Path, meta: Dict) -> Path:
    preferred = meta.get("bin_files", [])
    for rel in preferred:
        candidate = run_dir / rel
        if candidate.exists():
            return candidate
    bins = sorted(run_dir.glob("*.bin"))
    if not bins:
        raise SystemExit(f"No .bin files found in {run_dir}")
    return bins[0]


def extract_from_log(bin_path: Path) -> Dict:
    data = bytearray(bin_path.read_bytes())
    log = ppl.get_root(data)
    metadata = log.GraphMetadata()
    config = log.PartitionConfiguration()
    runtime = log.RunTime()
    memory = log.MemoryConsumption()
    metrics = log.Metrics()

    return {
        "graph_file": metadata.Filename() if metadata else "",
        "num_nodes": metadata.NumNodes() if metadata else 0,
        "num_edges": metadata.NumEdges() if metadata else 0,
        "k": config.K() if config else 0,
        "batch_size": config.BatchSize() if config else 0,
        "buffer_size": config.MaxBufferSize() if config else 0,
        "buffer_score": config.BufferScoreType() if config else -1,
        "haa_beta": config.HaaBeta() if config else 0.0,
        "haa_theta": config.HaaTheta() if config else 0.0,
        "ghost_enabled": config.GhostNeighborsEnabled() if config else False,
        "num_streams_passes": config.NumStreamsPasses() if config else 0,
        "buffer_neighbor_w": config.BufferNeighborWeight() if config else 0.0,
        "batch_frontier_w": config.BatchFrontierWeight() if config else 0.0,
        "collect_locality_metrics": config.CollectLocalityMetrics() if config else False,
        "io_time": runtime.IoTime() if runtime else 0.0,
        "partition_time": runtime.PartitionTime() if runtime else 0.0,
        "model_construction_time": runtime.ModelConstructionTime() if runtime else 0.0,
        "mapping_time": runtime.MappingTime() if runtime else 0.0,
        "total_time": runtime.TotalTime() if runtime else 0.0,
        "max_rss_kb": memory.MaxRSS() if memory else 0,
        "edge_cut": metrics.EdgeCut() if metrics else 0,
        "vertex_cut": metrics.VertexCut() if metrics else 0,
        "replicas": metrics.Replicas() if metrics else 0,
        "replication_factor": metrics.ReplicationFactor() if metrics else 0.0,
        "balance": metrics.Balance() if metrics else 0.0,
        "internal_edge_ratio": metrics.InternalEdgeRatio() if metrics else 0.0,
        "neighbor_coverage": metrics.NeighborCoverage() if metrics else 0.0,
        "conductance": metrics.Conductance() if metrics else 0.0,
    }


def format_buffer_score(score_type: int) -> str:
    reverse_map = {
        0: "cbs",
        1: "cbsq",
        2: "anr",
        3: "haa",
        4: "cms",
        5: "nss",
        6: "gts",
    }
    return reverse_map.get(score_type, str(score_type))


def collect_rows(run_dirs: Iterable[Path]) -> List[Dict]:
    rows: List[Dict] = []
    for run_dir in run_dirs:
        meta = load_meta(run_dir)
        log_file = pick_log_file(run_dir, meta)
        log_data = extract_from_log(log_file)

        row = {
            "experiment": meta.get("experiment", ""),
            "graph": meta.get("graph", ""),
            "variant": meta.get("variant", ""),
            "run_dir": str(run_dir),
            "status": meta.get("status", ""),
            "returncode": meta.get("returncode", 0),
            "duration": meta.get("duration_seconds", 0.0),
            **log_data,
        }
        row["buffer_score"] = format_buffer_score(row["buffer_score"])
        rows.append(row)
    return rows


def print_table(rows: List[Dict], output_csv: Path = None) -> None:
    if not rows:
        print("No runs found.")
        return

    columns = [
        "experiment",
        "graph",
        "variant",
        "status",
        "k",
        "batch_size",
        "buffer_size",
        "buffer_score",
        "ghost_enabled",
        "partition_time",
        "total_time",
        "max_rss_kb",
        "edge_cut",
        "replication_factor",
        "balance",
        "internal_edge_ratio",
        "neighbor_coverage",
        "conductance",
    ]

    header = ",".join(columns)
    lines = [header]
    for row in rows:
        values = [str(row.get(col, "")) for col in columns]
        lines.append(",".join(values))

    print("\n".join(lines))
    if output_csv:
        output_csv.write_text("\n".join(lines) + "\n")


def summarize(rows: List[Dict]) -> None:
    grouped: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
    grouped_by_k: Dict[Tuple[str, str, int], List[Dict]] = defaultdict(list)
    for row in rows:
        key = (row["experiment"], row["variant"])
        grouped[key].append(row)
        grouped_by_k[(row["experiment"], row["variant"], row.get("k", 0))].append(row)

    def avg(entries: List[Dict], field: str) -> float:
        if not entries:
            return 0.0
        return sum(e.get(field, 0.0) for e in entries) / len(entries)

    print("\n== Variant averages (all k combined) ==")
    print("experiment,variant,num_runs,partition_time,total_time,edge_cut,replication_factor,internal_edge_ratio,neighbor_coverage,conductance")
    for (experiment, variant), entries in grouped.items():
        print(
            f"{experiment},{variant},{len(entries)},"
            f"{avg(entries, 'partition_time'):.6f},"
            f"{avg(entries, 'total_time'):.6f},"
            f"{avg(entries, 'edge_cut'):.2f},"
            f"{avg(entries, 'replication_factor'):.6f},"
            f"{avg(entries, 'internal_edge_ratio'):.6f},"
            f"{avg(entries, 'neighbor_coverage'):.6f},"
            f"{avg(entries, 'conductance'):.6f}"
        )

    print("\n== Variant averages grouped by k ==")
    print("experiment,variant,k,num_runs,partition_time,total_time,edge_cut,replication_factor,internal_edge_ratio,neighbor_coverage,conductance")
    for (experiment, variant, k_val), entries in sorted(grouped_by_k.items(), key=lambda x: (x[0][0], x[0][1], x[0][2])):
        print(
            f"{experiment},{variant},{k_val},{len(entries)},"
            f"{avg(entries, 'partition_time'):.6f},"
            f"{avg(entries, 'total_time'):.6f},"
            f"{avg(entries, 'edge_cut'):.2f},"
            f"{avg(entries, 'replication_factor'):.6f},"
            f"{avg(entries, 'internal_edge_ratio'):.6f},"
            f"{avg(entries, 'neighbor_coverage'):.6f},"
            f"{avg(entries, 'conductance'):.6f}"
        )


def main() -> None:
    args = parse_cli()
    run_dirs = discover_run_dirs(args.paths)
    rows = collect_rows(run_dirs)
    print_table(rows, args.output_csv)
    summarize(rows)


if __name__ == "__main__":
    main()

# Experiment harness

This folder contains JSON configurations that describe repeatable BuffCut experiments.
Every config is executed by `tools/run_experiments.py`, which handles parameter sweeps,
organizes results, and records metadata needed for reproducible analysis.

## Running an experiment

1. Build a BuffCut binary (e.g., `cmake -S . -B build && cmake --build build -j`).
2. Launch the experiment runner sequentially:

   ```bash
   python tools/run_experiments.py --config experiments/configs/buffer_sweep.json
   ```

   Use `--dry-run` to inspect the commands without executing them, `--binary` to
   override the executable, and `--keep-going` to continue after failures.

Each run creates `experiments/results/<experiment>/<graph>/k_<k>/<variant>/`
containing:

- `run.log`: stdout/stderr from the BuffCut binary.
- `run_meta.json`: the merged argument map, command, duration, and generated `.bin` files.
- Flatbuffer outputs produced via `--write_log`, scoped to that directory.

## Config schema

```json
{
  "experiment_name": "buffer_sweep",
  "description": "Increase buffer size while keeping batch size fixed.",
  "base_command": "build/buffcut",
  "k_values": [4, 8, 16, 32, 64, 128, 256],
  "graphs": [
    {"name": "delaunay_n15", "path": "examples/delaunay_n15.graph"}
  ],
  "common_args": {
    "--batch_size": 16384,
    "--b_score": "haa",
    "--write_log": true,
    "--collect_locality_metrics": true
  },
  "variants": [
    {"name": "buffer_32k", "args": {"--buffer_size": 32768}},
    {"name": "buffer_64k", "args": {"--buffer_size": 65536}}
  ]
}
```

Fields:

- `experiment_name`: also used as the result folder name.
- `base_command`: path to the executable (overridden via `--binary`).
- `k_values`: optional list of k targets; if omitted, provide `--k` via `common_args`.
- `graphs`: list of `{name, path}` entries; optional per-graph overrides go in `"args"`.
- `common_args`: dictionary of CLI flags shared by all runs. Booleans translate to
  presence/absence of a flag, numbers/strings become `"--flag value"` pairs.
- `variants`: each entry provides a `name` (sub-folder) and `args` overriding or extending
  the shared arguments. Variants are executed for every listed graph.

## Included experiment templates

- `buffer_sweep.json`: varying buffer sizes, fixed batch size.
- `batch_sweep.json`: varying batch sizes, fixed buffer size.
- `scoring_compare.json`: compare PQ scoring functions.
- `ghost_toggle.json`: toggle the ghost-neighbor pipeline.

Create additional configs by copying one of these templates and editing the
`variants` list (or the `graphs` array) as needed.

## Analyzing results

`tools/analyze_experiments.py` scans run directories, parses the FlatBuffer logs, and
prints run-level CSV rows plus per-variant averages:

```bash
python tools/analyze_experiments.py experiments/results/buffer_sweep \
       --output-csv buffer_sweep_summary.csv
```

The analyzer relies on the metadata recorded by `run_experiments.py` and uses the
FlatBuffer reader from `tools/parse_partition_log.py`. It prints both aggregated averages
over all `k` values and a second table broken down per `k`, which makes it easy to
compare scalability across partition counts.

## Running with GNU parallel and memory limits

To fan out many runs in parallel (with a per-process RAM cap), use the planning mode plus
the helper script:

```bash
# Plan tasks and execute them with GNU parallel (4 jobs, 8GB per run)
tools/run_experiments_parallel.sh \
    experiments/configs/buffer_sweep.json \
    deploy/buffcut \
    4 \
    8192 \
    /tmp/buffcut_tasks
```

Internally this:
1. Calls `python tools/run_experiments.py --plan-output-dir ... --plan-only` to create
   per-run task JSON files describing the command, run directory, and metadata paths.
2. Invokes GNU parallel so each task is executed via
   `python tools/run_experiments.py --run-task task.json --memory-limit-mb <limit>`.
   The Python runner applies an OS-level address-space limit (RLIMIT_AS), reports if a
   run is killed by SIGKILL due to overuse, and still writes the standard logs/metadata.

You can also create the task directory manually and inspect the JSON files before
launching GNU parallel if you need to customize the invocation further.

# BuffCut


**BuffCut** is a **buffered streaming graph partitioner** developed as part of my master’s thesis.
It builds on the [HeiStream](https://github.com/KaHIP/HeiStream) framework and leverages **prioritized buffering** and **multilevel refinement** to achieve more **robust partition quality**, especially under **adversarial node orderings**.

---

## 🧩 Concept

BuffCut processes the input graph as a stream and keeps a configurable **buffer** of recently seen nodes.
Nodes are scored and selected for partitioning in small **batches** using a multilevel partitioning strategy.

- Increasing the **buffer size** and **batch size** improves partition quality.
- Larger values also increase **memory consumption**.
- Typically, a ratio of about **1:8** to **1:16** (batch size : buffer size) works well.
- The **ghost neighbor** mechanism can be enabled to improve partition quality.

---

## 🚀 Minimal Example

### Building BuffCut

To build the software, run
```shell
./compile.sh
```

After building the project, two executables are available in `deploy/`:
- `buffcut` – serial version
- `par_buffcut` – parallel version (faster, slightly higher memory usage)



### Basic Usage
To partition a graph in METIS format using the standard parameters of BuffCut, run

```bash
./par_buffcut <graph_file> --k=4
```


Adjusting buffer and batch size (serial):
```shell
./par_buffcut <graph file> --k=4 --batch_size=16384 --buffer_size=131072
```

Ghost neighbors enabled:
```shell
./par_buffcut <graph file> --k=4 --batch_size=16384 --buffer_size=131072 --enable_ghost
```


For a complete list of parameters alongside with descriptions, run:

```shell
./par_buffcut --help
```

For a description of the graph format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).


Note:
- The program stores the results of the executed command in a [flatbuffer](https://github.com/google/flatbuffers) `.bin`
  file identified by `graph_k_batchsize_buffersize.bin` if you pass the flag `write_log`.

---

## 📊 Experiment Automation

The `experiments/` folder now contains ready-to-use JSON templates that describe common
parameter sweeps (buffer size, batch size, scoring variants, and ghost-neighbor toggles).
Use `tools/run_experiments.py` to execute a config and automatically organize the FlatBuffer
logs per graph and variant:

```bash
python tools/run_experiments.py --config experiments/configs/buffer_sweep.json
```

Each run records `run_meta.json` (arguments, duration, generated `.bin` files) and `run.log`
alongside the FlatBuffer output under `experiments/results/<experiment>/<graph>/k_<k>/<variant>/`.
Configs may include `k_values: [4,8,16,32,64,128,256]` (or any list) to sweep over multiple
partition counts automatically; per-graph/variant runs are spawned for every listed `k`.
Point `tools/analyze_experiments.py` at any of these folders to obtain CSV summaries and
per-variant averages of quality, locality, runtime, and memory consumption, including a
breakdown by `k`.

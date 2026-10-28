# Benchmark results

Results reported in the project report, produced with `benchmark/lisp_benchmark.ipynb` on a Google Colab CPU runtime (Intel Xeon CPU @ 2.20 GHz, 2 cores, 12 GB RAM, Ubuntu 24.04.5 LTS, Python 3.13.16, R 4.6.1, localsp 0.2.0, gdverse 1.7, sdsfun 0.9.0).

| File | Content |
| --- | --- |
| `raw_runs.csv` | One row per run: language, size, runtime, memory |
| `summary_table.csv` | Median runtime and memory per language and size |
| `correctness.csv` | Python vs R output comparison at each size |
| `runtime_vs_size.png`, `memory_vs_size.png` | Plots used in the report |
| `subset_<n>.csv` | Input subsets (seed 42) |
| `py_<n>.csv`, `r_<n>.csv` | LISP output of each implementation at each size |

"""Summarise the benchmark: results table, correctness check per size, and plots.

Usage (from the repository root):
    python benchmark/analyse_benchmarks.py [results_folder]      (default: benchmark/results)
Reads <results_folder>/raw_runs.csv and writes into the same folder:
    summary_table.csv, correctness.csv, runtime_vs_size.png, memory_vs_size.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
raw = pd.read_csv(os.path.join(out, "raw_runs.csv"))
raw["lisp_mb"] = raw["peak_mb"] - raw["rss_before_mb"]

# ---- 1. Summary table: median over repeats --------------------------------
summary = (raw.groupby(["lang", "n"])
              .agg(runs=("seconds", "size"),
                   median_seconds=("seconds", "median"),
                   min_seconds=("seconds", "min"),
                   max_seconds=("seconds", "max"),
                   median_peak_mb=("peak_mb", "median"),
                   median_lisp_mb=("lisp_mb", "median"))
              .reset_index())
wide = summary.pivot(index="n", columns="lang", values="median_seconds")
if {"Python", "R"} <= set(wide.columns):
    summary = summary.merge((wide["Python"] / wide["R"]).rename("py_over_r_time").reset_index(), on="n")
summary.round(3).to_csv(os.path.join(out, "summary_table.csv"), index=False)
print("=== Median over repeats ===")
print(summary.round(3).to_string(index=False))

# ---- 2. Correctness: does Python match R at every size? -------------------
rows = []
for n in sorted(raw["n"].unique()):
    rp, pp = os.path.join(out, f"r_{n}.csv"), os.path.join(out, f"py_{n}.csv")
    if not (os.path.exists(rp) and os.path.exists(pp)):
        continue
    r = pd.read_csv(rp).sort_values("rid").reset_index(drop=True)
    p = pd.read_csv(pp).sort_values("rid").reset_index(drop=True)
    num = [c for c in r.columns if c != "rid" and pd.api.types.is_numeric_dtype(r[c])]
    lab = [c for c in r.columns if c.endswith("_Interaction")]
    rows.append({"n": n,
                 "same_columns": set(r.columns) == set(p.columns),
                 "numeric_cols": len(num),
                 "max_abs_diff": (r[num] - p[num]).abs().max().max(),
                 "label_cols": len(lab),
                 "label_mismatches": int(sum((r[c].astype(str) != p[c].astype(str)).sum() for c in lab))})
if rows:
    corr = pd.DataFrame(rows)
    corr.to_csv(os.path.join(out, "correctness.csv"), index=False)
    print("\n=== Python vs R output at each size ===")
    print(corr.to_string(index=False))

# ---- 3. Plots --------------------------------------------------------------
for col, ylabel, fname in [("median_seconds", "Runtime (s, median)", "runtime_vs_size.png"),
                           ("median_lisp_mb", "Memory used by LISP run (MB, median)", "memory_vs_size.png")]:
    fig, ax = plt.subplots(figsize=(6, 4))
    for lang, marker in [("R", "o"), ("Python", "s")]:
        d = summary[summary["lang"] == lang]
        if len(d):
            ax.plot(d["n"], d[col], marker=marker, label=lang)
    ax.set_xlabel("Number of observations")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname), dpi=200)
print(f"\nSaved table, correctness check and plots to {out}")

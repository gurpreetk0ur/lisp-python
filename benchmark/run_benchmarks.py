"""Run the R vs Python LISP benchmark.

Usage (from the repository root):
    python benchmark/run_benchmarks.py [--langs py,r] [--repeats 3] [--threshold 8] [--out benchmark/results]

Creates in the output folder:
    subset_<n>.csv            nested random subsets (seed 42) used by BOTH languages
    py_<n>.csv / r_<n>.csv    LISP outputs (from repeat 1), for the correctness check
    raw_runs.csv              one row per run (appended after every run, so nothing is lost)
Each run is a separate process, so memory and start-up do not carry over between runs.
The R runs need R with the localsp package installed (see README).
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SIZES = [100, 200, 400, 600, 908]

ap = argparse.ArgumentParser()
ap.add_argument("--langs", default="py,r")
ap.add_argument("--repeats", type=int, default=3)
ap.add_argument("--threshold", type=float, default=8.0)
ap.add_argument("--data", default=os.path.join(REPO, "data", "gtc.csv"))
ap.add_argument("--out", default=os.path.join(HERE, "results"))
a = ap.parse_args()

out = os.path.abspath(a.out)
os.makedirs(out, exist_ok=True)

# 1. Nested subsets: the first n rows of one fixed random permutation (seed 42),
#    kept in original file order. The same files are read by R and Python.
gtc = pd.read_csv(a.data)
perm = np.random.default_rng(42).permutation(len(gtc))
for n in SIZES:
    gtc.iloc[np.sort(perm[:n])].to_csv(os.path.join(out, f"subset_{n}.csv"), index=False)

# 2. Run every (language, size, repeat) in its own process.
cmds = {
    "py": lambda sub, res: [sys.executable, os.path.join(HERE, "bench_py.py"), sub, res, str(a.threshold)],
    "r":  lambda sub, res: ["Rscript", os.path.join(HERE, "bench_r.R"), sub, res, str(a.threshold)],
}
raw_path = os.path.join(out, "raw_runs.csv")
for lang in a.langs.split(","):
    for n in SIZES:
        for rep in range(1, a.repeats + 1):
            sub = os.path.join(out, f"subset_{n}.csv")
            res = os.path.join(out, f"{lang}_{n}.csv" if rep == 1 else f"{lang}_{n}_rep{rep}.csv")
            p = subprocess.run(cmds[lang](sub, res), capture_output=True, text=True)
            if p.returncode != 0:
                print(f"FAILED {lang} n={n} rep={rep}:\n{p.stderr[-2000:]}")
                continue
            row = json.loads(p.stdout.strip().splitlines()[-1])
            row.update(rep=rep, threshold=a.threshold)
            pd.DataFrame([row]).to_csv(raw_path, mode="a", index=False,
                                       header=not os.path.exists(raw_path))
            print(f"{row['lang']:6s} n={n:4d} rep={rep}  {row['seconds']:8.2f} s  "
                  f"peak {row['peak_mb']:7.1f} MB  (LISP used {row['peak_mb'] - row['rss_before_mb']:.1f} MB)")
            if rep > 1:
                os.remove(res)  # only repeat 1's output is kept for the correctness check
print(f"\nDone. Raw results in {raw_path}")

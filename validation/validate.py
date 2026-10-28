"""Validate the Python LISP implementation against output from the R package.

Runs lisp5.lisp() on data/gtc_subset_100.csv and compares the result with
data/r_ground_truth_100.csv (output of localsp::lisp() on the same rows and
settings). Numeric columns are compared by absolute difference and interaction
type columns by exact text match.

Usage (from the repository root):
    python validation/validate.py
    python validation/validate.py --r-output path/to/other_r_output.csv

Exits with code 0 if every numeric difference is below --tolerance and every
interaction label matches, and 1 otherwise.
"""
import argparse
import os
import sys
import time

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from lisp5 import lisp  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default=os.path.join(REPO, "data", "gtc_subset_100.csv"))
ap.add_argument("--r-output", default=os.path.join(REPO, "data", "r_ground_truth_100.csv"))
ap.add_argument("--threshold", type=float, default=8.0)
ap.add_argument("--tolerance", type=float, default=1e-6)
ap.add_argument("--save", default=None, help="optional path to save the Python output as CSV")
a = ap.parse_args()

# Settings must match those used to create the R output
df = pd.read_csv(a.data)
t0 = time.perf_counter()
py = lisp("GTC ~ .", data=df.drop(columns=["X", "Y"]), coords=df[["X", "Y"]].values,
          threshold=a.threshold, discnum=range(3, 6), discmethod=["quantile"])
seconds = time.perf_counter() - t0
if a.save:
    py.to_csv(a.save, index=False)

r = pd.read_csv(a.r_output).sort_values("rid").reset_index(drop=True)
py = py.sort_values("rid").reset_index(drop=True)

same_cols = set(r.columns) == set(py.columns)
num = [c for c in r.columns if c != "rid" and pd.api.types.is_numeric_dtype(r[c])]
lab = [c for c in r.columns if c.endswith("_Interaction")]
max_diff = (r[num] - py[num]).abs().max().max()
mismatches = int(sum((r[c].astype(str) != py[c].astype(str)).sum() for c in lab))

print(f"Rows: {len(py)}   Python runtime: {seconds:.1f} s")
print(f"Same columns as R:            {same_cols} ({len(py.columns)} columns)")
print(f"Numeric columns compared:     {len(num)}")
print(f"Maximum absolute difference:  {max_diff:.3e}")
print(f"Interaction label columns:    {len(lab)}")
print(f"Interaction label mismatches: {mismatches} of {len(lab) * len(r)}")

ok = same_cols and max_diff < a.tolerance and mismatches == 0
print("RESULT:", "PASS, Python output matches R" if ok else "FAIL, Python output differs from R")
sys.exit(0 if ok else 1)

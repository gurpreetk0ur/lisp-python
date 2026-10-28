"""Benchmark one run of the Python LISP (lisp5.py).

Usage: python benchmark/bench_py.py <subset.csv> <out_result.csv> <threshold>
Prints one line of JSON: {"lang", "n", "seconds", "rss_before_mb", "peak_mb"}

Memory is read from /proc/self/status (Linux), the same way bench_r.R does it:
  rss_before_mb = memory in use just before the timed block (interpreter + libraries + data)
  peak_mb       = highest memory the process ever reached (VmHWM)
so peak_mb - rss_before_mb is the extra memory used by the LISP run itself.
"""
import json
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lisp5 import lisp  # noqa: E402


def mem_mb(field):
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith(field + ":"):
                return int(line.split()[1]) / 1024  # kB -> MB
    return float("nan")


subset_csv, out_csv, threshold = sys.argv[1], sys.argv[2], float(sys.argv[3])

df = pd.read_csv(subset_csv)
coords = df[["X", "Y"]].values
data = df.drop(columns=["X", "Y"]).reset_index(drop=True)

rss_before = mem_mb("VmRSS")
t0 = time.perf_counter()
# lisp5 builds the distance matrix inside lisp(), so it is included in the timing
res = lisp("GTC ~ .", data=data, coords=coords, threshold=threshold,
           discnum=range(3, 6), discmethod=["quantile"])
seconds = time.perf_counter() - t0
peak = mem_mb("VmHWM")

res.to_csv(out_csv, index=False)
print(json.dumps({"lang": "Python", "n": len(df), "seconds": seconds,
                  "rss_before_mb": rss_before, "peak_mb": peak}))

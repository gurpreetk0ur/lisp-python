"""Minimal working example: run LISP on the 100-row glacier subset (a few seconds).

Usage (from the repository root):
    python example.py
"""
import pandas as pd

from lisp5 import lisp

gtc = pd.read_csv("data/gtc_subset_100.csv")
coords = gtc[["X", "Y"]].values                         # longitude, latitude
data = gtc.drop(columns=["X", "Y"])                    # GTC + 8 explanatory variables

result = lisp("GTC ~ .", data=data, coords=coords, threshold=8,
              discnum=range(3, 6), discmethod=["quantile"])

print(result.shape)                                    # (100, 101): one row per location
print(result[["rid", "pd_Pre", "sig_Pre", "Pre_Slope_Interaction", "Pre_Slope_pd"]].head().round(4))

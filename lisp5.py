"""
Local Indicator of Stratified Power (LISP): pure Python implementation.
No R / rpy2 dependency.

Copyright (C) 2026 Gurpreet Kour
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, version 3. It is distributed WITHOUT ANY WARRANTY. See the LICENSE
file in this repository or <https://www.gnu.org/licenses/gpl-3.0.html>.

ATTRIBUTION
This file is a Python translation (modified work, 2026) of code from the
following GPL-3 licensed R packages. Copyright in the original code belongs to
their authors:
    localsp  (J. Hu, W. Lv, Y. Song)        https://github.com/ausgis/localsp
        R/lisp.R           -> lisp(), _calcul_localq()
    gdverse  (W. Lv et al.)                 https://github.com/stscl/gdverse
        R/opgd.R           -> optimal_disc() (per-variable optimisation)
        R/geodetector.R    -> q_statistic() (factor_detector),
                              _classify_interaction() (interaction_detector)
    sdsfun   (W. Lv)                        https://github.com/stscl/sdsfun
        src/Disc.cpp       -> discretise() (_disc_equal/_disc_sd/
                              _disc_quantile/_disc_geometric/_disc_natural)
The method is described in: J. Hu et al., "A local indicator of stratified
power", Int. J. Geogr. Inf. Sci., 2024, doi:10.1080/13658816.2024.2437811.

Each function below is commented with the R/C++ source it mirrors. The package's
published vignette output predates the interaction detector, so validation was
done against output generated from the current R source instead (see README).

Install:
    pip install numpy pandas scipy
(neighbourhoods are computed directly from a full pairwise distance matrix,
mirroring R's `dm[rowindice,] <= bw[rowindice]` exactly; no spatial-weights
library is needed.)

Usage:
    import pandas as pd
    from lisp5 import lisp

    gtc = pd.read_csv("data/gtc_subset_100.csv")
    coords = gtc[["X", "Y"]].values
    data   = gtc.drop(columns=["X", "Y"]).reset_index(drop=True)

    result = lisp("GTC ~ .", data=data, coords=coords, threshold=8,
                  discnum=range(3, 6), discmethod=["quantile"])
    print(result)
"""

from __future__ import annotations

import warnings
import itertools
from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import ncf


# ── Discretisation ────────────────────────────────────────────────────────────
# Ported 1:1 from stscl/sdsfun src/Disc.cpp (sdDisc, equalDisc, quantileDisc,
# geometricDisc, ArmaJenksBreaks + naturalDisc). These are NOT the "textbook"
# versions of these binning schemes -- they replicate this package's specific
# (sometimes idiosyncratic) edge handling, since that's what the R pipeline
# actually calls via sdsfun::discretize_vector().

def _disc_equal(x: np.ndarray, k: int) -> np.ndarray:
    """Mirrors equalDisc(): index = ceil((x - min) / interval), clipped to [1,k]."""
    min_x, max_x = x.min(), x.max()
    interval = (max_x - min_x) / k
    if interval == 0:
        return np.ones(len(x), dtype=int)
    idx = np.ceil((x - min_x) / interval).astype(int)
    return np.clip(idx, 1, k)


def _disc_sd(x: np.ndarray, k: int) -> np.ndarray:
    """
    Mirrors sdDisc():
        diff  = (x - mean) / std        # std uses ddof=1, matching arma::stddev default
        index = floor((diff + n/2) / n * n)   # algebraically == floor(diff + n/2)
        clipped to [1, n]
    NOTE: n/2 uses true (float) division even for odd k, exactly as the C++ does
    (n is a double in the C++ signature).
    """
    mean_x = x.mean()
    std_x = x.std(ddof=1)
    if std_x == 0:
        return np.ones(len(x), dtype=int)
    diff = (x - mean_x) / std_x
    idx = np.floor(diff + k / 2.0).astype(int)
    return np.clip(idx, 1, k)


def _disc_quantile(x: np.ndarray, k: int) -> np.ndarray:
    """
    Mirrors quantileDisc():
        sorted_x = sort(x)
        quantiles[i] = sorted_x[ floor(i * (N-1) / n) ]   for i = 0..n
        (this is a *truncating* nearest-rank quantile, NOT linear-interpolated
        like np.percentile's default -- important for small local windows.)
        label = first j (1-indexed) such that x <= quantiles[j]
    """
    xs = np.sort(x)
    n_elem = len(xs)
    positions = [int((i * (n_elem - 1)) / k) for i in range(k + 1)]
    quantiles = xs[positions]

    idx = np.empty(len(x), dtype=int)
    for i, val in enumerate(x):
        label = k  # fallback: matches C++'s implicit "else last class" behaviour
        for j in range(k):
            if val <= quantiles[j + 1]:
                label = j + 1
                break
        idx[i] = label
    return idx


def _disc_geometric(x: np.ndarray, k: int) -> np.ndarray:
    """
    Mirrors geometricDisc():
        factor = (max/min) ** (1/n)
        index  = floor(log(x/min) / log(factor)) + 1, clipped to [1, n]
    NOTE: the R/C++ version assumes x > 0 and does NOT shift negative/zero
    values -- ported faithfully (will emit NaN/inf for non-positive x, same
    as the original would error/produce garbage in R).
    """
    min_x, max_x = x.min(), x.max()
    factor = (max_x / min_x) ** (1.0 / k)
    with np.errstate(divide="ignore", invalid="ignore"):
        idx = np.floor(np.log(x / min_x) / np.log(factor)).astype(int) + 1
    idx = np.nan_to_num(idx, nan=1, posinf=k, neginf=1)
    return np.clip(idx, 1, k)


def _jenks_breaks(xs_sorted: np.ndarray, k: int) -> np.ndarray:
    """
    Mirrors ArmaJenksBreaks(): exact Fisher-Jenks dynamic program.
    Returns the k-1 INTERIOR breakpoints (not including min/max), matching
    the C++ function's return shape exactly.
    """
    n = len(xs_sorted)
    # 1-indexed via padded arrays (index 0 row/col unused, mirrors C++ arma indices)
    lower = np.zeros((n + 1, k + 1), dtype=int)
    var_comb = np.full((n + 1, k + 1), np.inf)
    for i in range(1, k + 1):
        lower[1, i] = 1
        var_comb[1, i] = 0.0

    for l in range(2, n + 1):
        s1 = s2 = 0.0
        variance = 0.0
        for m in range(1, l + 1):
            lower_class_limit = l - m + 1
            val = xs_sorted[lower_class_limit - 1]
            s1 += val
            s2 += val * val
            variance = s2 - (s1 * s1) / m
            i4 = lower_class_limit - 1
            if i4 != 0:
                for j in range(2, k + 1):
                    cand = variance + var_comb[i4, j - 1]
                    if var_comb[l, j] >= cand:
                        lower[l, j] = lower_class_limit
                        var_comb[l, j] = cand
        lower[l, 1] = 1
        var_comb[l, 1] = variance

    kclass = np.zeros(k - 1)
    kk = n
    for count_num in range(k - 1, 0, -1):
        start = lower[kk, count_num + 1]
        kclass[count_num - 1] = xs_sorted[start - 2]  # 1-indexed start -> 0-indexed (start-2)
        kk = start - 1
    return kclass


def _disc_natural(x: np.ndarray, k: int) -> np.ndarray:
    """
    Mirrors naturalDisc(): exact Jenks (no sampling -- local windows here are
    always << 3000 points, so the C++ sampling branch never triggers).
    Assignment rule ported exactly: first j with x[i] < breaks[j] wins;
    otherwise (x[i] >= last break) falls into the top class.
    """
    if k < 2:
        return np.ones(len(x), dtype=int)
    xs_sorted = np.sort(x)
    breaks = _jenks_breaks(xs_sorted, k)
    bn = len(breaks)

    idx = np.zeros(len(x), dtype=int)
    for i, val in enumerate(x):
        label = 0
        for j in range(bn):
            if val < breaks[j]:
                label = j + 1
                break
        if bn == 0 or val >= breaks[bn - 1]:
            label = k
        idx[i] = label
    return idx


_DISC_METHODS = {
    "equal":     _disc_equal,
    "quantile":  _disc_quantile,
    "sd":        _disc_sd,
    "geometric": _disc_geometric,
    "natural":   _disc_natural,
}


def discretise(x: np.ndarray, k: int, method: str) -> np.ndarray:
    """Discretise continuous vector x into k classes using the given method."""
    if method not in _DISC_METHODS:
        raise ValueError(f"Unknown discretisation method '{method}'. "
                          f"Choose from {list(_DISC_METHODS)}")
    return _DISC_METHODS[method](np.asarray(x, dtype=float), k).astype(int)


# ── Q-statistic & significance test ──────────────────────────────────────────
# Ported 1:1 from stscl/gdverse R/geodetector.R :: factor_detector()

def q_statistic(y: np.ndarray, strata: np.ndarray) -> tuple[float, float]:
    """
    Faithful port of gdverse::factor_detector().

    R reference:
        gdf = tibble(x, y) |> group_by(x) |> filter(n() > 1) |> ungroup()
        x = gdf$x; y = gdf$y                      # <- y is REASSIGNED to the
                                                   #    filtered subset here
        rss = \\(y) (length(y) - 1) * var(y)       # pure sum-of-squares (SS)
        qv  = 1 - sum(tapply(y, x, rss)) / rss(y)  # SS_within / SS_total, both
                                                   # computed on the FILTERED y
        N = length(y); L = length(levels(x))
        Fv = ((N-L)*qv) / ((L-1)*(1-qv))
        lambda = (sum(hmean^2) - (sum(sqrt(Nh)*hmean))^2/N) / (var(y)*(N-1)/N)
                                                   # denominator is POPULATION
                                                   # variance of the filtered y
        pv = pf(Fv, L-1, N-L, ncp=lambda, lower.tail=FALSE)

    Returns (q, p_value). Returns (nan, nan) if computation is impossible.
    """
    y = np.asarray(y, dtype=float)
    strata = np.asarray(strata)

    labels, counts = np.unique(strata, return_counts=True)
    valid = labels[counts > 1]                       # drop singleton strata
    mask = np.isin(strata, valid)
    yf, sf = y[mask], strata[mask]                    # R reassigns x,y to this subset

    N, L = len(yf), len(valid)
    if N < 2 or L < 2:
        return np.nan, np.nan

    def ss(v: np.ndarray) -> float:
        """Pure sum-of-squares -- dof-invariant, matches R's rss()."""
        return float((len(v) - 1) * np.var(v, ddof=1)) if len(v) > 1 else 0.0

    total_ss = ss(yf)
    if total_ss == 0:
        return np.nan, np.nan

    within_ss = sum(ss(yf[sf == lab]) for lab in valid)
    q = float(np.clip(1.0 - within_ss / total_ss, 0.0, 1.0))

    df1, df2 = L - 1, N - L
    if df1 < 1 or df2 < 1 or (1 - q) == 0:
        return q, np.nan
    f_stat = (df2 * q) / (df1 * (1 - q))

    hmean = {lab: float(np.mean(yf[sf == lab])) for lab in valid}
    Nh = {lab: int(np.sum(sf == lab)) for lab in valid}
    v1 = sum(hmean[lab] ** 2 for lab in valid)
    v2 = (sum(np.sqrt(Nh[lab]) * hmean[lab] for lab in valid)) ** 2 / N
    pop_var = np.var(yf, ddof=1) * (N - 1) / N        # population variance, filtered y
    if pop_var == 0:
        return q, np.nan
    ncp = (v1 - v2) / pop_var

    p_val = 1.0 - ncf.cdf(f_stat, df1, df2, ncp)
    return q, float(p_val)


# ── Optimal discretisation (mirrors gdverse::gd_optunidisc per-variable logic) ─

def optimal_disc(
    y: np.ndarray,
    x: np.ndarray,
    discnum: list[int],
    discmethod: list[str],
) -> np.ndarray:
    """
    Mirrors gd_optunidisc(): for a single variable, try every (k, method) combo
    in the crossing grid and keep the strata that maximise the factor-detector Q.

    Tie-breaking note: R builds the grid with tidyr::crossing(x, k, method),
    which SORTS k ascending and method ALPHABETICALLY, then dplyr::slice_max(
    ..., with_ties = FALSE) keeps the first row encountered for ties. We
    replicate that iteration order (k ascending, method alphabetical) and use
    strict '>' so the first-encountered maximum wins, exactly like R's tie
    behaviour.
    """
    best_q = -np.inf
    best_strata = None

    for k in sorted(discnum):
        for method in sorted(discmethod):
            try:
                strata = discretise(x, k, method)
            except Exception:
                continue
            q, _ = q_statistic(y, strata)
            if np.isfinite(q) and q > best_q:
                best_q = q
                best_strata = strata

    return best_strata


# ── Interaction type classification (Wang et al. 2016) ───────────────────────
# Ported 1:1 from stscl/gdverse R/geodetector.R :: interaction_detector()

def _classify_interaction(q1: float, q2: float, q12: float) -> str:
    """
    Faithful port of the R if/elif chain (order and comparisons matter):
        if      qv12 <  min(q1,q2):                    "Weaken, nonlinear"
        else if qv12 >= min(q1,q2) & qv12 <= max(q1,q2): "Weaken, uni-"
        else if qv12 >  max(q1,q2) & qv12 <  q1+q2:      "Enhance, bi-"
        else if qv12 == q1+q2:                           "Independent"
        else:                                             "Enhance, nonlinear"
    """
    if any(not np.isfinite(v) for v in (q1, q2, q12)):
        return "NA"
    mn, mx = min(q1, q2), max(q1, q2)
    total = q1 + q2
    if q12 < mn:
        return "Weaken, nonlinear"
    elif q12 <= mx:
        return "Weaken, uni-"
    elif q12 < total:
        return "Enhance, bi-"
    elif q12 == total:
        return "Independent"
    else:
        return "Enhance, nonlinear"


# ── Per-row local computation ─────────────────────────────────────────────────

def _calcul_localq(
    rowindice: int,
    dep_var: str,
    indep_vars: list[str],
    df: pd.DataFrame,
    neighbour_idx: list[int],
    discvar: Optional[list[str]],
    discnum: list[int],
    discmethod: list[str],
) -> pd.DataFrame:
    """
    Compute local interaction Q-statistics for one focal observation.
    Mirrors calcul_localq() in R/lisp.R: discretise each variable ONCE
    (optimally), then reuse that same discretisation for BOTH the factor
    detector and every pairwise interaction detector.
    Returns a single-row DataFrame; NaN row on failure.
    """
    rid = rowindice + 1  # 1-based to match R output

    localdf = df.iloc[neighbour_idx].reset_index(drop=True)
    if len(localdf) < 3:
        warnings.warn(f"Row {rid}: neighbourhood too small ({len(localdf)}), skipping.",
                      RuntimeWarning, stacklevel=2)
        return pd.DataFrame({"rid": [rid]})

    y = localdf[dep_var].values

    to_disc = discvar if discvar is not None else indep_vars

    strata_map: dict[str, np.ndarray] = {}
    for var in indep_vars:
        if var in to_disc:
            st = optimal_disc(y, localdf[var].values, discnum, discmethod)
        else:
            st = localdf[var].values.astype(int)
        if st is None:
            st = np.zeros(len(localdf), dtype=int)
        strata_map[var] = st

    factor_results: dict[str, tuple[float, float]] = {}
    for var in indep_vars:
        factor_results[var] = q_statistic(y, strata_map[var])

    # R's opgd() internally routes variables through gd_optunidisc(), which
    # uses dplyr::group_by(x) %>% slice_max() %>% ungroup() -- dplyr's default
    # group_by SORTS GROUPS ALPHABETICALLY. So by the time combn() builds
    # interaction pairs, the variable list is already alphabetically ordered,
    # and every "var1_var2" column name follows alphabetical order, not the
    # original CSV column order. Match that here so column names line up.
    pairs = list(itertools.combinations(sorted(indep_vars), 2))
    interaction_results = {}
    for v1, v2 in pairs:
        s1, s2 = strata_map[v1], strata_map[v2]
        # Combine as a joint category, mirroring R's paste0(x1,'_',x2) exactly.
        # NOTE: an earlier version used `s1 * (s1.max()+1) + s2`, which silently
        # collides whenever s2's discretisation has more classes than s1's --
        # very possible since each variable is optimally discretised
        # independently. String-concat + factorize is collision-proof and
        # matches R's actual mechanism (paste, not arithmetic).
        combined_str = np.char.add(np.char.add(s1.astype(str), "_"), s2.astype(str))
        _, combined = np.unique(combined_str, return_inverse=True)
        q12, p12 = q_statistic(y, combined)
        q1, _ = factor_results[v1]
        q2, _ = factor_results[v2]
        interaction_type = _classify_interaction(q1, q2, q12)
        interaction_results[(v1, v2)] = (interaction_type, q12, p12)

    row: dict = {}
    for var in indep_vars:
        q, p = factor_results[var]
        row[f"pd_{var}"] = q
        row[f"sig_{var}"] = p

    for (v1, v2), (itype, q12, p12) in interaction_results.items():
        row[f"{v1}_{v2}_Interaction"] = itype
        row[f"{v1}_{v2}_pd"] = q12
        row[f"{v1}_{v2}_sig"] = p12

    row["rid"] = rid
    return pd.DataFrame([row])


# ── Formula parsing ───────────────────────────────────────────────────────────

def _parse_formula(formula: str, columns: list[str]) -> tuple[str, list[str]]:
    """Parse 'Y ~ .' or 'Y ~ X1 + X2' into (dep_var, [indep_vars])."""
    lhs, rhs = formula.split("~")
    dep = lhs.strip()
    rhs = rhs.strip()
    if rhs == ".":
        indep = [c for c in columns if c != dep]
    else:
        indep = [v.strip() for v in rhs.split("+")]
    return dep, indep


# ── Public API ────────────────────────────────────────────────────────────────

def lisp(
    formula: str,
    data: pd.DataFrame,
    coords: np.ndarray,
    threshold: Union[float, int],
    discvar: Optional[list] = None,
    discnum: Union[list, range] = range(3, 9),
    discmethod: Union[str, list] = ("sd", "equal", "geometric", "quantile", "natural"),
) -> pd.DataFrame:
    """
    Local Indicator of Stratified Power -- pure Python implementation.

    Parameters
    ----------
    formula : str
        Formula string, e.g. "GTC ~ ." or "GTC ~ X1 + X2".
    data : pd.DataFrame
        Observation data (no geometry/coordinate columns).
    coords : np.ndarray
        (n, 2) array of X, Y coordinates for neighbourhood construction.
    threshold : float or int
        Distance threshold for neighbourhood selection (mirrors R's
        `dm[rowindice,] <= bw[rowindice]`, i.e. a closed <= threshold).
    discvar : list of str, optional
        Variables to discretise. None -> all independent variables.
    discnum : list or range of int
        Candidate class counts. Default range(3, 9) (i.e. 3:8 in R).
    discmethod : str or list of str
        Discretisation method(s). Default: all five methods.

    Returns
    -------
    pd.DataFrame
        One row per observation, sorted by rid, in wide format.
    """
    data = data.reset_index(drop=True)
    discnum = list(discnum)
    discmethod = [discmethod] if isinstance(discmethod, str) else list(discmethod)

    dep_var, indep_vars = _parse_formula(formula, list(data.columns))

    coords = np.asarray(coords, dtype=float)
    dist_mat = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))

    results = []
    for i in range(len(data)):
        neighbours = list(np.where(dist_mat[i] <= threshold)[0])  # includes self (dist=0)
        row_result = _calcul_localq(
            rowindice=i,
            dep_var=dep_var,
            indep_vars=indep_vars,
            df=data,
            neighbour_idx=neighbours,
            discvar=discvar,
            discnum=discnum,
            discmethod=discmethod,
        )
        results.append(row_result)

    out = pd.concat(results, axis=0, ignore_index=True, join="outer")
    out = out.sort_values("rid").reset_index(drop=True)

    cols = ["rid"] + [c for c in out.columns if c != "rid"]
    return out[cols]

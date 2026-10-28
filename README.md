# lisp-python

A pure Python implementation of the **Local Indicator of Stratified Power (LISP)** [1], translated from the original R package [localsp](https://github.com/ausgis/localsp). It needs no R installation.

The implementation was validated against the R package: on a glacier dataset of up to 908 locations, all Q-statistics and p-values agreed with R to within 1 × 10⁻⁹ and all interaction types were identical. With the validated settings it ran 13 to 21 times faster than R.

This code accompanies the Honours report *Python Implementation and Validation of the Local Indicator of Spatial Processes (LISP)* (G. Kour, Curtin University, 2026).

## Attribution and licence

`lisp5.py` is a translation of code from three R packages, all licensed under GPL-3. Copyright in the original code belongs to their authors:

| Package | Authors | Translated code |
| --- | --- | --- |
| [localsp](https://github.com/ausgis/localsp) | J. Hu, W. Lv, Y. Song | `R/lisp.R` |
| [gdverse](https://github.com/stscl/gdverse) | W. Lv, Y. Lei, F. Liu, J. Yan, Y. Song, W. Zhao | `R/opgd.R`, `R/geodetector.R` |
| [sdsfun](https://github.com/stscl/sdsfun) | W. Lv | `src/Disc.cpp` |

As a derivative work, this repository is released under the **GNU General Public License v3.0** (see `LICENSE`). `data/gtc.csv` is copied from localsp. Each function in `lisp5.py` names the R or C++ function it translates.

## Repository contents

| Path | Content |
| --- | --- |
| `lisp5.py` | The Python implementation (`lisp()` and its components) |
| `example.py` | Minimal working example |
| `requirements.txt` | Python package versions used |
| `data/` | Glacier dataset, 100-row validation subset, R reference output and the row numbers of every sampled subset (see `data/README.md`) |
| `validation/validate.py` | Compares the Python output with the R reference output |
| `validation/lisp_python_replication.ipynb` | Development notebook: Python version of the localsp vignette, and the validation |
| `benchmark/` | Scripts and Colab notebook comparing R and Python runtime and memory |
| `benchmark/results/` | Benchmark results reported in the report |

## Requirements

- Python 3.13 (tested with 3.13.16). `lisp5.py` needs only numpy, pandas and scipy.
- Install the versions used for the results with:

  ```
  pip install -r requirements.txt
  ```

- R is needed only to regenerate the R reference output or run the R side of the benchmark. Versions used: R 4.6.1, localsp 0.2.0 (commit 532f618), gdverse 1.7 (commit 2d97584) and sdsfun 0.9.0 (commit 35b4521), installed from GitHub. Installation steps are in `benchmark/lisp_benchmark.ipynb`.

## Data

`data/gtc.csv` contains 908 glacier locations in High Mountain Asia, with glacier thickness change (GTC) as the response and eight explanatory variables. It is copied from the localsp package. It can also be obtained in R with `system.file("extdata/gtc.csv", package = "localsp")` or from the [localsp repository](https://github.com/ausgis/localsp/blob/main/inst/extdata/gtc.csv). Details of all data files are in `data/README.md`.

## Working example

From the repository root:

```python
import pandas as pd
from lisp5 import lisp

gtc = pd.read_csv("data/gtc_subset_100.csv")
coords = gtc[["X", "Y"]].values      # longitude, latitude
data = gtc.drop(columns=["X", "Y"])  # GTC + 8 explanatory variables

result = lisp("GTC ~ .", data=data, coords=coords, threshold=8,
              discnum=range(3, 6), discmethod=["quantile"])

print(result.shape)
print(result[["rid", "pd_Pre", "sig_Pre", "Pre_Slope_Interaction", "Pre_Slope_pd"]].head().round(4))
```

Output (takes about 3 seconds):

```
(100, 101)
   rid  pd_Pre  sig_Pre Pre_Slope_Interaction  Pre_Slope_pd
0    1  0.3798   0.0014    Enhance, nonlinear        0.6465
1    2  0.2763   0.0207          Enhance, bi-        0.5312
2    3  0.4180   0.0008    Enhance, nonlinear        0.7759
3    4  0.4433   0.0021          Enhance, bi-        0.7498
4    5  0.4433   0.0021          Enhance, bi-        0.7498
```

The same code is in `example.py`.

### Arguments of `lisp()`

| Argument | Meaning |
| --- | --- |
| `formula` | `"Y ~ ."` (all other columns) or `"Y ~ X1 + X2"` |
| `data` | DataFrame with the response and explanatory variables, without coordinates |
| `coords` | Array of shape (n, 2) with the point coordinates |
| `threshold` | Distance within which observations form each location's neighbourhood (same units as `coords`) |
| `discvar` | Variables to discretise (default: all explanatory variables) |
| `discnum` | Numbers of classes to try (default 3 to 8) |
| `discmethod` | Discretisation methods to try: `"sd"`, `"equal"`, `"geometric"`, `"quantile"`, `"natural"` (default: all five) |

### Output

One row per location with 101 columns for the glacier data:

- `rid`: location number, starting at 1
- `pd_<var>`, `sig_<var>`: Q-statistic and p-value for each explanatory variable
- `<a>_<b>_Interaction`, `<a>_<b>_pd`, `<a>_<b>_sig`: interaction type, Q-statistic and p-value for each pair of variables, with pairs named in alphabetical order

Column names match the output of `localsp::lisp()`.

## Reproducing the results

### Validation (Python only, about 3 seconds)

```
python validation/validate.py
```

This runs `lisp()` on `data/gtc_subset_100.csv` and compares the result with `data/r_ground_truth_100.csv`, the output of `localsp::lisp()` on the same rows with the same settings (threshold 8, 3 to 5 classes, quantile discretisation). Expected result:

```
Maximum absolute difference:  9.834e-10
Interaction label mismatches: 0 of 2800
RESULT: PASS, Python output matches R
```

To regenerate the R reference output, see `data/README.md`.

### Benchmark (R and Python)

Open `benchmark/lisp_benchmark.ipynb` in Google Colab and run the cells in order. The notebook clones this repository, installs the R packages, checks that the installed R code reproduces the reference output, validates the Python code, runs both implementations on 100, 200, 400, 600 and 908 locations, and produces the results table and plots.

The benchmark can also be run from the command line on Linux (memory is read from `/proc`):

```
python benchmark/run_benchmarks.py --langs py,r --repeats 3 --out my_results
python benchmark/analyse_benchmarks.py my_results
```

The R runs take about an hour in total. Absolute times depend on the machine.

## Results summary

Google Colab CPU runtime (Intel Xeon @ 2.20 GHz, 2 cores, 12 GB RAM), threshold 8, 3 to 5 classes, quantile discretisation, R on one core. Full results are in `benchmark/results/`.

| Locations | R runtime (s) | Python runtime (s) | Python speed-up | R extra memory (MB) | Python extra memory (MB) |
| --- | --- | --- | --- | --- | --- |
| 100 | 127.0 | 6.1 | 20.7× | 159.1 | 9.5 |
| 200 | 254.6 | 15.9 | 16.0× | 159.4 | 19.1 |
| 400 | 510.9 | 33.3 | 15.3× | 160.7 | 38.3 |
| 600 | 776.8 | 54.1 | 14.4× | 162.1 | 57.4 |
| 908 | 1,191.4 | 90.6 | 13.1× | 166.4 | 88.7 |

At every size the Python output matched the R output (maximum absolute difference below 1 × 10⁻⁹, no interaction label mismatches).

## Differences from the LISP paper

This implementation follows the R package, which differs from the method described in the paper [1]. In particular, the R package computes pairwise interactions with the OPGD interaction detector rather than the paper's GOZH multi-variable model, and the localsp vignette uses a threshold of 4 × the variogram range where the paper uses 2 ×. Variogram fitting is not part of `lisp()`.

## Limitations

- **Validated settings.** Only quantile discretisation with 3 to 5 classes was compared with R. The other four methods are translated from the sdsfun C++ source but have not been checked against R output.
- **Geometric method with negative values.** The geometric method is undefined for zero or negative values (for example aspect, summer and winter temperature in the glacier data). Such values are placed in the first class; whether R does the same was not checked.
- **Natural breaks is slow.** It is written in plain Python loops. With the default settings (3 to 8 classes, all five methods) the full dataset is estimated to take about 37 minutes.
- **Tested scale.** Tested on one dataset of up to 908 locations with one threshold. Larger datasets are untested.
- **Memory.** The full distance matrix is held in memory, which grows with the square of the number of locations.
- **Distances.** Distances are straight-line distances in the units of `coords` (degrees for the glacier data), not distances on the Earth's surface.
- **Small neighbourhoods.** A location with fewer than three observations in its neighbourhood returns an empty row with a warning, whereas the R code stops with an error.
- **No parallel processing.** Unlike localsp, `lisp()` has no option to use several cores.
- **Not a package.** The code is provided as a single module and is not installable with pip.

## How to cite

G. Kour, *lisp-python: Python implementation of the Local Indicator of Stratified Power*, version 1.0.1, 2026. https://github.com/gurpreetk0ur/lisp-python

Please also cite the original method and R packages:

[1] J. Hu et al., "A local indicator of stratified power," *International Journal of Geographical Information Science*, 2024, doi: [10.1080/13658816.2024.2437811](https://doi.org/10.1080/13658816.2024.2437811).

[2] J. Hu, W. Lv, and Y. Song, *localsp: Local Indicator of Stratified Power*, R package. https://github.com/ausgis/localsp

[3] W. Lv et al., *gdverse: Analysis of Spatial Stratified Heterogeneity*, R package. https://github.com/stscl/gdverse

[4] W. Lv, *sdsfun: Spatial Data Science Complementary Features*, R package. https://github.com/stscl/sdsfun

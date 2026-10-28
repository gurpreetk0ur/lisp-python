# Data

## gtc.csv

Glacier thickness change from 2000 to 2020 for 908 lake-terminating glaciers in the Greater Himalayas, as used in the LISP paper (Hu et al., 2024). The file is copied unchanged from the localsp R package (GPL-3), where it is distributed as `inst/extdata/gtc.csv`.

  | Column | Meaning |
  | --- | --- |
  | X, Y | Longitude and latitude (decimal degrees, EPSG:4326) |
  | GTC | Glacier thickness change, 2000 to 2020 (response variable) |
  | WinTem, SumTem, Pre, SurAlbedo | Linear trends over 2000 to 2020 in winter temperature, summer temperature, annual cumulative precipitation and summer surface albedo |
  | Elev, Slope, Aspect | Elevation, slope and aspect (NASADEM) |
  | LakeArea | Area of the proglacial lake |

Units are not reported in the LISP paper or the localsp package. The value ranges of Elev (42 to 318) and LakeArea (2,789 to 6,101) suggest these two columns may be labelled the wrong way round; this does not affect the validation, because both implementations use the same file.
To obtain it from the original source instead:

- In R, after installing localsp: `gtc = readr::read_csv(system.file("extdata/gtc.csv", package = "localsp"))`
- From GitHub: https://github.com/ausgis/localsp/blob/main/inst/extdata/gtc.csv

## gtc_subset_100.csv

A random subset of 100 rows of `gtc.csv`, used for validation. It was drawn in R following the example in the localsp documentation (`set.seed(42); gtc[sample.int(nrow(gtc), size = 100), ]`). The rows it contains are listed in `sampled_rows.csv`.

## r_ground_truth_100.csv

Output of the original R function `localsp::lisp()` on `gtc_subset_100.csv`, used as the reference for validating the Python implementation. Settings: formula `GTC ~ .`, distance threshold 8, `discnum = 3:5`, `discmethod = "quantile"`, `cores = 1`.

To regenerate it (requires R with localsp installed, see the main README):

```
Rscript benchmark/bench_r.R data/gtc_subset_100.csv r_output_100.csv 8
```

Running this with localsp 0.2.0 (commit 532f618), gdverse 1.7 (commit 2d97584) and sdsfun 0.9.0 (commit 35b4521), installed from GitHub, reproduced the file exactly (maximum difference 0).

## sampled_rows.csv

The rows of `gtc.csv` used in every subset, so the subsets can be rebuilt exactly. One row per row of `gtc.csv`:

| Column | Meaning |
| --- | --- |
| `row_id` | Row number in `gtc.csv` (1 = first data row) |
| `validation_order` | Position in `gtc_subset_100.csv` (empty if the row is not in the validation subset) |
| `benchmark_order` | Position in the random ordering used by the benchmark (`numpy.random.default_rng(42).permutation`, plus 1) |

The benchmark subset of size n contains the rows with `benchmark_order` ≤ n, kept in their original `gtc.csv` order. These rows rebuild `benchmark/results/subset_<n>.csv` exactly.

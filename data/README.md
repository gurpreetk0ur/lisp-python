# Data

## gtc.csv

Glacier thickness change data for 908 glacier locations in High Mountain Asia, as used in the LISP paper (Hu et al., 2024). The file is copied unchanged from the localsp R package (GPL-3), where it is distributed as `inst/extdata/gtc.csv`.

| Column | Meaning |
| --- | --- |
| X, Y | Longitude and latitude (decimal degrees, EPSG:4326) |
| GTC | Glacier thickness change (response variable) |
| Slope, Elev, LakeArea, WinTem, SumTem, Pre, Aspect, SurAlbedo | Slope, elevation, lake area, winter temperature, summer temperature, precipitation, aspect and surface albedo (explanatory variables) |

To obtain it from the original source instead:

- In R, after installing localsp: `gtc = readr::read_csv(system.file("extdata/gtc.csv", package = "localsp"))`
- From GitHub: https://github.com/ausgis/localsp/blob/main/inst/extdata/gtc.csv

## gtc_subset_100.csv

A random subset of 100 rows of `gtc.csv`, used for validation. It was drawn in R following the example in the localsp documentation (`set.seed(42); gtc[sample.int(nrow(gtc), size = 100), ]`).

## r_ground_truth_100.csv

Output of the original R function `localsp::lisp()` on `gtc_subset_100.csv`, used as the reference for validating the Python implementation. Settings: formula `GTC ~ .`, distance threshold 8, `discnum = 3:5`, `discmethod = "quantile"`, `cores = 1`.

To regenerate it (requires R with localsp installed, see the main README):

```
Rscript benchmark/bench_r.R data/gtc_subset_100.csv r_output_100.csv 8
```

Running this with localsp 0.2.0, gdverse 1.7 and sdsfun 0.9.0 reproduced the file exactly (maximum difference 0).

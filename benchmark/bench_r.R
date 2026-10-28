# Benchmark one run of the original R LISP (localsp::lisp).
#
# Usage: Rscript bench_r.R <subset.csv> <out_result.csv> <threshold>
# Prints one line of JSON, same fields as bench_py.py.
# Memory is read from /proc/self/status exactly as in the Python script.

suppressPackageStartupMessages(library(localsp))

mem_mb <- function(field) {
  lines <- readLines("/proc/self/status")
  line <- lines[startsWith(lines, paste0(field, ":"))]
  as.numeric(strsplit(trimws(sub(paste0(field, ":"), "", line)), "\\s+")[[1]][1]) / 1024
}

args <- commandArgs(trailingOnly = TRUE)
subset_csv <- args[1]; out_csv <- args[2]; threshold <- as.numeric(args[3])

df <- read.csv(subset_csv)
coords <- as.matrix(df[, c("X", "Y")])
data <- df[, setdiff(names(df), c("X", "Y"))]   # drop coordinates so GTC ~ . uses the 8 predictors

rss_before <- mem_mb("VmRSS")
t0 <- proc.time()[["elapsed"]]
# distance matrix is built INSIDE the timed block so R and Python time the same work
distmat <- as.matrix(dist(coords))
res <- localsp::lisp(GTC ~ ., data = data, threshold = threshold, distmat = distmat,
                     discnum = 3:5, discmethod = "quantile", cores = 1)
seconds <- proc.time()[["elapsed"]] - t0
peak <- mem_mb("VmHWM")

write.csv(as.data.frame(res), out_csv, row.names = FALSE)
cat(sprintf('{"lang": "R", "n": %d, "seconds": %.6f, "rss_before_mb": %.3f, "peak_mb": %.3f}\n',
            nrow(df), seconds, rss_before, peak))

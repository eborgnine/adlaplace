#!/usr/bin/env Rscript
# Knit one vignette with ADLAPLACE_VIGNETTE_NUM_THREADS set.
# Args: <vignette_dir> <vignette.Rmd basename> <output.html> <num_threads>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: knit-one-vignette.R <dir> <basename.Rmd> <out.html> <threads>", call. = FALSE)
}
vignette_dir <- args[[1L]]
basename_rmd <- args[[2L]]
out_html <- args[[3L]]
threads <- as.integer(args[[4L]])
if (is.na(threads) || threads < 1L) {
  stop("threads must be a positive integer", call. = FALSE)
}

Sys.setenv(
  NOT_CRAN = "true",
  ADLAPLACE_VIGNETTE_NUM_THREADS = as.character(threads)
)
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("rmarkdown is required", call. = FALSE)
}
setwd(vignette_dir)
rmarkdown::render(
  basename_rmd,
  output_file = out_html,
  quiet = TRUE,
  envir = new.env(parent = globalenv())
)

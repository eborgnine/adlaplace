#!/usr/bin/env Rscript
# Same-process adlaplace testthat run (Windows abort hunter).
#
# Full suite (default): test_check() + Location reporter — closest to R CMD check.
# Prefix/slice: one-process alphabetical test_file() loop for binary search.
#
# Usage:
#   Rscript run-adlaplace-testthat-suite.R
#   Rscript run-adlaplace-testthat-suite.R --from=1 --to=20
#   Rscript run-adlaplace-testthat-suite.R --max-files=20

args <- commandArgs(trailingOnly = TRUE)
parse_flag <- function(name, default = NA_integer_) {
  pref <- paste0("--", name, "=")
  hit <- args[startsWith(args, pref)]
  if (!length(hit)) {
    return(default)
  }
  as.integer(sub(pref, "", hit[[1L]]))
}

max_files <- parse_flag("max-files")
from <- parse_flag("from")
to <- parse_flag("to")
slice <- !is.na(from) || !is.na(to) || !is.na(max_files)

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%OS3"), paste(...), "\n", sep = " ")
  flush.console()
}

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required", call. = FALSE)
}
if (!requireNamespace("adlaplace", quietly = TRUE)) {
  stop("adlaplace is not installed", call. = FALSE)
}

say(
  "start", R.version.string,
  "adlaplace", as.character(utils::packageVersion("adlaplace")),
  "RCppAD", as.character(utils::packageVersion("RCppAD"))
)

library(testthat)
library(adlaplace)

reporter <- MultiReporter$new(list(
  LocationReporter$new(),
  FailReporter$new()
))

if (!slice) {
  say("mode test_check(adlaplace)")
  # Mirrors tests/testthat.R under R CMD check.
  test_check("adlaplace", reporter = reporter)
  say("SUITE DONE ok (test_check)")
  quit(save = "no", status = 0L)
}

root <- normalizePath(getwd())
test_dir <- file.path(root, "adlaplace", "tests", "testthat")
files <- sort(list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE))
n <- length(files)
if (!n) {
  stop("no test-*.R under ", test_dir, call. = FALSE)
}

if (is.na(from)) {
  from <- 1L
}
if (is.na(to)) {
  to <- n
}
if (!is.na(max_files)) {
  to <- min(to, from + max_files - 1L)
}
from <- max(1L, as.integer(from))
to <- min(n, as.integer(to))
files <- files[seq.int(from, to)]

say("mode test_file slice", sprintf("%d:%d/%d", from, to, n))

for (i in seq_along(files)) {
  f <- files[[i]]
  idx <- from + i - 1L
  say(sprintf("===== FILE %d/%d %s =====", idx, n, basename(f)))
  test_file(
    f,
    reporter = reporter,
    package = "adlaplace",
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
  say(sprintf("===== END FILE %d/%d %s =====", idx, n, basename(f)))
}

say("SUITE DONE ok", sprintf("slice %d:%d", from, to))

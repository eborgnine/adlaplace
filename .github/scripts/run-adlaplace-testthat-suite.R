#!/usr/bin/env Rscript
# Same-process adlaplace testthat run (Windows abort hunter).
#
# Full suite (default): test_dir() on source tests/ + installed package.
# (test_check() looks under the installed library, which has no tests/.)
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

quit_status <- function(status) {
  quit(save = "no", status = as.integer(status))
}

results_failed <- function(res) {
  df <- as.data.frame(res)
  if (!nrow(df)) {
    return(FALSE)
  }
  failed <- if ("failed" %in% names(df)) sum(df$failed, na.rm = TRUE) else 0
  errored <- if ("error" %in% names(df)) sum(df$error, na.rm = TRUE) else 0
  (failed + errored) > 0
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

root <- normalizePath(getwd())
test_path <- file.path(root, "adlaplace", "tests", "testthat")
if (!dir.exists(test_path)) {
  stop("missing test dir: ", test_path, call. = FALSE)
}
files_all <- sort(list.files(test_path, pattern = "^test-.*\\.R$", full.names = TRUE))
n <- length(files_all)
if (!n) {
  stop("no test-*.R under ", test_path, call. = FALSE)
}
say("test_dir", test_path, "n_files", n)

reporter <- MultiReporter$new(list(
  LocationReporter$new(),
  FailReporter$new()
))

if (!slice) {
  # Closest to R CMD check's testthat.R, but tests come from the source tree
  # because installed packages do not ship tests/.
  say("mode test_dir(source tests, load_package=installed)")
  res <- test_dir(
    test_path,
    package = "adlaplace",
    load_package = "installed",
    reporter = reporter,
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
  if (results_failed(res)) {
    say("SUITE DONE with failures")
    quit_status(1L)
  }
  say("SUITE DONE ok (test_dir)")
  quit_status(0L)
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
files <- files_all[seq.int(from, to)]

say("mode test_file slice", sprintf("%d:%d/%d", from, to, n))

any_fail <- FALSE
for (i in seq_along(files)) {
  f <- files[[i]]
  idx <- from + i - 1L
  say(sprintf("===== FILE %d/%d %s =====", idx, n, basename(f)))
  res <- test_file(
    f,
    reporter = reporter,
    package = "adlaplace",
    stop_on_failure = FALSE,
    stop_on_warning = FALSE
  )
  if (results_failed(res)) {
    any_fail <- TRUE
  }
  say(sprintf("===== END FILE %d/%d %s =====", idx, n, basename(f)))
}

if (any_fail) {
  say("SUITE DONE with failures", sprintf("slice %d:%d", from, to))
  quit_status(1L)
}
say("SUITE DONE ok", sprintf("slice %d:%d", from, to))
quit_status(0L)

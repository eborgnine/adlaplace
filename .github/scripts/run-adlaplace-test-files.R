#!/usr/bin/env Rscript
# Run one testthat file against the installed adlaplace package.
# Args: path to test-*.R (relative to repo root or absolute).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: run-adlaplace-test-files.R <test-file>", call. = FALSE)
}

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%OS3"), paste(...), "\n")
  flush.console()
}

path <- args[[1L]]
if (!file.exists(path)) {
  stop("missing test file: ", path, call. = FALSE)
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required", call. = FALSE)
}
if (!requireNamespace("adlaplace", quietly = TRUE)) {
  stop("adlaplace is not installed", call. = FALSE)
}

say("BEGIN", path)
# Location reporter prints file + test as they run; last line before abort
# identifies the crashing test under a hard Windows exit.
ok <- TRUE
tryCatch(
  testthat::test_file(
    path,
    reporter = c("location", "fail"),
    package = "adlaplace",
    stop_on_failure = TRUE,
    stop_on_warning = FALSE
  ),
  error = function(e) {
    say("ERROR", conditionMessage(e))
    ok <<- FALSE
  }
)
say("END", path, if (ok) "ok" else "FAILED")
quit(save = "no", status = if (ok) 0L else 1L)

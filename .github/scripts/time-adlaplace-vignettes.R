#!/usr/bin/env Rscript
# Knit every adlaplace vignette at two thread counts and print a speedup table.
#
# Env:
#   ADLAPLACE_VIGNETTE_DIR      path to vignettes/ (required)
#   ADLAPLACE_VIGNETTE_OUT      html output directory (required)
#   ADLAPLACE_VIGNETTE_RESULTS  pipe-table path (required)
#   ADLAPLACE_TIMING_THREADS_A  default 1
#   ADLAPLACE_TIMING_THREADS_B  default 4
#   NOT_CRAN                    should be "true" for full vignettes

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%OS3"), paste(...), "\n", sep = " ")
  flush.console()
}

env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v) || is.na(suppressWarnings(as.integer(v)))) {
    as.integer(default)
  } else {
    as.integer(v)
  }
}

vignette_dir <- Sys.getenv("ADLAPLACE_VIGNETTE_DIR", "")
out_dir <- Sys.getenv("ADLAPLACE_VIGNETTE_OUT", "")
results <- Sys.getenv("ADLAPLACE_VIGNETTE_RESULTS", "")
threads_a <- env_int("ADLAPLACE_TIMING_THREADS_A", 1L)
threads_b <- env_int("ADLAPLACE_TIMING_THREADS_B", 4L)

if (!nzchar(vignette_dir) || !dir.exists(vignette_dir)) {
  stop("ADLAPLACE_VIGNETTE_DIR missing or not a directory", call. = FALSE)
}
if (!nzchar(out_dir) || !nzchar(results)) {
  stop("ADLAPLACE_VIGNETTE_OUT and ADLAPLACE_VIGNETTE_RESULTS are required", call. = FALSE)
}
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("rmarkdown is required", call. = FALSE)
}
if (!requireNamespace("adlaplace", quietly = TRUE)) {
  stop("adlaplace is not installed", call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
vignettes <- sort(list.files(vignette_dir, pattern = "\\.Rmd$", full.names = TRUE))
if (!length(vignettes)) {
  stop("no .Rmd files in ", vignette_dir, call. = FALSE)
}

say(
  "timing config: threads_a=", threads_a, " threads_b=", threads_b,
  " n_vignettes=", length(vignettes),
  " has_openmp=", adlaplace::has_openmp(),
  " NOT_CRAN=", Sys.getenv("NOT_CRAN", "<unset>")
)

write_row <- function(vignette, threads, elapsed, exit_code) {
  cat(
    paste(vignette, threads, sprintf("%.3f", elapsed), exit_code, sep = "|"),
    "\n",
    file = results,
    append = TRUE,
    sep = ""
  )
}

knit_one <- function(vignette_path, threads) {
  base <- sub("\\.Rmd$", "", basename(vignette_path))
  out_html <- file.path(out_dir, sprintf("%s-t%d.html", base, threads))
  say("START", base, "threads=", threads)

  helper <- normalizePath(
    file.path(dirname(dirname(vignette_dir)), ".github", "scripts", "knit-one-vignette.R"),
    mustWork = TRUE
  )

  # Child process: fresh env + fresh OpenMP/CppAD high-water mark per thread count.
  tm <- system.time(
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      args = c(
        helper,
        vignette_dir,
        basename(vignette_path),
        normalizePath(out_html, mustWork = FALSE),
        as.character(threads)
      ),
      stdout = "",
      stderr = ""
    ),
    gcFirst = FALSE
  )
  if (is.na(status)) {
    status <- 1L
  }
  elapsed <- unname(tm[["elapsed"]])
  if (status != 0L) {
    say("FAIL", base, "threads=", threads, "status=", status)
  } else {
    say("OK", base, "threads=", threads, "elapsed=", sprintf("%.3f", elapsed))
  }
  write_row(base, threads, elapsed, status)
  list(vignette = base, threads = threads, elapsed = elapsed, status = as.integer(status))
}

rows <- list()
for (threads in c(threads_a, threads_b)) {
  say("===== pass threads=", threads, "=====")
  for (v in vignettes) {
    rows[[length(rows) + 1L]] <- knit_one(v, threads)
  }
}

df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
a <- df[df$threads == threads_a, c("vignette", "elapsed", "status")]
b <- df[df$threads == threads_b, c("vignette", "elapsed", "status")]
names(a) <- c("vignette", "sec_a", "ec_a")
names(b) <- c("vignette", "sec_b", "ec_b")
m <- merge(a, b, by = "vignette", all = TRUE)
m$speedup <- m$sec_a / m$sec_b
m$delta_sec <- m$sec_a - m$sec_b

say("===== 1 vs", threads_b, "comparison =====")
print(m, row.names = FALSE, digits = 3)
ok <- m[is.finite(m$sec_a) & is.finite(m$sec_b) & m$ec_a == 0 & m$ec_b == 0, ]
if (nrow(ok)) {
  cat(sprintf(
    "TOTAL sec_%s=%.3f sec_%s=%.3f speedup=%.3f\n",
    threads_a, sum(ok$sec_a), threads_b, sum(ok$sec_b),
    sum(ok$sec_a) / sum(ok$sec_b)
  ))
}

any_fail <- any(vapply(rows, function(r) r$status != 0L, logical(1)))
if (any_fail) {
  say("DONE with failures")
  quit(save = "no", status = 1L)
}
say("DONE ok")
quit(save = "no", status = 0L)

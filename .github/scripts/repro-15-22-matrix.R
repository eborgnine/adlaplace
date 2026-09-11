#!/usr/bin/env Rscript
# Standalone repro for the Windows "file 15 then file 22" abort.
#
# Replicates test-cppad-teardown.R (phase A: parallel log_lik_laplace deriv=TRUE
# loop) and test-fun-obj-fdfh.R (phase B: multi-thread fun_obj_fdfh) without
# testthat, so we can vary knobs via env vars and isolate the trigger.
#
# Knobs (env vars):
#   ADLAPLACE_REPRO_SKIP_A=1              skip phase A (control)
#   ADLAPLACE_REPRO_PHASE_A_THREADS=N     default 4
#   ADLAPLACE_REPRO_PHASE_A_ITERS=N       default 20
#   ADLAPLACE_REPRO_PHASE_B_THREADS=N     default 2
#   ADLAPLACE_REPRO_INTERLUDE=none|gc|warm|gcwarm   between A and B
#
# Prints timestamped progress to stderr. Exits 0 on ok, 1 on a caught
# non-finite/bad-values result, and 127 if the process hard-aborts (the shell
# observes 127; we never emit it ourselves).

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%OS3"), paste(...), "\n", sep = " ", file = stderr())
  flush(stderr())
}

env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "" || is.na(suppressWarnings(as.integer(v)))) default else as.integer(v)
}
env_str <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else v
}

skip_a <- env_int("ADLAPLACE_REPRO_SKIP_A", 0L) == 1L
a_threads <- env_int("ADLAPLACE_REPRO_PHASE_A_THREADS", 4L)
a_iters <- env_int("ADLAPLACE_REPRO_PHASE_A_ITERS", 20L)
b_threads <- env_int("ADLAPLACE_REPRO_PHASE_B_THREADS", 2L)
interlude <- env_str("ADLAPLACE_REPRO_INTERLUDE", "none")

suppressPackageStartupMessages({
  library(adlaplace)
  library(Matrix)
})

say("repro config: skip_a=", skip_a, " a_threads=", a_threads,
    " a_iters=", a_iters, " b_threads=", b_threads, " interlude=", interlude)
say("adlaplace", as.character(utils::packageVersion("adlaplace")),
    "RCppAD", as.character(utils::packageVersion("RCppAD")),
    "openmp=", adlaplace::has_openmp())

# ---- helpers mirroring tests/testthat/helper-ad-data.R ---------------------

test_ad_data <- function(y, A, X, config, theta_local_row = 0L) {
  adlaplace:::density_data_from_config_matrices(
    y = y, A = A, X = X, config = config, theta_local_row = theta_local_row
  )
}

test_random_shard <- function(data, config, gamma_ids, theta_id, Q,
                               density = "random_diagonal") {
  n_gamma_full <- length(config$gamma)
  n_term <- length(gamma_ids)
  gamma_map <- Matrix::sparseMatrix(
    i = as.integer(gamma_ids),
    j = seq.int(0L, length.out = n_term),
    dims = c(n_gamma_full, n_term), index1 = FALSE, repr = "C"
  )
  adlaplace:::density_data(
    beta_map = Matrix::Matrix(nrow = length(config$beta), ncol = 0L),
    gamma_map = gamma_map, theta_map = data@theta_map,
    ad_kind = "random", density = density, precision = Q
  )
}

as_shard <- function(data, ad_kind, density, precision = NULL) {
  data@ad_kind <- ad_kind
  data@density <- density
  data@precision <- precision
  data
}

# ---- phase A: cppad-teardown fixture + parallel log_lik_laplace loop --------

phase_a <- function(num_threads, iters) {
  say("PHASE A start threads=", num_threads, " iters=", iters)
  set.seed(42)
  Nobs <- 80L
  X <- Matrix::Matrix(cbind(1, stats::rbinom(Nobs, 1, 0.5)))
  Amat <- Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(4L, Nobs, replace = TRUE))
  config <- list(
    beta = rep(0, ncol(X)), theta = c(-1, -1), transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    obs_groups = adlaplace::obs_groups(Amat, num_shards = 20L),
    verbose = FALSE, package = "adlaplace"
  )
  model <- test_ad_data(y = rpois(Nobs, 2), A = Amat, X = X, config, theta_local_row = 1L)
  random_shard <- test_random_shard(
    data = model, config = config,
    gamma_ids = seq.int(0L, length.out = ncol(Amat)), theta_id = 0L, Q = rep(1, ncol(Amat))
  )
  ad_ptr <- do.call(c, list(
    adlaplace::ad_pack_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_pack_ptr(random_shard, config),
    adlaplace::ad_pack_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  ad_pack <- adlaplace::ad_pack(ad_ptr, num_threads = num_threads)
  args <- list(
    x = c(config$beta, config$theta), config = list(verbose = FALSE),
    gamma = config$gamma, ad_pack = ad_pack,
    control = list(maxit = 3L, report.level = 0, report.freq = 0), deriv = TRUE
  )
  for (i in seq_len(iters)) {
    ll <- do.call(adlaplace::log_lik_laplace, args)
    if (!is.finite(ll$log_lik)) stop("phase A iter ", i, " non-finite log_lik")
  }
  say("PHASE A done")
  invisible(ad_pack)
}

# ---- phase B: fun-obj-fdfh multi-thread -----------------------------------

phase_b <- function(num_threads) {
  say("PHASE B start threads=", num_threads)
  set.seed(0)
  Nobs <- 120L
  Nrandom1 <- 4L; Nrandom2 <- 5L
  X <- Matrix::Matrix(cbind(1, stats::rbinom(Nobs, 1, 0.5)))
  Amat <- cbind(
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom1, Nobs, replace = TRUE)),
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom2, Nobs, replace = TRUE))
  )
  config <- list(
    beta = rep(0, ncol(X)), theta = c(-1, -1, -1), transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    obs_groups = adlaplace::obs_groups(Amat, num_shards = 40L),
    verbose = FALSE, package = "adlaplace"
  )
  n_gamma <- length(config$gamma)
  model <- test_ad_data(y = rpois(Nobs, 2), A = Amat, X = X, config,
                        theta_local_row = length(config$theta) - 1L)
  random_shard <- test_random_shard(
    data = model, config = config,
    gamma_ids = seq.int(0L, length.out = ncol(Amat)), theta_id = 0L, Q = rep(1, ncol(Amat))
  )
  ad_ptr <- do.call(c, list(
    adlaplace::ad_pack_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_pack_ptr(random_shard, config),
    adlaplace::ad_pack_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  dens_ptr <- adlaplace::clone_ad_pack_ptr(ad_ptr)
  ad_pack <- adlaplace::ad_pack(ad_ptr, num_threads = num_threads)
  gamma <- stats::rnorm(n_gamma, sd = 0.1)
  parameters <- c(config$beta, config$theta)
  fo <- adlaplace::fun_obj_fdfh(ad_pack, parameters, gamma, inner = TRUE)
  ok <- is.finite(fo$f) && all(is.finite(fo$grad))
  say("PHASE B result=", if (ok) "ok" else "bad-values")
  if (!ok) quit(save = "no", status = 1L)
  invisible(ad_pack)
}

# ---- run ------------------------------------------------------------------

if (!skip_a) {
  phase_a(a_threads, a_iters)
}

if (interlude == "gc" || interlude == "gcwarm") {
  say("INTERLUDE gc()")
  gc(); gc()
  say("INTERLUDE gc done")
}
if (interlude == "warm" || interlude == "gcwarm") {
  say("INTERLUDE warm_openmp_runtime()")
  # warm_openmp_runtime() is @keywords internal (not in NAMESPACE), so use
  # triple-colon internal access from this standalone script.
  adlaplace:::warm_openmp_runtime()
  say("INTERLUDE warm done")
}

phase_b(b_threads)
say("repro finished ok")
quit(save = "no", status = 0L)

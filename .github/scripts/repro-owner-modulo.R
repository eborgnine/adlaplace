#!/usr/bin/env Rscript
# Minimal repro for Windows CI expectation failures after the team-size latch.
#
# After any parallel team >2 has been used in-process, Windows latch raises
# later ad_pack(..., num_threads = 2) to the high-water mark. Tests that assert
#   owners == (0:(n-1)) %% 2
# then fail (saved as _problems/test-log-lik-deriv-parallel-48.R etc.).
#
# Knobs (env vars):
#   ADLAPLACE_REPRO_SKIP_WARMUP=1     skip the 4-thread warmup (control)
#   ADLAPLACE_REPRO_WARMUP_THREADS=N  default 4
#   ADLAPLACE_REPRO_PROBE_THREADS=N   default 2
#   ADLAPLACE_CLAMP_TEAM_THREADS=0    disable Windows clamp-up
#
# Exit: 0 if owners match probe %% N; 1 on mismatch / bad values; 127 if abort.

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%OS3"), paste(...), "\n", sep = " ", file = stderr())
  flush(stderr())
}

env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "" || is.na(suppressWarnings(as.integer(v)))) default else as.integer(v)
}

skip_warmup <- env_int("ADLAPLACE_REPRO_SKIP_WARMUP", 0L) == 1L
warmup_threads <- env_int("ADLAPLACE_REPRO_WARMUP_THREADS", 4L)
probe_threads <- env_int("ADLAPLACE_REPRO_PROBE_THREADS", 2L)

suppressPackageStartupMessages({
  library(adlaplace)
  library(Matrix)
})

say(
  "owner-modulo config: skip_warmup=", skip_warmup,
  " warmup_threads=", warmup_threads,
  " probe_threads=", probe_threads,
  " clamp_env=", Sys.getenv("ADLAPLACE_CLAMP_TEAM_THREADS", "<unset>"),
  " openmp=", adlaplace:::has_openmp()
)

# Mirror test-reorder-shards / test-log-lik-deriv-parallel fixtures closely
# enough that n_groups > 2 and %%4 owners diverge from %%2.
make_ptr <- function(seed = 41L) {
  set.seed(seed)
  n_obs <- 60L
  X <- Matrix::Matrix(cbind(1, stats::rbinom(n_obs, 1, 0.5)))
  Amat <- Matrix::sparseMatrix(
    i = seq_len(n_obs),
    j = sample(5L, n_obs, replace = TRUE),
    x = 1
  )
  config <- list(
    beta = rep(0, ncol(X)),
    theta = -1,
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    obs_groups = adlaplace::obs_groups(Amat, num_shards = 12L),
    verbose = FALSE,
    package = "adlaplace"
  )
  model <- adlaplace:::density_data_from_config_matrices(
    y = stats::rpois(n_obs, 2),
    A = Amat,
    X = X,
    config = config,
    theta_local_row = 0L
  )
  as_shard <- function(data, ad_kind, density) {
    data@ad_kind <- ad_kind
    data@density <- density
    data
  }
  ptr <- do.call(c, list(
    adlaplace::ad_pack_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_pack_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  list(ptr = ptr, config = config)
}

run_pack <- function(num_threads, label) {
  built <- make_ptr()
  requested <- as.integer(num_threads)
  latched <- adlaplace:::latch_parallel_threads(requested)
  say(label, "requested=", requested, "latched=", latched)
  af <- adlaplace::ad_pack(
    built$ptr,
    num_threads = requested,
    reorder_shards = "none"
  )
  n <- adlaplace:::n_groups(af@ptr)
  owners <- vapply(seq_len(n) - 1L, function(g) {
    adlaplace:::get_thread_owner(af@ptr, g)
  }, integer(1))
  # Same assertion as the failing tests: modulo the *requested* thread count.
  expect_mod <- (seq_len(n) - 1L) %% requested
  n_unique <- length(unique(owners))
  ok <- identical(owners, expect_mod)
  say(
    label, "n_groups=", n, "n_unique_owners=", n_unique,
    "max_owner=", if (length(owners)) max(owners) else NA_integer_,
    "owners_match_%%", requested, "=", ok
  )
  if (!ok) {
    say(label, "owners=", paste(owners, collapse = ","))
    say(label, "expect=", paste(expect_mod, collapse = ","))
  }
  list(ok = ok, requested = requested, latched = latched, owners = owners)
}

if (!skip_warmup) {
  say("WARMUP start")
  wu <- run_pack(warmup_threads, "WARMUP")
  if (!wu$ok) {
    say("WARMUP unexpected owner mismatch")
    quit(save = "no", status = 1L)
  }
  say("WARMUP done")
}

probe <- run_pack(probe_threads, "PROBE")
if (!probe$ok) {
  say("PROBE owner mismatch (latch raised team above requested)")
  quit(save = "no", status = 1L)
}
say("owner-modulo repro finished ok")
quit(save = "no", status = 0L)

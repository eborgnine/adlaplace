# Debug script for macOS AD stale-parameter behavior (plan Phase 1-2).
# Run: Rscript inst/debug/macos_ad_debug.R
# Or source after vignette setup.

suppressPackageStartupMessages({
  library(Matrix)
  library(adlaplace)
  library(adlaplaceExample)
})

# Vignette data setup (matches adlaplace_test.qmd)
Nobs <- 1000L
Nrandom1 <- 10L
Nrandom2 <- 25L
set.seed(0L)

X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
Adf <- data.frame(
  r1 = sample(Nrandom1, Nobs, replace = TRUE),
  r2 = sample(Nrandom2, Nobs, replace = TRUE)
)
AmatList <- list(
  r1 = Matrix::sparseMatrix(i = seq_len(Nobs), j = Adf$r1),
  r2 = Matrix::sparseMatrix(i = seq_len(Nobs), j = Adf$r2)
)
Amat <- do.call(cbind, AmatList)

beta <- rep(1, ncol(X))
thetaOrig <- c(sd1 = 0.1, sd2 = 0.1, omega = 0.7, alpha = 0.8)
gamma <- rnorm(ncol(Amat), sd = rep(thetaOrig[1:2], c(Nrandom1, Nrandom2)))
eta_true <- as.vector(X %*% beta + Amat %*% gamma)

if (requireNamespace("sn", quietly = TRUE)) {
  y <- sn::rsn(Nobs, xi = eta_true, omega = thetaOrig["omega"], alpha = thetaOrig["alpha"])
} else {
  y <- stats::rnorm(Nobs, mean = eta_true, sd = thetaOrig["omega"])
}

build_model <- function(num_threads, num_shards) {
  config <- list(
    beta = beta,
    theta = c(log(thetaOrig[c("sd1", "sd2", "omega")]), thetaOrig["alpha"]),
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    shards = adlaplace::ad_shards(Amat, num_shards = num_shards),
    num_threads = as.integer(num_threads),
    verbose = FALSE,
    package = "adlaplaceExample"
  )
  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  n_theta <- length(config$theta)

  data_obs <- adlaplace::ad_data(
    y = y, A = Amat, X = X,
    beta_map = Matrix::Diagonal(n_beta),
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(c(
      which(names(thetaOrig) == "omega"),
      which(names(thetaOrig) == "alpha")
    ), n_theta),
    ad_kind = "observations",
    ad_fun = "skewnormal_obs",
    package = "adlaplaceExample"
  )
  data_extra <- adlaplace::ad_data(
    y = data_obs@y,
    beta_map = data_obs@beta_map,
    gamma_map = data_obs@gamma_map,
    theta_map = data_obs@theta_map,
    ad_kind = "parameters",
    ad_fun = "skewnormal_extra",
    package = "adlaplaceExample"
  )
  ad_fun_obs <- adlaplace::ad_fun_ptr(data = data_obs, config = config)
  ad_fun_extra <- adlaplace::ad_fun_ptr(data = data_extra, config = config)

  list(
    config = config, n_beta = n_beta, n_gamma = n_gamma, n_theta = n_theta,
    ad_fun_obs = ad_fun_obs, ad_fun_extra = ad_fun_extra,
    xx = c(beta, gamma, config$theta)
  )
}

omega_scan_cmp <- function(m) {
  xx <- m$xx
  omega_idx <- m$n_beta + m$n_gamma + which(names(thetaOrig) == "omega")
  log_omega_true <- m$config$theta[which(names(thetaOrig) == "omega")]
  seq_log_omega <- log_omega_true + seq(-0.5, 0.5, length.out = 7L)

  ad_obs <- if (!is.null(m$ad_fun_obs)) {
    vapply(seq_log_omega, function(lo) {
      x <- xx
      x[omega_idx] <- lo
      adlaplace::joint_log_dens(m$ad_fun_obs, x, negative = FALSE)
    }, numeric(1))
  } else {
    rep(NA_real_, length(seq_log_omega))
  }
  ad_extra <- vapply(seq_log_omega, function(lo) {
    x <- xx
    x[omega_idx] <- lo
    adlaplace::joint_log_dens(m$ad_fun_extra, x, negative = FALSE)
  }, numeric(1))

  data.frame(
    log_omega = seq_log_omega,
    ad_obs = ad_obs,
    ad_extra = ad_extra,
    obs_range = if (all(is.na(ad_obs))) NA_real_ else max(ad_obs) - min(ad_obs),
    extra_range = max(ad_extra) - min(ad_extra)
  )
}

fd_vs_grad <- function(m, shard = c("obs", "extra")) {
  shard <- match.arg(shard)
  ptr <- if (shard == "obs") m$ad_fun_obs else m$ad_fun_extra
  xx <- m$xx
  omega_idx <- m$n_beta + m$n_gamma + which(names(thetaOrig) == "omega")
  lo0 <- m$config$theta[which(names(thetaOrig) == "omega")]
  eps <- 1e-5

  f <- function(lo) {
    x <- xx
    x[omega_idx] <- lo
    adlaplace::joint_log_dens(ptr, x, negative = FALSE)
  }
  fd <- (f(lo0 + eps) - f(lo0 - eps)) / (2 * eps)
  x0 <- xx
  x0[omega_idx] <- lo0
  g <- adlaplace::grad(ptr, x0, negative = FALSE)
  list(shard = shard, fd = fd, grad_omega = g[omega_idx], f0 = f(lo0))
}

build_extra_only <- function(num_threads = 1L) {
  config <- list(
    beta = beta,
    theta = c(log(thetaOrig[c("sd1", "sd2", "omega")]), thetaOrig["alpha"]),
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    shards = adlaplace::ad_shards(Amat, num_shards = 1L),
    num_threads = as.integer(num_threads),
    verbose = FALSE,
    package = "adlaplaceExample"
  )
  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  data_extra <- adlaplace::ad_data(
    y = y,
    beta_map = Matrix::Diagonal(n_beta),
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(c(
      which(names(thetaOrig) == "omega"),
      which(names(thetaOrig) == "alpha")
    ), length(config$theta)),
    ad_kind = "parameters",
    ad_fun = "skewnormal_extra",
    package = "adlaplaceExample"
  )
  list(
    config = config,
    n_beta = n_beta,
    n_gamma = n_gamma,
    ad_fun_extra = adlaplace::ad_fun_ptr(data = data_extra, config = config),
    xx = c(beta, config$gamma, config$theta)
  )
}

cat("=== Phase 1.1: link OK (see shell otool) ===\n")
cat("=== Phase 1.2/2: scans and FD vs grad ===\n\n")

for (label in list(
  c(threads = 10L, shards = 100L),
  c(threads = 1L, shards = 100L),
  c(threads = 1L, shards = 1L)
)) {
  cat(sprintf("--- num_threads=%d, num_shards=%d ---\n", label[["threads"]], label[["shards"]]))
  m <- build_model(label[["threads"]], label[["shards"]])
  cmp <- omega_scan_cmp(m)
  print(cmp)
  cat("FD/grad obs:  "); print(fd_vs_grad(m, "obs"))
  cat("FD/grad extra:"); print(fd_vs_grad(m, "extra"))
  cat("\n")
}

cat("=== Extra-only build (no obs shards first) ===\n")
m_extra <- build_extra_only(1L)
cmp_extra <- omega_scan_cmp(m_extra)
print(cmp_extra)
cat("extra_range:", cmp_extra$extra_range[1], "\n")

cat("\n=== Mixed-DSO ad_fun_ptr shards: parallel deriv=TRUE ===\n")
config_parallel <- list(
  beta = beta,
  theta = c(log(thetaOrig[c("sd1", "sd2", "omega")]), thetaOrig["alpha"]),
  transform_theta = TRUE,
  gamma = rep(0, ncol(Amat)),
  shards = adlaplace::ad_shards(Amat, num_shards = 20L),
  num_threads = 4L,
  verbose = FALSE
)
n_beta_p <- length(beta)
n_gamma_p <- ncol(Amat)
n_theta_p <- length(config_parallel$theta)
obs_theta_idx <- (n_theta_p - 1L):n_theta_p
rand_theta_idx <- seq_len(n_theta_p - 2L)
obs_p <- adlaplace::ad_data(
  y = y, A = Amat, X = X,
  beta_map = Matrix::Diagonal(n_beta_p),
  gamma_map = Matrix::Diagonal(n_gamma_p),
  theta_map = list(obs_theta_idx, n_theta_p),
  ad_kind = "observations",
  ad_fun = "skewnormal_obs",
  package = "adlaplaceExample"
)
rand_p <- adlaplace::ad_data(
  beta_map = n_beta_p,
  gamma_map = Matrix::Diagonal(n_gamma_p),
  theta_map = list(rand_theta_idx, n_theta_p),
  ad_kind = "random",
  ad_fun = "random_diagonal",
  precision = rep(1, n_gamma_p)
)
extra_p <- adlaplace::ad_data(
  y = y,
  beta_map = n_beta_p,
  gamma_map = n_gamma_p,
  theta_map = list(obs_theta_idx, n_theta_p),
  ad_kind = "parameters",
  ad_fun = "skewnormal_extra",
  package = "adlaplaceExample"
)
ptrs_p <- c(
  adlaplace::ad_fun_ptr(obs_p, config_parallel),
  adlaplace::ad_fun_ptr(rand_p, config_parallel),
  adlaplace::ad_fun_ptr(extra_p, config_parallel)
)
ad_fun_parallel <- adlaplace::ad_fun(ptrs_p, num_threads = 4L)
ll_deriv <- adlaplace::log_lik_laplace(
  x = c(config_parallel$beta, config_parallel$theta),
  config = list(verbose = FALSE),
  gamma = config_parallel$gamma,
  ad_fun = ad_fun_parallel,
  control = list(maxit = 5L, report.level = 0, report.freq = 0),
  deriv = TRUE
)
cat(
  "log_lik finite:", is.finite(ll_deriv$log_lik),
  " trace3 all finite:", all(is.finite(ll_deriv$extra$trace3)),
  "\n"
)

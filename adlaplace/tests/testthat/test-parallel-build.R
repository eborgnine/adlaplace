# Parallel build_ad_fun_obs: taping + Hessian colouring/calibration over
# observation groups must reproduce the serial build exactly.

par_obs_fixture <- function(family, n = 240L, n_re = 24L, n_shards = 12L) {
  X <- Matrix::Matrix(cbind(1, rnorm(n), rbinom(n, 1, 0.4)))
  A <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = sample(n_re, n, replace = TRUE),
    x = 1,
    dims = c(n, n_re)
  )
  beta <- c(0.3, -0.2, 0.1)
  gamma <- rnorm(n_re, sd = 0.25)
  eta <- as.vector(X %*% beta + A %*% gamma)
  weights <- numeric(0)

  if (identical(family, "poisson_obs")) {
    y <- rpois(n, lambda = exp(pmin(eta, 4)))
    theta <- numeric(0)
    theta_map <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), dims = c(0L, 1L)
    )
  } else if (identical(family, "nbinom_obs")) {
    y <- rnbinom(n, size = 5, mu = exp(pmin(eta, 3)))
    theta <- log(1 / sqrt(5))
    theta_map <- Matrix::sparseMatrix(i = 1L, j = 1L, dims = c(1L, 1L))
  } else if (identical(family, "gamma_obs")) {
    group <- sample(c(0L, 1L), n, replace = TRUE)
    cv <- c(0.4, 0.6)
    shape <- 1 / cv[group + 1L]^2
    y <- rgamma(n, shape = shape, rate = shape / exp(eta - min(eta) + 1))
    theta <- log(cv)
    theta_map <- Matrix::Diagonal(2L)
    weights <- as.numeric(group)
  } else {
    stop("unknown family: ", family)
  }

  config <- list(
    beta = beta,
    gamma = gamma,
    theta = theta,
    transform_theta = TRUE,
    obs_groups = adlaplace::obs_groups(A, num_shards = n_shards),
    compact_tape = TRUE,
    verbose = FALSE
  )
  model <- adlaplace:::density_data(
    y = y,
    A = A,
    X = X,
    theta_map = theta_map,
    ad_kind = "observations",
    density = family,
    weights = weights
  )
  list(model = model, config = config, x = c(beta, gamma, theta))
}

expect_parallel_matches_serial <- function(fx, threads = 2L) {
  cfg_serial <- fx$config
  cfg_serial$num_threads <- 1L
  cfg_par <- fx$config
  cfg_par$num_threads <- as.integer(threads)

  ptr_serial <- adlaplace::ad_pack_ptr(fx$model, cfg_serial)
  ptr_par <- adlaplace::ad_pack_ptr(fx$model, cfg_par)

  n_g <- adlaplace:::n_groups(ptr_serial)
  expect_gt(n_g, 1L)
  expect_equal(adlaplace:::n_groups(ptr_par), n_g)
  for (g in seq_len(n_g) - 1L) {
    expect_identical(
      adlaplace:::get_sizes(ptr_par, g),
      adlaplace:::get_sizes(ptr_serial, g)
    )
    expect_identical(
      adlaplace:::get_sparse_pattern(ptr_par, g),
      adlaplace:::get_sparse_pattern(ptr_serial, g)
    )
  }

  x <- fx$x
  expect_identical(
    adlaplace::joint_log_dens(ptr_par, x, negative = FALSE),
    adlaplace::joint_log_dens(ptr_serial, x, negative = FALSE)
  )
  for (inner in c(FALSE, TRUE)) {
    expect_identical(
      adlaplace::grad(ptr_par, x, inner = inner, negative = FALSE),
      adlaplace::grad(ptr_serial, x, inner = inner, negative = FALSE)
    )
    expect_identical(
      as.matrix(adlaplace::hessian(
        ptr_par, x, inner = inner, negative = FALSE
      )),
      as.matrix(adlaplace::hessian(
        ptr_serial, x, inner = inner, negative = FALSE
      ))
    )
  }
  invisible(ptr_par)
}

test_that("parallel obs build matches serial (nbinom, 2 threads)", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(101)
  fx <- par_obs_fixture("nbinom_obs")
  expect_parallel_matches_serial(fx, threads = 2L)
})

test_that("parallel obs build matches serial (poisson, 2 threads)", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(102)
  fx <- par_obs_fixture("poisson_obs")
  expect_parallel_matches_serial(fx, threads = 2L)
})

test_that("parallel obs build matches serial (gamma, 2 threads)", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(103)
  fx <- par_obs_fixture("gamma_obs")
  expect_parallel_matches_serial(fx, threads = 2L)
})

test_that("parallel obs build matches serial (nbinom, 4 threads)", {
  skip_on_cran()
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(104)
  fx <- par_obs_fixture("nbinom_obs")
  expect_parallel_matches_serial(fx, threads = 4L)
})

test_that("error inside a parallel build worker surfaces as one R error", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(105)
  fx <- par_obs_fixture("gamma_obs")
  # gamma_obs rejects non-positive observations from inside the density,
  # i.e. on a worker thread in the parallel build.
  bad <- fx$model
  y <- bad@y
  y[1] <- 0
  bad@y <- y
  cfg <- fx$config
  cfg$num_threads <- 2L
  expect_error(
    adlaplace::ad_pack_ptr(bad, cfg),
    "gamma|parallel"
  )
  # The process survived the failed parallel build: a good parallel
  # build still works (and still uses its threads).
  cfg$num_threads <- 2L
  ptr <- adlaplace::ad_pack_ptr(fx$model, cfg)
  expect_gt(adlaplace:::n_groups(ptr), 1L)
  cfg_serial <- fx$config
  cfg_serial$num_threads <- 1L
  ptr_serial <- adlaplace::ad_pack_ptr(fx$model, cfg_serial)
  expect_identical(
    adlaplace::grad(ptr, fx$x, inner = FALSE, negative = FALSE),
    adlaplace::grad(ptr_serial, fx$x, inner = FALSE, negative = FALSE)
  )
})

test_that("ad_pack(num_threads) reaches the observation tape build", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(106)
  fx <- par_obs_fixture("nbinom_obs")
  fx$config$verbose <- TRUE
  out <- paste(
    capture.output(
      ap <- adlaplace::ad_pack(fx$model, fx$config, num_threads = 2L)
    ),
    collapse = "\n"
  )
  expect_match(out, "build_ad_fun_obs timing: threads 2", fixed = TRUE)
  expect_equal(
    as.integer(adlaplace:::get_configured_num_threads(ap@ptr)),
    2L
  )
})

test_that("parallel-built pack survives parallel eval, clone and teardown", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(107)
  fx <- par_obs_fixture("nbinom_obs")
  cfg <- fx$config
  cfg$num_threads <- 2L

  ptr_obs <- adlaplace::ad_pack_ptr(fx$model, cfg)

  clone <- adlaplace::clone_ad_pack_ptr(ptr_obs)
  expect_identical(
    adlaplace::grad(clone, fx$x, inner = FALSE, negative = FALSE),
    adlaplace::grad(ptr_obs, fx$x, inner = FALSE, negative = FALSE)
  )

  random_shard <- test_random_shard(
    data = fx$model,
    config = cfg,
    gamma_ids = seq.int(0L, length.out = length(cfg$gamma)),
    theta_id = 0L,
    Q = rep(1, length(cfg$gamma))
  )
  extra <- adlaplace::density_data(
    fx$model, ad_kind = "parameters", density = "nbinom_extra"
  )
  ptr <- do.call(c, list(
    ptr_obs,
    adlaplace::ad_pack_ptr(random_shard, cfg),
    adlaplace::ad_pack_ptr(extra, cfg)
  ))
  ap <- adlaplace::ad_pack(ptr, num_threads = 2L)

  args <- list(
    x = c(cfg$beta, cfg$theta),
    config = list(verbose = FALSE),
    gamma = cfg$gamma,
    ad_pack = ap,
    control = list(maxit = 4L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  ll1 <- do.call(adlaplace::log_lik_laplace, args)
  ll2 <- do.call(adlaplace::log_lik_laplace, args)
  expect_true(is.finite(ll1$log_lik))
  expect_true(is.finite(ll2$log_lik))
  expect_true(all(is.finite(ll2$deriv$d_neg_log_lik)))
  expect_true(all(is.finite(ll2$extra$trace3)))
})

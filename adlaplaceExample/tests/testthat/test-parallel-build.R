# Parallel build_ad_fun_obs in this backend DSO must match the serial
# build (the CppAD parallel guard is compiled per-DSO from adlaplace's
# installed headers).

test_that("parallel skewnormal_obs build matches serial", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  obs <- adlaplace::density_data(
    y = c(1, 2, 3, 0.5, 1.5),
    X = Matrix::Matrix(cbind(1, seq_len(5))),
    beta_map = Matrix::Diagonal(2L),
    gamma_map = Matrix::Matrix(nrow = 0L, ncol = 0L),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "observations",
    density = "skewnormal_obs",
    package = "adlaplaceExample"
  )

  make_config <- function(num_threads) {
    list(
      beta = c(0, 0),
      theta = c(0, log(0.5)),
      transform_theta = TRUE,
      gamma = numeric(0),
      num_threads = as.integer(num_threads),
      obs_groups = adlaplace::obs_groups(
        Matrix::sparseMatrix(
          i = seq_len(5), j = seq_len(5), x = 1, dims = c(5L, 5L)
        ),
        num_shards = 5L
      ),
      verbose = FALSE
    )
  }

  ptr_serial <- adlaplace::ad_pack_ptr(obs, make_config(1L))
  ptr_par <- adlaplace::ad_pack_ptr(obs, make_config(2L))

  expect_equal(
    adlaplace:::n_groups(ptr_par),
    adlaplace:::n_groups(ptr_serial)
  )
  expect_gt(adlaplace:::n_groups(ptr_par), 1L)
  for (g in seq_len(adlaplace:::n_groups(ptr_serial)) - 1L) {
    expect_identical(
      adlaplace:::get_sparse_pattern(ptr_par, g),
      adlaplace:::get_sparse_pattern(ptr_serial, g)
    )
  }

  x <- c(0.1, -0.2, 0.05, log(0.6))
  expect_identical(
    adlaplace::joint_log_dens(ptr_par, x, negative = FALSE),
    adlaplace::joint_log_dens(ptr_serial, x, negative = FALSE)
  )
  expect_identical(
    adlaplace::grad(ptr_par, x, inner = FALSE, negative = FALSE),
    adlaplace::grad(ptr_serial, x, inner = FALSE, negative = FALSE)
  )
  expect_identical(
    as.matrix(adlaplace::hessian(ptr_par, x, inner = FALSE, negative = FALSE)),
    as.matrix(adlaplace::hessian(ptr_serial, x, inner = FALSE, negative = FALSE))
  )
})

# Three FEM ssq modes:
#   analytic=TRUE,  tape=FALSE  (default)  closed-form, no tape recorded
#   analytic=TRUE,  tape=TRUE              closed-form over a calibrated tape
#   analytic=FALSE                         AD EvalShard (tape forced on)
# The AD shard is the reference oracle for values. Taped analytic is the
# reference for pattern indices, since those come from for_jac_sparsity.

ssq_pack <- function(alpha, transform, analytic, tape = !analytic,
                     seed_theta = c(1.4, 0.8)) {
  degree <- if (alpha == 2L) 2L else 3L
  nk <- if (alpha == 2L) 4L else 5L
  knots_list <- list(
    x = seq(0, 1, length.out = nk),
    y = seq(0, 1, length.out = nk)
  )
  g <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(g, knots_list, degree = degree)
  prec <- fem_precision_payload(fem, alpha = alpha)
  nr <- nrow(fem$C)
  theta <- if (transform) log(seed_theta) else seed_theta
  config <- list(
    gamma = rep(0, nr),
    theta = theta,
    transform_theta = transform,
    fem_analytic = analytic,
    fem_tape = tape
  )
  dd <- adlaplace::density_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    density = paste0("random_fem_ssq_", alpha),
    package = "adlaplaceFem",
    precision = prec
  )
  list(
    ptr = adlaplace::ad_pack_ptr(dd, config),
    n_gamma = nr,
    theta = theta
  )
}

ssq_grid <- expand.grid(
  alpha = c(2L, 3L),
  transform = c(TRUE, FALSE),
  KEEP.OUT.ATTRS = FALSE
)

expect_pattern_equal <- function(a, b, label) {
  pa <- adlaplace:::get_sparse_pattern(a$ptr, 0L)
  pb <- adlaplace:::get_sparse_pattern(b$ptr, 0L)
  sa <- adlaplace:::get_sizes(a$ptr, 0L)
  sb <- adlaplace:::get_sizes(b$ptr, 0L)
  expect_equal(sa, sb, info = label)
  expect_identical(pa$grad, pb$grad, info = label)
  expect_identical(pa$grad_inner, pb$grad_inner, info = label)
  expect_identical(pa$row_hess, pb$row_hess, info = label)
  expect_identical(pa$col_hess, pb$col_hess, info = label)
  expect_identical(pa$row_hess_inner, pb$row_hess_inner, info = label)
  expect_identical(pa$col_hess_inner, pb$col_hess_inner, info = label)
}

test_that("no-tape analytic patterns match the tape-discovered ones", {
  for (row in seq_len(nrow(ssq_grid))) {
    alpha <- ssq_grid$alpha[row]
    transform <- ssq_grid$transform[row]
    label <- sprintf("alpha=%d transform=%s", alpha, transform)
    notape <- ssq_pack(alpha, transform, analytic = TRUE, tape = FALSE)
    taped <- ssq_pack(alpha, transform, analytic = TRUE, tape = TRUE)
    expect_pattern_equal(notape, taped, label)

    tz <- adlaplace:::get_tape_sizes(notape$ptr, 0L)
    tt <- adlaplace:::get_tape_sizes(taped$ptr, 0L)
    expect_identical(tz$domain, 0L, info = label)
    expect_identical(tz$size_op, 0L, info = label)
    expect_gt(tt$domain, 0L)
    expect_gt(tt$size_op, 0L)
    expect_equal(tz$n_global, tt$n_global, info = label)
    expect_equal(tz$nnz_grad, tt$nnz_grad, info = label)
    expect_equal(tz$nnz_hes, tt$nnz_hes, info = label)
  }
})

test_that("all three ssq modes agree on value, gradient, Hessian", {
  for (row in seq_len(nrow(ssq_grid))) {
    alpha <- ssq_grid$alpha[row]
    transform <- ssq_grid$transform[row]
    label <- sprintf("alpha=%d transform=%s", alpha, transform)

    notape <- ssq_pack(alpha, transform, analytic = TRUE, tape = FALSE)
    taped <- ssq_pack(alpha, transform, analytic = TRUE, tape = TRUE)
    ad <- ssq_pack(alpha, transform, analytic = FALSE)

    set.seed(3)
    x <- c(rnorm(notape$n_gamma), notape$theta)

    f_ad <- adlaplace::joint_log_dens(ad$ptr, x, negative = FALSE)
    expect_equal(
      adlaplace::joint_log_dens(notape$ptr, x, negative = FALSE),
      f_ad,
      tolerance = 1e-10,
      info = label
    )
    expect_equal(
      adlaplace::joint_log_dens(taped$ptr, x, negative = FALSE),
      f_ad,
      tolerance = 1e-10,
      info = label
    )

    for (inner in c(FALSE, TRUE)) {
      g_ad <- as.numeric(adlaplace::grad(ad$ptr, x, inner = inner, negative = FALSE))
      expect_equal(
        as.numeric(adlaplace::grad(notape$ptr, x, inner = inner, negative = FALSE)),
        g_ad,
        tolerance = 1e-10,
        info = paste(label, "inner =", inner)
      )
      expect_equal(
        as.numeric(adlaplace::grad(taped$ptr, x, inner = inner, negative = FALSE)),
        g_ad,
        tolerance = 1e-10,
        info = paste(label, "inner =", inner)
      )
      H_ad <- as.matrix(adlaplace::hessian(ad$ptr, x, inner = inner, negative = FALSE))
      expect_equal(
        as.matrix(adlaplace::hessian(notape$ptr, x, inner = inner, negative = FALSE)),
        H_ad,
        tolerance = 1e-10,
        info = paste(label, "inner =", inner)
      )
      expect_equal(
        as.matrix(adlaplace::hessian(taped$ptr, x, inner = inner, negative = FALSE)),
        H_ad,
        tolerance = 1e-10,
        info = paste(label, "inner =", inner)
      )
    }
  }
})

test_that("analytic ssq Hessian gamma block is exactly -Q", {
  # d2/dgamma2 of -0.5 gamma' Q gamma is -Q and does not depend on gamma, so
  # the inner Hessian must reproduce the assembled precision directly.
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  g <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(g, knots_list, degree = 2L)
  nr <- nrow(fem$C)

  kappa <- 1.5
  range <- sqrt(8 * 1) / kappa
  sd <- 0.8
  tau <- 1 / (kappa * sd * sqrt(4 * pi))

  pk <- ssq_pack(2L, TRUE, analytic = TRUE, tape = FALSE, seed_theta = c(range, sd))
  set.seed(9)
  x <- c(rnorm(nr), log(range), log(sd))

  H <- as.matrix(adlaplace::hessian(pk$ptr, x, inner = TRUE, negative = FALSE))
  Q <- as.matrix(fem_precision(
    kappa = kappa, tau = tau, fem$C, fem$G, fem$G2, alpha = 2L
  ))
  expect_equal(H[seq_len(nr), seq_len(nr)], -Q, tolerance = 1e-10)
})

test_that("all three ssq modes agree on trace_hinv_t and vanish on gamma", {
  for (row in seq_len(nrow(ssq_grid))) {
    alpha <- ssq_grid$alpha[row]
    transform <- ssq_grid$transform[row]
    label <- sprintf("alpha=%d transform=%s", alpha, transform)

    notape <- ssq_pack(alpha, transform, analytic = TRUE, tape = FALSE)
    taped <- ssq_pack(alpha, transform, analytic = TRUE, tape = TRUE)
    ad <- ssq_pack(alpha, transform, analytic = FALSE)
    af_notape <- adlaplace::ad_pack(notape$ptr, num_threads = 1L)
    af_taped <- adlaplace::ad_pack(taped$ptr, num_threads = 1L)
    af_ad <- adlaplace::ad_pack(ad$ptr, num_threads = 1L)

    set.seed(11)
    x <- c(rnorm(notape$n_gamma), notape$theta)

    half <- af_notape@chol_inner_list$half_H_inv
    half <- methods::as(methods::as(half, "generalMatrix"), "CsparseMatrix")
    set.seed(5)
    half@x[] <- rnorm(length(half@x))
    cols <- af_notape@chol_inner_list$trace_columns

    t_notape <- adlaplace::trace_hinv_t(af_notape, x, half, cols)
    t_taped <- adlaplace::trace_hinv_t(af_taped, x, half, cols)
    t_ad <- adlaplace::trace_hinv_t(af_ad, x, half, cols)

    scale <- max(abs(t_ad))
    expect_lt(max(abs(t_notape - t_ad)) / scale, 1e-10)
    expect_lt(max(abs(t_taped - t_ad)) / scale, 1e-10)
    expect_identical(
      t_notape[seq_len(notape$n_gamma)],
      rep(0, notape$n_gamma),
      info = label
    )
    expect_lt(max(abs(t_ad[seq_len(notape$n_gamma)])) / scale, 1e-10)
  }
})

test_that("no-tape ssq pack clones and evaluates under OpenMP", {
  pk <- ssq_pack(2L, TRUE, analytic = TRUE, tape = FALSE)
  set.seed(4)
  gamma <- rnorm(pk$n_gamma)
  parameters <- pk$theta
  dens_ptr <- adlaplace::clone_ad_pack_ptr(pk$ptr)
  af <- adlaplace::ad_pack(pk$ptr, num_threads = 2L)
  fo <- adlaplace::fun_obj_fdfh(af, parameters, gamma, inner = TRUE)
  expect_true(is.finite(fo$f))
  expect_true(all(is.finite(fo$grad)))

  x <- c(gamma, parameters)
  expect_equal(
    fo$f,
    adlaplace::joint_log_dens(dens_ptr, x, negative = TRUE),
    tolerance = 1e-10
  )
})

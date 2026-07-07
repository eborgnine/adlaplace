simulate_fit_data <- function(n = 400L, n_re = 8L, seed = 1L) {
  set.seed(seed)
  g <- sample(n_re, n, replace = TRUE)
  re <- rnorm(n_re, sd = 0.3)
  x <- rbinom(n, 1, 0.5)
  eta <- 1 + 0.5 * x + re[g]
  nb_sd <- 0.3
  z <- rgamma(n, nb_sd^(-2), nb_sd^(-2))
  data.frame(
    y = rpois(n, exp(eta) * z),
    x = x,
    g = factor(g)
  )
}

test_that("adlaplace() fits a nbinom GLMM and methods work", {
  dat <- simulate_fit_data()
  fit <- adlaplace(
    nbinom(y, lower = 1e-9, init = 0.15) ~ x + iid(g, init = 0.3),
    data = dat,
    num_shards = 10L,
    control = list(maxit = 100L),
    verbose = FALSE
  )

  expect_s3_class(fit, "adlaplace_fit")
  expect_true(is.finite(fit$logLik))
  expect_identical(fit$optim$convergence, 0L)

  co <- coef(fit)
  expect_true(all(is.finite(co)))
  expect_equal(length(co), nrow(fit$par_info))
  # transformed parameters reported on the natural (positive) scale
  expect_true(all(co[fit$par_info$transform %in% TRUE] > 0))
  # rough recovery of the fixed effects (matched by label, order-free)
  beta_labels <- fit$par_info$label[
    seq_len(nrow(fit$model_data$data$info$beta))
  ]
  int_label <- grep("intercept", beta_labels, value = TRUE)[1]
  x_label <- setdiff(beta_labels, int_label)[1]
  expect_equal(unname(co[int_label]), 1, tolerance = 0.35)
  expect_equal(unname(co[x_label]), 0.5, tolerance = 0.35)

  vc <- vcov(fit)
  expect_true(isSymmetric(vc, tol = 1e-6))
  expect_true(all(diag(vc) > 0))

  ci <- confint(fit)
  expect_equal(dim(ci), c(length(co), 2L))
  expect_true(all(ci[, 1] < ci[, 2]))
  # estimates inside their own intervals
  expect_true(all(co > ci[, 1] & co < ci[, 2]))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "nobs"), nrow(dat))
  expect_equal(nobs(fit), nrow(dat))

  expect_equal(length(fitted(fit)), nrow(dat))
  expect_equal(length(fit$gamma), nlevels(dat$g))

  sm <- summary(fit)
  expect_s3_class(sm, "summary.adlaplace_fit")
  expect_true(all(is.finite(sm$fixed[, "Std. Error"])))
  expect_output(print(sm), "Fixed effects")
  expect_output(print(fit), "Coefficients")
})

test_that("adlaplace() with hessian = FALSE has no vcov", {
  dat <- simulate_fit_data(n = 200L, seed = 2L)
  fit <- adlaplace(
    nbinom(y, lower = 1e-9, init = 0.15) ~ x + iid(g, init = 0.3),
    data = dat,
    num_shards = 5L,
    hessian = FALSE,
    control = list(maxit = 50L)
  )
  expect_null(fit$vcov)
  expect_error(vcov(fit), "hessian = TRUE")
  expect_true(all(is.na(confint(fit))))
})

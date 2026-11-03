test_that("model_data captures y from plain formula LHS", {
  dat <- data.frame(y = 1:5, x1 = 1:5)
  expect_no_warning({
    md <- adlaplace::model_data(y ~ linear(x1), data = dat)
  })
  expect_equal(md$term_data$y, as.numeric(dat$y))
  # bare response defaults to a gaussian observation term
  expect_true(length(md$observations) == 1L)
  expect_identical(md$observations[[1]]@density, "gaussian_obs")
  expect_identical(names(md$parameters), "y_extra")
})

test_that("model_data captures y from explicit observations term", {
  skip_if_not_installed("mgcv")
  dat <- quiet_gamsim(6, n = 80, scale = 0.2, dist = "poisson", verbose = FALSE)

  md <- adlaplace::model_data(
    adlaplace::nbinom(y, lower = 1e-9) ~
      x1 +
      adlaplace::iwp(
        x2,
        p = 2,
        knots = seq(0, 1, len = 11),
        ref_value = 0.5
      ) +
      adlaplace::iid(fac),
    data = dat,
    verbose = FALSE
  )

  expect_equal(length(md$term_data$y), nrow(dat))
  expect_true(length(md$observations) >= 1L)
})

test_that("parameters label uses beta_label for fixed effects", {
  dat <- data.frame(
    y = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L),
    trt = factor(
      rep(c("placebo", "drug", "drug+"), each = 4L),
      levels = c("placebo", "drug", "drug+")
    ),
    wk = rep(c(0L, 1L), each = 6L),
    ID = factor(rep(paste0("X", 1:4), 3L))
  )
  setup <- adlaplace::term_data_setup(
    adlaplace::binom(y) ~ trt + wk + adlaplace::iid(ID),
    data = dat
  )
  info <- setup$info
  n_beta <- nrow(info$beta)
  n_theta <- nrow(info$theta)

  expect_identical(
    info$parameters$label[seq_len(n_beta)],
    info$beta$beta_label
  )
  expect_identical(anyDuplicated(info$parameters$label), 0L)
  expect_identical(
    info$parameters$label[seq(n_beta + 1L, length.out = n_theta)],
    info$theta$label
  )
  expect_true(all(info$beta$label[info$beta$term == "trt"] == "trt_linear"))
  expect_identical(
    info$beta$beta_label[info$beta$term == "trt"],
    c("trt_linear_trtdrug", "trt_linear_trtdrug+")
  )
})

test_that("term_data_setup warns when no response/observations term is present", {
  dat <- data.frame(y = 1:5, x1 = 1:5)
  expect_warning(
    adlaplace:::term_data_setup(list(adlaplace::linear("x1")), data = dat),
    "no response variable"
  )
})

test_that("model_data omits random shard when rpoly precision is NULL (sd = Inf)", {
  dat <- data.frame(y = rnorm(20), x = runif(20))
  md <- adlaplace::model_data(
    y ~ adlaplace::rpoly(x, p = 2, sd = Inf),
    data = dat,
    verbose = FALSE
  )
  expect_length(md$random, 0L)
})

test_that("model_data na_omit drops rows with NA covariates", {
  dat <- data.frame(y = c(1, 2, 3), x1 = c(1, NA, 3))
  md <- adlaplace::model_data(
    y ~ adlaplace::linear(x1),
    data = dat,
    na_omit = TRUE
  )
  expect_equal(nrow(md$term_data$data), 2L)
})

test_that("model_data na_omit=FALSE keeps all rows", {
  dat <- data.frame(y = 1:3, x1 = 1:3, unused = c(1, NA, 3))
  md <- adlaplace::model_data(y ~ adlaplace::linear(x1), data = dat, na_omit = FALSE)
  expect_equal(nrow(md$term_data$data), 3L)
  md2 <- adlaplace::model_data(y ~ adlaplace::linear(x1), data = dat, na_omit = TRUE)
  expect_equal(nrow(md2$term_data$data), 3L)
})

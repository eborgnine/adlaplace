test_that("hgp term constructors produce design and precision matrices", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(
    x = seq(0, 1, length.out = 20),
    mult = rnorm(20),
    g = rep(c("a", "b"), each = 10)
  )
  knots <- c(0, 0.25, 0.5, 0.75, 1)

  hiwp_terms <- hiwp(x = "x", knots = knots, by = "g", include_global = FALSE)
  hiwp_term <- hiwp_terms[[1]]
  expect_s4_class(hiwp_term, "hiwp")
  A_hiwp <- design(hiwp_term, data)
  P_hiwp <- precision(hiwp_term, data)
  expect_gt(ncol(A_hiwp), 0L)
  expect_gt(nrow(P_hiwp), 0L)

  rsiwp_terms <- rsiwp(x = "x", mult = "mult", knots = knots, include_linear = FALSE)
  rsiwp_term <- rsiwp_terms[[1]]
  expect_s4_class(rsiwp_term, "rsiwp")
  expect_equal(rsiwp_term@density, "random_diagonal")
  expect_equal(rsiwp_term@ad_kind, "random")
  expect_equal(rsiwp_term@label, "x_mult_rsiwp")
  A_rsiwp <- design(rsiwp_term, data)
  P_rsiwp <- precision(rsiwp_term, data)
  expect_s4_class(A_rsiwp, "dgCMatrix")
  expect_gt(ncol(A_rsiwp), 0L)
  expect_gt(nrow(P_rsiwp), 0L)

  gamma_info <- random_info(rsiwp_term, data)
  theta_row <- theta_info(rsiwp_term)
  expect_equal(unique(gamma_info$label), rsiwp_term@label)
  expect_equal(theta_row$label, rsiwp_term@label)

  rsiid_terms <- rsiid(x = "x", mult = "mult")
  rsiid_term <- rsiid_terms[[1]]
  expect_s4_class(rsiid_term, "rsiid")
  expect_equal(rsiid_term@density, "random_diagonal")
  expect_equal(rsiid_term@ad_kind, "random")
  expect_equal(rsiid_term@label, "x_mult_rsiid")
  A_rsiid <- design(rsiid_term, data)
  P_rsiid <- precision(rsiid_term, data)
  expect_gt(ncol(A_rsiid), 0L)
  expect_gt(nrow(P_rsiid), 0L)

  rsiid_gamma <- random_info(rsiid_term, data)
  rsiid_theta <- theta_info(rsiid_term)
  expect_equal(unique(rsiid_gamma$label), rsiid_term@label)
  expect_equal(rsiid_theta$label, rsiid_term@label)
  expect_equal(rsiid_gamma$gamma_label, colnames(A_rsiid))

  rsfpoly_term <- rsfpoly(x = "x", mult = "mult", p = 2)
  expect_s4_class(rsfpoly_term, "rsfpoly")
  X_rsfpoly <- design(rsfpoly_term, data)
  expect_gt(ncol(X_rsfpoly), 0L)
  expect_null(precision(rsfpoly_term, data))

  rsrpoly_term <- rsrpoly(x = "x", mult = "mult", p = 2, sd = 1)
  expect_s4_class(rsrpoly_term, "rsrpoly")
  A_rsrpoly <- design(rsrpoly_term, data)
  P_rsrpoly <- precision(rsrpoly_term, data)
  expect_gt(ncol(A_rsrpoly), 0L)
  expect_gt(nrow(P_rsrpoly), 0L)
})

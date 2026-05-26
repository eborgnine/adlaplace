test_that("f() and constructor formulas yield identical collect_terms", {
  skip_if_not_installed("mgcv")
  dat <- mgcv::gamSim(6, n = 80, scale = 0.2, dist = "poisson")

  with(dat, {
    f_style <-
      adlaplace::f(y, model = "nbinom", lower = 1e-9) ~
      x1 +
      adlaplace::f(
        x2,
        model = "iwp",
        p = 2,
        knots = seq(0, 1, len = 11),
        ref_value = 0.5
      ) +
      adlaplace::f(fac, model = "iid")

    ctor_style <-
      adlaplace::nbinom(y, lower = 1e-9) ~
      x1 +
      adlaplace::iwp(
        x2,
        p = 2,
        knots = seq(0, 1, len = 11),
        ref_value = 0.5
      ) +
      adlaplace::iid(fac)

    t1 <- adlaplace::collect_terms(f_style, verbose = FALSE)
    t2 <- adlaplace::collect_terms(ctor_style, verbose = FALSE)
    expect_identical(t1, t2)
    expect_false("response" %in% names(t1))
    expect_false("response" %in% names(t2))

    nbinom_idx <- which(vapply(t1, function(x) {
      methods::is(x, "nbinom")
    }, logical(1L)))
    expect_length(nbinom_idx, 1L)
    expect_equal(t1[[nbinom_idx]]@ad_fun, "nbinom_obs")
    expect_equal(t1[[nbinom_idx]]@lower, 1e-9)
  })
})

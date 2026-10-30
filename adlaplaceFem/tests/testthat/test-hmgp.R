test_that("hmgp constructor returns global terms by default", {
  knots <- seq(1, 5, length.out = 8)
  full <- hmgp("x", a = 2, knots = knots, by = "g", levels = c("A", "B"))
  expect_equal(
    names(full),
    c("x_mgp", "x_mgp_boundary", "x_hmgp", "x_hslope_mgp")
  )
  dev <- hmgp("x", a = 2, knots = knots, by = "g", include_global = FALSE)
  expect_equal(names(dev), c("x_hmgp", "x_hslope_mgp"))
  expect_s4_class(full$x_hmgp, "mgp")
  expect_equal(full$x_hmgp@by@levels, c("A", "B"))
  expect_error(hmgp("x", knots = 1, by = "g"), "knots")
  expect_error(hmgp("x", knots = knots), "by")
})

test_that("hmgp design and precision are the mgp blocks repeated by group", {
  knots <- seq(1, 5, length.out = 8)
  terms <- hmgp(
    "x", a = 2, knots = knots, ref_value = "left",
    by = "g", levels = c("A", "B"), include_global = TRUE
  )
  dat <- data.frame(
    x = c(1.2, 2, 3, 4.5),
    g = c("A", "A", "B", "B")
  )
  A <- design(terms$x_hmgp, dat)
  A_mgp <- design(terms$x_mgp, dat)
  n <- ncol(A_mgp)
  expect_equal(ncol(A), 2L * n)
  expect_equal(
    unname(as.matrix(A[1:2, seq_len(n)])),
    unname(as.matrix(A_mgp[1:2, , drop = FALSE]))
  )
  expect_equal(unname(as.matrix(A[1:2, n + seq_len(n)])), matrix(0, 2L, n))
  expect_equal(
    unname(as.matrix(A[3:4, n + seq_len(n)])),
    unname(as.matrix(A_mgp[3:4, , drop = FALSE]))
  )

  prec <- precision(terms$x_hmgp, dat)
  q0 <- precision(terms$x_mgp, dat)$Q
  expect_equal(as.matrix(prec$Q), as.matrix(Matrix::bdiag(q0, q0)))
  expect_equal(prec$log_det, 2 * precision(terms$x_mgp, dat)$log_det)
  expect_equal(prec$rank, 2L * precision(terms$x_mgp, dat)$rank)
  expect_identical(random_info(terms$x_hmgp, dat)$gamma_label, colnames(A))

  Ab <- design(terms$x_hslope_mgp, dat)
  mb <- as.numeric(design(terms$x_mgp_boundary, dat))
  expect_equal(as.numeric(Ab[, 1]), c(mb[1:2], 0, 0))
  expect_equal(as.numeric(Ab[, 2]), c(0, 0, mb[3:4]))
  expect_identical(colnames(precision(terms$x_hslope_mgp, dat)), colnames(Ab))
  expect_identical(random_info(terms$x_hslope_mgp, dat)$gamma_label, colnames(Ab))
  expect_equal(nrow(theta_info(terms$x_hmgp)), 1L)
  expect_equal(nrow(theta_info(terms$x_hslope_mgp)), 1L)
  expect_equal(theta_info(terms$x_hslope_mgp)$label, "x_hslope_mgp")
})

test_that("hmgp theta labels do not collide under grep", {
  knots <- seq(1, 5, length.out = 8)
  terms <- hmgp("pm", a = 2, knots = knots, by = "region", levels = "A")
  labs <- unlist(lapply(terms, function(t) {
    info <- theta_info(t)
    if (is.null(info)) character(0) else info$label
  }))
  for (t in terms) {
    info <- theta_info(t)
    if (is.null(info) || !nrow(info)) next
    expect_equal(length(grep(t@label, labs)), 1L)
  }
})

test_that("hmgp inherits psd and psd_prior targets the deviation term", {
  knots <- seq(1, 5, length.out = 8)
  terms <- hmgp("x", a = 2, knots = knots, by = "g", include_global = TRUE)
  expect_equal(
    psd(terms$x_hmgp, h = 0.4, x = 2),
    psd(terms$x_mgp, h = 0.4, x = 2)
  )
  pr <- psd_prior(terms$x_hmgp, u = 1, alpha = 0.5, h = 0.5, x = 2)
  expect_equal(pr@theta, "x_hmgp")
})

test_that("hmgp rejects unknown or missing groups and keeps levels order", {
  knots <- seq(1, 5, length.out = 8)
  term <- hmgp(
    "x", a = 2, knots = knots, by = "g", levels = c("B", "A"),
    include_global = FALSE
  )$x_hmgp
  dat <- data.frame(x = c(2, 3), g = c("A", "B"))
  cols <- colnames(design(term, dat))
  expect_true(grepl("_gB$", cols[1]))
  expect_error(design(term, data.frame(x = 2, g = "C")), "levels not in the term")
  expect_error(design(term, data.frame(x = 2, g = NA)), "missing values")
})

test_that("replace_vars updates hmgp covariate and grouping names", {
  knots <- seq(1, 5, length.out = 8)
  terms <- hmgp(
    "pm", a = 2, knots = knots, by = "region", levels = c("A", "B"),
    include_global = FALSE
  )
  out <- replace_vars(terms, c(pm = "pm2", region = "reg"))
  expect_equal(out$pm2_hmgp@name, "pm2")
  expect_equal(out$pm2_hmgp@by@term, "reg")
  expect_equal(out$pm2_hmgp@label, "pm2_hmgp")
  expect_equal(out$pm2_hslope_mgp@label, "pm2_hslope_mgp")
  dat <- data.frame(pm2 = c(1.5, 3), reg = c("A", "B"))
  expect_true(all(grepl("^pm2_hmgp_k", colnames(design(out$pm2_hmgp, dat)))))
  expect_true(all(grepl("_g(A|B)$", colnames(design(out$pm2_hmgp, dat)))))
})

test_that("four-region hmgp design has one block per region", {
  knots <- seq(8, 40, by = 4)
  regions <- c("BC", "Prairies", "Ontario", "Atlantic")
  term <- hmgp(
    "pm", a = 2, knots = knots, by = "region", levels = regions,
    include_global = FALSE
  )$pm_hmgp
  toy <- data.frame(
    pm = seq(8, 40, length.out = 12),
    region = rep(regions, length.out = 12)
  )
  A <- design(term, toy)
  expect_equal(ncol(A), 4L * term@fem$rank)
  expect_true(all(endsWith(colnames(A), paste0("_g", rep(regions, each = term@fem$rank)))))
})

test_that("two-group hmgp case-crossover fit is finite", {
  skip_on_cran()
  set.seed(3)
  knots <- seq(1, 5, length.out = 6)
  x <- seq(1.2, 4.8, length.out = 6)
  dat <- data.frame(
    y = rpois(12, 2),
    x = rep(x, 2),
    g = rep(c("A", "B"), each = 6),
    stratum = rep(1:4, each = 3)
  )
  fit <- adlaplace(
    dirichlet_multinom(y, by = "stratum", init = 0.05) ~ 0 +
      hmgp(x, a = 2, knots = knots, ref_value = "left", by = "g", levels = c("A", "B")),
    data = dat,
    config = list(num_threads = 1L, num_shards = 1L),
    control = list(maxit = 5L),
    verbose = FALSE,
    hessian = TRUE
  )
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  expect_true(all(c("x_mgp", "x_hmgp", "x_hslope_mgp") %in% fit$par_info$label))
  gnames <- names(fit$gamma)
  expect_true(any(grepl("^x_hmgp_k", gnames)))
  expect_true(any(grepl("^x_hslope_mgp_g", gnames)))
  expect_false(is.null(fit$details$hessian$half_H_inv))
})

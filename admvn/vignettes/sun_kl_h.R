## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
library(admvn)

## ----true-params--------------------------------------------------------------
par_true <- c(
  xi1 = 0.2, xi2 = -0.1, xi3 = 0.05,
  nu1 = 1.0, nu2 = 0.8, nu3 = 1.2,
  z21 = 0.15, z31 = 0.05, z32 = 0.10,
  z41 = atanh(0.5), z42 = 0, z43 = 0,
  z51 = 0, z52 = atanh(0.45), z53 = 0, z54 = 0,
  z61 = 0, z62 = 0, z63 = atanh(0.4), z64 = 0.1, z65 = 0.05
)
dp_true <- make_sun_hs_params(par_true)
str(dp_true, max.level = 1)
diag(dp_true$Gamma)
diag(dp_true$Delta)

## ----target, eval = requireNamespace("sn", quietly = TRUE)--------------------
library(sn)
set.seed(1)
n <- 600L
x <- sn::rsun(n, dp = dp_true)
w <- rep(1 / n, n)
logp_true <- as.numeric(sn::dsun(x, dp = dp_true, log = TRUE))
quad <- list(x = x, w = w, logp = logp_true)

kl_at_true <- evaluate_sun33_hs_slice(par_true, quad)
kl_at_true

## ----start, eval = requireNamespace("sn", quietly = TRUE)---------------------
par_start <- make_sun33_hs_start_from_normal(
  mu = colMeans(x),
  Sigma = stats::cov(x),
  skew_strength = 0.05
)
kl_at_start <- evaluate_sun33_hs_slice(par_start, quad)
kl_at_start

sel <- c(
  "xi1", "xi2", "xi3", "nu1", "nu2", "nu3",
  "z21", "z31", "z32", "z41", "z52", "z63", "z64", "z65"
)
rbind(true = par_true[sel], start = par_start[sel])

## ----fit, eval = requireNamespace("sn", quietly = TRUE)-----------------------
fit <- fit_sun33_hs_quad(
  quad = quad,
  start = par_start,
  mc.cores = 1L,
  obj_scale = 1e4,
  mvquad_opts = list(n_points = 64L, n_shifts = 2L),
  optim_opts = list(maxit = 250L, factr = 1e3, pgtol = 1e-8)
)
print_sun_kl_fit(fit, label = "SUN(3,3) hs KL match")

## ----recovery, eval = requireNamespace("sn", quietly = TRUE)------------------
par_hat <- fit$par
names(par_hat) <- names(par_true)
dp_hat <- fit$dp
kl_at_hat <- evaluate_sun33_hs_slice(par_hat, quad)

data.frame(
  kl_true = unname(kl_at_true[["kl"]]),
  kl_start = unname(kl_at_start[["kl"]]),
  kl_hat = unname(kl_at_hat[["kl"]])
)

# Hold-out mean log-density under true vs fitted dp
set.seed(2)
x_hold <- sn::rsun(400L, dp = dp_true)
ll_true <- mean(sn::dsun(x_hold, dp = dp_true, log = TRUE))
ll_hat <- mean(sn::dsun(x_hold, dp = dp_hat, log = TRUE))
ll_start <- mean(sn::dsun(x_hold, dp = make_sun_hs_params(par_start), log = TRUE))
c(holdout_ll_true = ll_true, holdout_ll_start = ll_start, holdout_ll_hat = ll_hat)

cmp <- data.frame(
  true = par_true,
  start = as.numeric(par_start),
  hat = as.numeric(par_hat),
  abs_err = abs(as.numeric(par_hat) - par_true)
)
print(cmp[sel, ], digits = 4)

list(
  xi = rbind(true = dp_true$xi, hat = dp_hat$xi),
  Omega_diag = rbind(true = diag(dp_true$Omega), hat = diag(dp_hat$Omega)),
  Delta_diag = rbind(true = diag(dp_true$Delta), hat = diag(dp_hat$Delta)),
  max_abs_Omega = max(abs(dp_hat$Omega - dp_true$Omega)),
  max_abs_Delta = max(abs(dp_hat$Delta - dp_true$Delta)),
  max_abs_Gamma = max(abs(dp_hat$Gamma - dp_true$Gamma)),
  max_abs_par = max(abs(par_hat - par_true))
)

## ----fig-recovery, eval = requireNamespace("sn", quietly = TRUE), fig.height = 4----
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

logq_hat <- as.numeric(sn::dsun(x, dp = dp_hat, log = TRUE))
plot(
  logp_true, logq_hat,
  xlab = "true log density", ylab = "fitted log density",
  main = "Density match on fitting sample"
)
abline(0, 1, col = "red", lty = 2)

plot(
  as.vector(dp_true$Omega), as.vector(dp_hat$Omega),
  xlab = "true Omega", ylab = "fitted Omega",
  main = expression(Omega ~ entries)
)
abline(0, 1, col = "red", lty = 2)

## ----kl-sun33-hs-profile-setup, eval = requireNamespace("sn", quietly = TRUE)----
par0 <- as.numeric(fit$par)
names(par0) <- names(par_true)
w_fit <- as.numeric(fit$quad$w)
w_fit <- w_fit / sum(w_fit)

tape33 <- dsun_hs_fun(
  x = as.matrix(fit$quad$x),
  par_seed = par0,
  weights = w_fit,
  n_points = 64L,
  n_shifts = 2L,
  n_threads = 1L
)

## dsun_hs_fun value = sum w log q; objective = -sum w log q (= KL + const)
obj_ad <- function(p) {
  v <- tape33$eval(p, deriv = 0L)$value
  if (!is.finite(v)) {
    stop("non-finite AD objective at profile point", call. = FALSE)
  }
  -as.numeric(v)
}
grad_ad <- function(p) {
  g <- as.numeric(tape33$eval(p, deriv = 1L)$gradient)
  if (any(!is.finite(g))) {
    stop("non-finite AD gradient at profile point", call. = FALSE)
  }
  -g
}

bnd <- sun33_hs_bounds()
col_kl <- "steelblue"

profile_seq <- function(j, n = 11L, half_width = 0.5) {
  center <- par0[j]
  lo <- bnd$lower[j]
  hi <- bnd$upper[j]
  hw <- half_width * max(abs(center), 1)
  lo_j <- center - hw
  hi_j <- center + hw
  if (is.finite(lo)) lo_j <- max(lo_j, lo)
  if (is.finite(hi)) hi_j <- min(hi_j, hi)
  if (!(is.finite(lo_j) && is.finite(hi_j) && lo_j < hi_j)) {
    stop(
      "Cannot build an interior profile sequence for parameter ", j,
      " (", names(par0)[j], ").",
      call. = FALSE
    )
  }
  seq(lo_j, hi_j, length.out = n)
}

## Plot obj + FD/AD gradient profiles for a block of parameter indices.
profile_par_block <- function(idx, n_col = 4L) {
  idx <- as.integer(idx)
  n_panels <- 2L * length(idx)
  n_row <- as.integer(ceiling(n_panels / n_col))
  op <- par(
    mfrow = c(n_row, n_col),
    mar = c(3, 3, 1.2, 0.4),
    mgp = c(1.6, 0.5, 0)
  )
  on.exit(par(op), add = TRUE)
  for (j in idx) {
    lab <- names(par0)[j]
    xj_seq <- profile_seq(j)
    obj_seq <- vapply(xj_seq, function(xj) {
      p <- par0
      p[j] <- xj
      obj_ad(p)
    }, numeric(1))
    ad_dobj <- vapply(xj_seq, function(xj) {
      p <- par0
      p[j] <- xj
      grad_ad(p)[j]
    }, numeric(1))
    ## Central FD at the same nodes as AD (Richardson on two step sizes).
    lo <- bnd$lower[j]
    hi <- bnd$upper[j]
    fd_at <- function(xj, h) {
      x_lo <- xj - h
      x_hi <- xj + h
      if (is.finite(lo)) x_lo <- max(x_lo, lo)
      if (is.finite(hi)) x_hi <- min(x_hi, hi)
      if (!(x_hi > x_lo)) {
        stop("central FD step collapsed for ", lab, " at ", format(xj))
      }
      p_lo <- p_hi <- par0
      p_lo[j] <- x_lo
      p_hi[j] <- x_hi
      (obj_ad(p_hi) - obj_ad(p_lo)) / (x_hi - x_lo)
    }
    fd_dobj <- vapply(xj_seq, function(xj) {
      h1 <- 1e-4 * max(abs(xj), 1)
      h2 <- 2 * h1
      fd1 <- fd_at(xj, h1)
      fd2 <- fd_at(xj, h2)
      (4 * fd1 - fd2) / 3
    }, numeric(1))

    obj_finite <- obj_seq[is.finite(obj_seq)]
    if (length(obj_finite)) {
      ylim_lo <- min(obj_finite)
      ylim_hi <- min(10, max(obj_finite))
      if (!(ylim_hi > ylim_lo)) {
        ylim_hi <- ylim_lo + 1e-8
      }
    } else {
      ylim_lo <- 0
      ylim_hi <- 10
    }
    plot(
      xj_seq, obj_seq,
      type = "b", pch = 16, col = col_kl, cex = 0.7,
      xlab = lab, ylab = expression(-sum(w * log(q))),
      ylim = c(ylim_lo, ylim_hi)
    )
    abline(v = par0[j], lty = 3, col = "grey60")
    plot(
      xj_seq, fd_dobj,
      type = "b", pch = 4, col = "grey30", cex = 1.4, lwd = 1.5,
      xlab = lab, ylab = paste0("d obj / d ", lab)
    )
    lines(xj_seq, ad_dobj, type = "b", pch = 3, col = col_kl, cex = 0.9, lwd = 1)
    abline(h = 0, lty = 3, col = "grey60")
    abline(v = par0[j], lty = 3, col = "grey60")
    legend(
      "topright",
      legend = c("FD", "AD"),
      col = c("grey30", col_kl), pch = c(4, 3),
      pt.cex = c(1.4, 0.9), bty = "n", cex = 0.65
    )
  }
}

## ----kl-sun33-hs-profile-1, eval = requireNamespace("sn", quietly = TRUE), fig.width = 7, fig.height = 9, fig.cap = "SUN(3,3) hs slice profiles for parameters 1-7 (cross-entropy objective; FD vs AD gradient)."----
profile_par_block(1:7)

## ----kl-sun33-hs-profile-2, eval = requireNamespace("sn", quietly = TRUE), fig.width = 7, fig.height = 9, fig.cap = "SUN(3,3) hs slice profiles for parameters 8-14."----
profile_par_block(8:14)

## ----kl-sun33-hs-profile-3, eval = requireNamespace("sn", quietly = TRUE), fig.width = 7, fig.height = 9, fig.cap = "SUN(3,3) hs slice profiles for parameters 15-21."----
profile_par_block(15:21)


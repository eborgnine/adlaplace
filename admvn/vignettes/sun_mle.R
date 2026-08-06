## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
library(admvn)

## ----params-------------------------------------------------------------------
true_params <- c(
  xi1 = 0, xi2 = 0, xi3 = 0,
  nu1 = 1, nu2 = 1, nu3 = 1,
  omega12 = 0.2, omega13 = 0.1, omega23 = 0.2,
  L11 = 0.4, L12 = 0.1, L13 = 0.1,
  L21 = 0.1, L22 = 0.4, L23 = 0.1,
  L31 = 0, L32 = 0.1, L33 = 0.35,
  gamma12 = 0.3, gamma13 = 0.2, gamma23 = 0.3
)
true_sun_params <- make_sun_params(true_params)
str(true_sun_params, max.level = 1)
# Sigma_vv = Gamma has unit diagonal by construction
diag(true_sun_params$Gamma)

## ----sim, eval = requireNamespace("sn", quietly = TRUE)-----------------------
library(sn)
set.seed(0)
n <- 250
ysim <- sn::rsun(n, dp = true_sun_params)
colnames(ysim) <- paste0("Y", 1:3)
head(ysim)

## ----tape, eval = requireNamespace("sn", quietly = TRUE)----------------------
# Unused for d<=3 specialized CDF values; kept for API compatibility.
qmc_n_points <- 64L
qmc_n_shifts <- 2L

build_time <- system.time({
  tape <- dsun_fun(
    ysim,
    true_params,
    n_points = qmc_n_points,
    n_shifts = qmc_n_shifts,
    n_threads = 2L
  )
})["elapsed"]

check <- tape$eval(true_params, log = TRUE, deriv = 1L)
sn_ll <- sum(sn::dsun(ysim, dp = true_sun_params, log = TRUE))
c(
  sn = sn_ll,
  admvn = check$value,
  build_sec = unname(build_time),
  n_threads = tape$n_threads
)

## ----mle-setup, eval = requireNamespace("sn", quietly = TRUE)-----------------
start_params <- true_params
start_params[] <- c(
  xi1 = 0.2, xi2 = -0.1, xi3 = 0.1,
  nu1 = 1.4, nu2 = 0.7, nu3 = 1.2,
  omega12 = 0.0, omega13 = 0.3, omega23 = 0.0,
  L11 = 0.7, L12 = 0.2, L13 = 0.1,
  L21 = 0.2, L22 = 0.8, L23 = 0.1,
  L31 = 0.1, L32 = 0.3, L33 = 0.9,
  gamma12 = 0.1, gamma13 = 0.0, gamma23 = 0.2
)

bnd <- sun33_bounds()
lower <- bnd$lower
upper <- bnd$upper
# Keep nu floors a bit above the package default for this demo.
lower[4:6] <- pmax(lower[4:6], 1e-3)
# Rescale free parameters so L-BFGS-B does not take wild Omega steps.
parscale <- pmax(abs(as.numeric(start_params)), 0.1)
optim_ctrl <- list(
  fnscale = -1,
  maxit = 40,
  parscale = parscale
)

obj <- tape$optim_fns()

## ----mle-admvn, eval = requireNamespace("sn", quietly = TRUE)-----------------
fit_ad_time <- system.time({
  fit_ad <- optim(
    start_params,
    fn = obj$fn,
    gr = obj$gr,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = optim_ctrl
  )
})

fit_ad$convergence
fit_ad$value
unname(fit_ad_time["elapsed"])

## ----mle-trust, eval = requireNamespace("sn", quietly = TRUE) && requireNamespace("trustOptim", quietly = TRUE)----
fit_trust_time <- system.time({
  fit_trust <- sun_mle(
    tape,
    start_params,
    control = list(maxit = 40L, report.freq = 0L, report.level = 0L)
  )
})
fit_trust$value
fit_trust$iterations
unname(fit_trust_time["elapsed"])

## ----compare, eval = requireNamespace("sn", quietly = TRUE)-------------------
sn_ll_at <- function(par) {
  sum(sn::dsun(ysim, dp = make_sun_params(par), log = TRUE))
}
ad_ll_at <- function(par) {
  tape$eval(par, log = TRUE, deriv = 0L)$value
}

compare <- data.frame(
  method = "ad_optim",
  sn_loglik = sn_ll_at(fit_ad$par),
  admvn_loglik = ad_ll_at(fit_ad$par),
  iterations = unname(fit_ad$counts[["gradient"]]),
  elapsed_sec = unname(fit_ad_time["elapsed"]),
  stringsAsFactors = FALSE
)
if (exists("fit_trust")) {
  compare <- rbind(
    data.frame(
      method = "ad_trustOptim",
      sn_loglik = sn_ll_at(fit_trust$par),
      admvn_loglik = ad_ll_at(fit_trust$par),
      iterations = as.numeric(fit_trust$iterations),
      elapsed_sec = unname(fit_trust_time["elapsed"]),
      stringsAsFactors = FALSE
    ),
    compare
  )
}
rownames(compare) <- NULL
compare

## ----grad-check, eval = requireNamespace("sn", quietly = TRUE), fig.width = 7, fig.height = 5----
par_hat <- as.numeric(if (exists("fit_trust")) fit_trust$par else fit_ad$par)
par_names <- names(true_params)
if (is.null(par_names)) {
  par_names <- paste0("p", seq_along(par_hat))
}
names(par_hat) <- par_names

idx <- match(c("xi1", "nu1", "L11", "gamma12"), par_names)

ll <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
ad_grad_j <- function(p, j) tape$eval(p, log = TRUE, deriv = 1L)$gradient[j]
fd_grad_j <- function(p, j, h) {
  p_up <- p
  p_dn <- p
  p_up[j] <- p[j] + h
  p_dn[j] <- p[j] - h
  (ll(p_up) - ll(p_dn)) / (2 * h)
}

op <- par(mfrow = c(2, 2), mar = c(3.2, 3.2, 2.2, 0.8), oma = c(0, 0, 1.5, 0))
on.exit(par(op), add = TRUE)

for (j in idx) {
  center <- par_hat[j]
  span <- pmax(0.15 * abs(center), 0.05)
  grid <- center + seq(-span, span, length.out = 9L)
  if (par_names[j] %in% c("nu1", "nu2", "nu3")) {
    grid <- grid[grid > 1e-3]
  }

  ad_vals <- numeric(length(grid))
  fd_vals <- numeric(length(grid))
  for (k in seq_along(grid)) {
    p <- par_hat
    p[j] <- grid[k]
    h <- pmax(1e-5, 1e-4 * abs(grid[k]))
    ad_vals[k] <- ad_grad_j(p, j)
    fd_vals[k] <- fd_grad_j(p, j, h)
  }

  ylim <- range(c(ad_vals, fd_vals), finite = TRUE)
  plot(
    grid, ad_vals,
    type = "l", lwd = 2, col = "#1b4f72",
    xlab = par_names[j], ylab = expression(partialdiff * ell / partialdiff * theta),
    main = par_names[j],
    ylim = ylim
  )
  lines(grid, fd_vals, lwd = 2, col = "#b03a2e", lty = 2)
  abline(v = center, col = "gray50", lty = 3)
}
legend(
  "topleft",
  legend = c("AD", "finite diff"),
  col = c("#1b4f72", "#b03a2e"),
  lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8
)


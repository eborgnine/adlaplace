# Rscript applications/synthetic/main-sim.R > applications/synthetic/log.txt 2>&1
#
# Simulation study with hierarchical PM smooth; development likelihood checks.

if (file.exists("adlaplace/DESCRIPTION")) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all("adlaplace", quiet = TRUE)
    devtools::load_all("adlaplaceHgp", quiet = TRUE)
    devtools::load_all("hpolcc", quiet = TRUE)
  }
}
library(adlaplace)
library(adlaplaceHgp)
library(hpolcc)
library(data.table)

# Data --------------------------------------------------------------------
genPM <- function(date) {
  l <- length(date) + 1
  date_int <- as.integer(date)
  date_int <- c(date_int[1] - 1, date_int)
  pm <- 10
  pm <- pm + 10 * sin(2 * pi * (date_int + rnorm(l, 0, 1)) / 365.25)
  pm <- pm + rnorm(l, 0, 1)
  pm <- (2 * pm[-l] + pm[-1]) / 3
  pmax(0, pm)
}

genHum <- function(date) {
  l <- length(date) + 1
  date_int <- as.integer(date)
  date_int <- c(date_int[1] - 1, date_int)
  hum <- 50
  hum <- hum + 5 * cos(2 * pi * (date_int + rnorm(l, 0, 1)) / 365.25)
  hum <- hum + rnorm(l, 0, 1)
  hum <- (2 * hum[-l] + hum[-1]) / 3
  hum <- pmax(0, hum)
  hum / max(hum)
}

genPmEffect <- function(pm, r1, r2, coseffect = 1) {
  0.1 * (r1 * pm / 10 + 0.2 * pm^(r2 / 2) - coseffect * cos(pm * pi / 25))
}

genCount <- function(hum, pm, r1 = 1, r2 = 1) {
  l <- length(pm)
  odSd <- 0.3
  odPrec <- odSd^(-2)
  od <- rgamma(l, shape = odPrec, rate = odPrec)
  rpois(l, exp(-log(1) + hum + genPmEffect(pm, r1, r2) + log(od)))
}

set.seed(0)
region_effect1 <- rnorm(10, 1, 0.2)
region_effect2 <- rnorm(10, 1, 0.2)

data <- lapply(1:10, \(i) {
  data <- data.table(date = as.Date(1:1000), region = i)
  data$hum <- genHum(data$date)
  data$pm <- genPM(data$date)
  data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i])
  data
})
data <- do.call(rbind, data)
data$monthDow <- format(data$date, "%Y-%m-%a")
dataOrig <- data

ref_values <- list(pm = 10)
knots_pm <- seq(0, 22, by = 2)

cc_formula <- dirichlet_multinom(
  count,
  by = c("region", "monthDow")
) ~ hum + iwp(
  pm,
  p = 2,
  ref_value = 10,
  knots = knots_pm
)

# Fit ---------------------------------------------------------------------
res <- hnlm(
  formula = cc_formula,
  data = dataOrig,
  config = list(
    transform_theta = TRUE,
    num_threads = 1L,
    num_shards = 50L,
    verbose = FALSE
  ),
  control = list(
    maxit = 200,
    trace = 1
  ),
  control_inner = list(
    maxit = 1000,
    report.level = 0
  )
)

cat("log-lik:", res$log_lik, " converged:", res$converged, "\n")

# Development bundle for likelihood probes --------------------------------
forres <- hnlm(
  formula = cc_formula,
  data = dataOrig,
  config = res$config,
  for_dev = TRUE
)

x_outer <- res$optim$par
gamma_mode <- res$cache$gamma

lik <- log_lik_laplace(
  x = x_outer,
  gamma = gamma_mode,
  ad_fun = forres$ad_fun,
  config = forres$config,
  control = res$control_inner,
  deriv = TRUE
)

cat("profile log-lik check:", lik$log_lik, "\n")

# 1-D likelihood slice over the first theta parameter ---------------------
sz <- forres$ad_fun@sizes
theta_pos <- as.integer(sz["beta"]) + 1L
bnd <- 0.1
SxL <- seq(-bnd, bnd, length.out = 6) + x_outer[theta_pos]

theL <- lapply(SxL, function(val) {
  x_try <- x_outer
  x_try[theta_pos] <- val
  log_lik_laplace(
    x = x_try,
    gamma = gamma_mode,
    ad_fun = forres$ad_fun,
    config = forres$config,
    control = res$control_inner,
    deriv = TRUE
  )
})

plot(
  SxL,
  vapply(theL, `[[`, numeric(1), "neg_log_lik"),
  type = "b",
  xlab = "theta perturbation",
  ylab = "neg log-lik"
)

if (length(theL) >= 3L && !is.null(theL[[3]]$deriv$d_log_lik)) {
  abline(
    h = -theL[[3]]$deriv$d_log_lik[1],
    col = "red",
    lty = 2
  )
}

# Conditional PM effect curves --------------------------------------------
xseq <- seq(knots_pm[1], rev(knots_pm)[1], by = 0.5)
newx <- list(pm = expand.grid(
  pm = xseq,
  hum = 0,
  region = unique(dataOrig$region)[1]
))

sim <- cond_sim_iwp(
  fit = lik,
  model_data = forres$model_data,
  newx = newx,
  n = 500,
  probs = c(0.005, 0.025, 0.1, 0.5, 0.9, 0.975, 0.995)
)

pm_q <- sim$quantiles$common$pm
Sprob <- c(0.005, 0.025, 0.1, 0.5, 0.9, 0.975, 0.995)
Slwd <- 4 * sqrt(pmin(Sprob, 1 - Sprob))
Scol <- c("grey", "black")[1 + (Sprob == 0.5)]
Sorder <- order(Slwd)

matplot(
  sim$x$pm,
  sim$sim$pm[, Sorder],
  type = "l",
  lty = 1,
  col = Scol[Sorder],
  lwd = Slwd[Sorder],
  xaxs = "i",
  yaxs = "i",
  xlab = "pm",
  ylab = "RR"
)
lines(
  xseq,
  exp(genPmEffect(xseq, 1, 1) - genPmEffect(ref_values$pm, 1, 1)),
  col = "yellow",
  lwd = 2
)

# Legacy third-derivative debugging (removed from hpolcc) -----------------
# thirdDeriv(), thirdStrata(), thirdOffDiagonals(), objectiveFunctionC(),
# get_ad_fun(), and wrappers_outer referred to the pre-adlaplace backend.
# Use log_lik_laplace(..., deriv = TRUE) and adlaplace::grad() / hessian()
# for derivative checks instead.

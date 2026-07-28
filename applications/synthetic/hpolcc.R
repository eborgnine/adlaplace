# Development check: evaluate Laplace likelihood at fixed parameters
# Requires simData.rds with components data, par (named beta/gamma/theta)

xx <- readRDS("simData.rds")
data <- xx$data

knots_pm <- seq(0, max(ceiling(max(data$pm) / 5) * 5, 20), 2.5)
data$monthDow <- format(data$date, "%Y-%m-%a")

library(adlaplace)
library(adlaplaceHgp)
library(hpolcc)

cc_formula <- dirichlet_multinom(
  count,
  by = c("region", "monthDow")
) ~ hum + iwp(
  pm,
  p = 2,
  ref_value = 10,
  knots = knots_pm
)

forres <- hnlm(
  formula = cc_formula,
  data = data,
  config = list(
    transform_theta = FALSE,
    num_threads = 4L,
    num_shards = 50L,
    verbose = FALSE
  ),
  for_dev = TRUE
)

beta_hat <- xx$par[grep("beta", names(xx$par))]
gamma_start <- xx$par[grep("gamma", names(xx$par))]
theta_hat <- xx$par[grep("theta", names(xx$par))]
x_outer <- c(beta_hat, theta_hat)

lik <- log_lik_laplace(
  x = x_outer,
  gamma = gamma_start,
  ad_pack = forres$ad_pack,
  config = forres$config,
  control = list(maxit = 300, report.level = 0),
  deriv = TRUE
)

cat("log-lik at fixed parameters:", lik$log_lik, "\n")
cat("inner gamma length:", length(lik$opt$solution), "\n")

if (requireNamespace("numDeriv", quietly = TRUE)) {
  bob <- function(x) {
    log_lik_laplace(
      x = c(x, theta_hat),
      gamma = gamma_start,
      ad_pack = forres$ad_pack,
      config = forres$config,
      control = list(maxit = 300, report.level = 0),
      deriv = FALSE
    )$neg_log_lik
  }
  H_num <- numDeriv::hessian(bob, beta_hat[seq_len(min(5L, length(beta_hat)))])
  print(H_num)
}

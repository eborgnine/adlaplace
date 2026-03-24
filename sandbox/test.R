 x = readRDS("~/research/healthcanada/pkg/cancc/x.rds")
 x$sqrt_pm = sqrt(pmin(x$pm2p5_lag_mean*1e9, 500))
x$sqrt_o3 = sqrt(x$gtco3_lag_mean*1e3)
x$temperature = x$t2m_lag_mean - 273.15
x$yearMonthDow = format(x$date, '%Y-%b-%a')

res <- hpolcc::hnlm(
  formula = y ~ sqrt_pm + sqrt_o3,
  data = x,
  cc_design = c("cdcode", "yearMonthDow"),
  config = list(
    transform_theta = TRUE,
    dirichlet = TRUE,
    boundary_is_random = TRUE,
    num_threads = parallel::detectCores(),
    verbose = 20
  ),
  control = list(
    maxit = 1000,
    report.level = 4,
    report.freq = 1
      )
)


log_pmf_sum <- function(y, tau, mu=c()) {

    if(!length(mu)) {
    mu = rep(1/length(y), length(y))
  }

  if (!is.numeric(y) || any(y < 0) || any(y != floor(y))) {
    stop("y must be a vector of nonnegative integers.")
  }
  if (!is.numeric(mu) || any(mu <= 0)) {
    stop("mu must be a positive numeric vector.")
  }
  if (length(y) != length(mu)) {
    stop("y and mu must have the same length.")
  }
  if (abs(sum(mu) - 1) > 1e-12) {
    stop("mu must sum to 1.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || tau < 0) {
    stop("tau must be a single nonnegative number.")
  }

  N <- sum(y)

  log_num <- 0
  for (k in seq_along(y)) {
    if (y[k] > 0) {
      log_num <- log_num + sum(log(mu[k] + (0:(y[k] - 1)) * tau))
    }
  }

  log_den <- if (N > 0) sum(log(1 + (0:(N - 1)) * tau)) else 0

  lgamma(N + 1) - sum(lgamma(y + 1)) + log_num - log_den
}


dirichlet_multinomial_density <- function(y, alpha, mu=c()) {
  N <- sum(y)

  if(!length(mu)) {
    mu = rep(1/length(y), length(y))
  }

  alpha_here <- alpha * mu

  lgamma(alpha) + lgamma(N + 1) - lgamma(N + alpha) +
    sum(lgamma(y + alpha_here)) -
    sum(lgamma(alpha_here)) -
    sum(lgamma(y + 1))
}

dirichlet_multinomial_cc = function(
  y, cc_matrix, overdisp, mu=c(), FUN=dirichlet_multinomial_density
  ) {
  result = rep(0, ncol(cc_matrix))
  for(D in 1:length(result)) {
    y_here = y[
      seq(cc_matrix@p[D], cc_matrix@p[D+1]-1L)+1L
    ]
    result[D] = FUN(y_here,overdisp, mu)
  }
  sum(result)
}

my_cc_matrix = res$objects$tmb_data$elgm_matrix[,1:5]
my_cc_matrix = my_cc_matrix[apply(my_cc_matrix,1,sum)>0,]



SlogSd = seq(-6, -1, len=41)
Ssd = exp(SlogSd)
Salpha = Ssd^(-2)
Slogalpha = log(Salpha)


Sll = mapply(
  dirichlet_multinomial_cc,
  overdisp = Salpha,
  MoreArgs = list(
    y = res$objects$tmb_data$y,
    cc_matrix = my_cc_matrix
    )
  )


Sll2 = mapply(
  dirichlet_multinomial_cc,
  overdisp = 1/Salpha,
  MoreArgs = list(
    y = res$objects$tmb_data$y,
    cc_matrix = my_cc_matrix,
    FUN = log_pmf_sum
    )
  )

matplot(Ssd, cbind(Sll, Sll2), log='x', type='l', lwd=c(3,2), lty=c(1,2))



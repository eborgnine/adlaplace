
# packages ----------------------------------------------------------------
library(TMB)
library(Matrix)
library(data.table)


# setup -------------------------------------------------------------------
# set.seed(123)

nn <- c(5,10,15)*300
nc <- length(nn)
np <- 2 # number of exposures
n <- sum(nn)
data <- data.table(
  id = 1:n,
  city = lapply(1:nc, \(k) rep(LETTERS[k], nn[k])) |> unlist(),
  dow = factor(weekdays(as.Date(1:n), T), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")),
  hum = runif(n),
  pol1 = rgamma(n, 20, .9),
  pol2 = rgamma(n, 50, 1.25)
)
for_dev = F

# order data by city and time
cities <- unique(data$city)
head_id <- lapply(cities, \(ct) which(data$city == ct)[1:5]) |> unlist()

source("R/hm.R")
source("R/hm_helpers.R")
source("R/case_crossover.R")
source("R/effects.R")
source("R/ziangs_function.R")
# knots_pol1 <- floor(min(data$pol1)):ceiling(max(data$pol1))
knots_pol1 <- seq(floor(min(data$pol1)/5)*5,ceiling(max(data$pol1)/5)*5,5)
# knots_pol2 <- floor(min(data$pol1)):ceiling(max(data$pol1))
formula <- hum ~ hum + od(id) +
  hiwp(pol1, p=2, ref_value = knots_pol1[round(length(knots_pol1)/2)], knots = knots_pol1, group_var = city)

cc_design <- ccDesign(time_var = "id", strat_vars = "city")
data <- data[sample(nrow(data)),]

list2env(hm(formula, data, cc_design = cc_design, for_dev = T), envir = environment())

# Log-precision parameters (2 parameters for 2 pollutants)
theta_id <- theta_info$id
theta_true <- log(rexp(length(unique(theta_id)), .0005))  # Log precision for each pollutant (will be exp(theta) to scale Q)
theta_true[1] <- 6
theta_true_long <- theta_true[theta_id]  # Log precision for each pollutant (will be exp(theta) to scale Q)

# Observed counts (Poisson responses)
beta_true <- .05 + rnorm(ncol(X), 0, .01)
Q <- .bdiag(Qs)
gamma_true <- lapply(seq_along(gamma_split), \(k){
  ii <- cumsum(c(0,gamma_split))[k + 0:1] + c(1,0)
  Q0 <- exp(theta_true_long[theta_id[k]+1]) * Q[ii[1]:ii[2],ii[1]:ii[2]]
  if(is.numeric(Q0) || Matrix::isDiagonal(Q0)){
    if(is.numeric(Q0)) Q0 <- as.matrix(Q0)
    rnorm(nrow(Q0), sd = 1/sqrt(diag(Q0)))
  }else{
    S <- solve(Q0)
    mvtnorm::rmvnorm(1, sigma = S |> as.matrix())
  }
}) |> unlist()
plot(gamma_true)

lambda_true <- exp(2 + X %*% beta_true + A %*% gamma_true) |> as.numeric()
y <- rpois(n, lambda_true)
hist(y, breaks=c(0:max(y)))

data$y <- integer(nrow(data))
data$y[new_order] <- y # NOTE HERE new_order mimicks how the data was permuted in hm.

formula <- paste0("y ~ ", as.character(formula)[3]) |> as.formula()
cc_design <- ccDesign(time_var = "id", strat_vars = "city", time_size = 4)
opt <- hm(formula, data, cc_design,
          tmb_parameters = list(beta = beta_true,
                                gamma = gamma_true,
                                theta = theta_true_long))

opt_par <- opt$env$last.par.best
beta_hat <- opt_par[names(opt_par) == "beta"]
gamma_hat <- opt_par[names(opt_par) == "gamma"]
theta_hat <- opt_par[names(opt_par) == "theta"]

plot(beta_true, col="green", ylim = range(c(beta_hat, beta_true))); points(beta_hat)
plot(gamma_true, col="green", ylim = range(c(gamma_hat, gamma_true))); points(gamma_hat)
plot(theta_true, col="green", ylim = range(c(theta_hat, theta_true))); points(theta_hat)

z_id <- which.max(gamma_split)
z1 <- cumsum(c(0,gamma_split))[z_id] + 1
z2 <- cumsum(c(0,gamma_split))[z_id+1]
plot(gamma_true[-(z1:z2)], col="green", ylim = range(c(gamma_hat, gamma_true))); points(gamma_hat[-(z1:z2)])
hist(gamma_true[z1:z2])
hist(gamma_hat[z1:z2])



parameters_list <- list(
  beta = beta_true,
  gamma = gamma_true,
  theta = theta_true
)

n_beta <- length(parameters_list$beta)
n_gamma <- length(parameters_list$gamma)
n_theta <- length(parameters_list$theta)
par <- unlist(parameters_list)

objective_function <- function(par) {
  
  # Extract the parameters (beta, gamma, theta) from the 'par' vector
  beta <- par[1:n_beta]
  gamma <- par[(n_beta + 1):(n_beta + n_gamma)]
  theta <- par[(n_beta + n_gamma + 1):length(par)]
  
  # Data inputs
  X <- data_list$X              # Design matrix for fixed effects
  A <- data_list$A              # Design matrix for random effects
  y <- data_list$y              # Observed counts (Poisson responses)
  Q <- data_list$Q              # Full precision matrix
  
  gamma_split <- data_list$gamma_split  # Length of gamma for each city
  theta_id <- data_list$theta_id        # Length of gamma for each city
  
  # Initialize negative log-likelihood
  eta_fixed <- X %*% beta
  eta_random <- A %*% gamma
  log_lambda <- (eta_fixed + eta_random) |> as.numeric()
  nll <- sum(dpois(y, exp(log_lambda), log = TRUE))
  
  # Random effects' contribution
  exp_theta <- rep(NA, sum(gamma_split))
  idx_start <- 1
  len_j <- 0
  
  # Assign the appropriate exp(theta) values based on theta_id
  for (j in 1:length(gamma_split)) {
    len_j <- gamma_split[j]
    exp_theta[idx_start:(idx_start + len_j - 1)] <- exp(theta[theta_id[j]])
    idx_start <- idx_start + len_j
  }
  
  # Apply the scaling of gamma
  gamma <- gamma * sqrt(exp_theta)
  
  # Loop over cities (or pollutants) to compute the random effects contribution
  d <- 0
  for (i in seq_along(theta_id)) if(theta_id[i] < np) d <- d + sum(gamma_split[i])
  nc <- length(gamma) / d  # number of cities (plus one for global effects)
  
  for (i in 1:nc) {
    gamma_i <- gamma[((i - 1) * d + 1):(i * d)]
    nll <- nll + 0.5 * sum(gamma_i * (Q %*% gamma_i))
  }
  
  # Return the negative log-likelihood
  return(nll)
}

int d = 0;
for (i in 1:theta_id.size()) if(theta_id(i) < np) d += sum(gamma_split(i))


# # Example of setting up the data and parameters for fitting:
# data_list <- list(
#   X = matrix(runif(100), ncol = 5),  # example fixed effect design matrix
#   A = matrix(runif(100), ncol = 5),  # example random effect design matrix
#   y = rpois(20, lambda = 5),         # example observed counts (Poisson responses)
#   Q = diag(5),                       # example precision matrix (identity)
#   gamma_split = c(2, 3),             # example split for gamma (cities and pollutants)
#   theta_id = c(1, 1)                 # example indices for log-precision
# )
# 
# # Initialize parameter vector with random values
# n_beta <- 5
# n_gamma <- sum(data$gamma_split)
# par <- c(runif(n_beta), runif(n_gamma), runif(length(data$theta_id)))
# 
# # Call the objective function for optimization
# result <- objective_function(par)
# 
# print(result)

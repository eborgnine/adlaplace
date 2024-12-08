library(Matrix)  # For sparse matrix operations

gamma_split = gamma_info$split
gamma_nreplicate = gamma_info$nreplicate # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
Q = Qs |> .bdiag()
theta_id = theta_info$id
theta_hyper = theta_info$hyper

beta = rep(0, ncol(X))
gamma = rep(0, ncol(A))
theta = theta_info$init



# Function to compute negative log-likelihood
# compute_nll <- function(X, A, y, Q, gamma_split, beta, gamma, theta) {
  # Initialize negative log-likelihood
  nll <- 0
  
  # Fixed effects contribution
  eta_fixed <- as.vector(X %*% beta)
  
  # Random effects contribution
  eta_random <- as.vector(A %*% gamma)
  
  # Compute log-lambda and Poisson log-likelihood
  log_lambda <- eta_fixed + eta_random
  nll <- -sum(dpois(y, exp(log_lambda), log = TRUE))
  
  # Contribution from random effects
  exp_theta <- rep(0, ncol(A))  # Placeholder for exp(theta)
  idx_start <- 1
  for (j in seq_along(gamma_split)) {
    len_j <- gamma_split[j]
    nll <- nll - 0.5 * len_j * theta[j]  # Log-determinant contribution
    exp_theta[idx_start:(idx_start + len_j - 1)] <- exp(theta[j])  # Compute exp(theta)
    idx_start <- idx_start + len_j
  }
  
  # Modify gamma for scaled precision matrix
  sqrt_exp_theta <- sqrt(exp_theta)
  gamma <- gamma * sqrt_exp_theta
  
  # Add precision matrix contribution
  nll <- nll + 0.5 * sum(gamma * (Q %*% gamma))
  
#   return(nll)
# }

# Example usage
set.seed(123)

# Example data
n <- 100
p <- 3
cities <- 5

X <- Matrix(rnorm(n * p), n, p, sparse = TRUE)
A <- Matrix(rnorm(n * cities), n, cities, sparse = TRUE)
y <- rpois(n, lambda = 10)
Q <- Diagonal(cities, x = runif(cities, 0.1, 1))  # Example sparse precision matrix
gamma_split <- rep(1, cities)  # Each city has one random effect

# Example parameters
beta <- rnorm(p)
gamma <- rnorm(cities)
theta <- log(runif(cities, 0.1, 1))

# Compute negative log-likelihood
nll <- compute_nll(X, A, y, Q, gamma_split, beta, gamma, theta)
print(nll)

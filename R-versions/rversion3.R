list2env(tmb_data, envir = environment())
list2env(tmb_parameters, envir = environment())

n_cc <- nrow(cc_matrix)
d_cc <- ncol(cc_matrix)


# Compute eta = log(lambda)
nll <- 0
eta <- rep(0, length(y))
if (ncol(X) > 0) eta <- eta + X %*% beta
if (ncol(A) > 0) eta <- eta + A %*% gamma

# Compute negative log likelihood (multinomial -- case crossover)
lsa <- -Inf
for (i in 1:n_cc) {
  lsa <- -Inf
  for (j in 1:d_cc) {
    if (cc_matrix[i, j] != 0) {
      lsa <- log(exp(lsa) + exp(eta[cc_matrix[i, j] - 1]))
    }
  }
  
  # - y_ij * log(p_ij), where p_ij = exp(eta_ij)/sum(exp(eta_i)) = exp(eta_ij - lsa)
  for (j in 1:d_cc) {
    if (cc_matrix[i, j] != 0) {
      nll <- nll - y[cc_matrix[i, j] - 1] * (eta[cc_matrix[i, j] - 1] - lsa)
    }
  }
}

# Random effects contribution
exp_theta <- numeric(ncol(A))
idx_start <- 1
len_j <- 0
for (j in 1:length(gamma_split)) {
  len_j <- gamma_split[j]
  # nll <- nll - 0.5 * (len_j * theta[j] + log_dets[j]) # log-determinant part
  nll <- nll - 0.5 * len_j * theta[j]
  for (i in idx_start:(idx_start + len_j - 1)) {
    exp_theta[i] <- exp(theta[j])
  }
  idx_start <- idx_start + len_j
}

# do sqrt(exp_theta) * gamma instead of exp_theta * Q
gamma <- gamma * sqrt(exp_theta)
nll <- nll + 0.5 * sum(gamma * (Q %*% gamma))

nll

list2env(data_list, envir = environment())
list2env(parameters_list, envir = environment())

eta_fixed = X %*% beta
eta_random = A %*% gamma
log_lambda = as.numeric(eta_fixed + eta_random)
nll = -sum(dpois(y, exp(log_lambda), TRUE))


exp_theta = rep(NA, ncol(A))
idx_start = 1
len_j = 0
for(j in 1:length(gamma_split)){
  len_j = gamma_split[j]
  for(i in idx_start:(idx_start+gamma_split[j]-1)){
    exp_theta[i] = exp(theta[theta_id[j]+1]);                     
  }
  nll = nll - 0.5 * exp(theta[theta_id[j]]) * log_dets[theta_id[j]]; 
  idx_start = idx_start + len_j;
}

gamma = gamma * sqrt(exp_theta)
d = sum(gamma_split[1:np])
# for(i in 1:np) d <- d+ gamma_split[i];
for (i in 0:nc) {
  gamma_i = gamma[i*d + 1:d];
  nll = nll + 0.5 * sum(gamma_i * (Q %*% gamma_i));
}

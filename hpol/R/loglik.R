library(hpolcc)
library(trustOptim)
library(Matrix)

get_gamma_hat_and_derivative <- function(beta, theta, gamma_start, fitD) {

  # Indices setup
  Sbeta <- seq_along(beta)
  Stheta <- seq_along(theta) + length(beta)
  Sgamma <- seq_along(gamma_start) + length(beta) + length(theta)
  
  # Complete parameter vector
  param_all <- c(beta, theta, gamma_start)

  # Optimize gamma keeping beta and theta fixed
  wrappers_gamma <- hpolcc::make_trustoptim_wrappers(
    data = fitD$tmb_data,
    config = c(fitD$config, list(
      beta = beta,
      theta = theta,
      transform_theta = TRUE,
      dirichelet = TRUE,
      num_threads = fitD$config$num_threads,
      sparsity = fitD$config$sparsity[
        fitD$config$sparsity$i %in% Sgamma & 
          fitD$config$sparsity$j %in% Sgamma, 
      ] - (length(beta) + length(theta))
    )),
    debug = FALSE
  )

  result_gamma <- trust.optim(
    x = gamma_start,
    fn = wrappers_gamma$fn,
    gr = wrappers_gamma$gr,
    hs = wrappers_gamma$hs,
    method = "Sparse",
    control = list(
      start.trust.radius = 50, maxit = 500, trust.iter = 5000,
      preconditioner = 1, report.level = 0
    )
  )

  gamma_hat <- result_gamma$solution

  # Hessian for gamma
  res_gamma <- hpolcc:::objectiveFunctionC(
    gamma_hat, fitD$tmb_data, wrappers_gamma$config
  )

  H_gamma <- sparseMatrix(
    i = res_gamma$hessian$i,
    j = res_gamma$hessian$j,
    x = res_gamma$hessian$x,
    symmetric = TRUE
  )

  # Full Hessian at (beta, theta, gamma_hat)
  res_full <- hpolcc:::objectiveFunctionC(
    c(beta, theta, gamma_hat), fitD$tmb_data, fitD$config
  )

  # Cross derivative: gamma vs (beta, theta)
  cross_hessian_entries <- res_full$hessian[
    res_full$hessian$i %in% Sgamma & res_full$hessian$j %in% c(Sbeta, Stheta),
  ]

  H_gamma_beta_theta <- sparseMatrix(
    i = match(cross_hessian_entries$i, Sgamma),
    j = match(cross_hessian_entries$j, c(Sbeta, Stheta)),
    x = cross_hessian_entries$x,
    dims = c(length(Sgamma), length(Sbeta) + length(Stheta))
  )

  # Compute derivative via implicit function theorem
  gamma_hat_derivative <- -solve(H_gamma, H_gamma_beta_theta)

  return(list(
    gamma_hat = gamma_hat,
    derivative = gamma_hat_derivative
  ))
}

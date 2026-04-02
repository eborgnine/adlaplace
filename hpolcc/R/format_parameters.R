#' parameters_info is a list with elements beta, gamma, theta
#' each needs columns var, name, theta needs log

format_parameters <- function(x) {
  full_parameters <- x$extra$full_parameters
  parameters_info <- x$objects$parameters_info

  Ntheta <- nrow(parameters_info$theta)
  Nbeta <- max(c(0, nrow(parameters_info$beta)))
  Ngamma <- nrow(parameters_info$gamma)

  parameters_info$beta$mle <- full_parameters[seq(from = 1, length.out = Nbeta)]
  parameters_info$theta$mle <- full_parameters[
    seq(to = length(full_parameters), length.out = Ntheta)
  ]
  if (x$objects$config$transform_theta) {
    parameters_info$theta$log_mle <- parameters_info$theta$mle
    not_logged = parameters_info$theta$model == "dispersion"
        parameters_info$theta$log_mle[not_logged] = NA
    parameters_info$theta$mle[!not_logged] <- exp(parameters_info$theta$log_mle[!not_logged])
  }
  parameters_info$gamma$mode <- full_parameters[seq(Nbeta + 1, length.out = Ngamma)]

  list(
    gamma = parameters_info$gamma,
    beta = parameters_info$beta,
    theta = parameters_info$theta
  )
}

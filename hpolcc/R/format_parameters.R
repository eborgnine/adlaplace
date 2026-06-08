#' Format Laplace parameter blocks into coefficient tables
#'
#' @param laplace Output of \code{\link[adlaplace]{log_lik_laplace}} with
#'   \code{full_parameters}.
#' @param parameters_info List with \code{beta}, \code{gamma}, and \code{theta}
#'   data frames from \code{model_data()$data$info}.
#' @param transform_theta Logical; back-transform logged theta columns when
#'   \code{TRUE}.
#' @keywords internal
format_parameters <- function(
  laplace,
  parameters_info,
  transform_theta = TRUE
) {
  if (is.list(laplace) && !is.null(laplace$objects)) {
    return(format_parameters(
      laplace = laplace$extra,
      parameters_info = laplace$objects$parameters_info,
      transform_theta = isTRUE(laplace$objects$config$transform_theta)
    ))
  }

  full_parameters <- laplace$full_parameters
  Ntheta <- nrow(parameters_info$theta)
  Nbeta <- max(c(0L, nrow(parameters_info$beta)))
  Ngamma <- nrow(parameters_info$gamma)

  parameters_info$beta$mle <- full_parameters[seq(from = 1L, length.out = Nbeta)]
  parameters_info$theta$mle <- full_parameters[
    seq(to = length(full_parameters), length.out = Ntheta)
  ]
  if (transform_theta) {
    parameters_info$theta$log_mle <- parameters_info$theta$mle
    not_logged <- parameters_info$theta$model == "overdispersion"
    parameters_info$theta$log_mle[not_logged] <- NA
    parameters_info$theta$mle[!not_logged] <- exp(
      parameters_info$theta$log_mle[!not_logged]
    )
  }
  parameters_info$gamma$mode <- full_parameters[seq(Nbeta + 1L, length.out = Ngamma)]

  list(
    gamma = parameters_info$gamma,
    beta = parameters_info$beta,
    theta = parameters_info$theta
  )
}

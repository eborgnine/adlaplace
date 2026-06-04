#' Build backend AD function handle
#'
#' Constructs skew-normal observation, random, and extra shards and returns an
#' \code{ad_fun} object with Hessian templates attached.
#'
#' @param data Model data list with \code{y}, \code{ATp}, \code{XTp}, and
#'   \code{Qdiag} (diagonal precision weights for random effects).
#' @param config Model configuration list (\code{beta}, \code{gamma}, \code{theta},
#'   optional \code{shards}, \code{transform_theta}, etc.).
#'
#' @return An \code{ad_fun} S4 object from \pkg{adlaplace}.
#'
#' @importFrom adlaplace ad_data ad_fun ad_fun_ptr
#' @seealso \code{\link[adlaplace]{ad_fun}}
#' @export
getAdFun_r <- function(data, config) {
  n_beta <- nrow(data$XTp)
  n_gamma <- nrow(data$ATp)
  n_theta <- length(config$theta)
  obs_theta_idx <- (n_theta - 1L):n_theta
  rand_theta_idx <- seq_len(n_theta - 2L)

  obs_shard <- adlaplace::ad_data(
    y = data$y,
    A = Matrix::t(data$ATp),
    X = Matrix::t(data$XTp),
    beta_map = Matrix::Diagonal(n_beta),
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(obs_theta_idx, n_theta),
    ad_kind = "observations",
    ad_fun = "skewnormal_obs"
  )

  random_shard <- adlaplace::ad_data(
    beta_map = n_beta,
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(rand_theta_idx, n_theta),
    ad_kind = "random",
    ad_fun = "random_diagonal",
    precision = data$Qdiag
  )

  extra_shard <- adlaplace::ad_data(
    y = data$y,
    beta_map = n_beta,
    gamma_map = n_gamma,
    theta_map = list(obs_theta_idx, n_theta),
    ad_kind = "parameters",
    ad_fun = "skewnormal_extra"
  )

  ptrs <- list(
    adlaplace::ad_fun_ptr(obs_shard, config),
    adlaplace::ad_fun_ptr(random_shard, config),
    adlaplace::ad_fun_ptr(extra_shard, config)
  )

  ad_ptr <- do.call(c, ptrs)
  num_threads <- if (!is.null(config$num_threads)) {
    as.integer(config$num_threads)[1L]
  } else {
    1L
  }
  adlaplace::ad_fun(ad_ptr, num_threads = num_threads)
}

#' Assemble model terms and per-shard \code{ad_data} objects
#'
#' Parses a formula, builds design matrices, and returns formula terms plus
#' \code{ad_data} shards for the observation density, parameter densities, and
#' each random-effect term.
#'
#' @param formula Model formula with \code{f()} or constructor terms.
#' @param data Data frame containing variables referenced in \code{formula}.
#' @param verbose Print extra information while parsing terms.
#' @return A list with:
#' \describe{
#'   \item{\code{terms}}{Named list of term objects from \code{\link{collect_terms}}.}
#'   \item{\code{shards}}{Named list of \code{ad_data} objects, including
#'     \code{parent} (full layout) and one entry per AD shard
#'     (observations, parameters, random effects).}
#'   \item{\code{info}}{List with \code{beta}, \code{gamma}, \code{theta},
#'     \code{parameters} data frames for optimization metadata.}
#'   \item{\code{precisions}}{Named list of precision lists for random shards
#'     (parallel names to random entries in \code{shards}).}
#' }
#' @export
#' @examples
#' \dontrun{
#' md <- model_data(
#'   y ~ x1 + f(x2, model = "iwp", p = 2, knots = seq(0, 1, len = 11)),
#'   data = dat
#' )
#' }
model_data <- function(formula, data, verbose = FALSE) {
  all_data <- data_setup(formula, data, verbose = verbose)
  the_terms <- all_data$terms

  n_beta <- nrow(all_data$info$beta)
  n_gamma <- nrow(all_data$info$gamma)
  n_theta <- nrow(all_data$info$theta)

  observations <- list()
  random <- list()

  for (term_here in the_terms) {
    if (!methods::is(term_here, "model")) {
      warning("term doesnt appear to be a model")
      next
    }
    ad_name <- term_here@ad_fun
    kind <- term_here@ad_kind

    if (identical(kind, "observations")) {
      observations[[term_here@term]] <- ad_data(
        y = all_data$y,
        A = all_data$A,
        X = all_data$X,
        beta_map = Matrix::Diagonal(n_beta),
        gamma_map = Matrix::Diagonal(n_gamma),
        theta_map = list(
          grep(term_here@label, all_data$info$theta$label),
          as.integer(n_theta)
        ),
        ad_fun = ad_name,
        ad_kind = "observations"
      )
    }

    if (identical(kind, "random")) {
      prec_mat <- precision(term_here, all_data$data)
      # will be NULL if sd = Inf, diffuse prior.
      if (!is.null(prec_mat)) {
        random[[term_here@term]] <- ad_data(
          beta_map = n_beta,
          gamma_map = list(
            which(all_data$info$gamma$label == term_here@label),
            n_gamma
          ),
          theta_map = list(
            which(all_data$info$theta$label == term_here@label),
            as.integer(n_theta)
          ),
          ad_fun = ad_name,
          ad_kind = "random",
          precision = as.numeric(Matrix::diag(prec_mat))
        )
      }
    }
  }

  parameters <- lapply(
    observations, function(xx) {
      ad_data(
        y = xx@y,
        theta_map = xx@theta_map,
        beta_map = nrow(xx@beta_map),
        gamma_map = nrow(xx@gamma_map),
        ad_fun = gsub("_obs$", "_extra", xx@ad_fun),
        ad_kind = "parameters"
      )
    }
  )
  names(parameters) <- paste0(names(parameters), "_extra")


  list(
    data = all_data,
    observations = observations,
    random = random,
    parameters = parameters,
    terms = the_terms
  )
}

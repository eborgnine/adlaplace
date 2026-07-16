#' Assemble model terms and per-shard \code{ad_data} objects
#'
#' Parses a formula, builds design matrices, and returns formula terms plus
#' \code{ad_data} shards for the observation density, parameter densities, and
#' each random-effect term.
#'
#' @param formula Model formula with constructor terms (e.g. \code{iwp()}, \code{iid()}).
#' @param data Data frame containing variables referenced in \code{formula}.
#' @param verbose Print extra information while parsing terms.
#' @param na_omit When \code{TRUE} (default), drop rows with \code{NA} in
#'   covariates, outcome, stratification (\code{by}), or random-slope
#'   \code{mult} variables; impute remaining outcome \code{NA}s as zero; and
#'   sort by observation-term \code{by} columns when present.
#' @return A list with:
#' \describe{
#'   \item{\code{terms}}{Named list of term objects from \code{\link{collect_terms}}.}
#'   \item{\code{data}}{Output of \code{\link{data_setup}}, including \code{elgm_matrix}
#'     when an observation term defines an \code{elgm_matrix} method.}
#'   \item{\code{observations}}{Named list of observation \code{ad_data} shards.}
#'   \item{\code{parameters}}{Named list of parameter \code{ad_data} shards.}
#'   \item{\code{random}}{Named list of random-effect \code{ad_data} shards.}
#' }
#' @export
#' @examples
#' \dontrun{
#' md <- model_data(
#'   y ~ x1 + iwp(x2, p = 2, knots = seq(0, 1, len = 11)),
#'   data = dat
#' )
#' }
model_data <- function(formula, data, verbose = FALSE, na_omit = TRUE) {
  formula_in <- formula
  the_terms <- parse_model_terms(formula, verbose = verbose)

  if (na_omit) {
    data <- prepare_model_rows(data, the_terms, verbose = verbose)
  }

  elgm_mats <- list()
  obs_terms <- Filter(function(t) {
    methods::is(t, "model") && identical(t@ad_kind, "observations")
  }, the_terms)

  for (term in obs_terms) {
    term_class <- class(term)[1L]
    if (methods::hasMethod("elgm_matrix", term_class)) {
      elgm_mats[[term@term]] <- elgm_matrix(term, data)
    }
  }

  if (length(elgm_mats) > 1L) {
    warning(
      "multiple observation terms with elgm_matrix methods; ",
      "using the first for data$elgm_matrix",
      call. = FALSE
    )
  }

  all_data <- data_setup(formula = formula_in, data = data, verbose = verbose)
  the_terms <- all_data$terms

  if (length(elgm_mats) == 1L) {
    all_data$elgm_matrix <- elgm_mats[[1L]]
  }

  n_beta <- max(c(0, nrow(all_data$info$beta))) # nrow might be NULL
  n_gamma <- max(c(0, nrow(all_data$info$gamma)))
  n_theta <- max(c(0, nrow(all_data$info$theta)))

  observations <- list()
  random <- list()
  random_parameters <- list()

  for (term_here in the_terms) {
    if (!methods::is(term_here, "model")) {
      warning("term doesnt appear to be a model")
      next
    }
    ad_name <- term_here@ad_fun
    kind <- term_here@ad_kind
    elgm_here <- elgm_mats[[term_here@term]]

    if (identical(kind, "observations")) {
      beta_map <- if (n_beta > 0L && ncol(all_data$X) > 0L) {
        list(
          match(colnames(all_data$X), all_data$info$beta$beta_label),
          as.integer(n_beta)
        )
      } else {
        n_beta
      }
      gamma_map <- if (n_gamma > 0L && ncol(all_data$A) > 0L) {
        list(
          match(colnames(all_data$A), all_data$info$gamma$gamma_label),
          as.integer(n_gamma)
        )
      } else {
        n_gamma
      }
      observations[[term_here@term]] <- ad_data(
        y = all_data$y,
        A = all_data$A,
        X = all_data$X,
        beta_map = beta_map,
        gamma_map = gamma_map,
        theta_map = list(
          grep(term_here@label, all_data$info$theta$label),
          as.integer(n_theta)
        ),
        elgm_matrix = elgm_here,
        ad_fun = ad_name,
        ad_kind = "observations",
        package = term_here@package,
        weights = {
          if (methods::is(term_here, "binomial") &&
              length(term_here@size) == 1L &&
              nzchar(term_here@size)) {
            size_col <- term_here@size
            if (!size_col %in% names(all_data$data)) {
              stop(
                "binomial size column '", size_col, "' not found in data",
                call. = FALSE
              )
            }
            as.numeric(all_data$data[[size_col]])
          } else {
            numeric(0)
          }
        }
      )
    }

    if (identical(kind, "random")) {
      prec_mat <- precision(term_here, all_data$data)
      if (!is.null(prec_mat)) {
        # Key by @label, not @term: multiple random terms can share a grouping
        # variable while still having distinct labels / gamma columns.
        random_name <- term_here@label
        if (random_name %in% names(random)) {
          stop(
            "duplicate random-term label '", random_name, "'; ",
            "each random term needs a unique @label",
            call. = FALSE
          )
        }
        prec_payload <- if (is.list(prec_mat)) {
          prec_mat
        } else {
          as.numeric(Matrix::diag(prec_mat))
        }
        theta_map_here <- list(
          grep(term_here@label, all_data$info$theta$label),
          as.integer(n_theta)
        )
        random[[random_name]] <- ad_data(
          beta_map = n_beta,
          gamma_map = list(
            which(all_data$info$gamma$label == term_here@label),
            n_gamma
          ),
          theta_map = theta_map_here,
          ad_fun = ad_name,
          ad_kind = "random",
          package = term_here@package,
          precision = prec_payload
        )
        extra_name <- extra_ad_fun(term_here)
        if (is.character(extra_name) && length(extra_name) == 1L &&
            nzchar(extra_name)) {
          random_parameters[[paste0(random_name, "_det")]] <- ad_data(
            beta_map = n_beta,
            gamma_map = n_gamma,
            theta_map = theta_map_here,
            ad_fun = extra_name,
            ad_kind = "parameters",
            package = term_here@package,
            precision = prec_payload
          )
        }
      }
    }
  }

  parameters <- Filter(
    Negate(is.null),
    lapply(
      observations,
      function(xx) {
        if (ncol(xx@theta_map) == 0L) {
          return(NULL)
        }
        ad_data(
          y = xx@y,
          theta_map = xx@theta_map,
          beta_map = nrow(xx@beta_map),
          gamma_map = nrow(xx@gamma_map),
          elgm_matrix = xx@elgm_matrix,
          ad_fun = parameter_ad_fun(xx@ad_fun),
          ad_kind = "parameters",
          package = xx@package,
          weights = xx@weights
        )
      }
    )
  )
  if (length(parameters) > 0L) {
    nm <- names(parameters)
    if (is.null(nm)) {
      nm <- rep("", length(parameters))
    }
    names(parameters) <- paste0(nm, "_extra")
  }
  if (length(random_parameters) > 0L) {
    parameters <- c(parameters, random_parameters)
  }

  list(
    data = all_data,
    observations = observations,
    random = random,
    parameters = parameters,
    terms = the_terms
  )
}

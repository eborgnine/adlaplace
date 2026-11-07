#' Bivariate iid random intercept
#'
#' One pair of coefficients per level of \code{x}, ordered station-major and
#' \code{by}-minor. The pair is normal with standard deviations \code{sd1},
#' \code{sd2} and correlation \code{rho}, independent across levels of \code{x}.
#' \code{rho} is uniform on its bounds (no extra prior).
#'
#' @name biv_iid-class
#' @aliases biv_iid
#' @docType class
#' @exportClass biv_iid
NULL

setClass(
  "biv_iid",
  slots = list(
    by = "by_group",
    rho_init = "numeric",
    rho_lower = "numeric",
    rho_upper = "numeric"
  ),
  contains = "model_term",
  prototype = prototype(
    ref_value = numeric(0),
    p.order = as.integer(0),
    knots = numeric(0),
    model_role = factor("random", levels = .model_role_levels),
    density = "biv_iid",
    ad_kind = "random",
    by = by_group(),
    rho_init = 0,
    rho_lower = -0.999,
    rho_upper = 0.999
  )
)

#' @param x Grouping variable (one pair of intercepts per level).
#' @param by Two-level variable selecting the pair. Pass a column name.
#' @param init,lower,upper,parscale Standard-deviation hyperparameter controls.
#'   Both standard deviations share these values and are optimized on the log scale.
#' @param rho_init,rho_lower,rho_upper Correlation parameter, on the natural scale.
#' @return A named list containing one \code{biv_iid} term.
#' @rdname biv_iid-class
#' @export
biv_iid <- function(x,
                    by,
                    init = .my_theta_init,
                    lower = .my_theta_lower,
                    upper = .my_theta_upper,
                    parscale = .my_theta_parscale,
                    rho_init = 0,
                    rho_lower = -0.999,
                    rho_upper = 0.999) {
  x <- strip_term_name(as.character(x))
  by <- strip_term_name(as.character(by))
  if (length(rho_init) != 1L || !is.finite(rho_init) ||
      rho_init <= rho_lower || rho_init >= rho_upper) {
    stop("rho_init must lie strictly inside (rho_lower, rho_upper)", call. = FALSE)
  }
  term <- methods::new(
    "biv_iid",
    name = x,
    label = paste(x, "biv", sep = "_"),
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    log = TRUE,
    by = by_group(term = by),
    rho_init = rho_init,
    rho_lower = rho_lower,
    rho_upper = rho_upper
  )
  result <- list(term)
  names(result) <- term@label
  result
}

.biv_pair_labels <- function(term, data) {
  if (length(term@by@levels) != 2L) {
    stop("biv_iid requires exactly two by levels", call. = FALSE)
  }
  station <- data[[term@name]]
  if (is.factor(station)) {
    s_levels <- levels(station)
  } else {
    s_levels <- as.character(sort(unique(station)))
  }
  if (!length(s_levels)) {
    stop("biv_iid grouping variable has no levels", call. = FALSE)
  }
  p_levels <- term@by@levels
  labels <- as.vector(vapply(
    s_levels,
    function(s) paste(term@label, s, p_levels, sep = "_"),
    character(length(p_levels))
  ))
  list(s_levels = s_levels, p_levels = p_levels, labels = labels)
}

#' @describeIn biv_iid-class One column per station-pollutant pair.
#' @param term A \code{biv_iid} term.
#' @param data Data frame containing the grouping columns.
#' @export
setMethod("design", "biv_iid", function(term, data) {
  info <- .biv_pair_labels(term, data)
  station <- data[[term@name]]
  if (!is.factor(station)) {
    station <- factor(as.character(station), levels = info$s_levels)
  }
  pol <- match(as.character(data[[term@by@term]]), info$p_levels)
  if (anyNA(pol)) {
    stop(
      "biv_iid by column has levels outside ",
      paste(info$p_levels, collapse = ", "),
      call. = FALSE
    )
  }
  s_idx <- as.integer(station)
  if (anyNA(s_idx)) {
    stop("biv_iid grouping variable has values outside its levels", call. = FALSE)
  }
  n_p <- length(info$p_levels)
  j <- (s_idx - 1L) * n_p + pol
  B <- Matrix::sparseMatrix(
    i = seq_len(nrow(data)),
    j = j,
    x = 1,
    dims = c(nrow(data), length(info$labels))
  )
  colnames(B) <- info$labels
  B
})

#' @describeIn biv_iid-class List payload so the 2x2 precision is built from theta.
#' @export
setMethod("precision", "biv_iid", function(term, data) {
  list(n_pollutant = 2L)
})

#' @describeIn biv_iid-class Two standard deviations and one correlation.
#' @export
setMethod("theta_info", "biv_iid", function(term) {
  if (length(term@by@levels) != 2L) {
    stop(
      "biv_iid by levels are empty; add_by_levels() must run before theta_info()",
      call. = FALSE
    )
  }
  sd_labels <- paste(term@label, "sd", term@by@levels, sep = "_")
  data.frame(
    term = term@name,
    model = "biv_iid",
    label = c(sd_labels, paste(term@label, "rho", sep = "_")),
    init = c(rep(term@init, 2L), term@rho_init),
    lower = c(rep(term@lower, 2L), term@rho_lower),
    upper = c(rep(term@upper, 2L), term@rho_upper),
    parscale = c(rep(term@parscale, 2L), 0.2),
    model_role = rep(term@model_role, 3L),
    log = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
})

#' @describeIn biv_iid-class One latent coefficient per station-pollutant pair.
#' @export
setMethod("random_info", "biv_iid", function(term, data) {
  info <- .biv_pair_labels(term, data)
  # expand.grid varies the first argument fastest, matching the design columns.
  grid <- expand.grid(
    pol = info$p_levels,
    station = info$s_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data.frame(
    term = term@name,
    model = "biv_iid",
    label = term@label,
    by = term@by@term,
    by_labels = grid$pol,
    basis = grid$station,
    order = NA,
    gamma_label = info$labels,
    stringsAsFactors = FALSE
  )
})

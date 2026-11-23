#' Cosine and sine columns for seasonal harmonics
#'
#' Column \eqn{2j-1} is \eqn{\cos(2\pi h_j u / \mathrm{period})} and column
#' \eqn{2j} is \eqn{\sin(2\pi h_j u / \mathrm{period})}, with
#' \eqn{u = x - \mathrm{ref\_value}} and \eqn{h_j} the \eqn{j}th entry of
#' \code{harmonics}. No column names and no exposure multiplier: those belong
#' to [sinusoids()].
#'
#' @param x Numeric covariate.
#' @param period Period of harmonic 1. Must be positive and finite.
#' @param harmonics Positive integers, in column order. Harmonic \code{h}
#'   has frequency \code{h / period}.
#' @param ref_value Value subtracted from \code{x} before the angle is formed.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{2 * length(harmonics)} columns.
#' @export
sinusoid_basis <- function(x, period, harmonics = 1L, ref_value = 0) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("x must be finite", call. = FALSE)
  }
  period <- .sinusoids_period(period)
  harmonics <- .sinusoids_harmonics(harmonics)
  ref_value <- .sinusoids_ref(ref_value, "ref_value")
  u <- x - ref_value
  n_h <- length(harmonics)
  D <- matrix(0, length(u), 2L * n_h)
  for (j in seq_len(n_h)) {
    angle <- 2 * pi * harmonics[[j]] * u / period
    D[, 2L * j - 1L] <- cos(angle)
    D[, 2L * j] <- sin(angle)
  }
  D
}

.sinusoids_period <- function(period) {
  period <- as.numeric(period)
  if (length(period) != 1L || !is.finite(period) || period <= 0) {
    stop("period must be a single positive finite value", call. = FALSE)
  }
  period
}

.sinusoids_harmonics <- function(harmonics) {
  harmonics <- as.numeric(harmonics)
  if (!length(harmonics) ||
    any(!is.finite(harmonics)) ||
    any(harmonics != as.integer(harmonics)) ||
    any(harmonics < 1)) {
    stop("harmonics must be unique positive integers", call. = FALSE)
  }
  harmonics <- as.integer(harmonics)
  if (anyDuplicated(harmonics)) {
    stop("harmonics must be unique positive integers", call. = FALSE)
  }
  harmonics
}

.sinusoids_ref <- function(value, what) {
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value)) {
    stop(what, " must be a single finite value", call. = FALSE)
  }
  value
}

#' Absent multiplier is \code{NA_character_}; a column name is length 1.
#' @noRd
.sinusoids_mult <- function(mult) {
  if (is.null(mult) || !length(mult) || anyNA(mult) || !nzchar(mult)) {
    return(NA_character_)
  }
  mult <- as.character(mult)
  if (length(mult) != 1L) {
    stop("mult must be a single column name", call. = FALSE)
  }
  mult
}

.sinusoids_has_mult <- function(mult) {
  length(mult) == 1L && !is.na(mult) && nzchar(mult)
}

.sinusoids_label <- function(name, mult) {
  if (.sinusoids_has_mult(mult)) {
    paste(c(name, mult, "sinusoids"), collapse = "_")
  } else {
    paste(name, "sinusoids", sep = "_")
  }
}

.sinusoids_colnames <- function(label, harmonics) {
  as.vector(rbind(
    paste0(label, "_c", harmonics),
    paste0(label, "_s", harmonics)
  ))
}

#' Sinusoid fixed-effect term
#'
#' Fixed effects for \eqn{\cos(2\pi h u/\mathrm{period})} and
#' \eqn{\sin(2\pi h u/\mathrm{period})}, \eqn{u = x - \mathrm{ref\_value}},
#' one pair per entry of \code{harmonics}. An optional \code{mult} column
#' scales every column by \code{mult - ref_mult}, so the linear predictor
#' contains terms \code{exposure * cos(2 * pi * h * x / period)}.
#'
#' @slot period Period of harmonic 1.
#' @slot harmonics Positive integers selecting the harmonics.
#' @slot mult Exposure column, or \code{NA} when the basis is not multiplied.
#' @slot ref_mult Value subtracted from \code{mult}.
#' @exportClass sinusoids
#' @rdname sinusoids-class
setClass(
  "sinusoids",
  slots = list(
    period = "numeric",
    harmonics = "integer",
    mult = "character",
    ref_mult = "numeric"
  ),
  contains = "model_term",
  prototype = prototype(
    period = 365.25,
    harmonics = 1L,
    mult = NA_character_,
    ref_mult = 0,
    model_role = factor("fixed", levels = .model_role_levels),
    density = NA_character_,
    ad_kind = NA_character_
  )
)

#' @param x Variable name (character), or a symbol in a formula.
#' @param period Period of harmonic 1. Default 365.25.
#' @param harmonics Positive integers. Default \code{1L} (one cosine and one sine).
#' @param mult Exposure column. \code{NA} (default) leaves the basis unscaled.
#' @param ref_value Value subtracted from \code{x}.
#' @param ref_mult Value subtracted from \code{mult}.
#' @param init,lower,upper,parscale Fixed-effect settings, recycled to two
#'   parameters per harmonic (cosine, then sine).
#' @return A \code{sinusoids} term.
#' @rdname sinusoids-class
#' @export
sinusoids <- function(x,
                      period = 365.25,
                      harmonics = 1L,
                      mult = NA_character_,
                      ref_value = 0,
                      ref_mult = 0,
                      init = .my_beta_init,
                      lower = .my_beta_lower,
                      upper = .my_beta_upper,
                      parscale = .my_beta_parscale) {
  if (is.symbol(x) || is.name(x)) {
    x <- as.character(x)
  } else {
    x <- as.character(x)
  }
  if (length(x) != 1L || !nzchar(x) || anyNA(x)) {
    stop("x must be a single variable name", call. = FALSE)
  }
  period <- .sinusoids_period(period)
  harmonics <- .sinusoids_harmonics(harmonics)
  mult <- .sinusoids_mult(mult)
  ref_value <- .sinusoids_ref(ref_value, "ref_value")
  ref_mult <- .sinusoids_ref(ref_mult, "ref_mult")
  n_beta <- 2L * length(harmonics)
  methods::new(
    "sinusoids",
    name = x,
    label = .sinusoids_label(x, mult),
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    ref_value = ref_value,
    period = period,
    harmonics = harmonics,
    mult = mult,
    ref_mult = ref_mult,
    init = rep_len(init, n_beta),
    lower = rep_len(lower, n_beta),
    upper = rep_len(upper, n_beta),
    parscale = rep_len(parscale, n_beta)
  )
}

#' @describeIn sinusoids-class Cosine and sine columns, scaled by
#'   \code{mult - ref_mult} when \code{mult} is set.
#' @param term A \code{sinusoids} term.
#' @param data Data frame containing \code{x} and, when set, \code{mult}.
#' @export
setMethod("design", "sinusoids", function(term, data) {
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  D <- sinusoid_basis(
    data[[term@name]],
    period = term@period,
    harmonics = term@harmonics,
    ref_value = term@ref_value
  )
  if (.sinusoids_has_mult(term@mult)) {
    if (!term@mult %in% names(data)) {
      stop("column '", term@mult, "' not found", call. = FALSE)
    }
    scale <- as.numeric(data[[term@mult]]) - term@ref_mult
    if (any(!is.finite(scale))) {
      stop(
        "sinusoids multiplier '", term@mult, "' has non-finite values",
        call. = FALSE
      )
    }
    D <- D * scale
  }
  colnames(D) <- .sinusoids_colnames(term@label, term@harmonics)
  D
})

#' @describeIn sinusoids-class One fixed-effect row per cosine and sine column.
#' @export
setMethod("beta_info", "sinusoids", function(term, data) {
  the_colnames <- colnames(design(term, data))
  n <- length(the_colnames)
  data.frame(
    term = term@name,
    model = "sinusoids",
    label = term@label,
    order = NA,
    beta_label = the_colnames,
    init = rep_len(term@init, n),
    lower = rep_len(term@lower, n),
    upper = rep_len(term@upper, n),
    parscale = rep_len(term@parscale, n)
  )
})

#' @rdname replace_vars
#' @export
setMethod("replace_vars", "sinusoids", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  x <- methods::callNextMethod(x, map)
  if (.sinusoids_has_mult(x@mult) && x@mult %in% names(map)) {
    x@mult <- as.character(map[[x@mult]])
  }
  x@label <- .sinusoids_label(x@name, x@mult)
  x
})

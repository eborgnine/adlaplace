#' Posterior mean, SD, and simulations of a matern field on a raster
#'
#' Evaluates the fitted Matern GRF on `eval_grid` under the Laplace
#' approximation for the random effects (outer parameters fixed at the MLE).
#' Returns a SpatRaster with layers `mean`, `sd`, and optionally `sim1` ... `simn`.
#'
#' @param fit An `"adlaplace_fit"` from [adlaplace::adlaplace()] that includes
#'   exactly one [matern()] term.
#' @param eval_grid A terra `SpatRaster` whose cell centers are evaluation sites.
#' @param n Number of posterior simulations of the field. Default `0` (mean and
#'   sd only).
#' @return A `SpatRaster` with `n + 2` layers: `mean`, `sd`, and `sim1`...`simn`
#'   when `n > 0`.
#' @export
matern_est <- function(fit, eval_grid, n = 0L) {
  if (!inherits(fit, "adlaplace_fit")) {
    stop("fit must be an adlaplace_fit object", call. = FALSE)
  }
  if (!inherits(eval_grid, "SpatRaster")) {
    stop("eval_grid must be a terra SpatRaster", call. = FALSE)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "terra is required for matern_est(); install.packages(\"terra\")",
      call. = FALSE
    )
  }
  n <- as.integer(n)[1L]
  if (is.na(n) || n < 0L) {
    stop("n must be a non-negative integer", call. = FALSE)
  }

  terms <- fit$model_data$terms
  mat_terms <- Filter(function(t) methods::is(t, "matern"), terms)
  if (length(mat_terms) != 1L) {
    stop("fit must contain exactly one matern term", call. = FALSE)
  }
  mat <- mat_terms[[1L]]

  ginfo <- fit$model_data$term_data$info$gamma
  idx <- which(ginfo$label == mat@label)
  if (!length(idx)) {
    stop("no random effects found for matern term '", mat@label, "'", call. = FALSE)
  }
  g_labels <- ginfo$gamma_label[idx]
  g_hat <- as.numeric(fit$gamma[g_labels])
  if (anyNA(g_hat)) {
    stop("missing matern coefficients in fit$gamma", call. = FALSE)
  }

  A <- fem_bspline(
    eval_grid,
    mat@fem$knots,
    degree = mat@fem$degree
  )$A
  mu <- as.numeric(A %*% g_hat)

  half <- adlaplace::laplace_half_H_inv(fit$details)
  B <- A %*% half[idx, , drop = FALSE]
  sdv <- sqrt(pmax(as.numeric(Matrix::rowSums(B * B)), 0))

  layers <- cbind(mean = mu, sd = sdv)
  if (n > 0L) {
    chol_inner <- fit$details$hessian$chol_inner
    if (is.null(chol_inner)) {
      stop(
        "fit$details must contain hessian$chol_inner ",
        "(from log_lik_laplace with deriv = TRUE)",
        call. = FALSE
      )
    }
    g_sims <- adlaplace::rmvnldl(
      n,
      mean = fit$details$inner_opt$solution,
      chol_prec = chol_inner
    )
    if (n == 1L) {
      g_sims <- matrix(g_sims, nrow = 1L)
    }
    u <- as.matrix(A %*% t(g_sims[, idx, drop = FALSE]))
    colnames(u) <- paste0("sim", seq_len(n))
    layers <- cbind(layers, u)
  }

  out <- terra::rast(eval_grid, nlyrs = ncol(layers))
  names(out) <- colnames(layers)
  terra::values(out) <- layers
  out
}

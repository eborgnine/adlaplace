#' Assemble FEM Matern precision Q2 or Q3 from Grams
#'
#' @param kappa,tau Positive SPDE parameters.
#' @param C,G,G2,G3 Sparse Grams from [fem_bspline()].
#' @param alpha `2` or `3`.
#' @return Sparse `dgCMatrix` precision.
#' @export
fem_precision <- function(kappa, tau, C, G, G2, G3 = NULL, alpha = 2L) {
  kappa <- as.numeric(kappa)
  tau <- as.numeric(tau)
  if (kappa <= 0 || tau <= 0) {
    stop("kappa and tau must be positive")
  }
  alpha <- as.integer(alpha)
  if (alpha == 2L) {
    Q <- tau^2 * (kappa^4 * C + 2 * kappa^2 * G + G2)
  } else if (alpha == 3L) {
    if (is.null(G3)) {
      stop("G3 is required for alpha = 3")
    }
    Q <- tau^2 * (kappa^6 * C + 3 * kappa^4 * G + 3 * kappa^2 * G2 + G3)
  } else {
    stop("alpha must be 2 or 3")
  }
  methods::as(
    methods::as(Matrix::drop0(Matrix::forceSymmetric(Q)), "generalMatrix"),
    "CsparseMatrix"
  )
}

#' Align Gram values onto an upper-triangle CSC pattern
#'
#' Both `M` and the target `(p, i)` pattern are CSC, so each nonzero gets a
#' column-major key and the two sets are matched in one pass. Reading
#' `M[row, col]` per pattern entry instead costs an S4 dispatch each time,
#' which runs into minutes for a few hundred thousand entries.
#' @keywords internal
#' @noRd
align_gram_to_pattern <- function(M, p, i, n) {
  M <- methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
  x <- numeric(length(i))
  if (!length(i) || !length(M@x)) {
    return(x)
  }
  # rep.int over diff() spreads column indices across nonzeros and handles
  # empty CSC columns (repeated pointers) without a special case.
  m_counts <- diff(M@p)
  m_col <- rep.int(seq_along(m_counts) - 1L, m_counts)
  upper <- M@i <= m_col
  if (!any(upper)) {
    return(x)
  }
  p_counts <- diff(p)
  p_col <- rep.int(seq_along(p_counts) - 1L, p_counts)
  # as.numeric keeps keys exact when n * n would overflow integer
  pos <- match(
    M@i[upper] + as.numeric(m_col[upper]) * n,
    i + as.numeric(p_col) * n
  )
  # A nonzero of M outside the pattern is ignored, matching the per-entry
  # read this replaces (it only ever visited pattern positions).
  found <- !is.na(pos)
  x[pos[found]] <- M@x[upper][found]
  x
}

#' Build upper-triangle CSC of a structural pattern
#' @keywords internal
#' @noRd
upper_csc_pattern <- function(S) {
  S <- methods::as(methods::as(Matrix::forceSymmetric(S), "generalMatrix"), "CsparseMatrix")
  n <- nrow(S)
  # Keep only upper triangle (row <= col) including diagonal
  U <- Matrix::triu(S)
  U <- methods::as(methods::as(U, "generalMatrix"), "CsparseMatrix")
  list(p = as.integer(U@p), i = as.integer(U@i), n = n)
}

#' Build `random_fem_ssq_*` / `random_fem_det_*` precision payload for adlaplace
#'
#' @param fem Result of [fem_bspline()] (or list with C, G, G2, optional G3).
#' @param alpha `2` or `3`.
#' @return List for `density_data@precision`: Grams, chol pattern, and Q CSC
#'   coefficients aligned for on-tape assembly.
#' @export
fem_precision_payload <- function(fem, alpha = 2L) {
  alpha <- as.integer(alpha)
  C <- fem$C
  G <- fem$G
  G2 <- fem$G2
  G3 <- if (alpha >= 3L) fem$G3 else NULL
  if (alpha >= 3L && is.null(G3)) {
    stop("fem$G3 is required for alpha = 3 (use degree >= 3 in fem_bspline)")
  }
  struct <- fem_Q_structure(C, G, G2, G3)
  if (any(diff(struct@p) == 0L)) {
    stop(
      "fem_precision_payload: Q structure has empty rows/columns",
      call. = FALSE
    )
  }
  for (nm in c("C", "G", "G2")) {
    Mx <- get(nm)
    if (any(!is.finite(Mx@x))) {
      stop(
        "fem_precision_payload: non-finite values in Gram matrix ", nm,
        call. = FALSE
      )
    }
  }
  if (!is.null(G3) && any(!is.finite(G3@x))) {
    stop(
      "fem_precision_payload: non-finite values in Gram matrix G3",
      call. = FALSE
    )
  }
  chol <- fem_chol_pattern(struct)
  pat <- upper_csc_pattern(struct)
  out <- list(
    C = C,
    G = G,
    G2 = G2,
    chol = chol,
    alpha = alpha,
    Q_p = pat$p,
    Q_i = pat$i,
    C_x = align_gram_to_pattern(C, pat$p, pat$i, pat$n),
    G_x = align_gram_to_pattern(G, pat$p, pat$i, pat$n),
    G2_x = align_gram_to_pattern(G2, pat$p, pat$i, pat$n)
  )
  if (!is.null(G3)) {
    out$G3 <- G3
    out$G3_x <- align_gram_to_pattern(G3, pat$p, pat$i, pat$n)
  } else {
    out$G3_x <- rep(0, length(pat$i))
  }
  out
}

#' Slots of a FEM quadratic form inside the joint inverse Hessian
#'
#' One 0-based index into \code{H_inv@x} per stored upper-triangle entry of
#' \code{Q}, in that CSC order. \code{-1L} marks a pair the symbolic inverse
#' does not store.
#'
#' @param shard A \code{density_data} shard.
#' @param H_inv Symbolic upper triangle of the joint inverse Hessian.
#' @return An integer vector, or \code{NULL} when \code{shard} is not a
#'   \code{random_fem_ssq_*} density.
#' @rdname hinv_trace_index
#' @export
setMethod("hinv_trace_index", "density_data", function(shard, H_inv) {
  if (length(shard@density) != 1L ||
    !grepl("^random_fem_ssq_", shard@density)) {
    return(NULL)
  }
  prec <- shard@precision
  if (is.null(prec) || is.null(prec$Q_p) || is.null(prec$Q_i)) {
    stop("random_fem_ssq_ shard is missing precision$Q_p / Q_i", call. = FALSE)
  }
  Q_p <- as.integer(prec$Q_p)
  Q_i <- as.integer(prec$Q_i)
  n_local <- length(Q_p) - 1L
  if (n_local < 1L || ncol(shard@gamma_map) != n_local) {
    stop(
      "gamma_map ncol must equal the Q dimension (", n_local, ")",
      call. = FALSE
    )
  }
  nnz_col <- diff(shard@gamma_map@p)
  if (any(nnz_col != 1L)) {
    stop("gamma_map must have exactly one nonzero per column", call. = FALSE)
  }
  # @i and @p are 0-based. One entry per column, in column order.
  g_row <- shard@gamma_map@i[shard@gamma_map@p[-length(shard@gamma_map@p)] + 1L]
  if (!length(Q_i)) {
    return(integer(0))
  }
  q_counts <- diff(Q_p)
  q_col <- rep.int(seq_along(q_counts) - 1L, q_counts)
  gr <- g_row[Q_i + 1L]
  gc <- g_row[q_col + 1L]
  lo <- pmin(gr, gc)
  hi <- pmax(gr, gc)
  H_inv <- methods::as(H_inv, "CsparseMatrix")
  n <- nrow(H_inv)
  h_counts <- diff(H_inv@p)
  h_col <- rep.int(seq_along(h_counts) - 1L, h_counts)
  h_key <- H_inv@i + as.numeric(h_col) * n
  q_key <- lo + as.numeric(hi) * n
  pos <- match(q_key, h_key)
  as.integer(ifelse(is.na(pos), -1L, pos - 1L))
})

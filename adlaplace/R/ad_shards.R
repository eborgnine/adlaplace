#' Partition observations into AD shards by sparsity pattern
#'
#' The function partitions columns into \code{num_shards} shards using
#' quantiles of the first right singular vector, optionally using \pkg{RSpectra}
#' for efficiency when available.
#'
#' The resulting grouping is returned as a sparse matrix whose columns
#' correspond to shards and whose entries are the singular-vector loadings.
#' Shards are ordered from most heterogeneous to most homogeneous.
#'
#' @param A Random-effects design matrix (\code{nrow(A)} = number of observations).
#' @param elgm_matrix A numeric matrix (or matrix-like object) for extended latent gaussian models.
#' @param num_shards Integer giving the maximum number of shards to construct.
#'   The actual number of shards may be smaller if fewer distinct loadings
#'   are present.
#' @param min_groups Integer giving the minimum number of shards.
#'
#' @details
#' If the \pkg{RSpectra} package is available, the leading singular vector
#' is computed using \code{RSpectra::svds}; otherwise, a full singular value
#' decomposition via \code{\link[base]{svd}} is used.
#'
#' Shard boundaries are defined by empirical quantiles of the loadings.
#' Shards are subsequently reordered so that shards with larger within-shard
#' variability appear first.
#'
#' @return
#' A sparse matrix of class \code{"dgCMatrix"} (from \pkg{Matrix}), with
#' one column per shard and one row per observation. Nonzero entries
#' correspond to singular-vector loadings.
#'
#' @examples
#' set.seed(1)
#' A <- matrix(rnorm(100), 20, 5)
#' G <- ad_shards(A, num_shards = 3)
#' G
#'
#' Default observation shard map (one column, all observations)
#'
#' @param n_obs Number of observations (\code{length(y)}).
#' @return A \code{dgCMatrix} with one shard column and structural ones mapping
#'   each observation row to that shard (0-based row indices).
#' @keywords internal
default_obs_shards <- function(n_obs) {
  n_obs <- as.integer(n_obs)[1L]
  if (is.na(n_obs) || n_obs < 1L) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, 1L),
      index1 = FALSE
    ))
  }
  Matrix::sparseMatrix(
    i = seq.int(0L, n_obs - 1L),
    j = rep(0L, n_obs),
    x = rep(1, n_obs),
    dims = c(n_obs, 1L),
    index1 = FALSE
  )
}

#' @export
ad_shards <- function(A, elgm_matrix, num_shards, min_groups = 0) {
  ATp <- Matrix::t(A)
  if (missing(num_shards)) {
    num_shards <- ncol(ATp)
  }

  # to do:  implement min_groups
  if (!missing(elgm_matrix)) {
    ATp_t <- methods::as(Matrix::Matrix(ATp), "TsparseMatrix")
    ATp_t <- data.frame(row = ATp_t@j, gamma = ATp_t@i)
    elgm_matrix_t <- methods::as(elgm_matrix, "TsparseMatrix")
    elgm_matrix_t <- data.frame(row = elgm_matrix_t@i, strata = elgm_matrix_t@j)

    data.table::setDT(ATp_t)
    data.table::setDT(elgm_matrix_t)

    A_elgm_merge <- merge(ATp_t, elgm_matrix_t)
    A_elgm_merge <- A_elgm_merge[, c("gamma", "strata")]
    A_elgm_merge <- A_elgm_merge[!duplicated(A_elgm_merge), ]
    A_elgm_merge <- A_elgm_merge[order(A_elgm_merge$strata, A_elgm_merge$gamma), ]

    ATp <- Matrix::sparseMatrix(
      i = A_elgm_merge$gamma,
      j = A_elgm_merge$strata,
      x = rep(1.0, nrow(A_elgm_merge)),
      index1 = FALSE, dims = c(nrow(ATp), ncol(elgm_matrix))
    )
  }

  if (inherits(ATp, "ngCMatrix")) {
    ATp <- methods::as(ATp, "dMatrix")
  }
  # RSpectra::svds expects numeric Matrix input in a compressed sparse form.
  ATp <- methods::as(methods::as(ATp, "dMatrix"), "CsparseMatrix")
  if (requireNamespace("RSpectra", quietly = TRUE) && (
    max(dim(ATp)) > 100
  )) {
    loadings <- RSpectra::svds(ATp, k = 1)$v[, 1]
  } else {
    if (max(dim(ATp) > 1e5)) {
      warning("ATp matrix is very large, consider installing the RSpectra package")
    }
    sv <- svd(ATp)
    loadings <- sv$v[, 1]
  }

  uniqueLoadings <- sort(unique(loadings))
  if (length(uniqueLoadings) < min_groups) {
    the_min_diff <- abs(diff(uniqueLoadings))
    if (length(the_min_diff) == 0) {
      the_min_diff <- 1
    } else {
      the_min_diff <- min(the_min_diff) / 2
    }
    loadings <- loadings + stats::runif(length(loadings), -the_min_diff, the_min_diff)
    groupCut <- stats::quantile(loadings, seq(0, 1, len = min_groups + 1))
    groupCut[1] <- groupCut[1] - 1
    groupCut[length(groupCut)] <- groupCut[length(groupCut)] + 1
  } else if (length(uniqueLoadings) <= num_shards) {
    num_shards <- length(uniqueLoadings)
    groupCut <- uniqueLoadings[-1] - diff(uniqueLoadings) / 2
    groupCut <- c(uniqueLoadings[1] - 1, groupCut, uniqueLoadings[length(uniqueLoadings)] + 1)
  } else {
    groupCut <- unique(stats::quantile(uniqueLoadings, seq(0, 1, len = num_shards + 1)))
    groupCut[1] <- groupCut[1] - 1
    groupCut[length(groupCut)] <- groupCut[length(groupCut)] + 1
  }
  loadingsCut <- cut(loadings, groupCut)
  loadingsCut2 <- factor(loadingsCut, names(sort(table(loadingsCut), decreasing = TRUE)))
  theJ <- as.integer(loadingsCut2) - 1L

  theSd <- pmax(tapply(loadings, theJ, stats::sd), 0, na.rm = TRUE)
  orderSd <- order(theSd, decreasing = TRUE) - 1L
  theJ2 <- match(theJ, orderSd) - 1L

  groupMat <- Matrix::sparseMatrix(
    i = seq(0, len = length(loadings)),
    j = theJ2,
    x = loadings,
    index1 = FALSE
  )
  groupMat
}

#' Per-shard layout for AD density evaluation
#'
#' @slot y Response vector.
#' @slot ATp Transpose of random-effects design matrix.
#' @slot XTp Transpose of fixed-effects design matrix.
#' @slot beta_map Beta parameter map (\code{ngCMatrix}).
#' @slot gamma_map Gamma parameter map (\code{ngCMatrix}); at most one
#'   structural nonzero per row and column.
#' @slot theta_map Theta parameter map (\code{ngCMatrix}).
#' @slot elgm_matrix Optional exposure-lag map (\code{ngCMatrix}; empty by default).
#' @slot ad_fun Registered AD density name for this shard.
#' @slot ad_kind Shard kind (\code{"observations"}, \code{"parameters"}, \code{"random"}).
#' @slot precision Optional precision payload (any R object).
#' @importClassesFrom Matrix Matrix ngCMatrix
#' @exportClass ad_data
setClass(
  "ad_data",
  slots = c(
    y = "numeric",
    ATp = "Matrix",
    XTp = "Matrix",
    beta_map = "ngCMatrix",
    gamma_map = "ngCMatrix",
    theta_map = "ngCMatrix",
    elgm_matrix = "ngCMatrix",
    ad_fun = "character",
    ad_kind = "character",
    precision = "ANY"
  )
)

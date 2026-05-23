#' Model object from \code{model_setup()} or \code{ad_model()}
#'
#' @slot data Original model \code{data.frame}.
#' @slot y Response vector.
#' @slot ATp Transpose of random-effects design matrix.
#' @slot XTp Transpose of fixed-effects design matrix.
#' @slot beta_map Beta parameter map (\code{ngCMatrix}).
#' @slot gamma_map Gamma parameter map (\code{ngCMatrix}); at most one
#'   structural nonzero per row and column.
#' @slot theta_map Theta parameter map (\code{ngCMatrix}).
#' @slot elgmMatrix Optional exposure–lag map (\code{ngCMatrix}; empty by default).
#' @slot terms Named list of model term objects.
#' @slot info List with \code{beta}, \code{gamma}, \code{theta}, \code{parameters}.
#' @importClassesFrom Matrix Matrix ngCMatrix
#' @exportClass ad_model
setClass(
  "ad_model",
  slots = c(
    data = "data.frame",
    y = "numeric",
    ATp = "Matrix",
    XTp = "Matrix",
    beta_map = "ngCMatrix",
    gamma_map = "ngCMatrix",
    theta_map = "ngCMatrix",
    elgmMatrix = "ngCMatrix",
    terms = "list",
    info = "list"
  )
)

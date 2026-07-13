#' @importFrom methods setGeneric setMethod setOldClass as
#' @importFrom Matrix Matrix drop0 kronecker tcrossprod Diagonal sparseMatrix
#' @importFrom splines splineDesign
NULL

# Allow S4 methods on data.frame (e.g. expand.grid output)
methods::setOldClass("data.frame")

#' AD function with Hessian templates attached
#'
#' @slot ptr Raw combined handle (\code{ad_fun_ptr}).
#' @slot group_sparsity Per-group inner gradient sparsity patterns.
#' @slot outer Outer Hessian template (\code{dgCMatrix}).
#' @slot inner Inner-gamma Hessian template (\code{dgCMatrix}).
#' @slot map_outer Outer Hessian shard map.
#' @slot map_inner Inner Hessian shard map.
#' @slot parallel_map Thread-affinity map (\code{ngCMatrix}); rows are shards,
#'   columns are threads, with one nonzero per shard row.
#' @slot chol_inner Symbolic LDL factor or empty sparse matrix.
#' @slot chol_inner_list Numeric LDL list for C++ (\code{L1}, \code{Linv},
#'   \code{perm}, \code{perm_inv}, \code{half_H_inv}, \code{H_inv}).
#' @slot sizes Named numeric vector \code{beta}/\code{gamma}/\code{theta}.
#' @slot info List of parameter metadata (\code{beta}, \code{gamma},
#'   \code{theta}, \code{parameters}); populated from \code{model_data()$data$info}
#'   when the handle is built from a model-data bundle, otherwise empty.
#' @importClassesFrom Matrix dgCMatrix ngCMatrix
#' @exportClass ad_fun
setClass(
  "ad_fun",
  slots = c(
    ptr = "ad_fun_ptr",
    group_sparsity = "list",
    outer = "dgCMatrix",
    inner = "dgCMatrix",
    map_outer = "list",
    map_inner = "list",
    parallel_map = "ngCMatrix",
    chol_inner = "ANY",
    chol_inner_list = "list",
    sizes = "numeric",
    info = "list"
  ),
  prototype = list(
    info = list()
  )
)

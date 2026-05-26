#' @keywords internal
new_ad_fun_from_ptr <- function(ptr) {
  if (!is(ptr, "ad_fun_ptr")) {
    stop("ptr must be an ad_fun_ptr external pointer")
  }
  sz <- get_sizes(ptr, 0L)
  n_beta <- sz$n_beta
  n_theta <- sz$n_theta
  n_gamma <- sz$n_outer - n_beta - n_theta
  sparsity <- lapply(
    seq_len(n_groups(ptr)) - 1L,
    function(g) c(get_sizes(ptr, g), get_sparse_pattern(ptr, g))
  )
  hessian_pack <- hessian_map(
    sparsity_list = sparsity,
    Nbeta = n_beta,
    Ngamma = n_gamma,
    Ntheta = n_theta
  )
  adlaplace_attach_hessian(ptr, hessian_pack)

  methods::new(
    "ad_fun",
    ptr = ptr,
    group_sparsity = lapply(sparsity, function(xx) xx$grad_inner),
    outer = hessian_pack$outer,
    inner = hessian_pack$inner,
    map_outer = hessian_pack$map_outer,
    map_inner = hessian_pack$map_inner,
    chol_inner = hessian_pack$chol_inner,
    chol_inner_list = hessian_pack$chol_inner_list,
    sizes = hessian_pack$sizes
  )
}

#' Build AD function with Hessian templates
#'
#' @param x One of: an \code{ad_fun_ptr}, an \code{ad_data} shard with
#'   \code{ad_kind}/\code{ad_fun} set, or a \code{model_data()} bundle (list
#'   with \code{observations}, \code{random}, \code{parameters}).
#' @param config Model configuration list. Required for \code{ad_data} and
#'   \code{model_data}; unused for \code{ad_fun_ptr}.
#' @param ... Unused.
#' @return Object of class \code{ad_fun}.
#' @export
setGeneric("ad_fun", function(x, config = NULL, ...) {
  standardGeneric("ad_fun")
})

#' @describeIn ad_fun Attach Hessian templates to a raw \code{ad_fun_ptr}.
#' @export
setMethod("ad_fun", signature = c(x = "ad_fun_ptr"), function(x, ...) {
  new_ad_fun_from_ptr(x)
})

#' @describeIn ad_fun Build pointer from one \code{ad_data} shard, then attach templates.
#' @export
setMethod("ad_fun", signature = c(x = "ad_data"), function(x, config, ...) {
  if (missing(config) || is.null(config)) {
    stop("config is required for ad_fun(ad_data, config)", call. = FALSE)
  }
  ad_fun(ad_fun_ptr(x, config))
})

#' @describeIn ad_fun Build pointers from a \code{model_data()} bundle, merge,
#'   and attach templates.
#' @export
setMethod("ad_fun", signature = c(x = "list"), function(x, config, ...) {
  if (!is_model_data_bundle(x)) {
    stop(
      "list `x` must be a bundle from model_data() with ",
      "`observations`, `random`, and `parameters`",
      call. = FALSE
    )
  }
  if (missing(config) || is.null(config)) {
    stop("config is required for ad_fun(model_data, config)", call. = FALSE)
  }
  shards <- unname(c(x$observations, x$random, x$parameters))
  ptrs <- lapply(shards, ad_fun_ptr, config = config)
  ad_fun(do.call(c, ptrs))
})

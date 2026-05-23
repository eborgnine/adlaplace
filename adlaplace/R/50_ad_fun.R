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
    function(g) {
      c(
        get_sizes(ptr, g),
        get_sparse_pattern(ptr, g)
      )
    }
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

#' @keywords internal
build_ad_fun_ptr <- function(model, config) {
  if (!is(model, "ad_model")) {
    stop("`model` must be an ad_model object from model_setup()")
  }
  validate_config_layout(model, config)

  shard_ptrs <- list()
  for (term in model@terms) {
    ad_name <- term@ad_fun
    if (is.na(ad_name) || !nzchar(ad_name)) {
      next
    }
    kind <- term@as_kind
    if (is.na(kind) || !nzchar(kind)) {
      warning("term ", term@term, " has ad_fun but no as_kind; skipping")
      next
    }

    ptr <- if (identical(kind, "observations")) {
      ad_fun_ptr(config, kind, ad_name, model)
    } else if (identical(kind, "parameters")) {
      ad_fun_ptr(config, kind, ad_name, model)
    } else if (identical(kind, "random")) {
      term_model <- random_term_model(term, model)
      precision <- random_shard_precision(term, model)
      if (is.null(term_model) || is.null(precision)) {
        next
      }
      ad_fun_ptr(config, kind, ad_name, term_model, precision = precision)
    } else {
      warning("unknown as_kind '", kind, "' for term ", term@term, "; skipping")
      next
    }
    shard_ptrs[[length(shard_ptrs) + 1L]] <- ptr
  }

  if (length(shard_ptrs) == 0L) {
    stop("no AD shards could be built from ad_model terms")
  }

  if (length(shard_ptrs) == 1L) {
    shard_ptrs[[1L]]
  } else {
    do.call(c, shard_ptrs)
  }
}

#' Build AD function with Hessian templates
#'
#' @param x An \code{ad_fun_ptr} or \code{ad_model}.
#' @param config Model configuration list. Required for \code{ad_model}; unused
#'   for \code{ad_fun_ptr} (layout sizes come from the handle).
#' @param ... Unused.
#' @return Object of class \code{ad_fun}.
#' @export
setGeneric("ad_fun", function(x, config = NULL, ...) {
  standardGeneric("ad_fun")
})

#' @describeIn ad_fun Attach templates to a combined raw handle
#' @export
setMethod(
  "ad_fun",
  signature = c(x = "ad_fun_ptr"),
  function(x, ...) {
    new_ad_fun_from_ptr(x)
  }
)

#' @describeIn ad_fun Build shards from ad_model, merge, attach templates
#' @export
setMethod(
  "ad_fun",
  signature = c(x = "ad_model"),
  function(x, config, ...) {
    if (missing(config) || is.null(config)) {
      stop("config is required for ad_fun(ad_model, config)")
    }
    ad_fun(build_ad_fun_ptr(x, config))
  }
)

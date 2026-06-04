#' @keywords internal
build_parallel_map <- function(n_shards, num_threads, owner_threads = NULL) {
  num_threads <- as.integer(num_threads)[1]
  if (is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be a positive integer", call. = FALSE)
  }
  if (n_shards < 1L) {
    stop("ad_fun pointer has no shards", call. = FALSE)
  }

  shard_ids <- seq_len(n_shards)
  if (is.null(owner_threads)) {
    thread_ids <- ((shard_ids - 1L) %% num_threads) + 1L
  } else {
    owner_threads <- as.integer(owner_threads)
    if (length(owner_threads) != n_shards) {
      stop("owner_threads must have length n_shards", call. = FALSE)
    }
    if (any(is.na(owner_threads) | owner_threads < 0L)) {
      stop("owner_threads must be non-negative integers", call. = FALSE)
    }
    max_owner <- max(owner_threads)
    if (max_owner + 1L > num_threads) {
      num_threads <- max_owner + 1L
    }
    thread_ids <- owner_threads + 1L
  }
  map <- Matrix::sparseMatrix(
    i = shard_ids,
    j = thread_ids,
    dims = c(n_shards, num_threads),
    index1 = TRUE,
    giveCsparse = TRUE
  )
  map
}

#' @keywords internal
new_ad_fun_from_ptr <- function(ptr, num_threads = 1L) {
  if (!is(ptr, "ad_fun_ptr")) {
    stop("ptr must be an ad_fun_ptr external pointer")
  }
  n_shards <- n_groups(ptr)
  owner_threads <- vapply(
    seq_len(n_shards) - 1L,
    function(g) get_thread_owner(ptr, g),
    integer(1)
  )
  parallel_map <- build_parallel_map(n_shards, num_threads, owner_threads)
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
    parallel_map = parallel_map,
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
#'   \code{model_data}; unused for \code{ad_fun_ptr}. For \code{model_data},
#'   \code{beta}, \code{gamma}, and \code{theta} are filled from \code{x$data$info}
#'   when omitted (only \code{shards}, \code{transform_theta}, etc. are needed).
#' @param num_threads Positive integer; number of thread columns in
#'   \code{parallel_map} used for shard thread affinity.
#' @param ... For \code{ad_fun_ptr} input, optional additional
#'   \code{ad_fun_ptr} shards to combine before attaching templates.
#' @return Object of class \code{ad_fun}.
#' @export
setGeneric("ad_fun", function(x, config = NULL, num_threads = 1L, ...) {
  standardGeneric("ad_fun")
})

#' @describeIn ad_fun Attach Hessian templates to a raw \code{ad_fun_ptr}.
#' @export
setMethod("ad_fun", signature = c(x = "ad_fun_ptr"), function(x, config = NULL, num_threads = 1L, ...) {
  dots <- list(...)
  extras <- dots
  if (!is.null(config)) {
    extras <- c(list(config), extras)
  }
  if (length(extras) > 0L) {
    if (!all(vapply(extras, function(ptr) is(ptr, "ad_fun_ptr"), logical(1)))) {
      stop("additional arguments must be ad_fun_ptr", call. = FALSE)
    }
    x <- do.call(c, c(list(x), extras))
  }
  new_ad_fun_from_ptr(x, num_threads = num_threads)
})

#' @describeIn ad_fun Build pointer from one \code{ad_data} shard, then attach templates.
#' @export
setMethod("ad_fun",
  signature = c(x = "ad_data"),
  function(x, config, num_threads = 1L, ...) {
    if (missing(config) || is.null(config)) {
      stop("config is required for ad_fun(ad_data, config)", call. = FALSE)
    }
    config_build <- modifyList(config, list(num_threads = as.integer(num_threads)[1]))
    ad_fun(ad_fun_ptr(x, config_build), num_threads = num_threads)
  }
)

#' @describeIn ad_fun Build pointers from a \code{model_data()} bundle, merge,
#'   and attach templates.
#' @export
setMethod("ad_fun", signature = c(x = "list"), function(x, config, num_threads = 1L, ...) {
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
  config_build <- modifyList(
    config,
    list(
      num_threads = as.integer(num_threads)[1],
      beta = x$data$info$beta$init,
      theta = x$data$info$theta$init,
      gamma = rep(0, nrow(x$data$info$gamma))
    )
  )
  if (identical(config_build$transform_theta, TRUE)) {
    config_build$theta <- log(config_build$theta)
  }
  ptrs <- lapply(shards, ad_fun_ptr, config = config_build)
  ad_fun(do.call(c, ptrs), num_threads = num_threads)
})

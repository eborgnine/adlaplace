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
new_ad_fun_from_ptr <- function(ptr, num_threads = 1L, info = list(), verbose = FALSE) {
  if (!is(ptr, "ad_fun_ptr")) {
    stop("ptr must be an ad_fun_ptr external pointer")
  }
  num_threads <- as.integer(num_threads)[1]
  if (is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be a positive integer", call. = FALSE)
  }
  n_shards <- n_groups(ptr)
  if (verbose) {
    cat(
      "ad_fun: attaching Hessian templates (",
      n_shards, " AD group(s), ",
      num_threads, " thread(s))...\n",
      sep = ""
    )
    flush.console()
  }
  if (verbose) {
    cat("  assigning owner threads...\n")
    flush.console()
  }
  assign_owner_threads(ptr, num_threads)
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
  if (verbose) {
    cat(
      "  collecting sparsity patterns (beta=", n_beta,
      ", gamma=", n_gamma, ", theta=", n_theta, ")...\n",
      sep = ""
    )
    flush.console()
  }
  shard_ids <- seq_len(n_shards) - 1L
  sparsity <- vector("list", n_shards)
  for (i in seq_along(shard_ids)) {
    g <- shard_ids[[i]]
    if (verbose && n_shards > 1L) {
      cat("    sparsity group ", i, "/", n_shards, "\n", sep = "")
      flush.console()
    }
    sparsity[[i]] <- c(get_sizes(ptr, g), get_sparse_pattern(ptr, g))
  }
  if (verbose) {
    cat("  building Hessian map...\n")
    flush.console()
  }
  hessian_pack <- hessian_map(
    sparsity_list = sparsity,
    Nbeta = n_beta,
    Ngamma = n_gamma,
    Ntheta = n_theta
  )
  chol_inner_list <- hessian_pack$chol_inner_list
  if (length(chol_inner_list) > 0L && !is.null(chol_inner_list$half_H_inv)) {
    if (verbose) {
      cat("  building trace column map...\n")
      flush.console()
    }
    hessian_pack$chol_inner_list$trace_columns <- trace_columns_from_pattern(
      group_sparsity = lapply(sparsity, function(xx) xx$grad_inner),
      n_beta = n_beta,
      n_gamma = n_gamma,
      half_H_inv_pat = chol_inner_list$half_H_inv
    )
    chol_inner_list <- hessian_pack$chol_inner_list
  }
  if (verbose) {
    cat("  attaching Hessian templates to C++ handle...\n")
    flush.console()
  }
  adlaplace_attach_hessian(ptr, hessian_pack)

  if (verbose) {
    cat("ad_fun: Hessian attach complete.\n")
    flush.console()
  }

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
    chol_inner_list = chol_inner_list,
    sizes = hessian_pack$sizes,
    info = info
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
#'   When \code{config$num_threads} is set, it overrides the \code{num_threads}
#'   argument for \code{ad_data} and \code{model_data} methods.
#' @param num_threads Positive integer; OpenMP thread count for \code{inner_opt}
#'   and parallel \code{trace_hinv_t}. Shards are assigned
#'   \code{owner_thread = shard_index \% num_threads} at attach time.
#'   Default \code{1L} (serial).
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
  verbose <- isTRUE(config[["verbose"]])
  if (!is.null(config)) {
    extras <- c(list(config), extras)
  }
  if (length(extras) > 0L) {
    if (!all(vapply(extras, function(ptr) is(ptr, "ad_fun_ptr"), logical(1)))) {
      stop("additional arguments must be ad_fun_ptr", call. = FALSE)
    }
    if (verbose) {
      cat("ad_fun: merging ", 1L + length(extras), " raw handle(s)...\n", sep = "")
      flush.console()
    }
    x <- do.call(c, c(list(x), extras))
  }
  new_ad_fun_from_ptr(x, num_threads = num_threads, verbose = verbose)
})

#' @describeIn ad_fun Build pointer from one \code{ad_data} shard, then attach templates.
#' @export
setMethod("ad_fun",
  signature = c(x = "ad_data"),
  function(x, config, num_threads = 1L, ...) {
    if (missing(config) || is.null(config)) {
      stop("config is required for ad_fun(ad_data, config)", call. = FALSE)
    }
    if (!is.null(config$num_threads)) {
      num_threads <- config$num_threads
    }
    ad_fun(
      ad_fun_ptr(x, config),
      num_threads = num_threads,
      config = config
    )
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
  if (!is.null(config$num_threads)) {
    num_threads <- config$num_threads
  }
  defaults <- list(
    beta = init_from_info_block(x$data$info$beta),
    theta = init_from_info_block(x$data$info$theta),
    gamma = rep(0, nrow(x$data$info$gamma))
  )
  config_build <- modifyList(defaults, config)
  if (identical(config_build$transform_theta, TRUE) && is.null(config$theta)) {
    config_build$theta <- apply_theta_log(
      x$data$info$theta,
      cols = "init",
      active = TRUE
    )$init
  }
  verbose <- isTRUE(config_build[["verbose"]])
  shard_list <- c(x$observations, x$random, x$parameters)
  shard_names <- names(shard_list)
  n_shards <- length(shard_list)
  if (verbose) {
    obs_groups <- if (!is.null(config_build$shards)) {
      ncol(config_build$shards)
    } else {
      NA_integer_
    }
    cat(
      "ad_fun: building ", n_shards, " density shard(s)",
      if (!is.na(obs_groups)) paste0(" (", obs_groups, " observation group(s) per obs shard)") else "",
      "...\n",
      sep = ""
    )
    flush.console()
  }
  ptrs <- vector("list", n_shards)
  for (i in seq_len(n_shards)) {
    shard <- shard_list[[i]]
    shard_name <- if (length(shard_names) >= i) shard_names[[i]] else NULL
    if (verbose) {
      cat(
        "  [", i, "/", n_shards, "] ",
        adlaplace:::ad_fun_shard_label(shard, shard_name),
        " (CppAD tape",
        if (identical(shard@ad_kind, "observations") && !is.null(config_build$shards)) {
          paste0(", ", ncol(config_build$shards), " groups")
        } else {
          ""
        },
        ")...\n",
        sep = ""
      )
      flush.console()
    }
    ptrs[[i]] <- ad_fun_ptr(shard, config = config_build)
    if (verbose) {
      cat("  [", i, "/", n_shards, "] done (", n_groups(ptrs[[i]]), " AD group(s)).\n", sep = "")
      flush.console()
    }
  }
  if (verbose) {
    cat("ad_fun: merging density handles...\n")
    flush.console()
  }
  new_ad_fun_from_ptr(
    do.call(c, ptrs),
    num_threads = num_threads,
    info = x$data$info,
    verbose = verbose
  )
})

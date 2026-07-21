#' @include RcppExports.R
#' @keywords internal
effective_num_threads <- function(num_threads) {
  num_threads <- as.integer(num_threads)[1]
  if (is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be a positive integer", call. = FALSE)
  }
  # Without OpenMP at compile time, ignore requested thread count.
  if (!isTRUE(has_openmp())) {
    return(1L)
  }
  num_threads
}

#' @keywords internal
build_parallel_map <- function(n_shards, num_threads, owner_threads = NULL) {
  num_threads <- effective_num_threads(num_threads)
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
    repr = "C"
  )
  map
}

#' @keywords internal
new_ad_fun_from_ptr <- function(ptr, num_threads = 1L, info = list(), verbose = FALSE) {
  if (!is(ptr, "ad_fun_ptr")) {
    stop("ptr must be an ad_fun_ptr external pointer")
  }
  num_threads <- effective_num_threads(num_threads)
  n_shards <- n_groups(ptr)
  if (verbose) {
    cat(
      "ad_fun: attaching Hessian templates (",
      n_shards, " AD group(s), ",
      num_threads, " thread(s))...\n",
      sep = ""
    )
    utils::flush.console()
  }
  if (verbose) {
    cat("  assigning owner threads...\n")
    utils::flush.console()
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
    utils::flush.console()
  }
  shard_ids <- seq_len(n_shards) - 1L
  sparsity <- vector("list", n_shards)
  for (i in seq_along(shard_ids)) {
    g <- shard_ids[[i]]
    if (verbose && n_shards > 1L) {
      cat("    sparsity group ", i, "/", n_shards, "\n", sep = "")
      utils::flush.console()
    }
    sparsity[[i]] <- c(get_sizes(ptr, g), get_sparse_pattern(ptr, g))
  }
  if (verbose) {
    cat("  building Hessian map...\n")
    utils::flush.console()
  }
  hessian_pack <- hessian_map(
    sparsity_list = sparsity,
    n_beta = n_beta,
    n_gamma = n_gamma,
    n_theta = n_theta
  )
  chol_inner_list <- hessian_pack$chol_inner_list
  if (length(chol_inner_list) > 0L && !is.null(chol_inner_list$half_H_inv)) {
    if (verbose) {
      cat("  building trace column map...\n")
      utils::flush.console()
    }
    hessian_pack$chol_inner_list$trace_columns <- trace_columns_from_pattern(
      group_sparsity = lapply(sparsity, function(shard) shard$grad_inner),
      n_beta = n_beta,
      n_gamma = n_gamma,
      half_H_inv_pat = chol_inner_list$half_H_inv
    )
    chol_inner_list <- hessian_pack$chol_inner_list
  }
  if (verbose) {
    cat("  attaching Hessian templates to C++ handle...\n")
    utils::flush.console()
  }
  adlaplace_attach_hessian(ptr, hessian_pack)

  if (verbose) {
    cat("ad_fun: Hessian attach complete.\n")
    utils::flush.console()
  }

  methods::new(
    "ad_fun",
    ptr = ptr,
    group_sparsity = lapply(sparsity, function(shard) shard$grad_inner),
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

#' @rdname ad_fun
#' @export
setMethod("ad_fun", signature = c(x = "ad_fun_ptr"), function(x, config = NULL, num_threads = 1L, ...) {
  extras <- list(...)
  verbose <- FALSE
  if (!is.null(config)) {
    if (is.list(config) && !is(config, "ad_fun_ptr")) {
      # a genuine config list: only verbose is used at attach time
      verbose <- isTRUE(config[["verbose"]])
    } else {
      # positional shorthand ad_fun(ptr1, ptr2, ...): treat as extra shard
      extras <- c(list(config), extras)
    }
  }
  if (length(extras) > 0L) {
    if (!all(vapply(extras, function(ptr) is(ptr, "ad_fun_ptr"), logical(1)))) {
      stop("additional arguments must be ad_fun_ptr", call. = FALSE)
    }
    if (verbose) {
      cat("ad_fun: merging ", 1L + length(extras), " raw handle(s)...\n", sep = "")
      utils::flush.console()
    }
    x <- do.call(c, c(list(x), extras))
  }
  new_ad_fun_from_ptr(x, num_threads = num_threads, verbose = verbose)
})

#' @rdname ad_fun
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

#' @rdname ad_fun
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
  config_build <- utils::modifyList(defaults, config)
  theta_info <- x$data$info$theta
  n_theta <- nrow(theta_info)
  force_no_log <- identical(config_build$transform_theta, FALSE)
  log_flags <- if (force_no_log) {
    rep(FALSE, n_theta)
  } else if (is.logical(config_build$transform_theta) &&
      length(config_build$transform_theta) == n_theta &&
      n_theta > 0L) {
    config_build$transform_theta
  } else if (n_theta > 0L && "log" %in% names(theta_info)) {
    flags <- theta_info$log
    flags[is.na(flags)] <- TRUE
    as.logical(flags)
  } else {
    rep(TRUE, n_theta)
  }
  config_build$transform_theta <- log_flags
  if (is.null(config$theta) && n_theta > 0L) {
    theta_for_log <- theta_info
    if (!"log" %in% names(theta_for_log)) {
      theta_for_log$log <- TRUE
    }
    theta_for_log$log <- log_flags
    config_build$theta <- apply_theta_log(
      theta_for_log,
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
    utils::flush.console()
  }
  ptrs <- vector("list", n_shards)
  for (i in seq_len(n_shards)) {
    shard <- shard_list[[i]]
    shard_name <- if (length(shard_names) >= i) shard_names[[i]] else NULL
    if (verbose) {
      cat(
        "  [", i, "/", n_shards, "] ",
        ad_fun_shard_label(shard, shard_name),
        " (CppAD tape",
        if (identical(shard@ad_kind, "observations") && !is.null(config_build$shards)) {
          paste0(", ", ncol(config_build$shards), " groups")
        } else {
          ""
        },
        ")...\n",
        sep = ""
      )
      utils::flush.console()
    }
    ptrs[[i]] <- ad_fun_ptr(shard, config = config_build)
    if (verbose) {
      cat("  [", i, "/", n_shards, "] done (", n_groups(ptrs[[i]]), " AD group(s)).\n", sep = "")
      utils::flush.console()
    }
  }
  if (verbose) {
    cat("ad_fun: merging density handles...\n")
    utils::flush.console()
  }
  new_ad_fun_from_ptr(
    do.call(c, ptrs),
    num_threads = num_threads,
    info = x$data$info,
    verbose = verbose
  )
})

#' C++ backend entry points
#'
#' Low-level entry points exposed to R. Accept an \code{ad_fun} S4 object
#' (via \code{@ptr}) or a raw \code{ad_fun_ptr}.
#'
#' @param ad_fun An \code{ad_fun} S4 object or \code{ad_fun_ptr} handle.
#' @param x Numeric parameter vector of length \code{Nparams}.
#' @param shards Optional integer vector of 0-based shard indices; \code{NULL} or
#'   \code{integer(0)} evaluates all shards.
#' @param negative Logical (default \code{TRUE}). If \code{TRUE}, return the
#'   **negative** log density \eqn{-\ell(x)} and its derivatives (minimization /
#'   \pkg{trustOptim} sign, consistent with \code{inner_opt()}). If \code{FALSE},
#'   return \eqn{\ell(x)}, \eqn{\nabla\ell}, and \eqn{\nabla^2\ell}.
#' @param inner Logical; if \code{TRUE}, evaluate inner-\eqn{\gamma}
#'   derivatives; if \code{FALSE}, evaluate outer derivatives at the full
#'   parameter vector.
#' @param verbose Logical; if \code{TRUE}, print threads, shards, and CppAD
#'   session phases.
#' @param parameters Numeric vector of length \code{Nbeta + Ntheta}.
#' @param gamma Numeric vector of length \code{Ngamma}.
#' @param LinvPt,LinvPtColumns Sparse factors for \code{trace_hinv_t()}.
#'
#' @section Sign convention:
#' With default \code{negative = TRUE}, \code{joint_log_dens()}, \code{grad()},
#' and \code{hessian()} match \code{inner_opt()} (negative log-density). Set
#' \code{negative = FALSE} for the joint log density and its derivatives at the
#' same \code{x}.
#'
#' @section Thread affinity:
#' These are **serial debug** APIs. They error if the handle was built with
#' \code{ad_fun(..., num_threads > 1)}. Use a plain \code{ad_fun_ptr} (or
#' \code{ad_fun(..., num_threads = 1)}), or \code{\link{clone_ad_fun_ptr}()}
#' before assigning multi-thread affinity. Parallel evaluation uses
#' \code{\link{fun_obj_fdfh}} / \code{\link{inner_opt}} /
#' \code{\link{log_lik_laplace}} on the affined handle.
#'
#' @name adlaplace_cpp
#' @return See individual functions.
NULL

#' @rdname adlaplace_cpp
#' @export
joint_log_dens <- function(ad_fun, x, shards = NULL, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_joint_log_dens`, ptr, x, shards, negative)
}

#' @rdname adlaplace_cpp
#' @export
grad <- function(ad_fun, x, shards = NULL, inner = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_grad`, ptr, x, shards, inner, negative)
}

#' @rdname adlaplace_cpp
#' @export
hessian <- function(ad_fun, x, shards = NULL, inner = FALSE, verbose = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_hessian`, ptr, x, shards, inner, verbose, negative)
}

#' @rdname adlaplace_cpp
#' @export
trace_hinv_t <- function(ad_fun, x, LinvPt, LinvPtColumns, verbose = FALSE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call("_adlaplace_trace_hinv_t", ptr, x, LinvPt, LinvPtColumns, verbose)
}

#' @rdname adlaplace_cpp
#' @export
fun_obj_fdfh <- function(ad_fun, parameters, gamma, inner = TRUE, verbose = FALSE) {
  .Call(`_adlaplace_fun_obj_fdfh`, parameters, gamma, ad_fun, inner, verbose)
}

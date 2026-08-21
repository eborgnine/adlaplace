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
match_reorder_shards <- function(reorder_shards) {
  reorder_shards <- as.character(reorder_shards)[1L]
  match.arg(reorder_shards, c("none", "gradient", "hessian", "third"))
}

#' Longest-processing-time thread assignment
#'
#' @param costs Numeric vector of per-shard costs (length n_shards).
#' @param num_threads Positive integer thread count.
#' @return Integer vector of length n_shards with 0-based owner thread ids.
#' @keywords internal
lpt_assign_owners <- function(costs, num_threads) {
  costs <- as.numeric(costs)
  num_threads <- as.integer(num_threads)[1L]
  n <- length(costs)
  if (n < 1L) {
    return(integer(0))
  }
  if (is.na(num_threads) || num_threads < 1L) {
    stop("num_threads must be a positive integer", call. = FALSE)
  }
  owners <- integer(n)
  loads <- numeric(num_threads)
  ord <- order(costs, decreasing = TRUE, na.last = TRUE)
  for (idx in ord) {
    t <- which.min(loads) - 1L
    owners[[idx]] <- t
    c_i <- costs[[idx]]
    if (is.na(c_i)) {
      c_i <- 0
    }
    loads[[t + 1L]] <- loads[[t + 1L]] + c_i
  }
  owners
}

#' @keywords internal
profile_shard_eval_times <- function(
  ptr,
  x,
  eval_fun,
  n_rep = 2L,
  n_warmup = 1L
) {
  n_rep <- as.integer(n_rep)[1L]
  n_warmup <- as.integer(n_warmup)[1L]
  if (is.na(n_rep) || n_rep < 1L) {
    stop("n_rep must be a positive integer", call. = FALSE)
  }
  if (is.na(n_warmup) || n_warmup < 0L) {
    stop("n_warmup must be non-negative", call. = FALSE)
  }
  # One GC before the timing loop (not per shard; avoids system.time gcFirst).
  invisible(gc())
  n_shards <- n_groups(ptr)
  times <- numeric(n_shards)
  for (s in seq_len(n_shards) - 1L) {
    if (n_warmup > 0L) {
      for (w in seq_len(n_warmup)) {
        eval_fun(ptr, x, s)
      }
    }
    elapsed <- numeric(n_rep)
    for (r in seq_len(n_rep)) {
      t0 <- proc.time()[["elapsed"]]
      eval_fun(ptr, x, s)
      elapsed[[r]] <- proc.time()[["elapsed"]] - t0
    }
    times[[s + 1L]] <- mean(elapsed)
  }
  times
}

#' @keywords internal
profile_shard_gradient_times <- function(ptr, x, n_rep = 2L, n_warmup = 0L) {
  profile_shard_eval_times(
    ptr,
    x,
    eval_fun = function(ptr, x, s) {
      grad(ptr, x, ad_shards = s, inner = FALSE)
    },
    n_rep = n_rep,
    n_warmup = n_warmup
  )
}

#' @keywords internal
profile_shard_hessian_times <- function(ptr, x, n_rep = 2L, n_warmup = 1L) {
  profile_shard_eval_times(
    ptr,
    x,
    eval_fun = function(ptr, x, s) {
      hessian(ptr, x, ad_shards = s, inner = FALSE)
    },
    n_rep = n_rep,
    n_warmup = n_warmup
  )
}

# Overrides RcppExports binding (Collate: ad_pack.R after RcppExports.R):
# one GC before C++ chrono timing, matching gradient/hessian helpers.
profile_shard_trace3_times <- function(
  handle,
  x,
  half_H_inv,
  trace_columns,
  n_rep = 4L,
  n_warmup = 1L
) {
  invisible(gc())
  .Call(
    `_adlaplace_profile_shard_trace3_times`,
    handle,
    x,
    half_H_inv,
    trace_columns,
    n_rep,
    n_warmup
  )
}

#' @keywords internal
build_parallel_map <- function(n_shards, num_threads, owner_threads = NULL) {
  num_threads <- effective_num_threads(num_threads)
  if (n_shards < 1L) {
    stop("ad_pack pointer has no shards", call. = FALSE)
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
new_ad_pack_from_ptr <- function(
  ptr,
  num_threads = 1L,
  info = list(),
  verbose = FALSE,
  reorder_shards = c("none", "gradient", "hessian", "third")
) {
  if (!is(ptr, "ad_pack_ptr")) {
    stop("ptr must be an ad_pack_ptr external pointer")
  }
  num_threads <- effective_num_threads(num_threads)
  reorder_shards <- match_reorder_shards(reorder_shards)
  n_shards <- n_groups(ptr)
  if (verbose) {
    cat(
      "ad_pack: attaching Hessian templates (",
      n_shards, " AD group(s), ",
      num_threads, " thread(s)",
      ", reorder_shards=", reorder_shards, ")...\n",
      sep = ""
    )
    utils::flush.console()
  }

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
    cat("ad_pack: Hessian attach complete.\n")
    utils::flush.console()
  }

  do_balance <- reorder_shards %in% c("gradient", "hessian", "third")
  if (num_threads <= 1L || n_shards <= 1L) {
    do_balance <- FALSE
  }

  owner_threads <- NULL
  if (do_balance) {
    n_params <- as.integer(sz$n_outer)
    x_profile <- rep(0, n_params)
    costs <- NULL
    if (identical(reorder_shards, "gradient")) {
      if (verbose) {
        cat("  profiling outer gradient times per shard...\n")
        utils::flush.console()
      }
      costs <- profile_shard_gradient_times(ptr, x_profile)
    } else if (identical(reorder_shards, "hessian")) {
      if (verbose) {
        cat("  profiling outer Hessian times per shard...\n")
        utils::flush.console()
      }
      costs <- profile_shard_hessian_times(ptr, x_profile)
    } else if (identical(reorder_shards, "third")) {
      if (
        !is.null(chol_inner_list$half_H_inv) &&
          !is.null(chol_inner_list$trace_columns)
      ) {
        if (verbose) {
          cat("  profiling Reverse3 (dummy direction) times per shard...\n")
          utils::flush.console()
        }
        half_pat <- methods::as(chol_inner_list$half_H_inv, "dgCMatrix")
        costs <- profile_shard_trace3_times(
          ptr,
          x_profile,
          half_pat,
          chol_inner_list$trace_columns,
          n_rep = 2L,
          n_warmup = 1L
        )
      } else if (verbose) {
        cat(
          "  reorder_shards='third' but no half_H_inv/trace_columns; ",
          "falling back to modulo owners.\n",
          sep = ""
        )
        utils::flush.console()
      }
    }
    if (!is.null(costs)) {
      owner_threads <- lpt_assign_owners(costs, num_threads)
      if (verbose) {
        cat(
          "  LPT owners: ",
          paste(owner_threads, collapse = ","),
          "\n",
          sep = ""
        )
        cat(
          "  shard times: ",
          paste(format(costs, digits = 3L), collapse = ","),
          "\n",
          sep = ""
        )
        utils::flush.console()
      }
    }
  }

  if (verbose) {
    cat("  assigning owner threads...\n")
    utils::flush.console()
  }
  if (is.null(owner_threads)) {
    assign_owner_threads(ptr, num_threads)
  } else {
    assign_owner_threads(ptr, num_threads, owners = owner_threads)
  }
  owner_threads <- vapply(
    seq_len(n_shards) - 1L,
    function(g) get_thread_owner(ptr, g),
    integer(1)
  )
  parallel_map <- build_parallel_map(n_shards, num_threads, owner_threads)

  methods::new(
    "ad_pack",
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
#' @param x One of: an \code{ad_pack_ptr}, an \code{density_data} shard with
#'   \code{ad_kind}/\code{ad_pack} set, or a \code{model_data()} bundle (list
#'   with \code{observations}, \code{random}, \code{parameters}).
#' @param config Model configuration list. Required for \code{density_data} and
#'   \code{model_data}; unused for \code{ad_pack_ptr}. For \code{model_data},
#'   \code{beta}, \code{gamma}, and \code{theta} are filled from \code{x$term_data$info}
#'   when omitted (only \code{shards}, \code{transform_theta}, etc. are needed).
#'   When \code{config$num_threads} is set, it overrides the \code{num_threads}
#'   argument for \code{density_data} and \code{model_data} methods.
#'   When \code{config$reorder_shards} is set, it overrides \code{reorder_shards}.
#' @param num_threads Positive integer; OpenMP thread count for \code{inner_opt}
#'   and parallel \code{trace_hinv_t}. Default \code{1L} (serial).
#' @param reorder_shards Character; how to assign OpenMP owner threads when
#'   \code{num_threads > 1} and there are multiple shards:
#'   \code{"none"} (default) uses \code{owner_thread = shard_index \% num_threads};
#'   \code{"gradient"} times each shard's outer
#'   \code{\link{grad}(..., inner = FALSE)} (mean of 2; no warmup), then
#'   longest-processing-time (LPT) assignment;
#'   \code{"hessian"} times each shard's outer
#'   \code{\link{hessian}(..., inner = FALSE)} (1 warmup + mean of 2), then LPT;
#'   \code{"third"} times a Reverse3 sweep with dummy directions matching the
#'   symbolic \code{half_H_inv} / \code{trace_columns} sparsity (no numeric
#'   Cholesky), then LPT.
#'   Physical shard order is unchanged; only \code{owner_thread} is rewritten.
#' @param ... For \code{ad_pack_ptr} input, optional additional
#'   \code{ad_pack_ptr} shards to combine before attaching templates.
#' @return Object of class \code{ad_pack}.
#' @export
setGeneric(
  "ad_pack",
  function(x, config = NULL, num_threads = 1L,
           reorder_shards = c("none", "gradient", "hessian", "third"), ...) {
    standardGeneric("ad_pack")
  }
)

#' @rdname ad_pack
#' @export
setMethod(
  "ad_pack",
  signature = c(x = "ad_pack_ptr"),
  function(x, config = NULL, num_threads = 1L,
           reorder_shards = c("none", "gradient", "hessian", "third"), ...) {
    extras <- list(...)
    verbose <- FALSE
    if (!is.null(config)) {
      if (is.list(config) && !is(config, "ad_pack_ptr")) {
        # a genuine config list: only verbose / reorder_shards used at attach
        verbose <- isTRUE(config[["verbose"]])
        if (!is.null(config[["reorder_shards"]])) {
          reorder_shards <- config[["reorder_shards"]]
        }
      } else {
        # positional shorthand ad_pack(ptr1, ptr2, ...): treat as extra shard
        extras <- c(list(config), extras)
      }
    }
    if (length(extras) > 0L) {
      if (!all(vapply(extras, function(ptr) is(ptr, "ad_pack_ptr"), logical(1)))) {
        stop("additional arguments must be ad_pack_ptr", call. = FALSE)
      }
      if (verbose) {
        cat("ad_pack: merging ", 1L + length(extras), " raw handle(s)...\n", sep = "")
        utils::flush.console()
      }
      x <- do.call(c, c(list(x), extras))
    }
    new_ad_pack_from_ptr(
      x,
      num_threads = num_threads,
      verbose = verbose,
      reorder_shards = reorder_shards
    )
  }
)

#' @rdname ad_pack
#' @export
setMethod(
  "ad_pack",
  signature = c(x = "density_data"),
  function(x, config, num_threads = 1L,
           reorder_shards = c("none", "gradient", "hessian", "third"), ...) {
    if (missing(config) || is.null(config)) {
      stop("config is required for ad_pack(density_data, config)", call. = FALSE)
    }
    if (!is.null(config$num_threads)) {
      num_threads <- config$num_threads
    }
    if (!is.null(config$reorder_shards)) {
      reorder_shards <- config$reorder_shards
    }
    ad_pack(
      ad_pack_ptr(x, config),
      num_threads = num_threads,
      reorder_shards = reorder_shards,
      config = config
    )
  }
)

#' @rdname ad_pack
#' @export
setMethod(
  "ad_pack",
  signature = c(x = "list"),
  function(x, config, num_threads = 1L,
           reorder_shards = c("none", "gradient", "hessian", "third"), ...) {
    if (!is_model_data_bundle(x)) {
      stop(
        "list `x` must be a bundle from model_data() with ",
        "`observations`, `random`, and `parameters`",
        call. = FALSE
      )
    }
    if (missing(config) || is.null(config)) {
      stop("config is required for ad_pack(model_data, config)", call. = FALSE)
    }
    if (!is.null(config$num_threads)) {
      num_threads <- config$num_threads
    }
    if (!is.null(config$reorder_shards)) {
      reorder_shards <- config$reorder_shards
    }
    defaults <- list(
      beta = init_from_info_block(x$term_data$info$beta),
      theta = init_from_info_block(x$term_data$info$theta),
      gamma = rep(0, nrow(x$term_data$info$gamma))
    )
    config_build <- utils::modifyList(defaults, config)
    theta_info <- x$term_data$info$theta
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
    config_build <- ensure_config_obs_groups(
      config_build,
      A = x$term_data$A,
      elgm_matrix = x$term_data$elgm_matrix,
      num_shards = config_build$num_shards,
      num_threads = num_threads
    )
    verbose <- isTRUE(config_build[["verbose"]])
    shard_list <- c(x$observations, x$random, x$parameters)
    shard_names <- names(shard_list)
    n_shards <- length(shard_list)
    if (verbose) {
      obs_groups <- if (!is.null(config_build$obs_groups)) {
        ncol(config_build$obs_groups)
      } else {
        NA_integer_
      }
      cat(
        "ad_pack: building ", n_shards, " density shard(s)",
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
          ad_shard_label(shard, shard_name),
          " (CppAD tape",
          if (identical(shard@ad_kind, "observations") && !is.null(config_build$obs_groups)) {
            paste0(", ", ncol(config_build$obs_groups), " groups")
          } else {
            ""
          },
          ")...\n",
          sep = ""
        )
        utils::flush.console()
      }
      ptrs[[i]] <- ad_pack_ptr(shard, config = config_build)
      if (verbose) {
        cat("  [", i, "/", n_shards, "] done (", n_groups(ptrs[[i]]), " AD group(s)).\n", sep = "")
        utils::flush.console()
      }
    }
    if (verbose) {
      cat("ad_pack: merging density handles...\n")
      utils::flush.console()
    }
    new_ad_pack_from_ptr(
      do.call(c, ptrs),
      num_threads = num_threads,
      info = x$term_data$info,
      verbose = verbose,
      reorder_shards = reorder_shards
    )
  }
)

#' C++ backend entry points
#'
#' Low-level entry points exposed to R. Accept an \code{ad_pack} S4 object
#' (via \code{@ptr}) or a raw \code{ad_pack_ptr}.
#'
#' @param ad_pack An \code{ad_pack} S4 object or \code{ad_pack_ptr} handle.
#' @param x Numeric parameter vector of length \code{Nparams}.
#' @param ad_shards Optional integer vector of 0-based shard indices; \code{NULL} or
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
#' \code{ad_pack(..., num_threads > 1)}. Use a plain \code{ad_pack_ptr} (or
#' \code{ad_pack(..., num_threads = 1)}), or \code{\link{clone_ad_pack_ptr}()}
#' before assigning multi-thread affinity. Parallel evaluation uses
#' \code{\link{fun_obj_fdfh}} / \code{\link{inner_opt}} /
#' \code{\link{log_lik_laplace}} on the affined handle.
#'
#' @name adlaplace_cpp
#' @return See individual functions.
NULL

#' @rdname adlaplace_cpp
#' @export
joint_log_dens <- function(ad_pack, x, ad_shards = NULL, negative = TRUE) {
  ptr <- if (isS4(ad_pack) && methods::.hasSlot(ad_pack, "ptr")) ad_pack@ptr else ad_pack
  .Call(`_adlaplace_joint_log_dens`, ptr, x, ad_shards, negative)
}

#' @rdname adlaplace_cpp
#' @export
grad <- function(ad_pack, x, ad_shards = NULL, inner = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_pack) && methods::.hasSlot(ad_pack, "ptr")) ad_pack@ptr else ad_pack
  .Call(`_adlaplace_grad`, ptr, x, ad_shards, inner, negative)
}

#' @rdname adlaplace_cpp
#' @export
hessian <- function(ad_pack, x, ad_shards = NULL, inner = FALSE, verbose = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_pack) && methods::.hasSlot(ad_pack, "ptr")) ad_pack@ptr else ad_pack
  .Call(`_adlaplace_hessian`, ptr, x, ad_shards, inner, verbose, negative)
}

#' @rdname adlaplace_cpp
#' @export
trace_hinv_t <- function(ad_pack, x, LinvPt, LinvPtColumns, verbose = FALSE) {
  ptr <- if (isS4(ad_pack) && methods::.hasSlot(ad_pack, "ptr")) ad_pack@ptr else ad_pack
  .Call("_adlaplace_trace_hinv_t", ptr, x, LinvPt, LinvPtColumns, verbose)
}

#' @rdname adlaplace_cpp
#' @export
fun_obj_fdfh <- function(ad_pack, parameters, gamma, inner = TRUE, verbose = FALSE) {
  .Call(`_adlaplace_fun_obj_fdfh`, parameters, gamma, ad_pack, inner, verbose)
}

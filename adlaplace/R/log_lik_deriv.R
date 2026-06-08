log_lik_deriv <- function(
  full_parameters,
  hessian_pack,
  grad,
  ad_fun,
  verbose = FALSE
) {
  if (!is(ad_fun, "ad_fun")) {
    stop("ad_fun must be an ad_fun object", call. = FALSE)
  }
  sz <- ad_fun@sizes
  if (length(sz) < 3L || any(is.na(sz[c("beta", "gamma", "theta")]))) {
    stop(
      "ad_fun@sizes must contain finite beta, gamma, and theta",
      call. = FALSE
    )
  }
  num_beta <- as.integer(sz["beta"])
  num_theta <- as.integer(sz["theta"])
  num_gamma <- as.integer(sz["gamma"])
  n_params <- num_beta + num_gamma + num_theta
  if (length(grad) != n_params) {
    stop(
      "length(grad) (", length(grad), ") must equal num_beta + num_gamma + ",
      "num_theta (", n_params, ")",
      call. = FALSE
    )
  }

  ldl <- as_ldl_list(hessian_pack$chol_inner)
  half_H_inv <- half_H_inv_from_ldl(ldl)
  Hstuff <- list(
    half_H_inv = half_H_inv,
    H_inv = Matrix::tcrossprod(half_H_inv)
  )

  seq_gamma1 <- seq.int(num_beta + 1L, length.out = num_gamma)
  seq_gamma0 <- seq_gamma1 - 1L

  sparsity <- ad_fun@group_sparsity
  which_columns_by_group1 <- lapply(
    sparsity,
    function(xx, refmat) {
      grad_inner_gamma <- match(xx, seq_gamma0)
      linv_here <- refmat[grad_inner_gamma, , drop = FALSE]
      which(diff(linv_here@p) > 0) - 1L
    },
    refmat = Hstuff$half_H_inv
  )

  which_columns_by_group <- Matrix::sparseMatrix(
    i = unlist(which_columns_by_group1),
    j = rep(
      seq(0, length.out = length(which_columns_by_group1)),
      unlist(lapply(which_columns_by_group1, length))
    ),
    index1 = FALSE,
    dims = c(num_gamma, length(which_columns_by_group1))
  )

  if (verbose) {
    message("log_lik_deriv: calling trace_hinv_t ...")
  }
  the_trace <- adlaplace::trace_hinv_t(
    ad_fun = ad_fun,
    x = full_parameters,
    LinvPt = Hstuff$half_H_inv,
    LinvPtColumns = which_columns_by_group,
    verbose = verbose
  )

  dU <- -Hstuff$H_inv %*% hessian_pack$outer[seq_gamma1, -seq_gamma1]

  result <- list(
    extra = c(
      list(dU = dU, trace3 = the_trace),
      Hstuff
    )
  )

  result$deriv <- data.frame(
    d_det_upart = -as.vector(the_trace[seq_gamma1] %*% dU),
    d_det_tpart = -the_trace[-seq_gamma1]
  )
  result$deriv$grad_theta <- -grad[-seq_gamma1]
  result$deriv$grad_u <- as.vector(-grad[seq_gamma1] %*% dU)
  result$deriv$d_det <- result$deriv$d_det_upart + result$deriv$d_det_tpart
  result$deriv$d_neg_log_lik <-
    result$grad <-
    -result$deriv$grad_theta +
    result$deriv$d_det - result$deriv$grad_u
  result$deriv$d_log_lik <- -result$deriv$d_neg_log_lik

  result
}

reformat_chol <- function(x) {
  ldl <- as_ldl_list(x)
  half_H_inv <- half_H_inv_from_ldl(ldl)
  H_inv <- Matrix::tcrossprod(half_H_inv)
  list(half_H_inv = half_H_inv, H_inv = H_inv)
}

log_lik_deriv <- function(
  full_parameters,
  hessian_pack,
  grad,
  config,
  ad_fun
) {
  Hstuff <- reformat_chol(hessian_pack$chol_inner)

  Sgamma1 <- seq.int(length(config$beta) + 1, length.out = length(config$gamma))
  Sgamma0 <- Sgamma1 - 1L

  whichColumnsByGroup1 <- lapply(
    ad_fun$group_sparsity, function(xx, refmat) {
      grad_inner_gamma <- match(xx, Sgamma0)
      linvHere <- refmat[grad_inner_gamma, , drop = FALSE]
      which(diff(linvHere@p) > 0) - 1L
    },
    refmat = Hstuff$half_H_inv
  )

  whichColumnsByGroup <- Matrix::sparseMatrix(
    i = unlist(whichColumnsByGroup1),
    j = rep(
      seq(0, length.out = length(whichColumnsByGroup1)),
      unlist(lapply(whichColumnsByGroup1, length))
    ),
    index1 = FALSE,
    dims = c(length(config$gamma), length(whichColumnsByGroup1))
  )

  theTrace <- adlaplace::traceHinvT(
    x = full_parameters,
    Hstuff$half_H_inv,
    whichColumnsByGroup,
    ad_fun = ad_fun,
    c(config$num_threads, 1L)[1]
  )

  dU <- -Hstuff$H_inv %*% hessian_pack$outer[Sgamma1, -Sgamma1]

  result <- list(
    extra = c(
      list(dU = dU, trace3 = theTrace),
      Hstuff
    )
  )

  result$deriv <- data.frame(
    d_det_upart = -as.vector(theTrace[Sgamma1] %*% dU),
    d_det_tpart = -theTrace[-Sgamma1]
  )
  result$deriv$grad_theta <- -grad[-Sgamma1]
  result$deriv$grad_u <- as.vector(-grad[Sgamma1] %*% dU)
  result$deriv$d_det <- result$deriv$d_det_upart + result$deriv$d_det_tpart
  result$deriv$d_neg_log_lik <-
    result$grad <-
    -result$deriv$grad_theta +
    result$deriv$d_det - result$deriv$grad_u
  result$deriv$d_log_lik <- -result$deriv$d_neg_log_lik

  result
}

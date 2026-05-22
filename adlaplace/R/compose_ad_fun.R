#' @keywords internal
is_adlaplace_handle <- function(x) {
  inherits(x, "adlaplace_handle_ptr")
}

#' Attach sparsity and hessian_map to an existing raw AD handle
#'
#' @keywords internal
attach_ad_fun_artifacts <- function(ad_ptr, config) {
  sparsity <- lapply(
    seq_len(n_groups(ad_ptr)) - 1L,
    function(g) {
      c(
        get_sparse_sizes(ad_ptr, g),
        get_sparse_pattern(ad_ptr, g)
      )
    }
  )
  hessian_pack <- hessian_map(
    sparsity_list = sparsity,
    Nbeta = length(config$beta),
    Ngamma = length(config$gamma),
    Ntheta = length(config$theta)
  )
  adlaplace_attach_hessian(ad_ptr, hessian_pack)

  out <- c(
    list(
      ad_fun = ad_ptr,
      group_sparsity = lapply(sparsity, function(xx) xx$grad_inner)
    ),
    hessian_pack
  )
  attr(out, "adlaplace.backend") <- "adlaplace"
  out
}

#' Build backend AD function handle
#'
#' Constructs grouped CppAD tapes and attaches Hessian sparsity templates
#' expected by \pkg{adlaplace}.
#'
#' @param data Model data list.
#' @param config Model configuration list.
#'
#' @return List with \code{ad_fun}, \code{group_sparsity}, and components from
#'   \code{\link[adlaplace]{hessian_map}}.
#'
#' @importFrom adlaplace n_groups get_sparse_sizes get_sparse_pattern hessian_map
#' @seealso \code{\link[adlaplace]{get_ad_fun}}
#' @export
getAdFun_r <- function(data, config) {
  ad_ptr <- get_ad_fun_raw_hpolcc(data, config)
  sparsity <- lapply(
    seq_len(adlaplace::n_groups(ad_ptr)) - 1L,
    function(g) {
      c(
        adlaplace::get_sparse_sizes(ad_ptr, g),
        adlaplace::get_sparse_pattern(ad_ptr, g)
      )
    }
  )
  hessian_pack <- adlaplace::hessian_map(
    sparsity_list = sparsity,
    Nbeta = length(config$beta),
    Ngamma = length(config$gamma),
    Ntheta = length(config$theta)
  )
  c(
    list(
      ad_fun = ad_ptr,
      group_sparsity = lapply(sparsity, function(xx) xx$grad_inner)
    ),
    hessian_pack
  )
}

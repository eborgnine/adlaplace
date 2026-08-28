# FEM covariance validation with hierarchical basis
# Append to fem-run chunk in adlaplace.Rmd after fem_corr_df is built

fem_cov_hierarchical <- function(sites_eval, knots_coarse, knots_refine) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(NULL)
  }
  coarse_r <- terra::rast(terra::ext(-0.2, 1.2, -0.2, 1.2), resolution = 0.2)
  fine_r <- terra::rast(terra::ext(0.2, 0.8, 0.2, 0.8), resolution = 0.1)
  fem <- adlaplaceFem::fem_bspline(
    sites_eval,
    list(coarse_r, fine_r),
    degree = 2L
  )
  Q <- adlaplaceFem::fem_precision(
    kappa, tau, fem$C, fem$G, fem$G2, alpha = 2L
  )
  A <- fem$A
  X <- Matrix::solve(Q, Matrix::t(A))
  list(
    cov = as.matrix(A %*% X),
    dof = ncol(A),
    label = "hierarchical"
  )
}

# Example usage in fem-run:
# fem_h <- fem_cov_hierarchical(sites_eval, knots_coarse, NULL)
# fem_corr_df <- rbind(fem_corr_df, data.frame(
#   dist = dist_vec, corr = fem_h$corr, grid = "hierarchical"
# ))

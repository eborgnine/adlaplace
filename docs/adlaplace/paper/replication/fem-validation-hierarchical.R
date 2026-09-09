# FEM covariance validation with hierarchical basis
# Append to fem-run chunk in adlaplace.Rmd after fem_corr_df is built

fem_cov_hierarchical <- function(sites_eval, knots_coarse, knots_refine) {
  kn <- adlaplaceFem::hb_knots(
    outer = list(
      xmin = -0.2, xmax = 1.2, ymin = -0.2, ymax = 1.2, resolution = 0.2
    ),
    inner = c(0.2, 0.8, 0.2, 0.8),
    fact = 2
  )
  fem <- adlaplaceFem::fem_bspline(
    sites_eval,
    kn,
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

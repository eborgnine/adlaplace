# Loaloa replication with hierarchical FEM knots
# Sourced from adlaplace.Rmd chunk loaloa-run

loaloa_ok <- FALSE

if (requireNamespace("adlaplaceFem", quietly = TRUE) &&
    requireNamespace("adlaplace", quietly = TRUE) &&
    requireNamespace("geostatsp", quietly = TRUE) &&
    requireNamespace("terra", quietly = TRUE)) {

  library(adlaplace)
  library(adlaplaceFem)
  library(geostatsp)
  library(terra)

  data("loaloa", package = "geostatsp")
  loaloa <- terra::unwrap(loaloa)
  elev_r <- terra::unwrap(elevationLoa)
  evi_r <- terra::unwrap(eviLoa)
  elev_s_r <- (elev_r - 500) / 100
  evi_s_r <- (evi_r - mean(terra::values(evi_r), na.rm = TRUE)) /
    sd(terra::values(evi_r), na.rm = TRUE)

  loaloa_for_fit <- geostatsp:::gm.dataSpatial(
    formula = y ~ elev + evi,
    data = loaloa,
    covariates = list(elev = elev_s_r, evi = evi_s_r),
    grid = geostatsp::squareRaster(loaloa, cells = 200, buffer = 1e4)
  )

  coarse <- geostatsp::squareRaster(loaloa, cells = 6, buffer = 300 * 1000)
  west <- terra::rast(
    terra::ext(c(3.9e5, 8.1e5, 3.7e5, 7.7e5)),
    resolution = 75 * 1000
  )
  east <- terra::rast(
    terra::ext(c(9.9e5, 1.19e6, 3.7e5, 7.7e5)),
    resolution = 75 * 1000
  )
  knots_grid <- list(coarse, list(west, east))
  hb_dof <- hb_summary(hb_knots(knots_grid, degree = 3L))

  loaloa_df <- as.data.frame(loaloa_for_fit$data, geom = "WKT")
  loaloa_model <- adlaplace::binomial(y, size = N) ~
    elev + evi +
    adlaplace::iid(villageID, init = 0.2) +
    adlaplaceFem::matern(
      geometry,
      knots = knots_grid,
      shape = 2L,
      init = c(100 * 1000, 0.2),
      lower = c(1 * 1000, 0.001)
    )

  fit_loaloa <- adlaplace::adlaplace(
    loaloa_model,
    data = loaloa_df,
    control = list(maxit = 200, trace = 0, REPORT = 0),
    config = list(num_threads = 2L, num_shards = 20L),
    verbose = FALSE
  )

  cf <- coef(fit_loaloa)
  ad_pars <- c(
    intercept = unname(cf[["intercept"]]),
    elev = unname(cf[["elev_linear"]]),
    evi = unname(cf[["evi_linear"]]),
    village_sd = unname(cf[["villageID_iid"]]),
    range = unname(cf[["geometry_matern_log_range"]]),
    sd = unname(cf[["geometry_matern_log_sd"]])
  )

  comparison_loaloa <- data.frame(
    parameter = names(ad_pars),
    adlaplace = unname(ad_pars),
    row.names = NULL
  )

  if (requireNamespace("INLA", quietly = TRUE)) {
    loaFit <- geostatsp::glgm(
      y ~ elev_s + evi_s +
        f(villageID, prior = "pc.prec", param = c(1, 0.5), model = "iid"),
      loaloa,
      grid = 100,
      covariates = list(elev_s = elev_s_r, evi_s = evi_s_r),
      family = "binomial",
      Ntrials = loaloa$N,
      shape = 2,
      buffer = 25000,
      prior = list(sd = 1, range = 500 * 1000),
      control.inla = list(strategy = "gaussian")
    )
    if (length(loaFit$parameters)) {
      glgm_sum <- loaFit$parameters$summary
      # append glgm columns when mapping is available
    }
  }

  plot_loaloa_maps <- function() {
    if (!exists("loaFit") || !length(loaFit$parameters)) {
      my_sim <- matern_est(fit_loaloa, loaloa_for_fit$grid, n = 0L)
      plot(my_sim[["mean"]], main = "adlaplace (hierarchical FEM)")
      return(invisible(NULL))
    }
    my_sim <- matern_est(fit_loaloa, loaloa_for_fit$grid, n = 0L)
    op <- par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
    on.exit(par(op), add = TRUE)
    plot(loaFit$raster[["random.mean"]], main = "geostatsp::glgm")
    plot(my_sim[["mean"]], main = "adlaplace (hierarchical FEM)")
    invisible(NULL)
  }

  plot_loaloa_knots <- function() {
    op <- par(mar = c(3, 3, 2, 1))
    on.exit(par(op), add = TRUE)
    plot(coarse, col = c(NA, "grey90"), legend = FALSE, main = NULL)
    terra::plot(west, add = TRUE, col = c(NA, "coral"), legend = FALSE)
    terra::plot(east, add = TRUE, col = c(NA, "coral"), legend = FALSE)
    points(terra::crds(loaloa), pch = 16, cex = 0.5)
    legend("topright",
      legend = c("coarse", "clusters"),
      fill = c("grey90", "coral"),
      bty = "n"
    )
    invisible(NULL)
  }

  loaloa_ok <- TRUE
}

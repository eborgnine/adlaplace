## ----setup, include=FALSE---------------------------------------------------------------------------------------------------------------------------
library("Matrix")
library("adlaplace")


## ----gamm-------------------------------------------------------------------------------------------------------------------------------------------
if (requireNamespace("mgcv", quietly = TRUE)) {
  dat <- mgcv::gamSim(6, n = 200, scale = .2, dist = "poisson")
  b2 <- mgcv::gamm(y ~ s(x2),
    family = poisson,
    data = dat, random = list(fac = ~1)
  )
  plot(b2$gam, pages = 1)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------
if (requireNamespace("mgcv", quietly = TRUE)) {
  md <- adlaplace::model_data(
    data = dat,
    formula = y ~ x1 +
      adlaplace::iwp(x2, p = 2, knots = seq(0, 1, len = 11)) +
      adlaplace::iid(fac) +
      adlaplace::overdispersion(lower = 1e-9)
  )

  config <- list(
    beta = md$info$beta$init,
    gamma = rep(0, nrow(md$info$gamma)),
    theta = adlaplace::apply_theta_log(md$info$theta, cols = "init")$init,
    transform_theta = TRUE,
    shards = adlaplace::ad_shards(
      Matrix::t(md$shards$parent@ATp), num_shards = 100
    ),
    num_threads = 2L,
    verbose = TRUE,
    package = "adlaplace"
  )


  ad_fun <- adlaplace::ad_fun(md$shards$parent, config)
  res <- adlaplace::log_lik_laplace(
    x = c(config$beta, config$theta),
    ad_fun = ad_fun,
    config = c(config, list(verbose = FALSE)),
    deriv = TRUE
  )
  res$parameters
  res$fval
  res$grad
}


## ----trustOptimOuterWrappers------------------------------------------------------------------------------------------------------------------------
if (requireNamespace("mgcv", quietly = TRUE)) {
  theta_for_opt <- adlaplace::apply_theta_log(
    md$info$theta,
    cols = c("init", "lower", "upper")
  )

  to_keep <- c("init", "lower", "upper", "parscale")

  config$opt <- as.list(
    rbind(
      md$info$beta[, to_keep],
      theta_for_opt[, to_keep]
    )
  )

  cache <- new.env(parent = emptyenv())
  cache$gamma <- rep(0, nrow(md$info$gamma))

  x0 <- config$opt$init

  adlaplace::outer_fn(x = x0, cache = cache, config = config, ad_fun = ad_fun)
  adlaplace::outer_gr(x = x0, cache = cache, config = config, ad_fun = ad_fun)

  outer_fit <- stats::optim(
    par = config$opt$init,
    fn = adlaplace::outer_fn,
    gr = adlaplace::outer_gr,
    method = "L-BFGS-B",
    control = list(
      maxit = 1000,
      trace = 3,
      REPORT = 1
    ),
    lower = config$opt$lower,
    upper = config$opt$upper,
    config = config,
    ad_fun = ad_fun,
    cache = cache
  )

  outer_fit$par
  outer_fit$value
}


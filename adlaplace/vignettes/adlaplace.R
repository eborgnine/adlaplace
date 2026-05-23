## ----theData----------------------------------------------------------------------------------------------------------------------------------------
Nobs <- 1000
Nrandom1 <- 10
Nrandom2 <- 25

set.seed(0)
X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
AmatList <- list(
  Matrix::sparseMatrix(
    i = 1:Nobs,
    j = sample(Nrandom1, Nobs, replace = TRUE)
  ),
  Matrix::sparseMatrix(
    i = 1:Nobs,
    j = sample(Nrandom2, Nobs, replace = TRUE)
  )
)
Amat <- do.call(cbind, AmatList)
beta <- rep(1, ncol(X))
thetaOrig <- c(0.1, 0.1, 0.1)
theta <- log(thetaOrig)
gamma <- rnorm(ncol(Amat), sd = rep(thetaOrig[1:2], c(Nrandom1, Nrandom2)))
eta <- as.vector(X %*% beta + Amat %*% gamma)
etaT <- exp(eta)
Z <- rgamma(length(etaT), 0.1^(-2), 0.1^(-2))
c(mean(Z), sd(Z))
y <- rpois(length(etaT), etaT * Z)

NperEffect <- unlist(lapply(AmatList, ncol))
# map = rep(seq(0, len=length(NperEffect)), NperEffect)
map <- Matrix::sparseMatrix(
  i = seq(0L, len = ncol(Amat)),
  j = rep(seq(0L, len = length(NperEffect)), NperEffect),
  x = 1L,
  index1 = FALSE,
  dims = c(ncol(Amat), length(thetaOrig))
)

data <- list(
  y = y,
  ATp = as(Matrix::t(Amat), "dMatrix"),
  XTp = as(Matrix::t(X), "CsparseMatrix"),
  map = map,
  Qdiag = rep(1, ncol(Amat))
)


config <- list(
  beta = rep(1, nrow(data$XTp)),
  theta = c(rep(-1, length(NperEffect)), -1),
  transform_theta = TRUE,
  gamma = rep(0, nrow(data$ATp)),
  shards = adlaplace::adFun_groups(data$ATp, Ngroups = 1000),
  num_threads = 2L,
  verbose = TRUE, package = "adlaplace"
)


## ----testLogLik-------------------------------------------------------------------------------------------------------------------------------------
model <- adlaplace:::ad_model_from_config_matrices(
  y = data$y,
  ATp = data$ATp,
  XTp = data$XTp,
  config = config,
  theta_local_row = length(config$theta) - 1L
)
ad_fun <- adlaplace::ad_fun(model, config)

res <- adlaplace::log_lik_laplace(
  x = c(config$beta, config$theta),
  ad_fun = ad_fun,
  config = modifyList(config, list(verbose = FALSE)),
  deriv = TRUE
)
res$full_parameters[1:5]
res$neg_log_lik
res$opt$fval
res$grad


## ----trustOptimOuterWrappers------------------------------------------------------------------------------------------------------------------------
model <- adlaplace:::ad_model_from_config_matrices(
  y = data$y,
  ATp = data$ATp,
  XTp = data$XTp,
  config = config,
  theta_local_row = length(config$theta) - 1L
)
ad_fun <- adlaplace::ad_fun(model, config)

cache <- new.env(parent = emptyenv())
cache$gamma <- config$gamma

x0 <- c(config$beta, config$theta)

adlaplace::outer_fn(x=x0, cache=cache, config=config, ad_fun = ad_fun)
adlaplace::outer_gr(x=x0, cache=cache, config=config, ad_fun = ad_fun)

outer_fit <- stats::optim(
  par = x0,
  fn = adlaplace::outer_fn,
  gr = adlaplace::outer_gr,
  method = "L-BFGS-B",
  lower = rep(-5, length(x0)),
  upper = rep(5, length(x0)),
  control = list(
    maxit = 1000,
    trace = 3,
    REPORT = 1
  ),
  config = config,
  ad_fun = ad_fun,
  cache = cache,
  control_inner = list(
    maxit = 100,
    report.level = 0,
    report.freq = 0
  )
)

outer_fit$par
outer_fit$value


## ----testLogLgrad-----------------------------------------------------------------------------------------------------------------------------------
model <- adlaplace:::ad_model_from_config_matrices(
  y = data$y,
  ATp = data$ATp,
  XTp = data$XTp,
  config = config,
  theta_local_row = length(config$theta) - 1L
)
ad_fun <- adlaplace::ad_fun(model, config)
x <- c(config$beta, config$theta)
Npar <- 13
Dpar <- length(x)

parDf <- data.frame(x)[, rep(1, Npar)]
Sx <- seq(-0.1, 0.1, len = 13) + parDf[Dpar, 1]
SxD <- Sx[-1] - diff(Sx) / 2
parDf[Dpar, ] <- Sx

res <- mapply(
  adlaplace::log_lik_laplace,
  x = as.list(parDf),
  MoreArgs = list(ad_fun = ad_fun, config = config, deriv = TRUE),
  SIMPLIFY = FALSE
)
Slik <- unlist(lapply(res, "[[", "log_lik"))
SnegLik <- unlist(lapply(res, "[[", "neg_log_lik"))
Sdet <- unlist(lapply(res, function(xx) xx$hessian$half_log_det))
dU <- do.call(
  abind::abind,
  c(lapply(res, function(xx) as.matrix(xx$extra$dU)), along = 3)
)
uHat <- do.call(rbind, lapply(res, function(xx) xx$opt$solution))
extraDf <- do.call(abind::abind, c(lapply(res, "[[", "deriv"), along = 3))
gradMat <- do.call(rbind, lapply(res, "[[", "grad"))
gradVec <- unlist(lapply(res, function(xx) xx$deriv[Dpar, "d_log_lik"]))


plot(Sx, SnegLik)

plot(Sx, gradMat[, Dpar], type = "l")
points(SxD, diff(SnegLik) / diff(Sx))

Du <- 1
plot(Sx, uHat[, Du])

plot(Sx, dU[Du, Dpar, ])
lines(SxD, diff(uHat[, Du]) / diff(Sx))

plot(Sx, Sdet)

plot(Sx, extraDf[Dpar, "d_det", ])
lines(SxD, diff(Sdet) / diff(Sx))


## ----testDeriv--------------------------------------------------------------------------------------------------------------------------------------
model <- adlaplace:::ad_model_from_config_matrices(
  y = data$y,
  ATp = data$ATp,
  XTp = data$XTp,
  config = config,
  theta_local_row = length(config$theta) - 1L
)
ad_fun <- adlaplace::ad_fun(model, modifyList(config, list(verbose = TRUE)))

x = c(config$beta, rep(0.1, length(config$gamma)), config$theta)
adlaplace::joint_log_dens(ad_fun, x, negative = FALSE)


str(adlaplace::grad(ad_fun, x, inner = FALSE, negative = FALSE))
str(adlaplace::grad(ad_fun, x, inner = TRUE, negative = FALSE))


Dgroup = ncol(config$shards)+1
str(h1 <- adlaplace::hessian(ad_fun, x, inner = FALSE, shards = Dgroup, negative = FALSE))
ad_fun@group_sparsity[[1L + Dgroup]]
(Shere = 1 + seq(
  ad_fun@map_outer[[Dgroup + 1L]]$p[Dgroup + 1L],
  ad_fun@map_outer[[Dgroup + 1L]]$p[Dgroup + 2L] - 1L
))
ad_fun@map_outer[[Dgroup + 1L]]$local[Shere]
ad_fun@map_outer[[Dgroup + 1L]]$global[Shere]
(hseq = which(h1@x != 0))
h1@x[hseq]


str(h2 <- adlaplace::hessian(ad_fun, x, inner = TRUE, negative = FALSE))


## ----trustOptimInterface, eval=FALSE----------------------------------------------------------------------------------------------------------------
# model <- adlaplace:::ad_model_from_config_matrices(
#   y = data$y,
#   ATp = data$ATp,
#   XTp = data$XTp,
#   config = config,
#   theta_local_row = length(config$theta) - 1L
# )
# ad_fun <- adlaplace::ad_fun(model, config)
# inner_res <- adlaplace::inner_opt(
#   parameters = c(config$beta, config$theta),
#   gamma = config$gamma,
#   control = list(
#     maxit = 50,
#     report.level = 0,
#     report.freq = 0
#   ),
#   config = config,
#   ad_fun = ad_fun
# )
# 
# str(inner_res)
# 
# quantile(inner_res$gradient)


## ----derivJointDens, eval=TRUE----------------------------------------------------------------------------------------------------------------------
config$gamma <- rep(1, length(config$gamma))
x <- c(config$beta, config$gamma, config$theta)
model <- adlaplace:::ad_model_from_config_matrices(
  y = data$y,
  ATp = data$ATp,
  XTp = data$XTp,
  config = config,
  theta_local_row = length(config$theta) - 1L
)
ad_fun <- adlaplace::ad_fun(model, modifyList(config, list(verbose = TRUE)))

inner <- FALSE
type <- c("outer", "inner")[1 + inner]

Sgamma1 <- seq.int(length(config$beta) + 1, length.out = length(config$gamma))

Npar <- 13
Dpar <- 3
if (inner) {
  Dpar <- match(Dpar, Sgamma1)
}
Dpar0 <- Dpar - 1L

h_outer <- adlaplace::hessian(ad_fun, x, inner = FALSE, negative = FALSE)
bob <- as(h_outer, "TsparseMatrix")
(whichIndex <- which(bob@i == Dpar0 & bob@j == Dpar0) - 1L)
map_outer <- ad_fun@map_outer[[1L]]
(whichGlobal <- which(map_outer$global == whichIndex) - 1L)

(whichP <- mapply(
  function(xx) min(which(map_outer$p >= xx)),
  xx = whichGlobal
) - 1L)

do.call(rbind, ad_fun@group_sparsity[[whichP[1]]][c("row_hess", "col_hess")])
do.call(rbind, ad_fun@group_sparsity[[whichP[1]]][c("row_hess_inner", "col_hess_inner")])

(shards <- whichP - 1L)
shards <- seq.int(from = 0, length.out = adlaplace::n_groups(ad_fun@ptr))

parDf <- data.frame(x)[, rep(1, Npar)]
Sx <- seq(-0.1, 0.1, len = 13) + parDf[Dpar, 1]
SxD <- Sx[-1] - diff(Sx) / 2
parDf[Dpar, ] <- Sx


dens1 <- mapply(
  adlaplace::joint_log_dens,
  x = as.list(parDf),
  MoreArgs = list(ad_fun = ad_fun, shards = shards, negative = FALSE),
  SIMPLIFY = TRUE
)

plot(Sx, dens1)


grad1 <- mapply(
  adlaplace::grad,
  x = as.list(parDf),
  MoreArgs = list(ad_fun = ad_fun, inner = inner, shards = shards),
  SIMPLIFY = TRUE
)
plot(Sx, grad1[Dpar, ])
lines(SxD, diff(dens1) / diff(Sx))


hes1 <- mapply(
  adlaplace::hessian,
  x = as.list(parDf),
  MoreArgs = list(ad_fun = ad_fun, inner = inner, shards = shards),
  SIMPLIFY = FALSE
)

hes2 <- do.call(abind::abind, c(
  lapply(hes1, as.matrix),
  along = 3
))

gradD <- apply(grad1, 1, diff) / mean(diff(Sx))

Dpar2 <- dim(hes2)[1] - 2
plot(Sx, hes2[Dpar, Dpar2, ],
  ylab = paste(Dpar, Dpar2),
  ylim = range(hes2[Dpar, Dpar2, ], na.rm = TRUE)
)
lines(SxD, gradD[, Dpar2])


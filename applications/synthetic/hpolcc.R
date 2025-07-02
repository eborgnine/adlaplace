
xx = readRDS("simData.rds")
data = xx$data 

betaHat  = xx$par[grep("beta", names(xx$par))]
logPrecHat = xx$par[grep("theta", names(xx$par))]
sigmaHat = exp(logPrecHat)^(-1/2)

gamma_start = xx$par[grep("gamma", names(xx$par))]


  knots_pm <- seq(0, max(ceiling(max(data$pm)/5)*5, 20), 2.5)

  formula <- as.formula(
    sprintf("count ~ f(date, model = 'iid') + hum + f(pm, model = 'hiwp', p = 2, ref_value = 10, knots = c(%s), group_var = region)",
            paste(knots_pm, collapse = ", "))
  ) 

library(hpolcc)

  cc_design <- ccDesign(time_var = "date", strat_vars = "region")


fit <- hnlm(formula, data, cc_design = cc_design, 
 dirichelet= FALSE, for_dev=TRUE)



fit$config$num_threads = 8
fit$config$transform_theta = FALSE


	# don't remaximize
res1 = hpolcc:::objectiveFunctionC( 
 c(betaHat, gamma_start, sigmaHat),
 fit$tmb_data, 
 fit$config
  )
res1$h2 = Matrix::forceSymmetric(res1$hessian)
res1$hessian[1:5,1:5]
xx$h[1:5,1:5]

allPar = c(betaHat, gamma_start, sigmaHat)
bob = function(x) {
	par = c(x,  allPar[-(1:length(x))])
	configHere = fit$config
	configHere$maxDeriv = 0
	hpolcc:::objectiveFunctionC( 
 		c(betaHat, gamma_start, sigmaHat),
 		fit$tmb_data, configHere
  )$value
}
bob(allPar)

res1$cholHessian = Matrix::chol(res1$h2)

# remaximize
res = loglik(c(betaHat, sigmaHat), 
             gamma_start, 
             data=fit$tmb_data, 
    config = fit$config,
    control = list(maxit=300, start.trust.radius = 100,
      report.level=4, report.freq=1))



cholVariance = Matrix::solve(res$cholHessian)


Nsim = 200
gammaSim =  res$gamma_hat + 
  cholVariance %*% matrix(
  rnorm(Nsim*nrow(cholVariance)),
  nrow(cholVariance), Nsim)


# only plot the global effects
rowsToKeep = grep("pm.*global", rownames(fitd$tmb_data$ATp))

Stemp = seq(0,max(data$pm),by=2.5)
SobsTemp = na.omit(match(Stemp, round(fit$data$tempClean, 1)))

par(mar=c(3,3,0,0), bty='n')
matplot(
  Stemp,
  crossprod(
    fit$tmb_data$XTp[, SobsTemp],
    parameters[1:Nbeta]
    ) + 
  crossprod(
    fit$tmb_data$ATp[rowsToKeep, SobsTemp], 
    gammaSim[rowsToKeep,]),
  type='l', col='#00000020', lty=1, xaxs='i',
  yaxs='i'
)




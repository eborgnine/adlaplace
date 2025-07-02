
data = readRDS("simData.rds")
betaHat  = c(1.0858353, 0.2385639)
logPrecHat = c( 5.6024098, 4.0603066, 5.5502284, 6.9152824)
sigmaHat = exp(logPrecHat)^(-1/2)

  knots_pm <- seq(0, max(ceiling(max(data$pm)/5)*5, 20), 2.5)

  formula <- as.formula(
    sprintf("count ~ f(date, model = 'iid') + hum + f(pm, model = 'hiwp', p = 2, ref_value = 10, knots = c(%s), group_var = region)",
            paste(knots_pm, collapse = ", "))
  ) 


  cc_design <- ccDesign(time_var = "date", strat_vars = "region")

library(hpolcc)

fit <- hnlm(formula, data, cc_design = cc_design, 
 dirichelet= FALSE, for_dev=TRUE)


gamma_start = rep(0, nrow(fit$tmb_data$ATp))

res = loglik(MLE, 
             gamma_start, 
             data=fit$tmb_data, 
    config = fit$config,
    control = list(maxit=300, start.trust.radius = 100,
      report.level=4, report.freq=1))



cholVariance = Matrix::solve(res$cholHessian)
quantile(sqrt(Matrix::diag(cholVariance)))

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


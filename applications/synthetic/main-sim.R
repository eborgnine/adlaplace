# Rscript applications/synthetic/main-sim.R   > applications/synthetic/log.txt 2>&1


# Packages ----------------------------------------------------------------
library(hpolcc)
library(data.table)



# Data --------------------------------------------------------------------
genPM <- function(date){
  l <- length(date)+1
  date_int <- as.integer(date)
  date_int <- c(date_int[1]-1, date_int)
  pm <- 10
  pm <- pm + 10*sin(2*pi * (date_int + rnorm(l,0,1))/365.25)
  pm <- pm + rnorm(l,0,1)
  pm <- (2*pm[-l] + pm[-1])/3
  pmax(0,pm)
}

genHum <- function(date){
  l <- length(date)+1
  date_int <- as.integer(date)
  date_int <- c(date_int[1]-1, date_int)
  hum <- 50
  hum <- hum + 5*cos(2*pi * (date_int + rnorm(l,0,1))/365.25)
  hum <- hum + rnorm(l,0,1)
  hum <- (2*hum[-l] + hum[-1])/3
  hum <- pmax(0,hum)
  hum/max(hum)
}

genPmEffect <- function(pm, r1, r2, coseffect=1){
  0.1*(r1*pm/10 + 0.2*pm^(r2/2) - coseffect*cos(pm*pi/25))
}
genCount <- function(hum, pm, r1=1, r2=1, od = FALSE){
  l <- length(pm)
#  od <- rnorm(l, 0, ifelse(od, .1, 0))
  odSd = 0.3;odPrec = odSd^(-2)
  od <- rgamma(l, shape = odPrec, rate = odPrec )
#  print(sd(od)/odSd)
  rpois(l, exp(-log(1) + hum + genPmEffect(pm, r1, r2) + log(od  )))
}


# sim loop ----------------------------------------------------------------

#res <- mclapply(1:500, \(dummy){
set.seed(0)
region_effect1 <- rnorm(10,1,.2)
region_effect2 <- rnorm(10,1,.2)
  data <- lapply(1:10, \(i){
    data <- data.table(date = as.Date(1:1000), region = i)
    data$hum <- genHum(data$date)
    data$pm <- genPM(data$date)
    data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i], TRUE) 
    data
  }) 
  data = do.call(rbind, data)

  data$monthDow = format(data$date, '%Y-%m-%a')

  dataOrig = data

  cc_design <- ccDesign(time_var = "date", strat_vars = c("region",'monthDow'))

  ref_values <- list("pm" = 10)
  knots_pm <- seq(0, 22, by=2)

  res  <- hnlm(
    formula = count ~  hum + 
    f(pm, model = 'hiwp', p = 2, ref_value = 10, 
      knots = knots_pm, group_var = region), 
    dataOrig, 
    cc_design = cc_design, 
    for_dev=TRUE, 
    verbose=FALSE, 
    dirichlet=TRUE,
    control_inner=list(maxit=1000, start.trust.radius = 1, prec=1e-7, stop.trust.radius = 1e-9,
        cg.tol = 1e-6, report.level=0),
    control=list(maxit=1000, start.trust.radius = 1, prec=1e-6, stop.trust.radius = 1e-9,
        cg.tol = 1e-6, report.level=4, report.freq=1, report.header.freq=10, report.precision=7),
    config = list(num_threads = parallel::detectCores(), strataPerIter=100, transform_theta = TRUE)
  )

    image(seq(0, len=nrow(res$groups$firstDeriv)), seq(0, len=ncol(res$groups$firstDeriv)), 
        as.matrix(res$groups$firstDeriv[, 1+res$groups$groups$i])
    )
    abline(h=res$groups$groups$p, col='blue', lty=1)


library(hpolcc)


adfun = res$adfun()
grad2 = grad(res$start_gamma, data=res$tmb_data, config=configHere, adfun=adfun)
hes2 = hessian(res$start_gamma, data=res$tmb_data, config=configHere, adfun=adfun)

innerRes = trustOptim::trust.optim(
    x = res$start_gamma, 
    fn = jointLogDens,
    gr = grad,
    hs = hessian,
    method = 'Sparse',
    control = list(start.trust.radius = 0.1, report.freq=1, report.level=10),
    data=res$tmb_data, config = res$config, adfun = res$adfun()
  )


innerResB <- BB::spg(par = innerRes$solution, 
    fn = jointLogDens,
    gr = grad,
  data=res$tmb_data, config = configHere)  

res2  = loglik(
  res$parameters,
  res$start_gamma,
  res$tmb_data, config=res$config,
  deriv=0
)


str(res$groups$group_sparsity[[1]]$second)

plot(as.matrix(hesSF), hesDF)

hesD[1:8,1:8]
hes[1:8,1:8]

plot(as.matrix(hes) ,hesD);abline(0,1)


hesN[1:8,1:8]



hes[1:5,1:5]
hesN[Sgamma1,Sgamma1][1:5,1:5]
hesD[Sgamma1,Sgamma1][1:5,1:5]


res$config$sparsity = grouped$sparsity
res$config$dense=FALSE
res$config$maxDeriv = 3
dQ = logLikOnlyQ(res$parameters_for_sparsity, res$tmb_data, res$config)


res$config$maxDeriv = 2
res$config$dense=FALSE
res$config$num_threads = 1
dL = logLikNoQStrata(res$parameters_for_sparsity, res$tmb_data, res$config,
    #grouped$groups, 
    c(grouped$groups[c('i','j')], list(p=grouped$groups$p) ), 
    grouped$group_sparsity)
res$config$dense=TRUE
dL2 = logLikNoQStrata(res$parameters_for_sparsity, res$tmb_data, res$config,
    #grouped$groups, 
    c(grouped$groups[c('i','j')], list(p=grouped$groups$p) ), 
    grouped$group_sparsity)


dL$hessianMat = Matrix::sparseMatrix(
  i=res$config$sparsity$second$full$i, 
  p=res$config$sparsity$second$full$p, 
  x=drop(dL$hessian), index1=0,
  symmetric=TRUE)

quantile(as.matrix(dL$hessianMat) - dL2$hessian)


bob = numDeriv::hessian(
  function(xx)  logLikNoQStrata(
    c(xx, res$parameters_for_sparsity[-(1:length(xx))]), 
    res$tmb_data, c(res$config[setdiff(names(res$config), 'maxDeriv')], list(maxDeriv=0)),
    grouped$groups, grouped$group_sparsity
  )$value, 
  res$parameters_for_sparsity[1:5])

bob
dL$hessianMat[1:5,1:5]
dL2$hessian[1:5,1:5]


#' To Do:
#' - combine third diag and offdiag.
#' - sum over groups in Third, add Q to third
#' - rename likelihood funcitons
#' - wrappers
#' - loglik
#' create ADfuns separately, pass instead of data.  createADfun(x, data,  config, groups)  and list of patterns and list of work
#' functions grad, hessiansparse, hessiandense, thirdsparse, thirddense which take x, adfunlist, config, sparsity as arguments
#' R interface 1) pass data, creates adfuns  2) pass external pointer stuff
#' save ADfun as EXTPTRSXP (and pattern, srowout, work?)



config2 = res$config
config2$verbose=FALSE
config2$dense=TRUE
config2$num_threads = 10
bob = hpolcc::thirdOffDiagonals(xx, res$tmb_data, config2)
bobA = thirdDeriv(xx, res$tmb_data, config2)



cache = new.env()
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)
assign("gamma_start", res$start_gamma, cache)
assign("file", "simbfgs.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))


mleX =  trustOptim::trust.optim(
    x = res$parameters, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
    method = 'SR1',
    control = list(start.trust.radius = 0.1, report.freq=1, report.level=10),
    data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
  )
mleX$solution

mleB <- BB::spg(par = mleX$solution, 
  fn = wrappers_outer$fn,
  gr = wrappers_outer$gr,
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner,
  control = list(maxit = 1e2, ftol=1e-12, gtol = 1e-9, M = 15, trace=TRUE, checkGrad=TRUE))  


mleN <- optimx::optimx(
  par = mleX$solution, 
  fn = wrappers_outer$fn,
  gr = wrappers_outer$gr,
  method = 'Rvmmin',
  control = list(trace = 10, maxit=1e3),
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
)

xx=c(mleB$value,mleX$fval,mleX3$fval, mleN$value);xx - min(xx)

mle = loglik(
  parameters=mleX$solution,
  gamma_start=get("gamma_start", cache),
  data=res$tmb_data, config=res$config, control=res$control_inner,
  check=TRUE
)


  Sigma <- mle$extra$invHessianRandom
  mu = mle$solution
  samp <- mvtnorm::rmvnorm(1e4, mu, as.matrix(Sigma))

  xseq = seq(knots_pm[1],rev(knots_pm)[1],by=0.5)
  df = expand.grid(
    pm = xseq,
    region =  c(-1, unique(res$data$region)), 
    hum=0)
  XA = hpolcc:::getNewXA(res$terms, df)
  df$fixed = as.vector(XA$X %*% mle$parameters[1:nrow(res$tmb_data$XTp)])

  samp2 = tcrossprod(XA$A, samp)
  Sprob = c(0.005, 0.025, 0.1, 0.5, 0.9, 0.975, 0.995)
  sampq = t(apply(samp2, 1, quantile, prob = Sprob))
  colnames(sampq) = paste0('q',Sprob)
  sampq = sampq + df$fixed


  Slwd = 4*sqrt(pmin(Sprob, 1-Sprob))
  Scol = c('grey','black')[1+(Sprob == 0.5)]
  Sorder = order(Slwd)

  toPlot = which(df$region == -1)
  matplot(df[toPlot, 'pm'], 
    sampq[toPlot,Sorder], type='l', 
  lty=1, col=Scol[Sorder], lwd=Slwd[Sorder], xaxs='i', yaxs='i')
  lines(xseq, genPmEffect(xseq, 1, 1) - genPmEffect(ref_values$pm, 1, 1), col='yellow')

  toPlot = which(df$region == 1)
  matplot(df[toPlot, 'pm'], 
    sampq[toPlot,Sorder], type='l', 
  lty=1, col=Scol[Sorder], lwd=Slwd[Sorder], xaxs='i', yaxs='i')





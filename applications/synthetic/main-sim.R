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
        cg.tol = 1e-6, report.level=4, report.freq=10, report.header.freq=10, report.precision=7),
    config = list(num_threads = 1, #parallel::detectCores(), 
      Nclusters = 1, transform_theta = TRUE)
  )
#save(res, file='.RData')
    image(seq(0, len=nrow(res$groups$firstDeriv)), seq(0, len=ncol(res$groups$firstDeriv)), 
        as.matrix(res$groups$firstDeriv[, 1+res$groups$groups$i])
    )
    abline(h=res$groups$groups$p, col='blue', lty=1)


res$tmb_data$Qdiag = rep(1, 0)
library(hpolcc)
adFun = getAdFun(res$gamma_start, res$tmb_data, res$config)
adFunFull = getAdFun(res$parameters_and_gamma, res$tmb_data, res$config)

res$config$dense = TRUE
bob =  thirdStrata(
    res$parameters_and_gamma, data=res$tmb_data, config=res$config, adFunFull
  ) 
bob$diagMat = matrix(bob$diag, length(res$parameters_and_gamma))
bob$Tmat = matrix(bob$Tijk[[1]], ncol=nrow(res$config$group_sparsity[[1]]$third$pairs))
bob$H = matrix(bob$second, length(res$parameters_and_gamma))

rawT = res$config$group_sparsity[[1]]$third$pairs[
  rep(
    1:nrow(res$config$group_sparsity[[1]]$third$pairs), 
    each = 137), c('i','j')]
rawT$k = rep(seq(0,len=137), nrow(res$config$group_sparsity[[1]]$third$pairs)) 
rawT$x = round(bob$Tijk[[1]],4)


Dpar = 2
Sx = seq(-0.001, 0.001, len=6) + res$parameters[Dpar]
Sx1 = Sx[-1] - diff(Sx)/2
parMat = as.data.frame(matrix(res$parameters_and_gamma, ncol=length(Sx), nrow=length(res$parameters_and_gamma), byrow=FALSE))
parMat[Dpar,] = Sx

theH = mapply(
  hessian,
  parameters = parMat,
  MoreArgs = list(
    data = res$tmb_data, config=res$config, adFun = adFunFull
  )
)
theH2 = lapply(theH, Matrix::as.matrix)
theH3 = do.call(abind::abind, c(theH2, list(along=3)))
theT = apply(theH3, 1:2, diff)/mean(diff(Sx))


round(diag(theT[3,,])[1:10],4)
round(bob$diagMat[Dpar, 1:10], 4)*2



Dpar1 = 3;Dpar2=8


Tiik = bob$diagMat

round(xx <- c(theT[3,Dpar1, Dpar2],Tiik[Dpar, Dpar1], Tiik[Dpar, Dpar2], Tiik[Dpar1, Dpar], Tiik[Dpar2, Dpar], 
  Tiik[Dpar1, Dpar2], Tiik[Dpar2, Dpar1]), 4)

 -1888.4776*c(0.5,1,2,4)

sum(xx[1:3]*c(-1,1,1))
sum(xx[c(1,4,5)]*c(-1,1,1))
sum(xx[1:3]*c(1,0.5, 0.5))
sum(xx[c(1,4,5)]*c(1,0.5, 0.5))
sum(xx[1:3]*c(-1,0.5, 0.5))
sum(xx[c(1,4,5)]*c(-1,0.5, 0.5))



round(theT[3,1:8,1:10],4)
rawT[136+1:9,]
round(bob$diagMat[1:8, 1:10], 4)

cbind(res$config$group_sparsity[[1]]$third$ijk, round(bob$raw$Tijk[[1]], 5))[1:5]

bob$third[apply(bob$third[,1:3], 1,function(xx) all((xx+1) %in% c(Dpar, Dpar1, Dpar2))),]

bob$raw$diagMat[1:8,1:8] # need to double!
bob$raw$Tmat[1:8,1:8] # need to double!



res$config$group_sparsity[[1]]$third$pairs[1:2,]
cbind(i=0, j=7, k=seq(0, len=length(res$parameters_and_gamma)),rawMat[,1])[1:5,]


plot(
  Sx, theH3[Dpar1, Dpar2,],
)


theT[3, 3:8,3:8]
bob$thirdList[[Dpar]][3:8,3:8]
bob$raw$diagMat = matrix(bob$raw$diag, nrow(parMat), nrow(parMat))
bob$raw$diagMat[1:6,1:6]


sum(bob$raw$diagMat[c(Dpar1, Dpar2),Dpar])/2
sum(bob$raw$diagMat[Dpar, c(Dpar1, Dpar2)])/2



plot(Sx1, theT[,Dpar1, Dpar2])
abline(h=bob$raw$diagMat[Dpar, Dpar1])# - bob$fullHessian[Dpar, Dpar1])
abline(v=res$parameters_and_gamma[Dpar])


abline(h=bob$thirdList[[Dpar]][Dpar1, Dpar2])
abline(h=bob$raw$diagMat[Dpar, Dpar1])

res$config$dense = TRUE
res2 = loglik(res$parameters, gamma_start = res$gamma_start, data = res$tmb_data, config=res$config, 
  adFun = adFun, adFunFull = adFunFull, control = res$control_inner)


wrappers_outer$fn(
  res$parameters, data = res$tmb_data, config=res$config, 
  adFun = adFun, adFunFull = adFunFull, control = res$control_inner,
  cache = res$cache
)

wrappers_outer$gr(
  res$parameters, data = res$tmb_data, config=res$config, 
  adFun = adFun, adFunFull = adFunFull, control = res$control_inner,
  cache = res$cache
)




mle <- trustOptim::trust.optim(
    x = res$parameters,
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
    method = "SR1",
    control = res$control,
    data=res$tmb_data, config = res$config, 
    control_inner = res$control_inner,
    adFun=adFun, adFunFull = adFunFull, 
    cache = res$cache
  )

mleB <- BB::spg(
  par = res$parameters, #mle$solution, 
  fn = wrappers_outer$fn,
  gr = wrappers_outer$gr,
  data=res$tmb_data, config = res$config, 
  control_inner = res$control_inner,
  adFun=adFun, adFunFull = adFunFull, 
  cache = res$cache,
  control = list(maxit = 1e2, ftol=1e-12, gtol = 1e-9, M = 15, trace=TRUE, checkGrad=TRUE))  




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





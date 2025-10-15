# Rscript applications/synthetic/main-sim.R   > applications/synthetic/log.txt 2>&1


# Packages ----------------------------------------------------------------
library(hpolcc)
library(data.table)


constructEffect <- function (fit, exposure_var, group_var, group, values, ref_values, pars = NULL, probs = c(.1,.5,.9)){
  vars <- c(as.character(fit$formula)[2], unique(sapply(fit$terms, 
                                                        "[[", "var")), unique(unlist(sapply(fit$terms, "[[", 
                                                                                            "group_var"))), fit$cc_design$time_var, fit$cc_design$strat_var)
  group_var[is.na(group_var)] <- "__NONE__"
  df <- data.frame(row.names = 1:(length(values) * length(group)))
  for (v in vars) df[[v]] <- ifelse(is.null(ref_values[[v]]), 
                                    0, ref_values[[v]])
  df[[group_var]] <- unlist(rep(group, each = length(values)))
  df[[exposure_var]] <- rep(values, times = max(1, length(group)))
  list2env(hpoltest:::getNewXA(fit$terms, df), envir = environment())
  if(is.null(pars)) pars <- fit$obj$env$last.par.best
  if(!is.matrix(pars)) pars <- as.matrix(pars)
  
  beta <- pars[,colnames(pars) == "beta",drop=F]
  gamma <- pars[,colnames(pars) == "gamma",drop=F]
  
  P <- X %*% t(beta) + A %*% t(gamma)
  Q <- t(apply(P,1,quantile,probs=probs))
  colnames(Q) <- paste0("q_",probs)
  
  df <- data.frame(variable = exposure_var, var_value = values, group = df[[group_var]])
  cbind(df, Q)
}



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

genPmEffect <- function(pm, r1, r2, coseffect=5){
  0.5*(r1*pm/10 + pm^(r2/2) - coseffect*cos(pm*pi/25))
}
genCount <- function(hum, pm, r1=1, r2=1, od = FALSE){
  l <- length(pm)
#  od <- rnorm(l, 0, ifelse(od, .1, 0))
  odSd = 0.3;odPrec = odSd^(-2)
  od <- rgamma(l, shape = odPrec, rate = odPrec )
#  print(sd(od)/odSd)
  rpois(l, exp(-log(0.1) + hum + genPmEffect(pm, r1, r2) + log(od  )))
}


# sim loop ----------------------------------------------------------------

#res <- mclapply(1:500, \(dummy){
set.seed(0)
region_effect1 <- rnorm(10,1,.1)
region_effect2 <- rnorm(10,1,.1)
  data <- lapply(1:10, \(i){
    data <- data.table(date = as.Date(1:1000), region = i)
    data$hum <- genHum(data$date)
    data$pm <- genPM(data$date)
    data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i], TRUE) 
    data
  }) 
  data = do.call(rbind, data)

  data$monthDow = format(data$date, '%Y-%m-%a')

  cc_design <- ccDesign(time_var = "date", strat_vars = c('monthDow',"region"))

  ref_values <- list("pm" = 10)
  knots_pm <- seq(0, 22, by=2)

  res  <- hnlm(count ~  hum + 
    f(pm, model = 'hiwp', p = 2, ref_value = 10, 
      knots = knots_pm, group_var = region), 
    data, 
    cc_design = cc_design, 
    for_dev=TRUE, 
    verbose=FALSE, 
    dirichlet=TRUE,
    control_inner=list(maxit=1000, start.trust.radius = 1, prec=1e-7, stop.trust.radius = 1e-9,
        cg.tol = 1e-6, report.level=0),
    control=list(maxit=1000, start.trust.radius = 1, prec=1e-6, stop.trust.radius = 1e-9,
        cg.tol = 1e-6, report.level=4, report.freq=1, report.header.freq=10, report.precision=7),
    config = list(num_threads = 10, strataPerIter=100, transform_theta = TRUE)
  )



cache = new.env()
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)
assign("gamma_start", res$start_gamma, cache)
assign("file", "sim.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))

mle <- BB::spg(par = res$parameters, #mle$solution, 
  fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
        data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner,
           control = list(maxit = 1e4, M = 10, trace=TRUE, checkGrad=TRUE))  # M = nonmonotone history

# element 4 has a problem





basePar = c(0.952726603,  0.416559358, -2.938788, -4.834099, -3.586874, -1.951997)
Dpar = 3
Spar = seq(-0.1, 0.1, len=5) + basePar[Dpar]
parMat = matrix(basePar, ncol=length(Spar), nrow=length(basePar))
parMat[Dpar,] = Spar

Nbeta = nrow(res$tmb_data$XTp)
Sgamma = seq(Nbeta+1, len=nrow(res$tmb_data$ATp))

bob = Matrix::sparseMatrix(i=res$config$sparsity$second$full$i, j=res$config$sparsity$second$full$j, symmetric=TRUE, index1=FALSE)
bob[1:25, -Sgamma]


  bobL =  mapply(loglik, 
   parameters = as.list(as.data.frame(parMat)), 
   MoreArgs = list(data=res$tmb_data, config = res$config, gamma_start=get("gamma_start", cache), control = res$control_inner),
   SIMPLIFY=FALSE)

bobL[[1]]$extra$fullHessian[1:25,-Sgamma]


Slik = unlist(lapply(bobL, function(xx) xx$minusLogLik))
Sdl = unlist(lapply(bobL, function(xx) xx$deriv[Dpar,'dL']))
Sddet = unlist(lapply(bobL, function(xx) xx$deriv[Dpar,'det']))

Sdet = unlist(lapply(bobL, function(xx) xx$extra$halfLogDet))


numD = diff(Slik)/diff(Spar)
plot(Spar, Sdl, type='o', col='red', ylim = range(c(Sdl, numD)))
points(Spar[-1] - diff(Spar)/2, numD)
abline(h=0);abline(v=basePar[Dpar])

numD = diff(Sdet)/diff(Spar)
plot(Spar, Sddet, type='o', col='red', ylim = range(c(numD, Sddet)))
points(Spar[-1] - diff(Spar)/2, numD)

Sh = lapply(bobL, function(xx) xx$hessian)
SdhAd = lapply(bobL, function(xx) xx$extra$dH[[Dpar]])
SdhAd2 = Sdh = list()
for(D in seq(1, length(Sh)-1)) {
  Sdh[[D]] = (Sh[[D+1]] - Sh[[D]])/(diff(Spar)[D])
  SdhAd2[[D]] = (SdhAd[[D+1]] + SdhAd[[D]])/2
}

bob=(SdhAd2[[1]] - Sdh[[1]])

mleX =  trustOptim::trust.optim(
    x = res$parameters, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
    method = 'BFGS',
    control = res$control,
    data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
  )
mle$fval


assign("Nfun", 0, cache)
assign("Ngr", 0, cache)
assign("file", "simbb.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))




assign("file", "sim5.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)

Nbeta = nrow(res$tmb_data$XTp)

mle5 = optim(
  par = res$parameters, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
  method = "L-BFGS-B",
  upper = c(rep(2, Nbeta), rep(0, length(res$parameters)-Nbeta)),
  lower = c(rep(-2, Nbeta), rep(-10, length(res$parameters)-Nbeta)),
  control= list(trace=5, REPORT=20, parscale = rep(c(1e-1, 1),  c(Nbeta, length(res$parameters)-Nbeta))),
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
)

assign("file", "sim2.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)

mle2 = optim(
  par = mle$solution, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
  fn = wrappers_outer$fn,
  method = "Nelder-Mead",
  control= list(trace=5, REPORT=20),
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
)

modGrad = function(...) {
  sum(wrappers_outer$gr(...)^2)
}

wrappers_outer$fn(mle$solution, data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner)
sum(wrappers_outer$gr(mle$solution, data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner)^2)
modGrad(mle$solution, data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner)


assign("file", "sim3.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)
mle3 = optim(
  par = mle2$par, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
  fn = modGrad,
  method = "Nelder-Mead",
  control= list(trace=5, REPORT=20, temp =  2, tmax = 5),
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
)



assign("file", "sim4.txt", cache)
if(file.exists(get("file", cache))) file.remove(get("file", cache))
assign("Nfun", 0, cache)
assign("Ngr", 0, cache)
mle4 = optim(
  par = mle2$par, #c(1, 0.3, log(c(0.05, 0.001, 0.05, 0.001))),
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
  method = "BFGS",
  control= list(trace=5, REPORT=20),
  data=res$tmb_data, config = res$config, cache =  cache, controlInner = res$control_inner
)


Nbeta = nrow(res$tmb_data$XTp)
Sgamma = seq(Nbeta+1, len=nrow(res$tmb_data$ATp))
parameters = res$start_parameters[-Sgamma]

config2 = res$config
config2$beta = parameters[1:Nbeta]
config2$theta = parameters[-(1:Nbeta)]
config2$strataPerIter = 10



wrappers_outer$fn(x=res$parameters, data=res$tmb_data, config=res$config, control=res$control_inner, cache=cache)
cache$gamma_start[1:5]
config2$beta = res$parameters[1:Nbeta]
config2$theta = res$parameters[-(1:Nbeta)]
wrappers_gamma$gr(
  cache$gamma_start,
  res$tmb_data, config2
)[1:5]
bob = loglik(res$parameters, cache$gamma_start, res$tmb_data, config2)
bob$deriv


wrappers_outer$fn(x=res$parameters+1, data=res$tmb_data, config=res$config, control=res$control_inner, cache=cache)
cache$gamma_start[1:5]
loglik(parameters=res$parameters, gamma_start = cache$gamma_start, data=res$tmb_data, 
  config=res$config, control=res$control_innern, deriv=0)

innerOpt = trustOptim::trust.optim(
  x=cache$gamma_start,
  fn = wrappers_gamma$fn,
  gr=wrappers_gamma$gr,
  hs = wrappers_gamma$hs,
  method = 'Sparse',
  control = res$control_inner,
  config=config2, data=res$tmb_data
)


mle$extra = loglik(
  mle$solution, 
  get("gamma_start", cache), 
  res$tmb_data, res$config, res$control_inner)



stuff = trustOptim::trust.optim(
  cache$gamma_start,
  fn = wrappers_gamma$fn,
  gr= wrappers_gamma$gr,
  hs = wrappers_gamma$hs,
  method = 'Sparse',
  dataList = res$tmb_data, 
  configList = config2
 )

allPar = mle$solution
Dpar = 1
Spar = seq(-0.1, 0.1, len=5)


wrappers_outer$gr(

)



  Sigma <- solve(sdr$jointPrecision)
  samp <- mvtnorm::rmvnorm(1e4, mu, as.matrix(Sigma))
  res1 <- constructEffect(fit, exposure_var = "pm", 
                          group_var = "region", group = unique(data$region), 
                          values = seq(knots_pm[1],rev(knots_pm)[1],.01),
                          ref_values = ref_values,
                          pars = samp, probs = c(.1,.5,.9))
  res1sub = res1[res1$group == 1, ]
  matplot(res1sub$var_value, res1sub[,grep("q", colnames(res1sub))], type='l', 
  lty=1, col=c('grey','black','grey'), lwd=2)

#  for(i in unique(res1$group)){
#    x <- genPmEffect(res1$var_value[res1$group == i], region_effect1[i], region_effect2[i])
#    k <- which(res1$var_value[res1$group == i] == ref_values$pm)
#    res1$true[res1$group == i] <- x - x[k]
#  }
#  (res1$true >= res1$q_0.1) & (res1$true <= res1$q_0.9)
  
#}, mc.cores=10)

#saveRDS(res, "applications/synthetic/coverage.rds")



x = seq(0, 1, len=21)[-1]

nu = x = 0.1


thel = unlist(mapply(function(xx) {
  objectiveFunctionC(
    c(rep(0, length(res$start_gamma) + length(res$parameters)-1), log(xx)),
    res$tmb_data, res$config)
  }, 
xx= x))

sumPerStrata
data2

cbind(r=-lR, c=thel)


matplot(x, cbind(-lR, thel), type='l')







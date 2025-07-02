# Rscript applications/synthetic/main-sim.R   > applications/synthetic/log.txt 2>&1


# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(parallel)
library(RhpcBLASctl)
blas_set_num_threads(1)

library(ggplot2)

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

region_effect1 <- rnorm(10,1,.1)
region_effect2 <- rnorm(10,1,.1)
genPmEffect <- function(pm, r1, r2){
  0.25*(r1*pm/10 + pm^(r2/2))
}
genCount <- function(hum, pm, r1=1, r2=1, od = FALSE){
  l <- length(pm)
  od <- rnorm(l, 0, ifelse(od, .1, 0))
  rpois(l, exp(-log(10) + hum + genPmEffect(pm, r1, r2) + od))
}




# model -------------------------------------------------------------------
# knots_pm <- seq(0, 40, 2.5)
ref_values <- list("pm" = 10)
# formula <- count ~ hum + #f(date, model = "iid") +
#   f(pm, model = "hiwp", p=2, ref_value = 10, knots = knots_pm, group_var = region)
cc_design <- ccDesign(time_var = "date", strat_vars = "region")



# sim loop ----------------------------------------------------------------

#res <- mclapply(1:500, \(dummy){
set.seed(0)
  data <- lapply(1:10, \(i){
    data <- data.table(date = as.Date(1:1000), region = i)
    data$hum <- genHum(data$date)
    data$pm <- genPM(data$date)
    data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i], TRUE) 
    data
  }) |> rbindlist()
  
  # knots_pm <- seq(0, max(ceiling(max(data$pm)/5)*5, 20), 2.5)
  knots_pm <- seq(0, max(ceiling(max(data$pm)/5)*5, 20), 2.5)
  # formula <- as.formula(
  #   sprintf("count ~ hum + f(pm, model = 'hiwp', p = 2, ref_value = 10, knots = c(%s), group_var = region)",
  #           paste(knots_pm, collapse = ", "))
  # )  
  formula <- as.formula(
    sprintf("count ~ f(date, model = 'iid') + hum + f(pm, model = 'hiwp', p = 2, ref_value = 10, knots = c(%s), group_var = region)",
            paste(knots_pm, collapse = ", "))
  ) 
  fit <- hnlm(formula, data, cc_design = cc_design)
  sdr <- sdreport(fit$obj, getJointPrecision=TRUE)
  mu <- fit$obj$env$last.par.best

mu[grep("gamma", names(mu), invert=TRUE)]

fit$obj$env$f(mu)
bob = function(xx) fit$obj$env$f(c(xx, mu[-(1:length(xx))]))
numDeriv::hessian(bob, mu[1:5])
sdr$jointPrecision[1:5,1:5]


  saveRDS(list(data=data, par = fit$obj$env$last.par.best, h=sdr$jointPrecision), file='simData.rds') 





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











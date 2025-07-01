# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(5)

library(ggplot2)



# Data --------------------------------------------------------------------
genPM <- function(date){
  l <- length(date)+1
  date_int <- as.integer(date)
  date_int <- c(date_int[1]-1, date_int)
  pm <- 10
  pm <- pm + 5*sin(2*pi * (date_int + rnorm(l,0,1))/365.25)
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
  r1*pm/10 + pm^(r2/3)
}
# region_effect <- rnorm(10,1,.5)
# genPmEffect <- function(pm, r){
#   (r-1) + (r-1)*pm/10
# }
genCount <- function(hum, pm, r1=1, r2=1){
  l <- length(pm)
  rpois(l, exp(hum + genPmEffect(pm, r1, r2)))
}



data <- lapply(1:2, \(i){
  data <- data.table(date = as.Date(1:5000), region = i)
  data$hum <- genHum(data$date)
  data$pm <- genPM(data$date)
  
  knots_pm <- seq(0, ceiling(max(data$pm)/5)*5, 2.5)
  ss <- knots_pm[1]:rev(knots_pm)[1]
  ref_value = 10
  pmE <- genPmEffect(ss, region_effect1[i], region_effect2[i])
  plot(pmE - pmE[which(ss == ref_value)])

  data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i]) 
  data
}) |> rbindlist()



# Variables of interest  --------------------------------------------------

knots_pm <- seq(0, ceiling(max(data$pm)/5)*5, 2.5)
formula <- count ~ hum + #f(date, model = "iid") +
  f(pm, model = "hiwp", p=2, ref_value = 10, knots = knots_pm, group_var = region)
# formula <- count ~ hum + #f(date, model = "iid") +
#   f(pm, model = "iwp", p=2, ref_value = 25, knots = knots_pm)
# formula <- count ~ hum + #f(date, model = "iid") +
#   f(pm, model = "hrpoly", p=1, ref_value = 25, group_var = region)
cc_design <- ccDesign(time_var = "date", strat_vars = "region")
fit <- hnlm(formula, data, cc_design = cc_design)


fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]


# results per group
ref_values <- list("pm" = 10)

res1 <- getEffect(fit, exposure_var = "pm", 
                  group_var = "region", group = unique(data$region), 
                  values = knots_pm[1]:rev(knots_pm)[1],
                  ref_values = ref_values)
res2 <- getEffect(fit, exposure_var = "pm", 
                  group_var = "region", group = NA, 
                  values = knots_pm[1]:rev(knots_pm)[1],
                  ref_values = ref_values)

library(ggplot2)
ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_line(data=res2[,-4], linetype=2) +
  facet_grid(group~variable, scales = "free_x")

ggs <- lapply(c("pm"), \(x){
  ggplot(res1[res1$variable == x,], aes(x=var_value, y=effect_value)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, col="gray75") +
    geom_line() +
    geom_line(data=res2[res2$variable == x,][,-4], linetype=2) +
    facet_grid(group~variable, scales = "free_x")
})
ggs[[1]]


ggs <- lapply(c("pm"), \(x){
  ggplot(res1[res1$variable == x,], aes(x=var_value, y=effect_value, colour=as.factor(group))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, col="gray75") +
    geom_line() +
    geom_line(data=res2[res2$variable == x,][,-4], linetype=2, colour="black") +
    facet_grid(~variable, scales = "free_x")
})
ggs[[1]]




sdr <- sdreport(fit$obj, getJointPrecision=TRUE)
c(sdr$cov.fixed, H2$diag.cov.random)
diag(solve(H2$jointPrecision))

keep1 <- colnames(sdr$cov.fixed) != "theta"
keep2 <- colnames(sdr$jointPrecision) != "theta"
plot(c(diag(sdr$cov.fixed)[keep1], sdr$diag.cov.random))
points(diag(solve(sdr$jointPrecision))[keep2], col=2)

mu <- fit$obj$env$last.par.best
Sigma <- solve(H2$jointPrecision)
samp <- mvtnorm::rmvnorm(1e4, mu, as.matrix(Sigma))

res1 <- constructEffect(fit, exposure_var = "pm", 
                        group_var = "region", group = unique(data$region), 
                        values = knots_pm[1]:rev(knots_pm)[1],
                        ref_values = ref_values,
                        pars = samp, probs = c(.1,.5,.9))
res2 <- constructEffect(fit, exposure_var = "pm", 
                        group_var = "region", group = NA, 
                        values = knots_pm[1]:rev(knots_pm)[1],
                        ref_values = ref_values,
                        pars = samp, probs = c(.1,.5,.9))


ggplot(res1[res1$variable == "pm",], aes(x=var_value, y=q_0.5, colour=as.factor(group))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_ribbon(aes(ymin = q_0.1, ymax = q_0.9, fill = as.factor(group)), col=NA, alpha=.3) + 
  geom_line(data=res2[res2$variable == "pm",][,-4], linetype=2, colour="black") +
  facet_grid(~variable, scales = "free_x")




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


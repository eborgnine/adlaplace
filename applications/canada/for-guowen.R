# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(RhpcBLASctl)
library(timeDate)
blas_set_num_threads(5)

library(ggplot2)

app_path <- "applications/canada/"
store_path <- "../../Data" # hum..


# Data --------------------------------------------------------------------
dataset_path <- file.path(store_path, "sam_data.RData")
dataset_path = '/home/patrick/Downloads/sam_data.RData'
attach(dataset_path)
ls(); search()
data <- as.data.table(sam_data)
detach(pos = grep(paste0("file:",dataset_path), search()))

names(data)
apply(data, 2, \(x) sum(is.na(x)))
data <- data[, which(apply(data, 2, \(x) sum(is.na(x))) == 0), with=F]

names(data)[grepl("Date", names(data))] <- "date"
names(data)[grepl("dailymorb", names(data))] <- "count"
names(data)[grepl("HCtemp", names(data))] <- "temp"
names(data)[grepl("pm25", names(data))] <- "pm25"
names(data)[grepl("no2", names(data))] <- "no2"
names(data)[grepl("o3", names(data))] <- "o3"
names(data)

data$spm25 <- sqrt(data$pm25)
data$sno2 <- sqrt(data$no2)
data$so3 <- sqrt(data$o3)
data <- hpoltest:::removeHolidays(data, "rm_all")


# Fit models --------------------------------------------------------------
knots_pm25 <- seq(floor(min(data$pm25)/5)*5,ceiling(max(data$pm25)/5)*5,5)
knots_spm25 <- sqrt(knots_pm25)
range(data$pm25)

knots_no2 <- seq(floor(min(data$no2)/5)*5,ceiling(max(data$no2)/5)*5,5)
knots_sno2 <- sqrt(knots_no2)
range(data$no2)

knots_o3 <- seq(floor(min(data$o3)/5)*5,ceiling(max(data$o3)/5)*5,5)
knots_so3 <- sqrt(knots_o3)
range(data$o3)

knots_temp <- seq(floor(min(data$temp)/5)*5,ceiling(max(data$temp)/5)*5,5)
range(data$temp)



formula <- count ~ f(date, model = "iid") +
  f(spm25, model = "hiwp", p=2, ref_value = sqrt(15), knots = knots_spm25, group_var = cdcode) +
  # f(sno2, model = "hiwp", p=2, ref_value = 10, knots = knots_sno2, group_var = cdcode) +
  # f(so3, model = "hiwp", p=2, ref_value = 20, knots = knots_so3, group_var = cdcode) +
  f(temp, model = "hiwp", p=2, ref_value = 10, knots = knots_temp, group_var = cdcode)

cc_design <- ccDesign(time_var = "date", strat_vars = "cdcode")
fit <- hnlm(formula, data, cc_design = cc_design)





# results -----------------------------------------------------------------
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]

max(fit$theta_info$id)
split(fit$theta_info$var, fit$theta_info$id)


# new function
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



# results per group
ref_values <- list("spm25" = sqrt(15), "temp" = 10, "no2" = 10, "o3" = 20)
knots <- list("spm25" = knots_spm25, "temp" = knots_temp, "no2" = knots_no2, "o3" = knots_o3)
sdr <- sdreport(fit$obj, getJointPrecision=TRUE)
mu <- fit$obj$env$last.par.best
Sigma <- solve(sdr$jointPrecision)
samp <- mvtnorm::rmvnorm(1e4, mu, as.matrix(Sigma))


# HERE. WHAT TO DO WITH OVERDISPERSION TERMS!! THERE 
all_vars <- sapply(fit$terms, "[[", "var") |> unique()
ress <- lapply(c("temp", "spm25", "sno2", "so3"), \(x){
  
  if(!(x %in% all_vars)) return(NULL)
  
  
  res1 <- constructEffect(fit, exposure_var = x, 
                          group_var = "cdcode", group = unique(data$cdcode), 
                          values = seq(knots[[x]][1],rev(knots[[x]])[1],.1),
                          ref_values = ref_values,
                          pars = samp, probs = c(.1,.5,.9))
  res2 <- constructEffect(fit, exposure_var = x, 
                          group_var = "cdcode", group = NA, 
                          values = seq(knots[[x]][1],rev(knots[[x]])[1],.1),
                          ref_values = ref_values,
                          pars = samp, probs = c(.1,.5,.9))

  if(x != "temp"){
    res1$var_value <- res1$var_value^2
    res2$var_value <- res2$var_value^2
    res1$variable <- res2$variable <- substr(x, 2, 100)
  }
  list(res1 = res1, res2 = res2)
})

res1 <- lapply(ress, "[[", "res1") |> do.call(what = "rbind")
res2 <- lapply(ress, "[[", "res2") |> do.call(what = "rbind")


ggplot(res1, aes(x=var_value, y=q_0.5, group=as.factor(group))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line(alpha=.4) +
  # geom_ribbon(aes(ymin = q_0.1, ymax = q_0.9, fill = as.factor(group)), col=NA, alpha=.05) + 
  geom_ribbon(aes(ymin = q_0.1, ymax = q_0.9), col=NA, alpha=.025) + 
  geom_line(data=res2, linetype=2, colour="black") +
  facet_wrap(~variable, scales = "free")

fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$theta_info
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]


ggplot(res2, aes(x=var_value, y=q_0.5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_ribbon(aes(ymin = q_0.1, ymax = q_0.9, fill = as.factor(group)), col=NA, alpha=.3) + 
  facet_wrap(~variable, scales = "free")

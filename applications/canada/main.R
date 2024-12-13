# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(5)

library(ggplot2)
source("R/removeHolidays.R") # requires lubridate and timeDate

app_path <- "applications/canada/"
store_path <- "/store/samuel/simulations-air-pollution-1" # hum..


# Data --------------------------------------------------------------------
dataset_path <- file.path(store_path, "data/AllDatabase1.RData")
attach(dataset_path)
ls(); search()
data <- as.data.table(AllDatabase)
detach(pos = grep(paste0("file:",dataset_path), search()))


# Variables of interest  --------------------------------------------------

cause0 <- "all"
group0 <- "morb"
data <- data[cause == cause0 & group == group0]
data$spm25lag <- sqrt(data$pm25lag)
data <- data[, .(count = sum(count)), .(group, cause, city, date, dow, temp4lag, pm25lag, spm25lag, hum_mean)]
data <- data[!apply(is.na(data[,c("hum_mean", "pm25lag", "temp4lag")]), 1, any),]
gc()


# Fit models --------------------------------------------------------------

# data <- data[city == "Toronto"]

# data (do we do that?)
# data <- removeHolidays(data, "rm_all")


knots_pm25lag <- seq(floor(min(data$pm25lag)/5)*5,ceiling(max(data$pm25lag)/5)*5,5)
knots_spm25lag <- sqrt(knots_pm25lag)
knots_temp4lag <- seq(floor(min(data$temp4lag)/5)*5,ceiling(max(data$temp4lag)/5)*5,5)
range(data$temp4lag)
formula <- count ~ hum_mean + #od(date) +
  hiwp(spm25lag, p=2, ref_value = sqrt(15), knots = knots_spm25lag, group_var = city) +
  hiwp(temp4lag, p=2, ref_value = 10, knots = knots_temp4lag, group_var = city)
  
cc_design <- ccDesign(time_var = "date", strat_vars = "city")
fit <- hnlm(formula, data, cc_design = cc_design)


fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]


# results per group
ref_values <- list("spm25lag" = sqrt(15), "temp4lag" = 10)

res1_temp <- getEffect(fit, exposure_var = "temp4lag", 
                    group_var = "city", group = unique(data$city), 
                    values = knots_temp4lag[1]:rev(knots_temp4lag)[1],
                    ref_values = ref_values)
res2_temp <- getEffect(fit, exposure_var = "temp4lag", 
                       group_var = "city", group = NA, # for global effect
                       values = knots_temp4lag[1]:rev(knots_temp4lag)[1],
                       ref_values = ref_values)

res1_pm <- getEffect(fit, exposure_var = "spm25lag", 
                     group_var = "city", group = unique(data$city), 
                     values = knots_spm25lag[1]:rev(knots_spm25lag)[1],
                     ref_values = ref_values)
res2_pm <- getEffect(fit, exposure_var = "spm25lag", 
                     group_var = "city", group = NA, 
                     values = knots_spm25lag[1]:rev(knots_spm25lag)[1],
                     ref_values = ref_values)
res1_pm$var_value <- res1_pm$var_value^2
res1_pm$variable <- "pm25lag"
res2_pm$var_value <- res2_pm$var_value^2
res2_pm$variable <- "pm25lag"

res1 <- rbind(res1_pm, res1_temp)
res2 <- rbind(res2_pm, res2_temp)

library(ggplot2)
ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_line(data=res2[,-4], linetype=2) +
  facet_grid(group~variable, scales = "free_x")

ggs <- lapply(c("pm25lag", "temp4lag"), \(x){
  ggplot(res1[res1$variable == x,], aes(x=var_value, y=effect_value)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, col="gray75") +
    geom_line() +
    geom_line(data=res2[res2$variable == x,][,-4], linetype=2) +
    facet_grid(group~variable, scales = "free_x")
})
ggs[[1]]
ggs[[2]]

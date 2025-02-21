# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(5)

library(ggplot2)
source("R/removeHolidays.R") # requires lubridate and timeDate

app_path <- "applications/canada/"
store_path <- "../../Data" # hum..


# Data --------------------------------------------------------------------
dataset_path <- file.path(store_path, "sam_data.RData")
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
data <- removeHolidays(data, "rm_all")


# Fit models --------------------------------------------------------------
knots_pm25 <- seq(floor(min(data$pm25)/5)*5,ceiling(max(data$pm25)/5)*5,5)
knots_spm25 <- sqrt(knots_pm25)
knots_temp <- seq(floor(min(data$temp)/5)*5,ceiling(max(data$temp)/5)*5,5)
range(data$temp)
formula <- count ~ f(date, model = "iid") +
  f(spm25, model = "hiwp", p=2, ref_value = sqrt(15), knots = knots_spm25, group_var = cdcode) +
  f(temp, model = "hiwp", p=2, ref_value = 10, knots = knots_temp, group_var = cdcode)

cc_design <- ccDesign(time_var = "date", strat_vars = "cdcode")
fit <- hnlm(formula, data, cc_design = cc_design)


fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]

max(fit$theta_info$map)
split(fit$theta_info$var, fit$theta_info$map)

# results per group
ref_values <- list("spm25" = sqrt(15), "temp" = 10)

# HERE. WHAT TO DO WITH OVERDISPERSION TERMS!! THERE 
res1_temp <- getEffect(fit, exposure_var = "temp", 
                       group_var = "cdcode", group = unique(data$cdcode), 
                       values = knots_temp[1]:rev(knots_temp)[1],
                       ref_values = ref_values)
res2_temp <- getEffect(fit, exposure_var = "temp", 
                       group_var = "cdcode", group = NA, # for global effect
                       values = knots_temp[1]:rev(knots_temp)[1],
                       ref_values = ref_values)

res1_pm <- getEffect(fit, exposure_var = "spm25", 
                     group_var = "cdcode", group = unique(data$cdcode), 
                     values = knots_spm25[1]:rev(knots_spm25)[1],
                     ref_values = ref_values)
res2_pm <- getEffect(fit, exposure_var = "spm25", 
                     group_var = "cdcode", group = NA, 
                     values = knots_spm25[1]:rev(knots_spm25)[1],
                     ref_values = ref_values)
res1_pm$var_value <- res1_pm$var_value^2
res1_pm$variable <- "pm25"
res2_pm$var_value <- res2_pm$var_value^2
res2_pm$variable <- "pm25"

res1 <- rbind(res1_pm, res1_temp)
res2 <- rbind(res2_pm, res2_temp)

ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_line(data=res2[,-4], linetype=2) +
  facet_grid(group~variable, scales = "free_x")


ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line(aes(col=group)) +
  geom_line(data=res2[,-4], linetype=2) +
  facet_wrap(~variable, scales = "free")



ggs <- lapply(c("pm25", "temp"), \(x){
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


ggs <- lapply(c("pm25", "temp"), \(x){
  ggplot(res1[res1$variable == x,], aes(x=var_value, y=effect_value, colour=group)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, col="gray75") +
    geom_line() +
    geom_line(data=res2[res2$variable == x,][,-4], linetype=2, colour="black") +
    facet_grid(~variable, scales = "free_x")
})
ggs[[1]]
ggs[[2]]

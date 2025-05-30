# Packages ----------------------------------------------------------------
library(hpoltest)
library(data.table)
library(RhpcBLASctl)
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


fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]

max(fit$theta_info$map)
split(fit$theta_info$var, fit$theta_info$map)

# results per group
ref_values <- list("spm25" = sqrt(15), "temp" = 10, "no2" = 10, "o3" = 20)
knots <- list("spm25" = knots_spm25, "temp" = knots_temp, "no2" = knots_no2, "o3" = knots_o3)


# HERE. WHAT TO DO WITH OVERDISPERSION TERMS!! THERE 
all_vars <- sapply(fit$terms, "[[", "var") |> unique()
ress <- lapply(c("temp", "spm25", "sno2", "so3"), \(x){
  
  if(!(x %in% all_vars)) return(NULL)
  
  res1 <- getEffect(fit,
                    exposure_var = x, 
                    group_var = "cdcode", group = unique(data$cdcode), 
                    values = knots[[x]][1]:rev(knots[[x]])[1],
                    ref_values = ref_values)
  res2 <- getEffect(fit, 
                    exposure_var = x, 
                    group_var = "cdcode", group = NA, # for global effect
                    values = knots[[x]][1]:rev(knots[[x]])[1],
                    ref_values = ref_values)
  
  if(x != "temp"){
    res1$var_value <- res1$var_value^2
    res2$var_value <- res2_pm$var_value^2
    res1$variable <- res2$variable <- substr(x, 2, 100)
  }
  list(res1 = res1, res2 = res2)
})

res1 <- lapply(ress, "[[", "res1") |> do.call(what = "rbind")
res2 <- lapply(ress, "[[", "res2") |> do.call(what = "rbind")




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

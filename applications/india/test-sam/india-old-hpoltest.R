
# Test : Trying to run india example with "old" package


# packages ----------------------------------------------------------------
# install.packages("~/Git/hpoltest_0.0.0.9000.tar.gz")
library(hpoltest)
library(data.table)





# setup -------------------------------------------------------------------
# from india.Rmd

SstoreDir = c('/store/patrick/copernicus/agg', '~/Downloads')
storeDir = SstoreDir[1]
if(!dir.exists(storeDir)) storeDir = SstoreDir[2]

deaths = as.data.table(readRDS( file.path(storeDir, "indiaCC.rds")))
indiaEcoRegions = terra::unwrap(readRDS(file.path(storeDir, "indiaEcoRegions.rds")))

# transform variables
deaths$sqrtPm = sqrt(pmin(deaths$pm,1000))
deaths$timeIid = as.numeric(deaths$date)
deaths$tempClean = pmax(deaths$temp, -5)
deaths$dom = as.numeric(format(deaths$date, '%d'))
deaths[deaths$dom == Hmisc::monthDays(deaths$date), 'dom'] = 31
deaths[format(deaths$date, '%m%d') == '0101', 'dom'] = 0
deaths$dom = relevel(factor(deaths$dom), '1')
deaths$month = format(deaths$date, '%m')

dataMedSimple = deaths[date >= as.Date('2006/1/1') & format(date, '%b%d') != 'Jan01' & tempClean >=0]
dataMedSimpler = deaths[sample(1e3)]




# fit model ---------------------------------------------------------------

# frm = medical ~ f(date, model='iid') +
#   f(tempClean, model='hiwp', ref_value = 20, p=3, knots = seq(0, 45, by=2.5), group_var = ecoRegion)
frm = medical ~ f(tempClean, model='hiwp', ref_value = 20, p=3, knots = seq(0, 45, by=5), group_var = ecoRegion)

my_data <- dataMedSimple
# my_data <- dataMedSimple[ecoRegion <= 7]
my_data <- dataMedSimple[ecoRegion <= 3]
table(dataMedSimple$ecoRegion)

dataMedSimple[, .(mean(medical)), .(ecoRegion)]
dataMedSimple[, .(mean(total)), .(ecoRegion)]


fit <- hnlm(formula = frm,
            data = my_data,
            cc_design = ccDesign(strat_vars = c("cell","ecoRegion", "yearWeek3Dow")),
            # tmb_parameters = list(beta = fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"],
            #                       gamma = fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"],
            #                       theta = c(4, 4, 4, 4, 4, 4)),
            # for_dev = TRUE)
            for_dev = FALSE)

saveRDS(fit, "applications/india/test-sam/fit_old_hpoltest.rds")
# H_tmb_old_hpoltest <- fit$obj$env$spHess()
# saveRDS(H_tmb_old_hpoltest, "applications/india/test-sam/H_tmb_old_hpoltest.rds")

image(H_tmb_old_hpoltest)




fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "gamma"]
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]
fit$theta_info

# max(fit$theta_info$map)
# split(fit$theta_info$var, fit$theta_info$map)

# results per group
ref_values <- list("tempClean" = 20)
# knots <- list("tempClean" = seq(0, 45, by=2.5))
knots <- list("tempClean" = seq(0, 45, by=5))


# HERE. WHAT TO DO WITH OVERDISPERSION TERMS!! THERE 
all_vars <- sapply(fit$terms, "[[", "var") |> unique()
ress <- lapply(c("tempClean"), \(x){
  
  if(!(x %in% all_vars)) return(NULL)
  res1 <- getEffect(fit,
                    exposure_var = x, 
                    group_var = c("ecoRegion"), group = unique(my_data$ecoRegion), 
                    values = knots[[x]][1]:rev(knots[[x]])[1],
                    ref_values = ref_values)
  res2 <- getEffect(fit, 
                    exposure_var = x, 
                    group_var = c("ecoRegion"), group = NA, 
                    values = knots[[x]][1]:rev(knots[[x]])[1],
                    ref_values = ref_values)
  
  list(res1 = res1, res2 = res2)
})

res1 <- lapply(ress, "[[", "res1") |> do.call(what = "rbind")
res2 <- lapply(ress, "[[", "res2") |> do.call(what = "rbind")

library(ggplot2)
ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_line(data=res2[,-4], linetype=2) +
  facet_grid(group~variable, scales = "free_x")

ggplot(res1, aes(x=var_value, y=effect_value)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line(aes(group=group, col=as.factor(group))) +
  geom_line(data=res2[,-4], linetype=2) +
  facet_grid(~variable, scales = "free_x")





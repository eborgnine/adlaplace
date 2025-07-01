
# Rscript applications/india/test-sam/tmb-old-hpoltest.R > applications/india/test-sam/output_tmb.txt 2>&1


# packages ----------------------------------------------------------------
# Should be using install.packages("~/Git/hpoltest_sam.tar.gz")
library(hpoltest)
library(data.table)


fit <- readRDS("applications/india/test-sam/fit_old_hpoltest.rds")
fit$obj$env$last.par.best

H_nD_old_hpoltest <- numDeriv::hessian(fit$obj$fn, fit$obj$env$last.par.best)
# my_fn <- \(x) fit$obj$fn(c(x, fit$obj$env$last.par.best[-seq_along(x)]))
# H_nD_old_hpoltest <- numDeriv::hessian(my_fn, fit$obj$env$last.par.best[1:10])
saveRDS(H_nD_old_hpoltest, "applications/india/test-sam/H_nD_old_hpoltest.rds")

H_tmb_old_hpoltest <- TMB::sdreport(fit$obj, getJointPrecision = TRUE)
saveRDS(H_tmb_old_hpoltest, "applications/india/test-sam/H_tmb_old_hpoltest.rds")



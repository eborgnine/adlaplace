

# Rscript applications/india/test-sam/numDeriv-old-hpoltest.R > applications/india/test-sam/output.txt 2>&1

# packages ----------------------------------------------------------------
# Should be using install.packages("~/Git/hpoltest_patrick.tar.gz")
# (which I got from R CMD BUILD hpoltest_orig)
library(hpoltest)
library(data.table)


fit <- readRDS("applications/india/test-sam/fit_new_hpoltest.rds")
my_fn <- \(x) fit$obj$fn(c(x, fit$obj$env$last.par.best[-seq_along(x)]))

# H_nD_old_hpoltest <- numDeriv::hessian(fit$obj$fn, fit$obj$env$last.par.best)
H_nD_new_hpoltest <- numDeriv::hessian(my_fn, fit$obj$env$last.par.best[1:4])
saveRDS(H_nD_new_hpoltest, "applications/india/test-sam/H_nD_new_hpoltest.rds")


x <- fit$obj$env$last.par.best[1:4]
fit$obj$fn(fit$obj$env$last.par.best[-seq_along(x)])

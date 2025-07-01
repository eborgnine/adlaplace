# oldies
H_tmb_old_hpoltest <- readRDS("applications/india/test-sam/H_tmb_old_hpoltest.rds")
H_nD_old_hpoltest <- readRDS("applications/india/test-sam/H_nD_old_hpoltest.rds")

diag(Matrix::solve(H_tmb_old_hpoltest$jointPrecision))[1:10]
diag(H_nD_old_hpoltest)

1/diag(Matrix::solve(H_tmb_old_hpoltest$jointPrecision))[1:10]
round(diag(H_nD_old_hpoltest), 5)

# newies
# H_tmb_new_hpoltest <- readRDS("applications/india/test-sam/H_tmb_new_hpoltest.rds")
# H_nD_new_hpoltest <- readRDS("applications/india/test-sam/H_nD_new_hpoltest.rds")

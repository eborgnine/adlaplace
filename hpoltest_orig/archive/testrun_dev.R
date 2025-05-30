compile("src/hpoltest_dev.cpp")
dyn.load(dynlib("src/hpoltest_dev"))

tmb_data <- list(
  X = X, A = A, y = y,
  # gamma_nreplicate = gamma_info$nreplicate, # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
  Q = Qs |> .bdiag(),
  # theta_id = theta_info$id
  gamma_split = gamma_info$split
)

tmb_parameters <- list(
  beta = rep(0, ncol(X)),
  gamma = rep(0, ncol(A)),
  theta = theta_info$init
)

theta_info$id[theta_info$id == 0] <- max(theta_info$id) + 1:sum(theta_info$id == 0)
map <- list(theta = factor(theta_info$id))

obj <- MakeADFun(data = tmb_data, 
                 parameters = tmb_parameters,
                 random = c("gamma"),
                 map = map,
                 DLL = "hpoltest_dev")


obj$fn(obj$par)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
              control = list(eval.max=2000, iter.max=2000))

subpar <- tmb_parameters[names(tmb_parameters)!="gamma"] |> unlist()
opt$par - (subpar)

plot(opt$par)
points(subpar, col=2)



opt$par - c(beta_true, theta_true)


# packages ----------------------------------------------------------------
library(TMB)
library(Matrix)
library(data.table)

n0 <- 200
nc <- 3
np <- 3
n <- nc*n0

db <- 5
X <- matrix(rnorm(n*db), n, db) |> abs() |> as("dgTMatrix")

dg0 <- 20
dg <- np*(nc+1)*dg0
A <- replicate(nc, matrix(rnorm(n0*dg0*np), n0, dg0*np) |> abs(), simplify = F)
A <- cbind(do.call("rbind", A), .bdiag(A))
dim(A)
image(A)
gamma_split <- rep(dg0, (nc+1)*np)

theta_id <- c(1:np, replicate(nc, np + 1:np)) - 1
theta_true <- log(rexp(2*np, .0001))  # Log precision for each pollutant (will be exp(theta) to scale Q)

# Qs <- lapply(1:np, \(j){
#   (diag(gamma_split[j])*100) |> as("dgTMatrix") 
# })
createD <- function(d, p) {
  if (p == 0) return(Diagonal(d, 1))
  D <- Matrix::bandSparse(d, k = c(0, 1), diagonals = list(rep(-1, d), rep(1, d - 1)))[-d, ]
  if (p == 1) return(D)
  return(createD(d, p - 1)[-1, -1] %*% D)
}

Qs <- lapply(1:np, \(j){
  Q <- Matrix::crossprod(createD(gamma_split[j], 2))
  kk <- round(ncol(Q)/2) + c(0,1)
  diag(Q) <- diag(Q) + .01
  Q
})


# Observed counts (Poisson responses)
beta_true <- rexp(ncol(X), 5)
gamma_true <- lapply(seq_along(gamma_split), \(k){
  S <- solve(exp(theta_true[theta_id[k]+1]) * Qs[[(theta_id[k] %% np) + 1]])
  mvtnorm::rmvnorm(1, sigma = S |> as.matrix())
}) |> unlist()
# plot(gamma_true)

lambda_true <- exp(X %*% beta_true + A %*% gamma_true) |> as.numeric()
y <- rpois(n, lambda_true)


# Test run ----------------------------------------------------------------

# Create a list of inputs
data_list <- list(
  X = X, A = A, y = y,
  Q = Qs |> .bdiag(), 
  log_dets = sapply(Qs, \(q) Matrix::determinant(q, logarithm=T)$modulus),
  gamma_split = gamma_split,
  theta_id = theta_id,
  nc = nc, np = np
)

# parameters_list <- list(
#   beta = beta_true,
#   # gamma = gamma_true + .001,
#   gamma = rep(0, length(gamma_true)),
#   theta = theta_true
# )
parameters_list <- list(
  beta = rep(0, length(beta_true)),
  gamma = rep(0, length(gamma_true)),
  theta = rep(6, length(theta_true))
)


# Compile the model
compile("src/hpoltest.cpp")
dyn.load(dynlib("src/hpoltest"))

# Create the object
map <- list(
  beta = seq_along(parameters_list$beta) |> factor(),
  # beta = rep(NA, length(parameters_list$beta)) |> factor(),
  gamma = seq_along(parameters_list$gamma) |> factor(),
  # gamma = rep(NA, length(parameters_list$gamma)) |> factor(),
  # gamma = c(rep(NA, sum(gamma_split[1:np])),
  #           seq_along(parameters_list$gamma)[-(1:sum(gamma_split[1:np]))]) |> factor(),
  theta = seq_along(parameters_list$theta) |> factor()
  # theta = rep(NA, length(parameters_list$theta)) |> factor()
)

obj <- MakeADFun(data = data_list, 
                 parameters = parameters_list,
                 DLL = "hpoltest",
                 random = c("gamma"),
                 # map = map,
                 checkParameterOrder = TRUE)

obj$fn()
obj$env$last.par.best - c(beta_true, gamma_true, theta_true)
plot(gamma_true)
points(obj$env$last.par.best[names(obj$env$last.par.best) == "gamma"], col=2)
# plot(gamma_true - obj$env$last.par.best[names(obj$env$last.par.best) == "gamma"], col=2)
# 
# opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
#               control = list(eval.max=200, iter.max=200))
# beta_id <- which(names(obj$par) == "beta")
# theta_temp <- theta_true + rnorm(length(theta_true),0,.2)
# opt <- nlminb(start = obj$par[beta_id], 
#               objective = \(x) obj$fn(c(x, theta_temp)), 
#               gradient =  \(x) obj$gr(c(x, theta_temp))[beta_id], 
#               control = list(eval.max=200, iter.max=200))


opt <- nlminb(start = obj$par, objective = \(x){print(x); obj$fn(x)}, gradient = obj$gr,
              control = list(eval.max=200, iter.max=200))



# plot(opt$par - (parameters_list |> unlist())[!is.na(unlist(map))])

# par_true <- c(beta_true, gamma_true, theta_true)[!is.na(unlist(map))]
par_true <- c(beta_true, theta_true)
plot(opt$par, ylim=range(c(opt$par, par_true)))
points(par_true, col=2)

obj$env$last.par.best - c(beta_true, gamma_true, theta_true)
plot(gamma_true)
points(obj$env$last.par.best[names(obj$env$last.par.best) == "gamma"], col=2)


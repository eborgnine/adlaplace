
# packages ----------------------------------------------------------------
library(TMB)
library(Matrix)
library(data.table)


# setup -------------------------------------------------------------------
set.seed(123)

nn <- c(50,100,150)*10
nc <- length(nn)
np <- 2 # number of exposures
n <- sum(nn)
data <- data.table(
  id = 1:n,
  city = lapply(1:nc, \(k) rep(LETTERS[k], nn[k])) |> unlist(),
  dow = factor(weekdays(as.Date(1:n), T), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")),
  hum = runif(n),
  pol1 = rgamma(n, 20, .9),
  pol2 = rgamma(n, 50, 1.25)
)


# order data by city and time
cities <- unique(data$city)
head_id <- lapply(cities, \(ct) which(data$city == ct)[1:5]) |> unlist()

# Fixed effects design matrix X
betai_names <- sparse.model.matrix(~ 0 + dow + hum + pol1 + pol2, data[1]) |> colnames()
beta_names <- lapply(cities, \(ct) paste0(substr(ct, 1, 3), ":", betai_names)) |> unlist()
X <- lapply(cities, \(ct){
  sparse.model.matrix(~ 0 + dow + hum + pol1 + pol2, data[city == ct],
                      drop.unused.levels = F)
}) |> .bdiag()
colnames(X) <- beta_names
image(X)
image(X[head_id,])

# Random effects design matrix A
pol1_bins <- round(data$pol1/2) |> range()
pol2_bins <- round(data$pol2/2) |> range()

data[, pol1_bin := factor(round(pol1/2), levels = pol1_bins[1]:pol1_bins[2])]
data[, pol2_bin := factor(round(pol2/2), levels = pol2_bins[1]:pol2_bins[2])]

gammai_names <- sparse.model.matrix(~ 0 + pol1_bin + pol2_bin, data[1]) |> colnames()
gamma_names <- lapply(c("global", cities), \(ct) paste0(substr(ct, 1, 3), ":", gammai_names)) |> unlist()
A <- cbind(
  sparse.model.matrix(~ 0 + pol1_bin + pol2_bin, data, drop.unused.levels = F),
  lapply(cities, \(ct){
    sparse.model.matrix(~ 0 + pol1_bin + pol2_bin, data[city == ct], drop.unused.levels = F)
  }) |> .bdiag()
)
colnames(A) <- gamma_names
image(A)
image(A[head_id,])

# Log-precision parameters (2 parameters for 2 pollutants)
theta_id <- c(1:np, replicate(nc, np + 1:np)) - 1
theta_true <- log(rexp(2*np, .1))  # Log precision for each pollutant (will be exp(theta) to scale Q)

gamma_split = c("pol1", "pol2") |> sapply(\(x) length(grep(x, gammai_names))) |> rep(nc+1)
createD <- function(d, p) {
  if (p == 0) return(Diagonal(d, 1))
  D <- Matrix::bandSparse(d, k = c(0, 1), diagonals = list(rep(-1, d), rep(1, d - 1)))[-d, ]
  if (p == 1) return(D)
  return(createD(d, p - 1)[-1, -1] %*% D)
}

Qs <- lapply(1:np, \(j){
  Q <- Matrix::crossprod(createD(gamma_split[j], 2))
  kk <- round(ncol(Q)/2) + c(0,1)
  Q[,kk] <- Q[kk,] <- 0
  Q[kk[1],kk[1]] <- Q[kk[2],kk[2]] <- 1
  Q
})


# Observed counts (Poisson responses)
beta_true <- rexp(ncol(X), 40)
gamma_true <- lapply(seq_along(gamma_split), \(k){
  S <- solve(exp(theta_true[theta_id[k]+1]) * Qs[[(theta_id[k] %% np) + 1]])
  mvtnorm::rmvnorm(1, sigma = S |> as.matrix())
}) |> unlist()
plot(gamma_true)

lambda_true <- exp(X %*% beta_true + A %*% gamma_true) |> as.numeric()
y <- rpois(n, lambda_true)

# precision matrix for random effects
# Q <- sparseMatrix(i = 1:length(gammai_names), j = 1:length(gammai_names),
#                   repr = "T") |> as("dgTMatrix")







# Test run ----------------------------------------------------------------
# Compile the model
compile("src/hpoltest.cpp")
dyn.load(dynlib("src/hpoltest"))

# Create a list of inputs
data_list <- list(
  X = X, A = A, y = y,
  gamma_split = gamma_split,
  Q = Qs |> .bdiag(), 
  log_dets = sapply(Qs, \(q) Matrix::determinant(q, logarithm=T)$modulus),
  # log_dets = rep(1, 2*np),
  theta_id = theta_id,
  nc = nc, np = np
)

parameters_list <- list(
  beta = beta_true,
  gamma = gamma_true,
  theta = theta_true
)
# parameters_list <- list(
#   beta = rep(.001, length(beta_true)),
#   gamma = rep(.001, length(gamma_true)),
#   theta = rep(-1, length(theta_true))
# )
parameters_list <- list(
  beta = beta_true + rnorm(length(beta_true), sd=.1),
  gamma = gamma_true + rnorm(length(gamma_true), sd=.1),
  theta = theta_true
)

# Create the object
map <- list(
  beta = seq_along(parameters_list$beta) |> factor(),
  # gamma = seq_along(parameters_list$gamma) |> factor(),
  gamma = rep(NA, length(parameters_list$gamma)) |> factor(),
  theta = seq_along(parameters_list$theta) |> factor()
  # theta = rep(NA, length(parameters_list$theta)) |> factor()
)
obj <- MakeADFun(data = data_list, 
                 parameters = parameters_list,
                 # random = c("gamma"),
                 map = map,
                 DLL = "hpoltest")

# Run optimization
obj$fn(obj$par)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
              control = list(eval.max=2000, iter.max=2000))

opt$par - (parameters_list |> unlist())

plot(opt$par)
points(unlist(parameters_list), col=2)

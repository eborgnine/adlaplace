Rcpp::compileAttributes()
devtools::load_all()

objectiveFunction = function(parameters, data, 
    config=list(dirichelet=TRUE)) {
	data = formatHpolData(data)
	result = objectiveFunctionC(parameters, data, config)
	try(result$hessian <- do.call(Matrix::sparseMatrix, resut$hessian))
	result
}

wrappers <- make_trustoptim_wrappers(formatHpolData(data))

library('trustOptim')

start_val <- c(rep(0, ncol(data$X)), rep(0, ncol(data$A)), rep(1, length(unique(data$map))+1))


result <- trust.optim(
  x = start_val,
  fn = wrappers$fn,
  gr = wrappers$gr,
  hs = wrappers$hs,
  method = "Sparse"
)


# Example inputs (adjust to your actual data)
beta <- rnorm(3)      # Fixed effects coefficients
gamma <- rnorm(5)     # Random effects
theta <- c(1.0, 0.5)  # Variance parameters
 # Design matrix (fixed effects)
A <- matrix(rnorm(30), ncol = 5)  # Design matrix (random effects)
y <- rpois(10, 1)     # Observed counts
Q <- diag(5)          # Precision matrix
gamma_split <- c(2, 3) # Grouping for random effects
psd_scale <- c(1.0, 1.0)  # Scaling factors
cc_matrix <- matrix(1:10, ncol = 2)  # Case-control matrix
dirichlet <- 0         # 0=Multinomial, 1=Dirichlet-Multinomial

# Call the function
result <- evaluate_objective(
    beta, gamma, theta,
    list(X))
#, A, y, Q,
#    gamma_split, psd_scale,
#    cc_matrix, dirichlet
#)

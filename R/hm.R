#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param x A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param y A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details Provide any additional details about the function, such as edge cases, assumptions, or implementation notes.
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
hm <- function(formula, data, cc_design = ccDesign(), weight_var, for_dev = F) {
  
  data <- as.data.frame(data)
  
  # # Check inputs
  # if (!is(formula, "formula")) stop("formula must be a formula.")
  # if (!is.data.frame(data)) stop("data must be a data.frame.")
  # if (!missing(weight_var) && !is.character(weight_var) && !(weight_var %in% colnames(data)))
  #   stop("weight_var must be a character vector.")

  # Order the rows of data appropriately.
  if(is.null(cc_design$strat_vars) & is.null(cc_design$time_var)) stop("Provide a valid stratification (or time) variable.")
  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)
  new_order <- eval(str2lang(paste0("order(", paste0("data[['", strat_time_vars, "']]", collapse = ", "), ")")))
  data <- data[new_order,]
  
  # setup the data for case-crossover
  cc_matrix <- setStrata(cc_design = cc_design, data = data)
  
  # setup of the design matrices and other parameters
  # terms carries all the information throughout
  terms <- collectTerms(formula)
  
  # design matrices
  X <- NULL # fixed effects
  A <- NULL # random effects
  Qs <- list() # precisions for random effects
  
  beta_info <- list()
  gamma_info <- list()
  theta_info <- list()

  # loop
  k <- 1
  while(k <= length(terms)){
    term <- terms[[k]] |> getExtra(data=data)
    term$id <- k
    if(!is.factor(data[[term$var]][1]) && is.null(term$range)) range <- term$range <- range(data[[term$var]])
    
    if(term$run_as_is){
      Xsub <- sparse.model.matrix(term$f, data)
      beta_info$names <- c(beta_info$names, colnames(Xsub))
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }
    

    if(term$type %in% "fpoly"){
      Xsub <- poly(data[[term$var]] - term$ref_value, raw = T, simple = T) |> as("dgTMatrix")
      beta_info$names <- c(beta_info$names, colnames(Xsub))
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }    
    
    # below takes care of random effects
    
    # design matrix
    Asub <- getDesign(term, data)
    A <- cbind(A, Asub)
    gamma_info$names <- c(gamma_info$names, rep(paste0(term$var, "_", term$id), ncol(Asub)))
    gamma_info$split <- c(gamma_info$split, getGammaSplit(term))
    # if(!is.null(term$group_var)) gamma_info$nrep <- c(gamma_info$nrep, length(split(1:nrow(data), data[[term$group_var]])))
    # Note: for iwp 1 knot removed for constraints
    
    # Add fized and random polynomial effects
    terms <- c(terms, addFPoly(term), addRPoly(term))
    
    
    # precision matrix
    Qs[[length(Qs) + 1]] <- getPrecision(term)
    
    # theta parameters
    theta_setup <- getThetaSetup(theta_info, term)
    
    theta_info$id <- c(theta_info$id, theta_setup$id)
    theta_info$init <- c(theta_info$init, theta_setup$init)

    # update term with new elements
    terms[[k]] <- term
    k <- k+1
  }
  
  if(for_dev) return(list(X = X, A = A, gamma_split = gamma_info$split, Qs = Qs, theta_info=theta_info, new_order = new_order))

  y <- data[[all.vars(formula)[1]]]
  tmb_data <- list(
    X = X, A = A, y = y,
    # gamma_nreplicate = gamma_info$nreplicate, # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
    Q = Qs |> .bdiag(), 
    gamma_split = gamma_info$split,
    # theta_id = theta_info$id
    cc_matrix = cc_matrix
  )
  
  tmb_parameters <- list(
    beta = rep(0, ncol(X)),
    gamma = rep(0, ncol(A)),
    theta = theta_info$init
  )

  theta_info$id[theta_info$id == 0] <- max(theta_info$id) + 1:sum(theta_info$id == 0)
  map <- list(theta = factor(theta_info$id))

  dyn.load(dynlib("src/hpoltest"))
  obj <- MakeADFun(data = tmb_data, 
                   parameters = tmb_parameters,
                   random = c("gamma"),
                   map = map,
                   DLL = "hpoltest")
  
  # Run optimization
  obj$fn(obj$par)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
                control = list(eval.max=2000, iter.max=2000))
  
  # Return the result
  return(obj)
}

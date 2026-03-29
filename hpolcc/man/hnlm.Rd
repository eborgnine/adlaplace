\name{hnlm}
\alias{hnlm}
\title{Fit Hierarchical Non-Linear Models}
\description{
  This function fits a hierarchical model using the specified formula and data, 
  incorporating case-crossover designs, fixed and random effects, and precision 
  matrices for random effects.
}
\usage{
  hnlm(formula, data, cc_design = c(), 
       config = list(dirichlet = TRUE, boundary_is_random = FALSE, 
                     transform_theta = TRUE), 
       control = list(fnscale = -1, trace = 3, REPORT = 1, maxit = 1000, 
                      parscale = NULL), 
       control_inner = list(report.level = 0), 
       for_dev = FALSE, ...)
}
\arguments{
  \item{formula}{A formula object specifying the model to be fitted.}
  \item{data}{A data frame containing the variables specified in the formula.}
  \item{cc_design}{An object specifying the case-crossover design.}
  \item{config}{A list of configuration options.}
  \item{control}{A list of control options for optimization.}
  \item{control_inner}{A list of control options for inner optimization.}
  \item{for_dev}{Logical; if TRUE, returns intermediate objects for development.}
  \item{...}{Additional arguments passed to methods.}
}
\details{
  The function handles fixed effects, random effects, and their associated 
  precision matrices. It also optimizes the model using TMB with options for 
  additional preprocessing and handling specific random effect structures.
}
\value{
  A list containing the fitted model object and related information.
}
\seealso{
  \code{\link{f}} for specifying model terms
}
\examples{
  # See package vignette for basic usage examples
}

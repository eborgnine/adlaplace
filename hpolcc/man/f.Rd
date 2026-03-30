\name{f}
\alias{f}
\title{Functions to Specify Terms of the Formula}
\description{
  This suite of functions implements various models (fixed and random effects), 
  polynomials, integrated Wiener processes, and hierarchical extensions. The 
  primary interface for these models is the f() function, which allows 
  specification of model terms via the model argument.
}
\usage{
  f(x, model = "iid", ...)
}
\arguments{
  \item{x}{A numeric vector representing the data to which the model is applied.}
  \item{model}{A character string specifying the model of model term. Options include:
    \describe{
      \item{"fpoly"}{Fixed polynomial effect terms.}
      \item{"rpoly"}{Random polynomial effect terms.}
      \item{"hrpoly"}{Random hierarchical polynomial effect terms.}
      \item{"iwp"}{Integrated Wiener process terms (random).}
      \item{"hiwp"}{Hierarchical integrated Wiener process terms (random).}
    }
  }
  \item{p}{A numeric value specifying the polynomial degree (default: varies by function).}
  \item{ref_value}{A numeric value specifying the reference point for the basis functions.}
  \item{knots}{A numeric vector specifying the knot locations for spline-based models.}
  \item{range}{A numeric vector of length 2 specifying the range of the data (optional).}
  \item{by}{A factor or grouping variable used to define hierarchical structures (optional).}
  \item{include_global}{A logical value indicating whether to include a global component in hierarchical models (default: TRUE).}
  \item{...}{Additional arguments passed to other methods or functions.}
}
\details{
  The f() function serves as the primary entry point for specifying model terms. 
  Internally, it dispatches to specific functions based on the model argument.
}
\seealso{
  \code{\link{hnlm}} for the main modeling function
}


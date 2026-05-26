% Generated manually
\docType{data}
\name{gamm}
\alias{gamm}
\title{Example GAMM data and fitted model}
\format{
This dataset loads two objects into the workspace:
\describe{
  \item{\code{dat}}{A data frame from \code{mgcv::gamSim(6, n = 200, scale = 0.2, dist = "poisson")}.}
  \item{\code{b2}}{A fitted \code{mgcv::gamm} object from \code{y ~ s(x2)} with random intercept \code{fac}.}
}
}
\source{
Generated with \code{mgcv::gamSim()} and \code{mgcv::gamm()}.
}
\examples{
data("gamm", package = "adlaplace")
exists("dat")
exists("b2")
}
\keyword{datasets}

# admvn

R package for multivariate normal cumulative probabilities with derivatives
with respect to the upper integration limits, computed via CppAD automatic
differentiation of Genz's quasi-Monte Carlo algorithm.

## Requirements

- CppAD headers (e.g. `brew install cppad`)

## Build

```sh
make admvn
R CMD INSTALL admvn_0.1.0.tar.gz
```

## Usage

```r
library(admvn)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
out <- pmvn(c(1, 1), lower = c(-Inf, -Inf), sigma = sigma)
out$value
out$gradient
out$hessian
```

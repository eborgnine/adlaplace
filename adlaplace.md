---
title: "Laplace Approximations with adlaplace"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Laplace Approximations with adlaplace}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


Neg Binoma, $E(y) = \mu$, $\theta$ is size, $1/\sqrt(\theta)$ is sd of Gamma random effect

$$
\begin{aligned}
f(y;\mu,\theta) =  & \log\Gamma(y + r) - \log\Gamma(r) - \log\Gamma(y + 1) 
	+ r \log(\frac{r}{r + \mu})
	+ y \log(\frac{\mu}{r + \mu}) \\
	= &  \log\Gamma(y + r) - \log\Gamma(r) - \log\Gamma(y + 1)  + r \log r - r \log(r + \mu) + y \log(\mu) - y\log(r + \mu)
\end{aligned}
$$

To Do:

- non-diagonal precision matrix
- getAdFun in parallel
- precision matrix has non-multiplicative parameters, AD cholesky


``` r
Nobs = 1000;Nrandom1=10;Nrandom2=25

set.seed(0)
X = Matrix::Matrix(cbind(1, rbinom(Nobs,1, prob= 0.5)))
AmatList = list(
	Matrix::sparseMatrix(
		i = 1:Nobs,
		j= sample(Nrandom1, Nobs, replace=TRUE)
	),
	Matrix::sparseMatrix(
		i = 1:Nobs,
		j= sample(Nrandom2, Nobs, replace=TRUE)
	)
)
Amat = do.call(cbind, AmatList)
beta = rep(1, ncol(X))
thetaOrig = c(0.1, 0.1, 0.1);theta = log(thetaOrig)
gamma = rnorm(ncol(Amat), sd=rep(thetaOrig[1:2], c(Nrandom1, Nrandom2)))
eta = as.vector(X %*% beta + Amat %*% gamma)
etaT = exp(eta)
Z = rgamma(length(etaT), 0.1^(-2), 0.1^(-2))
c(mean(Z), sd(Z))
```

```
## [1] 1.00976352 0.09797471
```

``` r
y = rpois(length(etaT), etaT*Z)

NperEffect = unlist(lapply(AmatList, ncol))   
#map = rep(seq(0, len=length(NperEffect)), NperEffect)
map = Matrix::sparseMatrix(
	i = seq(0L, len=ncol(Amat)),
	j = rep(seq(0L, len=length(NperEffect)), NperEffect),
	x = 1L,
	index1=FALSE,
	dims = c(ncol(Amat), length(thetaOrig))
	)

data = list(y = y, 
	ATp = as(Matrix::t(Amat), 'dMatrix'),
	XTp = as(Matrix::t(X), 'CsparseMatrix'),
	map = map, 
	Qdiag = rep(1, ncol(Amat)))


config = list(beta = rep(1, nrow(data$XTp)), 
	theta = c(rep(-1, length(NperEffect)), -1), 
	transform_theta=TRUE, 
	gamma = rep(0, nrow(data$ATp)),
	groups = adlaplace::adFun_groups(data$ATp, 1000),
	num_threads=parallel::detectCores(),
	verbose=TRUE, package='adlaplace')	
```


model is defined in the file `src/objectiveFunction.cpp`




``` r
adFun <- adlaplace::getAdFun(data, config)
```

```
## outer, groups 245 Nbeta 2 Ntheta 3 Ngamma 35 Nparams 40
## q ngamma  35 map.ncol 3 map@i.size 35 ntheta 3 exp theta map 0 0.367879 gamma0 0.
## logDensRandom 2.83715 qpart 0 det -2.83715
## extra .. done.
## 12..3456
```

``` r
bob = adlaplace::hessianMap(adFun$sparsity, 
	length(config$beta), length(config$gamma), length(config$theta))
identical(bob$map, adFun$hessians$map)
```

```
## [1] TRUE
```

``` r
identical(bob$hessian, adFun$hessians$hessian)
```

```
## [1] TRUE
```

``` r
res = adlaplace::logLikLaplace(
	x=c(config$beta, config$theta),
	adFun = adFun,
	config=modifyList(config, list(verbose=FALSE)),
	deriv=TRUE)
```

```
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2190.840694    2.147475                 Continuing    1.000000     10          Reached tolerance
##   2   2190.832081    0.001210                 Continuing    1.000000      9          Reached tolerance
##   3   2190.832081    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2190.832081    0.000000                    Success
```

``` r
res$parameters
```

```
## [1]  1  1 -1 -1 -1
```

``` r
res$fval
```

```
## [1] 2245.111
```

``` r
res$grad	
```

```
## [1]   0.4894629  -9.3982090   8.9584896  20.4065180 131.7467515
```




``` r
adFun <- adlaplace::getAdFun(data, config)

cache <- new.env(parent = emptyenv())
cache$gamma <- config$gamma

x0 <- c(config$beta, config$theta)

adlaplace::outer_fn(x=x0, cache=cache, config=config, adFun = adFun)
adlaplace::outer_gr(x=x0, cache=cache, config=config, adFun = adFun)

outer_fit <- trustOptim::trust.optim(
  x = x0,
  fn = adlaplace::outer_fn,
  gr = adlaplace::outer_gr,
  method = "SR1",
  config = config,
  adFun = adFun,
  cache = cache,
  control = list(
    maxit = 1000,
    report.level = 4,
    report.freq =1 
  ),
  control_inner = list(
    maxit = 100,
    report.level = 0,
    report.freq = 0
  )
)


outer_fit$solution
outer_fit$fval
```



``` r
adFun = adlaplace::getAdFun(data,config)
```

```
## outer, groups 245 Nbeta 2 Ntheta 3 Ngamma 35 Nparams 40
## q ngamma  35 map.ncol 3 map@i.size 35 ntheta 3 exp theta map 0 0.367879 gamma0 0.
## logDensRandom 2.83715 qpart 0 det -2.83715
## extra .. done.
## 12..3456
```

``` r
x = c(config$beta,config$theta)
Npar=13;Dpar = length(x)

parDf = data.frame(x)[,rep(1, Npar)]
Sx = seq(-0.1, 0.1, len=13) + parDf[Dpar,1]
SxD = Sx[-1] - diff(Sx)/2
parDf[Dpar,] = Sx

res = mapply(
  adlaplace::logLikLaplace,
  x = as.list(parDf),
  MoreArgs = list(adFun = adFun, config = config, deriv=TRUE),
  SIMPLIFY=FALSE
  )
```

```
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2177.986915     3.293335                 Continuing    1.000000     10          Reached tolerance
##   2   2177.969061     0.002895                 Continuing    1.000000      9          Reached tolerance
##   3   2177.969061     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2177.969061     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2179.911147     3.093387                 Continuing    1.000000     10          Reached tolerance
##   2   2179.895105     0.002538                 Continuing    1.000000      9          Reached tolerance
##   3   2179.895105     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2179.895105     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2181.919320     2.896781                 Continuing    1.000000     10          Reached tolerance
##   2   2181.904981     0.002214                 Continuing    1.000000      9          Reached tolerance
##   3   2181.904981     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2181.904981     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2184.013888     2.703696                 Continuing    1.000000     10          Reached tolerance
##   2   2184.001141     0.001921                 Continuing    1.000000      9          Reached tolerance
##   3   2184.001141     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2184.001141     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2186.197317     2.514319                 Continuing    1.000000     10          Reached tolerance
##   2   2186.186055     0.001657                 Continuing    1.000000      9          Reached tolerance
##   3   2186.186055     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2186.186055     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f         nrm_gr                     status         radCG iter                  CG result
##   1   2188.472089     2.328844                 Continuing    1.000000     10          Reached tolerance
##   2   2188.462204     0.001420                 Continuing    1.000000      9          Reached tolerance
##   3   2188.462204     0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2188.462204     0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2190.840694    2.147475                 Continuing    1.000000     10          Reached tolerance
##   2   2190.832081    0.001210                 Continuing    1.000000      9          Reached tolerance
##   3   2190.832081    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2190.832081    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2193.305634    1.970437                 Continuing    1.000000     10          Reached tolerance
##   2   2193.298188    0.001024                 Continuing    1.000000      9          Reached tolerance
##   3   2193.298188    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2193.298188    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2195.869413    1.797982                 Continuing    1.000000     10          Reached tolerance
##   2   2195.863033    0.000860                 Continuing    1.000000      9          Reached tolerance
##   3   2195.863033    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2195.863033    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2198.534536    1.628729                 Continuing    1.000000      9          Reached tolerance
##   2   2198.529126    0.000718                 Continuing    1.000000      9          Reached tolerance
##   3   2198.529126    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2198.529126    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2201.303525    1.466415                 Continuing    1.000000      9          Reached tolerance
##   2   2201.298981    0.000594                 Continuing    1.000000      9          Reached tolerance
##   3   2201.298981    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2201.298981    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2204.178885    1.309800                 Continuing    1.000000      9          Reached tolerance
##   2   2204.175111    0.000486                 Continuing    1.000000      9          Reached tolerance
##   3   2204.175111    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2204.175111    0.000000                    Success
## 
## done inner opt
## logLikLaplace using package  adlaplace for objective funcion
## Beginning optimization
## 
## iter           f        nrm_gr                     status         radCG iter                  CG result
##   1   2207.163120    1.159498                 Continuing    1.000000      9          Reached tolerance
##   2   2207.160023    0.000394                 Continuing    1.000000      9          Reached tolerance
##   3   2207.160023    0.000000                 Continuing    1.000000      9          Reached tolerance
## 
## Iteration has terminated
##   3   2207.160023    0.000000                    Success
## 
## done inner opt
```

``` r
Slik = unlist(lapply(res, "[[", 'logLik'))
Sdet = unlist(lapply(res, function(xx) xx$opt$halfLogDet))
dU = do.call(abind::abind, c(lapply(res, function(xx) as.matrix(xx$extra$dU)), along=3))
uHat = do.call(rbind, lapply(res, function(xx) xx$opt$solution))
extraDf = do.call(abind::abind, c(lapply(res, '[[', "deriv"), along=3))
gradMat = do.call(rbind, lapply(res, '[[', 'grad'))
gradVec = unlist(lapply(res, function(xx) xx$deriv[Dpar,'dL']))


plot(Sx, Slik)
```

![plot of chunk testLogLgrad](figure/testLogLgrad-1.png)

``` r
plot(Sx, gradVec)
lines(SxD, diff(Slik)/diff(Sx), lwd=2)
```

![plot of chunk testLogLgrad](figure/testLogLgrad-2.png)

``` r
Du = 1
plot(Sx, uHat[,Du])
```

![plot of chunk testLogLgrad](figure/testLogLgrad-3.png)

``` r
plot(Sx, dU[Du,Dpar,])
lines(SxD, diff(uHat[,Du])/diff(Sx))
```

![plot of chunk testLogLgrad](figure/testLogLgrad-4.png)

``` r
plot(Sx, Sdet)
```

![plot of chunk testLogLgrad](figure/testLogLgrad-5.png)

``` r
plot(Sx, extraDf[Dpar,'dDet',])
```

![plot of chunk testLogLgrad](figure/testLogLgrad-6.png)

``` r
lines(SxD, diff(Sdet)/diff(Sx))
```

```
## Error in `xy.coords()`:
## ! 'x' and 'y' lengths differ
```






``` r
adFun = adlaplace::getAdFun(data, modifyList(config, list(verbose=TRUE)))
```

```
## outer, groups 245 Nbeta 2 Ntheta 3 Ngamma 35 Nparams 40
## q ngamma  35 map.ncol 3 map@i.size 35 ntheta 3 exp theta map 0 0.367879 gamma0 0.
## logDensRandom 2.83715 qpart 0 det -2.83715
## extra .. done.
## 12..3456
```

``` r
x = c(config$beta, rep(0.1, length(config$gamma)), config$theta)
adlaplace::jointLogDens(x, adFun)
```

```
## [1] -2274.308
```

``` r
str(adlaplace::grad(x, adFun, FALSE))
```

```
##  num [1:40] -558.5 -337.6 -40 -88.7 -67.4 ...
```

``` r
str(adlaplace::grad(x, adFun, TRUE))
```

```
##  num [1:40] 0 0 -40 -88.7 -67.4 ...
```

``` r
Dgroup = ncol(config$groups)+1
str(h1 <- adlaplace::hess(x, adFun, FALSE, Sgroups = Dgroup))
```

```
## Formal class 'dsCMatrix' [package "Matrix"] with 7 slots
##   ..@ i       : int [1:814] 0 0 1 0 1 2 0 1 2 3 ...
##   ..@ p       : int [1:41] 0 1 3 6 10 15 21 28 36 45 ...
##   ..@ Dim     : int [1:2] 40 40
##   ..@ Dimnames:List of 2
##   .. ..$ : NULL
##   .. ..$ : NULL
##   ..@ x       : num [1:814] 0 0 0 0 0 0 0 0 0 0 ...
##   ..@ uplo    : chr "U"
##   ..@ factors : list()
```

``` r
adFun$sparsity[1+Dgroup]
```

```
## [[1]]
## [[1]]$grad
## [1] 39
## 
## [[1]]$grad_inner
## integer(0)
## 
## [[1]]$row_hess
## [1] 39
## 
## [[1]]$col_hess
## [1] 39
## 
## [[1]]$row_hess_inner
## integer(0)
## 
## [[1]]$col_hess_inner
## integer(0)
```

``` r
(Shere = 1+seq(
	adFun$hessians$map$outer$p[Dgroup+1],
	adFun$hessians$map$outer$p[Dgroup+2]-1))
```

```
## [1] 4179
```

``` r
adFun$hessians$map$outer$local[Shere]
```

```
## [1] 0
```

``` r
adFun$hessians$map$outer$global[Shere]
```

```
## [1] 813
```

``` r
(hseq = which(h1@x != 0))
```

```
## [1] 814
```

``` r
h1@x[hseq]
```

```
## [1] 121340.2
```

``` r
str(h2 <- adlaplace::hess(x, adFun, TRUE))
```

```
## Formal class 'dsCMatrix' [package "Matrix"] with 7 slots
##   ..@ i       : int [1:630] 0 0 1 0 1 2 0 1 2 3 ...
##   ..@ p       : int [1:36] 0 1 3 6 10 15 21 28 36 45 ...
##   ..@ Dim     : int [1:2] 35 35
##   ..@ Dimnames:List of 2
##   .. ..$ : NULL
##   .. ..$ : NULL
##   ..@ x       : num [1:630] -224 0 -250 0 0 ...
##   ..@ uplo    : chr "U"
##   ..@ factors : list()
```


``` r
adFun = adlaplace::getAdFun(data, config)
inner_res <- adlaplace::inner_opt(
  parameters = c(config$beta, config$theta),
  gamma = config$gamma,
  control = list(
    maxit = 50,
    report.level = 0,
    report.freq = 0),
  config = config,
  adFun = adFun)

str(inner_res)

quantile(inner_res$gradient)
```





``` r
config$gamma = rep(1, length(config$gamma))
x = c(config$beta, config$gamma, config$theta)
adFun = adlaplace::getAdFun(data, modifyList(config, list(verbose=TRUE)))

inner = FALSE
type = c('outer','inner')[1+inner]


Sgamma1 = seq.int(length(config$beta)+1, length.out = length(config$gamma))


Npar=13;Dpar = 3
if(inner) {
	Dpar = match(Dpar, Sgamma1)
}
Dpar0=Dpar-1L

bob = as(adFun$hessians$hessian[[type]], 'TsparseMatrix')
(whichIndex = which(bob@i == Dpar0 & bob@j == Dpar0)-1L)
(whichGlobal = which(adFun$hessians$map[[type]]$global==whichIndex)-1L)

(whichP = mapply(
	function(xx) min(which(adFun$hessians$map[[type]]$p >= xx)),
	xx = whichGlobal) -1L)



do.call(rbind, adFun$sparsity[[
	whichP[1] 
	]][c('row_hess','col_hess')])
do.call(rbind, adFun$sparsity[[
	whichP[1] 
	]][c('row_hess_inner', 'col_hess_inner')])
	

(Sgroups = whichP-1L)
Sgroups = seq.int(from=0, length.out = length(adFun$sparsity))

#Sgroups = c()


parDf = data.frame(x)[,rep(1, Npar)]
Sx = seq(-0.1, 0.1, len=13) + parDf[Dpar,1]
SxD = Sx[-1] - diff(Sx)/2
parDf[Dpar,] = Sx


dens1 = mapply(
  adlaplace::jointLogDens,
  x = as.list(parDf),
  MoreArgs = list(backendContext = adFun, Sgroups = Sgroups),
  SIMPLIFY=TRUE
  )

plot(Sx, dens1)


grad1 = mapply(
  adlaplace::grad,
  x = as.list(parDf),
  MoreArgs = list(backendContext = adFun, inner=inner, Sgroups = Sgroups),
  SIMPLIFY=TRUE
  )
plot(Sx, grad1[Dpar,])
lines(SxD, diff(dens1)/diff(Sx))


hes1 = mapply(
  adlaplace::hess,
  x = as.list(parDf),
  MoreArgs = list(backendContext = adFun, inner=inner, Sgroups = Sgroups),
  SIMPLIFY=FALSE)

hes2 = do.call(abind::abind, c(
	lapply(hes1, as.matrix), along=3))

gradD = apply(grad1, 1, diff) / mean(diff(Sx))

Dpar2 = dim(hes2)[1]-2
plot(Sx, hes2[Dpar,Dpar2, ], ylab = paste(Dpar, Dpar2),
	ylim = range(hes2[Dpar, Dpar2,], na.rm=TRUE))
lines(SxD, gradD[,Dpar2])
```



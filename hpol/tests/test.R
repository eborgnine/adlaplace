#Rcpp::compileAttributes("hpol")

library('hpolcc')

x = c(10,10,1e6,0)
u = c(1,0,0,0)
w = c(0,0,1)

res1 = hpolcc:::test3(x, u, w)
res1$taylor2
res = matrix(res1$taylor3, ncol=length(w), byrow=TRUE)
colnames(res) = c('one','two','three')
res



f = function(x) sum(x[1:2]^3)

f1 = round(numDeriv::grad(f, x),5)

f2 = round(numDeriv::hessian(f, x) %*% u, 5)


f3f = function(x, u) {
    tt = array(0, rep(length(x), 3))
    for(D in 1:length(x))
        tt[D,D,D] = 6
    result = rep(NA, length(x))
    for(D in 1:length(x))
        result[D] = drop(u %*% tt[D,,]%*% u)
    round(result, 5)
}

f3 = f3f(x, u)

cbind(
    cbind(
f1 + f2 + (f3)/2,
f1 + f2 ,
    f1
), NA, res)


cbind(
# third should be
three=f3(x, u)/2,
# hes should be
hes=round(f2 %*% u,5),
# grad should be
grad=f1) 
# %*% cbind(c(1,1,1), c(0,1,1), c(0,0,1))

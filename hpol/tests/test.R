#Rcpp::compileAttributes("hpol")


# f_j is deriv in j direction, f_jk is 2nd deriv in jk, f_jkl third

#| returned component         | closed form (scalar case)                                                                                           |
#| -------------------------- | ------------------------------------------------------------------------------------------------------------------- |
#| **order 0**<br>`dw[j*3+0]` | $w^{(0)}f_j + w^{(1)}\sum_k f_{jk}\,U_k +\dfrac{w^{(2)}}{2}\Bigl(f_{j k}V_k+ \sum_{k\ell}f_{jk\ell}U_kU_\ell\Bigr)$ |
#| **order 1**<br>`dw[j*3+1]` | $w^{(1)}f_j + w^{(2)}\sum_k f_{jk}\,U_k$                                                                            |
#| **order 2**<br>`dw[j*3+2]` | $\dfrac{w^{(2)}}{2}\,f_j$                                                                                           |


# first column is third deriv
# always specify u and v, and three-digit w
# see hessian3.md in inst

library('hpolcc')

x = c(5,1,-1,0)
u = c(0,1,0,0)
v = c(0,1,0,0)
w = 1#c(0,0,1)

res1 = hpolcc:::test3(x, u,v, w)
#res1$taylor2
res = matrix(res1$taylor3, ncol=3, byrow=TRUE)
colnames(res) = c('one','two','three')
res



f = function(x) sum(x*(seq(1,len=length(x))))^4

f1 = round(numDeriv::grad(f, x),5)


hh = numDeriv::hessian(f, x)
f2u = round(hh %*% u, 5)
f2v = round(hh %*% v, 5)





f3f = function(x, u) {
    eps = 0.01
    tt = array(0, rep(length(x), 3))
    for(D in 1:length(x)) {
        xplus = x
        xplus[D] = x[D] + eps
        hhplus = numDeriv::hessian(f, xplus)
        hdiff = (hhplus -hh)/eps

        tt[D, , ] = hdiff
    }
    result = rep(NA, length(x))
    for(D in 1:length(x))
        result[D] = drop(u %*% tt[D,,]%*% u)
    round(result, 5)
}

f3 = f3f(x, u)

cbind(
    cbind(
 f2v + (f3)/2,
 f2u ,
    f1
), NA, res)



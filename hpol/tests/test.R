#Rcpp::compileAttributes("hpol")

library('hpolcc')

x = c(0,0,0,0)
u = c(0,1,0,0)
w = c(0,0,1)


res = matrix(hpolcc:::test3(x, u, w)$taylor3, ncol=length(w), byrow=TRUE)
colnames(res) = c('grad','hes','three')
res

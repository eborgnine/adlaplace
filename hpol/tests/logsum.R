library(hpolcc)

hpolcc:::logspaceadd_forward_deriv(c(0,0,0)+10, 1L)    


hpolcc:::testAddGammaR(1:3)

numDeriv::grad(function(xx) hpolcc:::testAddGammaR(xx)$value, 1:3)
numDeriv::hessian(function(xx) hpolcc:::testAddGammaR(xx)$value, 1:3)


# A few sample input vectors
inputs <- list(
  c(1,10),
  c(1,1),
  c(0.5, 1, 2),
  c(1, 2, 3, 4),
  c(-2, 0, 1, 2),
  c(-5, 0, 5)
)

for(D in 1:length(inputs)) {
  print(c(
    hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 0L),
    log(sum(exp(inputs[[D]])))
  ))
}


for(D in 1:length(inputs)) {
  print(rbind(
    hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 1L),
    hpolcc:::logspaceadd_forward_deriv(inputs[[D]], 1L)[,1],
    exp(inputs[[D]]) / sum(exp(inputs[[D]]))
  ))
}



for(D in 1:length(inputs)) {
  print(inputs[[D]])
  print(
    round(cbind(
      hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 2L), NA,
      hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], -2L), NA,
      hpolcc:::logspaceadd_forward_deriv(inputs[[D]], 2L)
  ), 4))
}

D=3
(hes =  hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 2L))
(xx =  hpolcc:::logspaceadd_forward_deriv(inputs[[D]], 2L) )
hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], -2L)


hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 3L)

D=3
numDeriv::hessian(
 function(x) log(sum(exp(x))),
 inputs[[D]] 
)
numDeriv::grad(
 function(x) log(sum(exp(x))),
 inputs[[D]] 
)

D = 1
thethird = array(
  hpolcc:::logspaceadd_inbuilt_deriv(inputs[[D]], 3L),
  rep(length(inputs[[D]]), 3))

delta = 0.1
deriv2 = list()
d2 = hpolcc:::logspaceadd_inbuilt_deriv(
    inputs[[D]], 2L)
for(Dderiv in 1:length(inputs[[D]])) {
  deriv2[[Dderiv]] = (hpolcc:::logspaceadd_inbuilt_deriv(
    inputs[[D]] + delta * (1:length(inputs[[D]]) == Dderiv), 
    2L) - d2)/delta
}
inputs[[D]]
do.call(abind::abind, c(deriv2, list(along=length(inputs[[D]])+1)))
thethird

#Grad[log[Sum[exp\(40)Subscript[x,i]\(41),{i,1,4}]] , {Subscript[x,1], Subscript[x,2], Subscript[x,3], Subscript[x,4]}]
#Partial[log\(91)Sum[exp\(40)Subscript[x,i]\(41),{i,1,4}]\(93) ,Subscript[x,1],Subscript[x,2]]


x1 = seq(-1,1,len=101);x2=0
y = log(exp(x1) + exp(x2))
y0 = mapply(
  function(x1, x2) hpolcc:::logspaceadd_inbuilt_deriv(c(x1, x2), 0L),
  x1 = x1, MoreArgs = list(x2=x2))
plot(x1, y0)
lines(x1, y, col='green')

y1 = mapply(
  function(x1, x2) hpolcc:::logspaceadd_inbuilt_deriv(c(x1, x2), 1L)[1],
  x1 = x1, MoreArgs = list(x2=x2))

plot(x1, y1)
lines(x1[-1], diff(y)/diff(x1), col='green')

y2 = mapply(
  function(x1, x2) hpolcc:::logspaceadd_inbuilt_deriv(c(x1, x2), 2L)[1,1],
  x1 = x1, MoreArgs = list(x2=x2))

y2f = mapply(
  function(x1, x2) hpolcc:::logspaceadd_inbuilt_deriv(c(x1, x2), -2L)[1,1],
  x1 = x1, MoreArgs = list(x2=x2))

matplot(x1, cbind(y2, y2f), pch=c(1,2))
lines(x1[-c(1, length(x1))], diff(diff(y))/diff(x1)[-1]^2, col='green')

y3 = mapply(
  function(x1, x2) hpolcc:::logspaceadd_inbuilt_deriv(c(x1, x2), 3L)[1],
  x1 = x1, MoreArgs = list(x2=x2))
plot(x1, y3)
plot(x1[-c(1,2, length(x1))], diff(diff(diff(y)))/mean(diff(x1))^3)







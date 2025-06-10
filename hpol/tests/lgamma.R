

xseq = c(0.5, 1,2)
Sorder = 1:4

xx = expand.grid(x =xseq, order=Sorder)

library('hpolcc')
xx$forward = mapply(
	hpolcc:::lgamma_forward_deriv,
	x=xx$x, order=xx$order
)
xx$reverse = mapply(
	hpolcc:::lgamma_reverse_deriv,
	x=xx$x, order=xx$order
)

xx$R = mapply(
	psigamma,
	x = xx$x, deriv = xx$order-1
)

dd = 2
matplot(
	xx[xx$order == dd, 'x'],
	xx[xx$order == dd,c('forward','reverse','R')],
	ylim = c(-11,11)
)

xx


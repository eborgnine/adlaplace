

firsTwoElements = function(xx, k) sort(c(xx[xx!=k], rep(k,2))[1:2] )

thirdTensor = function(k, third, N) {
	thirdHere = third[apply(third[,c('i','j','k')] == k, 1, any), ]
	if(!nrow(thirdHere)) return(NULL)
		newxy = t(apply(thirdHere[,c('i','j','k')], 1,  firsTwoElements, k=k))
	try(Matrix::sparseMatrix(i=newxy[,1], j=newxy[,2], x=thirdHere[,'x'],
		dims = rep(N, 2), symmetric=TRUE, index1=FALSE))
}

sumTrace = function(Hinv, Tuux, Sgamma1) {
#      sum(Matrix::diag(Hinv %*% Tuux[Sgamma1,Sgamma1]))
	sum(as(Hinv * Tuux[Sgamma1,Sgamma1], "generalMatrix")@x)  
}    


lengthUnique = function(xx) length(unique(xx))


setdiffSingle = function(x, y) {
	result = setdiff(x, y)
	if(!length(result)) {
		tableX = sort(table(x), decreasing=TRUE)
		result = type.convert(names(tableX), as.is=TRUE)
	}
	result[1]
}

nIn = function(ref, check, n) (sum(check %in% ref)>=n) & all(ref %in% check)

Tijdot = function(x, pair) {
	pair = unlist(pair)
  	# get all elements of the tensor data frame that contain both elements in y
	x = as.data.frame(x[apply(as.matrix(x[,c('i','j','k')]), 1, nIn, ref = pair, n=2), ])
	newk = apply(as.matrix(x[,c('i','j','k')]), 1, setdiffSingle, y=pair)
	cbind(x[, setdiff(names(x), c('i','j','k')), drop=FALSE], k=newk)
}

getTijdotDu = function(pair, third, Sgamma1, dUhat, Nparameters) {
	thirdHere = Tijdot(third, pair)
	thirdHere = Matrix::sparseVector(thirdHere$x, thirdHere$k+1, length = Nparameters)
	as.vector(crossprod(dUhatDtheta, thirdHere[Sgamma1]) + thirdHere[-Sgamma1])
}


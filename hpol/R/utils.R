rowSortP = function(x) {

if(requireNamespace("Rfast", quietly=TRUE)) {      
      return(Rfast::rowSort(x))
    } else {
      return(t(apply(x, 1, sort)))
    }
    x
}

lengthUnique = function(xx) {
	length(unique(xx))
}


firstTwoElements = function(xx, k) sort(c(xx[xx!=k], rep(k,2))[1:2] )

thirdTensorSparse = function(third, sparsity) {
    ijkHere = as.matrix(sparsity$third$ijk[,c('i','j','k')])
    ijkHere = rowSortP(ijkHere)
    colnames(ijkHere) = c('i','j','k')
    ijkHere = as.data.frame(ijkHere)
    ijkHere$Tijk = third
    ijkHere
  }

thirdTensorDense = function(group_sparsity, Tijk, Nparameters)
      {


        Npairs = nrow(group_sparsity$third$pairs)
        matHere = as(Matrix::Matrix(Tijk, ncol=Npairs),'TsparseMatrix')
        matHerePairK = cbind(k = matHere@i, pair=matHere@j)

        toApply = cbind(
        	as.matrix(as.data.frame(
        		group_sparsity$third$pairs[c('i','j')]
        	))[1+matHerePairK[,'pair'],],
        	matHerePairK[,'k', drop=FALSE]
        )

        resultIJK = rowSortP(toApply)


        colnames(resultIJK) = c('i','j','k')
        duplicatedHere = duplicated(resultIJK)
        Nunique = apply(resultIJK, 1, lengthUnique)#Rfast::sort_unique.length)

        toKeep = (!duplicatedHere) & (Nunique == 3)

				resultIJK = as.data.frame(resultIJK[toKeep,])          
        resultIJK$Tijk = matHere@x[toKeep]
        resultIJK
      }

thirdTensor = function(k, third, N) {
	thirdHere = third[apply(third[,c('i','j','k')] == k, 1, any), ]
	if(!nrow(thirdHere)) return(NULL)
		newxy = t(apply(thirdHere[,c('i','j','k')], 1,  firstTwoElements, k=k))
	notDup = !duplicated(newxy)
	try(Matrix::sparseMatrix(i=newxy[notDup,1], j=newxy[notDup,2], x=thirdHere[notDup,'x'],
		dims = rep(N, 2), symmetric=TRUE, index1=FALSE))
}

sumTrace = function(Hinv, Tuux, Sgamma1) {
#      sum(Matrix::diag(Hinv %*% Tuux[Sgamma1,Sgamma1]))
	sum(as(Hinv * Tuux[Sgamma1,Sgamma1], "generalMatrix")@x)  
}



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

getDh = function(pair, third, Sgamma1, dUhat, Nparameters) {
	thirdHere = Tijdot(third, pair)
	thirdHere = thirdHere[order(thirdHere$k), ]
	thirdHere = thirdHere[!duplicated(thirdHere$k), ]
	thirdHere = Matrix::sparseVector(thirdHere$x, thirdHere$k+1, length = Nparameters)
	as.vector(crossprod(dUhat, thirdHere[Sgamma1]) + thirdHere[-Sgamma1])
}

forDhList = function(Tijk, dUp, Sgamma1) {
      if(is.null(Tijk)) return(matrix(nrow=0, ncol=2+length(dUp)))
      There =try( as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix'))
      cbind(j=There@i, i=There@j, outer(There@x, dUp))
    }

forTijpAdd =    function(Tijk, Sgamma1) {
  if(is.null(Tijk)) return(matrix(nrow=0, ncol=3))
  There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
  cbind(i=There@i, j=There@j, x=There@x)
  }


forDh = function(TU, ij, Tijp, dims) {
    toAgg = rbind(cbind(ij, x=TU), Tijp)
    theAgg = aggregate(toAgg[,'x', drop=FALSE], toAgg[,c('i','j')], sum, na.rm=TRUE)
    Matrix::sparseMatrix(
      i=pmax(theAgg$j, theAgg$i), 
      j=pmin(theAgg$i, theAgg$j), x=theAgg$x, 
      index1=FALSE, dims=dims, symmetric=TRUE)
  }


traceProd = function(Dhp, Hinv) {
# sum(Matrix::diag(Hinv %*% Dhp))
  sum((Hinv * Dhp))
}
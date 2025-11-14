rowSortP = function(x) {

  if(requireNamespace("Rfast", quietly=TRUE)) {      
    return(Rfast::rowSort(as.matrix(x)))
  } else {
    return(t(apply(x, 1, sort)))
  }
  x
}

lengthUnique = function(xx) {
	length(unique(xx))
}

sparsityByBlockFirst = function(strata, config, parameters, data, Sgamma0) {
  configHere = modifyList(config, list(groups = strata, dense=TRUE, num_threads=1, verbose=FALSE))
  adFunHere = getAdFun(parameters, data, configHere)
  whichAll = which(abs(grad(parameters, data, config=configHere, adFun=adFunHere)) > .Machine$double.eps)-1L
  whichRandom = as.vector(na.omit(match(whichAll, Sgamma0))-1L)
  list(full = whichAll, random=whichRandom)
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

indexThirdTensor = function(i, ijk) {

  idxI <- which(ijk[,'i'] == i) -1L
  idxJ <- which(ijk[,'j'] == i) -1L
  idxK <- which(ijk[,'k'] == i) -1L

  out = data.frame(
    index=c(idxI, idxJ, idxK),
    isI = rep(c(TRUE, FALSE), c(length(idxI), length(idxJ) + length(idxK))),
    isJ = rep(c(FALSE, TRUE, FALSE), c(length(idxI), length(idxJ), length(idxK)))
  )
  if(any(duplicated(out$index))) {
    warning("duplicated i j k")
  }

  out
}

indexThirdTensorDiag = function(i, nonSymmetric) {
  idx <- which(nonSymmetric$i == i | nonSymmetric$j == i)
  isTheSingle = nonSymmetric$i[idx] == i
  data.frame(index = idx - 1L, isTheSingle = isTheSingle)
}

thirdTensorFromIndex = function(
  Tijk, Tiij, thirdIndexOffDiag, thirdIndexDiag, N
) {

  TijkHere = Tijk[thirdIndexOffDiag$index+1L,,drop=FALSE]

  TijkHere[thirdIndexOffDiag$isI,'i'] = TijkHere[thirdIndexOffDiag$isI,'j']
  TijkHere[thirdIndexOffDiag$isI,'j'] = TijkHere[thirdIndexOffDiag$isI,'k']
  TijkHere[thirdIndexOffDiag$isJ,'j'] = TijkHere[thirdIndexOffDiag$isJ,'k']

  # i=nsJ, j=nsJ, k=nsI
  TiijHere = Tiij[thirdIndexDiag$index+1L,c('i','k','Tijk')] # is is the double (j==i)  
  colnames(TiijHere) = gsub("^k$", "j", colnames(TiijHere))
  # if it matches i, keep i,k
  # otherwise matches k, set k=i
  TiijHere[thirdIndexDiag$isTheSingle, 'j'] = TiijHere[thirdIndexDiag$isTheSingle, 'i']
  toSwitch = which(TiijHere$i > TiijHere$j)
  if(length(toSwitch)) {
    newI = TiijHere[toSwitch,'i'] 
    TiijHere[toSwitch,'i'] = TiijHere[toSwitch,'j']
    TiijHere[toSwitch,'j'] = newI
  }
  res = rbind(TiijHere, TijkHere[,colnames(TiijHere)])
  try(Matrix::sparseMatrix(i=res$i, j=res$j, x=res$Tijk, index1=FALSE, dims=c(N,N), symmetric=TRUE))
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
  if(is.null(Tijk)) {
    return(matrix(nrow=0, ncol=2+length(dUp)))
  }
  There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
  return(cbind(data.frame(j=There@i, i=There@j), outer(There@x, dUp)))
}

forTijpAdd =    function(Tijk, Sgamma1) {
  if(is.null(Tijk)) return(matrix(nrow=0, ncol=3))
    There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
  cbind(i=There@i, j=There@j, x=There@x)
}



forDh = function(Dpar, dHagg, thirdList, dims, Sgamma1) {

  DparC = paste0('p', Dpar)
  DparP1 = Dpar + 1L
  There = as(thirdList[[DparP1]][Sgamma1,Sgamma1, drop=FALSE], 'TsparseMatrix')
  Tijp = data.table::data.table(i=There@i, j=There@j, x=There@x)

  baseDT <- dHagg[, .(i, j, x=get(DparC))]
  allDT = data.table::rbindlist(list(baseDT, Tijp), use.names = TRUE)
  allDT[, `:=`(ii = pmax(i, j), jj = pmin(i, j))]
  aggDT <- allDT[, .(x = sum(x, na.rm = TRUE)), by = .(i = ii, j = jj)]

  Matrix::sparseMatrix(i = aggDT$i, j = aggDT$j, x = aggDT$x,
   dims = dims, index1 = FALSE, symmetric = TRUE)
}

traceProd = function(Dhp, Hinv) {
# sum(Matrix::diag(Hinv %*% Dhp))
  sum((Hinv * Dhp))
}
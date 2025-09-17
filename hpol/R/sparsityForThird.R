#' @export
sparsityForThird = function(x, data, config=list()) {


  if(is.vector(x)) {
    # x is a vector of parameters
    if(any(config$verbose > 0)) {
      message("input is parameters, computing sparse hessian")
    }
    configForSparse = config[setdiff(names(config), "sparsity")]
    configForSparse$maxDeriv = 2
    hessian1 = objectiveFunctionC(
      parameters = x, data=data, 
      config=configForSparse)
    hessian = hessian1$hessian
  } else {
    # assume it's the hessian
    hessian = x
    if(is.matrix(x)) hessian = Matrix::Matrix(x)
  }

if(is.data.frame(hessian)) {
  hessian = Matrix::sparseMatrix(
    i=hessian$i, j=hessian$j, 
    index1=FALSE, symmetric=TRUE)
} 


Nbeta = nrow(data$XTp)
Ngamma = nrow(data$ATp)
Sgamma = seq(from=Nbeta+1, len=Ngamma)
Sgamma0 = Sgamma -1



if(nrow(hessian) <= Ngamma ) 
  warning("hessian must be from the full set of parameters and random effects")

Ntotal = nrow(hessian)
Sparams = setdiff(1:Ntotal, Sgamma)
Sparams0 = Sparams-1
Nparams = length(Sparams)

hessianT = as(hessian,'dsTMatrix')
hessianC = as(hessian, 'CsparseMatrix')

  # store upper and lower triangle
hessianC2 = as(hessianC, "generalMatrix")
hessianT2 = as(hessianT, "generalMatrix")

hessianTdf = data.frame(
  i=hessianT2@i, 
  j=hessianT2@j)

#  hessianS = as(Matrix::forceSymmetric(hessian), "TsparseMatrix")


hessianRandom = hessianT[Sgamma,Sgamma]

parGammaDiag = hessianC
#parGammaDiag[Sparams,] = 0

parGamma = list(
  i=parGammaDiag@i,
  j = rep(seq(0, len=nrow(hessian)), diff(parGammaDiag@p)), 
  p=parGammaDiag@p,
  Sparams = Sparams0)

pastedParGamma = paste(parGamma$i, parGamma$j, sep='_')
indexHessian1 = match(pastedParGamma,
  paste(hessianTdf$i, hessianTdf$j, sep='_'))
indexHessian2 = match(pastedParGamma,
  paste(hessianTdf$j, hessianTdf$j, sep='_'))
parGamma$indexHessian = pmax(indexHessian1, indexHessian2, na.rm=TRUE)


# entries of third deriv which are non-zero
# for each T_ijk, H_ik H_jk H_ij nonzero
# for each ij pair, find H[,j] > 0 & H[,j] >0

forTriples = cbind(hessianT@i, hessianT@j)
forTriples = forTriples[forTriples[,1] !=  forTriples[,2], ]
# dont want both as parameters


forTriples1= forTriples+1

theK = apply(forTriples1, 1, function(xx, hessianC) {
  possibleK = which( (hessianC[,xx[1]] != 0) & (hessianC[,xx[2]] != 0))
  possibleK [ (possibleK < min(xx))]
}, hessianC = hessianC)
NperRow = unlist(lapply(theK, length))

Snk = rep(1:nrow(forTriples), NperRow)

ijk1 = cbind(
  forTriples1[Snk, ],
  unlist(theK)
) 
colnames(ijk1) = c('i','j','k')

ijk = ijk1-1

ijk = ijk[order(ijk[,'i'],ijk[,'j'],ijk[,'k']),]


if(FALSE) {
  table(duplicated(ijk[,c('i','j')]))
  table(duplicated(ijk[,c('i','k')]))
  table(duplicated(ijk[,c('j','k')]))
}

indexKii1 = match(
  apply(ijk[,c('i','k')], 1, paste, collapse='_'),
  paste(parGamma$i, parGamma$j, sep='_')
)
indexKii2 = match(
  apply(ijk[,c('i','k')], 1, paste, collapse='_'),
  paste(parGamma$j, parGamma$i, sep='_')
)
indexKii = pmin(indexKii1, indexKii2, na.rm=TRUE)

if(FALSE) {
  cbind(ijk, indexKii)[seq(-2,2) + min(which(is.na(indexKii))), ]
}

indexKjj1 = match(
  apply(ijk[,c('j','k')], 1, paste, collapse='_'),
  paste(parGamma$i, parGamma$j, sep='_')
)
indexKjj2 = match(
  apply(ijk[,c('j','k')], 1, paste, collapse='_'),
  paste(parGamma$j, parGamma$i, sep='_')
)
indexKjj = pmin(indexKjj1, indexKjj2, na.rm=TRUE)

if(FALSE) {
  cbind(ijk, indexKjj)[1:12,][seq(-2,2) + min(which(is.na(indexKjj))), ]
}


ijkp = cbind(
  ijk,
  p=seq(from=0, len=nrow(ijk)),
  indexKii = indexKii-1,
  indexKjj = indexKjj-1 
)

# apply(ijkp, 2, function(xx) sum(is.na(xx)))

ij = ijkp[!duplicated(ijkp[,c('i','j')]), c('i','j','p')]

#jk = ijkp[!duplicated(ijkp[,c('j','k')]), setdiff(colnames(ijkp), 'i')]


pEnd = c(ij[-1,'p'], nrow(ijk))
ij=cbind(
  ij[,c('i','j','p')],
  pEnd = pEnd,
  n = pEnd - ij[,'p'])

#quantile(ij[,'n'])  





storage.mode(ij) = 'integer'
storage.mode(ijkp) = 'integer'




sparsity = list(
  second = list(
    full = list(i=hessianT@i, j=hessianT@j, 
      p=hessianC@p),
    random = list(i=hessianRandom@i, j=hessianRandom@j, 
      p=as(hessianRandom, "CsparseMatrix")@p),
    parGamma = parGamma
  ),
  third = list(
    ijk = as.data.frame(ijkp),
    ij = as.data.frame(ij)
  )
)


return(sparsity)
}



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

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Sgamma = seq(from=Nbeta+1, len=Ngamma)


  if(is.data.frame(hessian)) {
    hessian = Matrix::sparseMatrix(
      i=hessian$i, j=hessian$j, 
      index1=FALSE, symmetric=TRUE)
  } 

  if(nrow(hessian) <= Ngamma ) 
    warning("hessian must be from the full set of parameters and random effects")

  Sparams = setdiff(1:nrow(hessian), Sgamma)
  Nparams = length(Sparams)

  if(all(class(hessian) == 'dsCMatrix')) {
    hessian = as(hessian,'dsTMatrix')
  }

  hessianTdf = data.frame(
      i=c(hessian@i,hessian@j), 
      j=c(hessian@j, hessian@i))
  hessianTdf = hessianTdf[!duplicated(hessianTdf), ]

# entries of third deriv which are non-zero
# each pair musth have non-zero 2nd deriv
  parametersGamma1 = hessianTdf[hessianTdf$j %in% (Sparams-1), ]
  parametersGamma2 = mapply(
    function(ifrom1, jfrom1, hessian) {
      res1 = apply(hessian[,c(ifrom1,jfrom1)]!=0, 2, which, simplify=FALSE)
      cbind(i=ifrom1, j=jfrom1, k=do.call(intersect, res1))-1
    },
    ifrom1 = hessianTdf$i+1, jfrom1 = hessianTdf$j+1, 
    MoreArgs = list(hessian=hessian)
  )
  parametersGamma3 = do.call(rbind, parametersGamma2)
  parametersGamma4 = t(apply(parametersGamma3, 1, sort, decreasing=TRUE))
  parametersGamma5 = parametersGamma4[!duplicated(parametersGamma4), ]


# need all three indices to be different
  parametersGamma = parametersGamma5[
  apply(parametersGamma5, 1, function(xx) length(unique(xx)))==3,]

  colnames(parametersGamma) = c('i','j','k')

  hessianS = Matrix::forceSymmetric(hessian)
  hessianRandom = hessianS[Sgamma,Sgamma]

  ijk = cbind(p=seq(from=0, len=nrow(parametersGamma)), parametersGamma)

  jk = ijk[!duplicated(ijk[,c('j','k')]), c('p','j','k') ]
  pEnd = c(jk[-1,'p'], nrow(jk)+1)
  jk=cbind(
    pEnd = pEnd,
    n = pEnd - jk[,'p'],
    jk)

  storage.mode(jk) = 'integer'
  storage.mode(ijk) = 'integer'


  sparsity = list(
    second = list(
      full = list(i=hessianS@i, j=hessianS@j, 
        p=as(hessianS, "CsparseMatrix")@p),
      random = list(i=hessianRandom@i, j=hessianRandom@j, 
        p=as(hessianRandom, "CsparseMatrix")@p)
    ),
    third = list(
      ijk = as.data.frame(ijk),
      jk = as.data.frame(jk)
    )
  )


  return(sparsity)
}



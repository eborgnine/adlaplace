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
    hessianC = hessian
    hessian = as(hessian,'dsTMatrix')
  } else {
    hessianC = as(hessian, 'CsparseMatrix')
  }
  # store upper and lower triangle
  hessianC2 = as(hessianC, "generalMatrix")
  hessianT2 = as(hessianC2, "TsparseMatrix")

  hessianTdf = data.frame(
      i=hessianT2@i, 
      j=hessianT2@j)

# entries of third deriv which are non-zero
# each pair musth have non-zero 2nd deriv
# i, j are gammas, k parameters  
  parametersGamma1 = hessianTdf[
#  hessianTdf$j %in% (Sparams-1) &
    hessianTdf$i %in% (Sgamma-1), ]
  parametersGamma2 = mapply(
    function(ifrom1, jfrom1, hessian) {
      res1 = apply(hessian[, c(ifrom1,jfrom1)]!=0, 2, 
        which, simplify=FALSE)
      inboth = do.call(intersect, res1)
      newj = intersect(inboth, Sgamma[Sgamma < ifrom1])
      cbind(i=rep(ifrom1,length(newj)), j=newj, 
        k=rep(jfrom1,length(newj)))-1
    },
    ifrom1 = parametersGamma1$i+1, jfrom1 = parametersGamma1$j+1, 
    MoreArgs = list(hessian=hessian)
  )
  parametersGamma = do.call(rbind, parametersGamma2)

  hessianS = Matrix::forceSymmetric(hessian)
  hessianRandom = hessianS[Sgamma,Sgamma]

  ijk = cbind(p=seq(from=0, len=nrow(parametersGamma)), 
    parametersGamma)

  #jk = ijk[!duplicated(ijk[,c('j','k')]), c('p','i','j','k') ]
  ij = ijk[!duplicated(ijk[,c('i','j')]), c('p','i','j','k') ]

if(FALSE) {
  table(duplicated(ijk[,c('i','j')]))
  table(duplicated(ijk[,c('i','k')]))
  table(duplicated(ijk[,c('j','k')]))
}

  pEnd = c(ij[-1,'p'], nrow(ij)+1)
  ij=cbind(
    pEnd = pEnd,
    n = pEnd - ij[,'p'],
    ij)

  parGammaDiag = hessianC#[,Sparams]
  parGammaDiag[Sparams,] = 0

  storage.mode(ij) = 'integer'
  storage.mode(ijk) = 'integer'


  sparsity = list(
    second = list(
      full = list(i=hessianS@i, j=hessianS@j, 
        p=as(hessianS, "CsparseMatrix")@p),
      random = list(i=hessianRandom@i, j=hessianRandom@j, 
        p=as(hessianRandom, "CsparseMatrix")@p),
      parGamma = list(
        i=parGammaDiag@i,
        j = rep(seq(0, len=nrow(hessian)), diff(parGammaDiag@p)), 
        p=parGammaDiag@p,
        Sparams = Sparams-1)
    ),
    third = list(
      ijk = as.data.frame(ijk),
      ij = as.data.frame(ij)
    )
  )


  return(sparsity)
}



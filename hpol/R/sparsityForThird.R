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

hessianC = as(hessian, 'CsparseMatrix')
hessianT = as(hessianC,'TsparseMatrix')

  # store upper and lower triangle
hessianC2 = as(hessianC, "generalMatrix")
hessianT2 = as(hessianC2, "TsparseMatrix")

hessianTdf = data.frame(
  i=hessianT2@i, 
  j=hessianT2@j)

#  hessianS = as(Matrix::forceSymmetric(hessian), "TsparseMatrix")


hessianRandom = hessianT[Sgamma,Sgamma]





# entries of third deriv which are non-zero
# for each T_ijk, H_ik H_jk H_ij nonzero
# for each ij pair, find H[,j] > 0 & H[,j] >0

forTriples = cbind(hessianT@i, hessianT@j)
forTriples = forTriples[forTriples[,1] !=  forTriples[,2], ]

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

# for debugging, all hessian entries
ijk1 = cbind(
  forTriples[rep(1:nrow(forTriples1), nrow(hessianT)), ],
  rep(1:nrow(hessianT), each=nrow(forTriples1))
)
colnames(ijk1) = c('i','j','k')

# get rid of tripples where all are parameters
NinParams = apply(ijk1[,c('i','j','k')],1, function(xx) sum(xx %in% (Sparams-1)))
ijk1 = ijk1[NinParams < 3, ]

# get rid of tripples there is a pair
Nunique = apply(ijk1[,c('i','j','k')],1, function(xx) length(unique(xx)))
ijk1 = ijk1[Nunique ==  3, ]


ijk = ijk1-1

ijk = ijk[order(ijk[,'i'],ijk[,'j'],ijk[,'k']),]


if(FALSE) {
  table(duplicated(ijk[,c('i','j')]))
  table(duplicated(ijk[,c('i','k')]))
  table(duplicated(ijk[,c('j','k')]))
}


nonsymmetric = list(
  i=hessianT2@i,
  j=hessianT2@j, 
  p=hessianC2@p
)

nsJi = paste(nonsymmetric$j, nonsymmetric$i, sep="_")

ijkp = cbind(
  ijk,
  p=seq(from=0, len=nrow(ijk)),
  diagIk1 = match(apply(ijk[,c('i','k')], 1, paste, collapse='_'),nsJi),
  diagJk1 = match(apply(ijk[,c('j','k')], 1, paste, collapse='_'),nsJi)  
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
    nonSymmetric = nonsymmetric
  ),
  third = list(
    ijk = as.data.frame(ijkp),
    ij = as.data.frame(ij)
  )
)


return(sparsity)
}



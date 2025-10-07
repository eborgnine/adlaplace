#' @export
sparsity_pattern = function(x, data, config=list()) {

#library('hpolcc')

  if(identical(config$dense, "TRUE")) {
    type = 'compute'
  } else {
    type = 'hessian'
  }

  configForDiag = config[setdiff(names(config), c("dense","sparsity"))]
  configForDiag$dense=TRUE
  resThirdDiag = thirdDiagonals(
    x, data, configForDiag
  ) 

  
  hessian = as(
    Matrix::forceSymmetric(Matrix::Matrix(resThirdDiag$second, sparse=TRUE)),
    'TsparseMatrix')  
  hessianUL = as(hessian, 'generalMatrix')

  hessianIJ = data.frame(i = hessianUL@i, j=hessianUL@j)
  hessianIJ = hessianIJ[order(hessianIJ$j, hessianIJ$i), ]

    pairs = hessianIJ[!duplicated(hessianIJ[,c('i','j')]), ]
    pairs = pairs[order(pairs$j, pairs$i), ]
    pairs = as.matrix(pairs)


  if(type[1] == 'hessian') {
    # find sparsity pattern based on the hessian

    Sk = apply(hessianIJ, 1, function(xx, ref) {
      intersect(ref[ref$i == xx['i'], 'j'], ref[ref$j == xx['j'], 'i'])
    }, ref = hessianIJ)
    ijk1 = hessianIJ[rep(1:nrow(hessianIJ), unlist(lapply(Sk, length))), ]
    ijk1$k = unlist(Sk)

  } else {

    # compute and check all entries

    resThirdSparse = thirdNonDiagonalsSparsity(
      x, data, config, pairs
    )
    ijk2 = as(resThirdSparse, "TsparseMatrix")
    ijk1 = data.frame(
      i=pairs[ijk2@i+1, 'i'], 
      j=pairs[ijk2@i+1, 'j'], 
      k=ijk2@j)
  } # else dense

ijk1 = t(apply(ijk1, 1, sort))
ijk1 = ijk1[!duplicated(ijk1), ]

Nunique = apply(ijk1,1, function(xx) length(unique(xx)))
ijk = as.data.frame(ijk1[Nunique == 3, ])
colnames(ijk) = c('i','j','k')
ijk$ij = apply(ijk[,c('i','j')], 1, paste, collapse='_')
ijk$jk = apply(ijk[,c('j','k')], 1, paste, collapse='_')
ijk$ik = apply(ijk[,c('i','k')], 1, paste, collapse='_')

Sdim = c('ij','jk','ik')
Spairs = c("pairs1","pairs2","pairs3")
ijk$pair = ijk$type = NA
Nna=1;Niter = 0
while(Nna > 0 & Niter < nrow(ijk)) {
  theNA = which(is.na(ijk$pair))
  theTable = table(unlist(ijk[theNA, Sdim]))
  biggestPair = names(theTable)[which.max(theTable)]
  newIJ =intersect(which(ijk$ij==biggestPair), theNA)
  newIK =intersect(which(ijk$ik==biggestPair), theNA)
  newJK =intersect(which(ijk$jk==biggestPair), theNA)


  if(length(newIJ)) {
    ijk[newIJ,'pair'] = biggestPair
    ijk[newIJ,'type'] = 'ij'
  }
  if(length(newJK)) {
    ijk[newJK,'pair'] = biggestPair
    ijk[newJK,'type'] = 'jk'
  }
  if(length(newIK)) {
    ijk[newIK,'pair'] = biggestPair
    ijk[newIK,'type'] = 'ik'
  }
  Niter = Niter +1
  Nna = sum(is.na(ijk$pair))
} # while


ijk$pairFac = factor(ijk$pair, levels=unique(apply(pairs, 1, function(xx) paste(sort(xx), collapse='_'))))
ijk[,Spairs] = as.integer(NA)
theIj = which(ijk$type == 'ij')
ijk[theIj,Spairs] = ijk[theIj,c('i','j','k')]
theJk = which(ijk$type == 'jk')
ijk[theJk,Spairs] = ijk[theJk,c('j','k','i')]
theIk = which(ijk$type == 'ik')
ijk[theIk,Spairs] = ijk[theIk,c('i','k','j')]



ijkp = cbind(
  ijk[order(ijk[,'pairFac'],ijk[,'pairs3']),],
  p=seq(from=0, len=nrow(ijk))
)
ijkp = ijkp[,c("pairs1","pairs2","pairs3","p")]

names(ijkp) = c('i','j','k', 'p')

ijkp[] <- lapply(ijkp, as.integer)
rownames(ijkp) = NULL

pairs = ijkp[!duplicated(ijkp[,c('i','j')]), setdiff(names(ijkp), 'k')]
pairs$pEnd = c(pairs[-1,'p'], nrow(ijkp))
pairs$n =pairs$pEnd - pairs$p

Nbeta = nrow(data$XTp)
Ngamma = nrow(data$ATp)
Ntotal = length(x)
Nparams = Ntotal - Ngamma  

Sgamma1 = seq(from=Nbeta+1, len=Ngamma)


hessianC = as(hessian, 'CsparseMatrix')
hessianT = as(hessianC,'TsparseMatrix')
hessianRandom = hessianT[Sgamma1, Sgamma1]

  # store upper and lower triangle
hessianC2 = as(hessianC, "generalMatrix")
hessianT2 = as(hessianC2, "TsparseMatrix")

nonsymmetric = list(
  i=hessianT2@i,
  j=hessianT2@j, 
  p=hessianC2@p
)
nonsymmetric$ij = paste(nonsymmetric$i, nonsymmetric$j, sep="_")


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
    pairs = as.data.frame(pairs)
  )
)



sparsity

}





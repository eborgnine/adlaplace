#' @export
sparsityForThird = function(x, data, config=list()) {

  configForDiag = config[setdiff(names(config), c("dense","sparsity"))]
  configForDiag$dense = TRUE
    resThirdDiag = thirdDiagonals(
      x, data, configForDiag
    ) 
  hessian = Matrix::forceSymmetric(Matrix::Matrix(resThirdDiag$second, sparse=TRUE))  
 
  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Ntotal = length(x)
  Nparams = Ntotal - Ngamma  

  Stotal0 = seq(0, len=Ntotal)

  Sgamma1 = seq(from=Nbeta+1, len=Ngamma)
  Sgamma0 = Sgamma1 -1

  Sparams1 = setdiff(1:Ntotal, Sgamma1)
  Sparams0 = Sparams1-1

    # do dense to find sparsity pattern
    lowertri = expand.grid(
      i=as.integer(Stotal0), j=as.integer(Stotal0)
    )
    lowertri = lowertri[lowertri$i > lowertri$j, ]
    ijk = cbind(
      lowertri[rep(1:nrow(lowertri), length(Sparams0)), ],
      k=rep(Sparams0, each=nrow(lowertri))
    )
    Nunique = apply(ijk[,c('i','j','k')],1, function(xx) length(unique(xx)))
    ijk = ijk[Nunique == 3, ]

    pairs = ijk[!duplicated(ijk[,c('i','j')]), setdiff(names(ijk), 'k')]
    pairs = pairs[order(pairs$j, pairs$i), ]
    pairs = as.matrix(pairs)

    resThirdSparse = thirdNonDiagonalsSparsity(
      x, data, config, pairs
    )

    resList = apply(resThirdSparse, 1, which)

    triplets = cbind(
      pairs[rep(1:nrow(pairs), unlist(lapply(resList, length))),],
      k = unlist(resList)-1
    )
    triplets2 = t(apply(triplets, 1, sort))
    triplets2 = triplets2[!duplicated(triplets2), ]

    NN = apply(triplets2, 1, function(xx) length(unique(xx)))
    ijk = as.data.frame(triplets2[NN==3,])
    names(ijk) = c('i','j','k')
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

ijk$pairFac = factor(ijk$pair, levels=apply(lowertri, 1, function(xx) paste(sort(xx), collapse='_')))
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

#sparsityType = c('full','pairs','allpairs')


sparsityForThirdOld = function(x, data, config=list()) {

  Ntotal = length(x)
  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Nparams = Ntotal - Ngamma  
  Sgamma1 = seq(from=Nbeta+1, len=Ngamma)
  Sgamma0 = Sgamma1-1

  Sparams1 = setdiff(1:Ntotal, Sgamma1)
  Sparams0 = Sparams1-1
  if(!length(config$sparsityType)) {
    config$sparsityType = 'pairs'
  }

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



if(nrow(hessian) <= Ngamma) 
  warning("hessian must be from the full set of parameters and random effects")

hessianC = as(hessian, 'CsparseMatrix')
hessianT = as(hessianC,'TsparseMatrix')

  # store upper and lower triangle
hessianC2 = as(hessianC, "generalMatrix")
hessianT2 = as(hessianC2, "TsparseMatrix")

nonsymmetric = list(
  i=hessianT2@i,
  j=hessianT2@j, 
  p=hessianC2@p
)
nonsymmetric$ij = paste(nonsymmetric$i, nonsymmetric$j, sep="_")


if(config$sparsityType == 'full') {
  lowertri = expand.grid(
    i=Sgamma0, j=Sgamma0
  )
  lowertri = lowertri[lowertri$i > lowertri$j, ]
  lowertri$ij = apply(lowertri, 1, paste, collapse='_')
  ijk = cbind(
    lowertri[rep(1:nrow(lowertri), length(Sparams0)), ],
    k=rep(Sparams0, each=nrow(lowertri))
  )
  Nunique = apply(ijk[,c('i','j','k')],1, function(xx) length(unique(xx)))
  ijk = ijk[Nunique == 3, ]
  ijk$pairFac = factor(ijk$ij, levels=lowertri$ij)

  ijkp = cbind(
    ijk[order(ijk[,'pairFac'],ijk[,'k']),],
    p=seq(from=0, len=nrow(ijk))
  )
  ijkp = ijkp[,c('i','j','k', 'p')]
  ijkp[] = lapply(ijkp, as.integer)


  pairs = ijkp[!duplicated(ijkp[,c('i','j')]), setdiff(names(ijkp), 'k')]
  pairs$pEnd = c(pairs[-1,'p'], nrow(ijkp))
  pairs$n =pairs$pEnd - pairs$p
  # ijkpFull = ijkp;pairsFull = pairs

} else { # not full

lowertri = data.frame(
  i=hessianT@i,j=hessianT@j
)
lowertri = lowertri[lowertri$j < lowertri$i,]
lowertri = lowertri[order(lowertri$j, lowertri$i),]
lowertri$ij = apply(lowertri, 1, paste, collapse='_')

if(config$sparsityType == 'pairs') {
  ijkList = mapply(
    function(k, lowertri) {
      lowertri$k = k
      lowertri
    },
    k = seq(0, len=Ntotal), 
    MoreArgs = list(lowertri=lowertri), # 466 593
    SIMPLIFY=FALSE
  )
} else if(config$sparsityType == 'allpairs') {
  ijkList = mapply(
    function(k, lowertri) {
      lowerK = lowertri[lowertri$j > k,c('i','j')]
      haveBoth = 
      (
        paste(lowerK$i, k, sep='_') %in% lowertri$ij
      ) & (
        paste(lowerK$j, k, sep='_') %in% lowertri$ij
      )
      lowerK = lowerK[haveBoth, ,drop=FALSE]
      lowerK$k=rep(k, sum(haveBoth))
      lowerK
    },
    k = seq(0, len=Ntotal), 
    MoreArgs = list(lowertri=lowertri), # 466 593
    SIMPLIFY=FALSE
  )
} else {
  warning("config$sparsityType should be full, pairs, or allpairs")
}

ijk = do.call(rbind, ijkList)
ijk$ij = apply(ijk[,c('i','j')], 1, paste, collapse='_')
ijk$jk = apply(ijk[,c('j','k')], 1, paste, collapse='_')
ijk$ik = apply(ijk[,c('i','k')], 1, paste, collapse='_')
ijk$index = 1:nrow(ijk)

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

ijk$pairFac = factor(ijk$pair, levels=lowertri$ij)
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

pairs = ijkp[!duplicated(ijkp[,c('i','j')]), setdiff(names(ijkp), 'k')]
pairs$pEnd = c(pairs[-1,'p'], nrow(ijkp))
pairs$n =pairs$pEnd - pairs$p
} # sparsity not full


hessianRandom = hessianT[Sgamma1,Sgamma1]

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


return(sparsity)
}



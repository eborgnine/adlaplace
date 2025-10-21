pairsStrataFun = function(Dstrata, x, config, data) {
# get third tensor for only one strata
config$sparsity$third$strata = rep(1, nrow(config$sparsity$third$pairs))
from3 = hpolcc::thirdOffDiagonals(x, data, config)
ijk= as.data.frame(config$sparsity$third$ijk)
ijk$x = abs(drop(from3))>0
ijk2 = aggregate(ijk[,'x', drop=FALSE], ijk[,c('i','j'), drop=FALSE], sum)
ijk2 = as.data.frame(
  ijk2[ijk2$x != 0, c('i','j'), drop=FALSE])
ijk2 = ijk2[order(ijk2$j, ijk2$i), ]
ijk2$strata = Dstrata
ijk2
}

#' @export
sparsity_pattern = function(x, data, config=list()) {

#library('hpolcc')

  if(identical(config$dense, "TRUE")) {
    type = 'compute'
  } else {
    type = 'hessian'
  }

  configForDiag = config[setdiff(names(config), c("dense","sparsity","beta","theta"))]
  configForDiag$dense=TRUE
  configForDiag$maxDeriv = 2

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Ntotal = length(x)
  Nparams = Ntotal - Ngamma  
  Sgamma1 = seq(from=Nbeta+1, len=Ngamma)
  Sparams = setdiff(seq(0, len=Ntotal), seq(from=Nbeta, len=Ngamma))

  # sometimes zero parameter values lead to zeros in the hessian.
  x[x==0] = 1e-2

  denseHessian = objectiveFunctionC(
    x, data, configForDiag
  )$denseHessian


  hessian = as(
    Matrix::forceSymmetric(Matrix::Matrix(denseHessian, sparse=TRUE)),
    'TsparseMatrix')  
  if(identical(config$verbose, TRUE)) {
    cat("done hessian ", nrow(hessian), " ", length(hessian@x), "\n")
  }

  if(any(is.na(hessian@x))) warning("NA's in hessian")

  hessianUL = as(hessian, 'generalMatrix')

  hessianIJ = data.frame(i = hessianUL@i, j=hessianUL@j)
  hessianIJ = hessianIJ[order(hessianIJ$j, hessianIJ$i), ]

  pairs = hessianIJ[!duplicated(hessianIJ[,c('i','j')]), ]
  pairs = pairs[order(pairs$j, pairs$i), ]
  pairs = as.matrix(pairs)

    if(identical(config$verbose, TRUE)) {
      cat("third non-zeros\n")
    }
  if(type[1] == 'hessian') {
    # find sparsity pattern based on the hessian
    if(identical(config$verbose, TRUE)) {
      cat("third non-zeros from hessian\n")
    }
    Sk =  apply(hessianIJ, 1, 
      function(xx, ref) {
        intersect(ref[ref$i == xx['i'], 'j'], ref[ref$j == xx['j'], 'i'])
      }, 
      ref = hessianIJ)
    ijk1 = hessianIJ[rep(1:nrow(hessianIJ), unlist(lapply(Sk, length))), ]
    ijk1$k = unlist(Sk)

  } else {
    if(identical(config$verbose, TRUE)) {
      cat("computing offdiag\n")
    }
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
  ijk2 = as.data.frame(ijk1[Nunique == 3, ])

  # only one parameter per trio, dont need Tijk for two parameters
  Nparams = apply(ijk2, 1, function(xx) sum(xx %in% Sparams))
  ijk = ijk2[Nparams < 2,,drop=FALSE]

  colnames(ijk) = c('i','j','k')
  ijk$ij = apply(ijk[,c('i','j')], 1, paste, collapse='_')
  ijk$jk = apply(ijk[,c('j','k')], 1, paste, collapse='_')
  ijk$ik = apply(ijk[,c('i','k')], 1, paste, collapse='_')

  Sdim = c('ij','jk','ik')
  Spairs = c("pairs1","pairs2","pairs3")
  ijk$pair = ijk$type = NA

  if(identical(config$verbose, TRUE)) {
    cat("finding pairs\n")
  }
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

  if(identical(config$verbose, TRUE)) {
    cat("done pairs ", Niter, "\n")
  }


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

    # find optimal indieces for evaluating hessian
  pairsRandom = data.frame(i=hessianRandom@i, j=hessianRandom@j)
  pairsHessian = pairsRandom[pairsRandom$j > pairsRandom$i, ]
  pairsHessian$pair = pairsHessian$other = NA
  Nna=1;Niter = 0
  Sdim = c('i','j')
  while(Nna > 0 & Niter < nrow(pairsHessian)) {
    theNA = which(is.na(pairsHessian$pair))
    theTable = table(unlist(pairsHessian[theNA, Sdim]))
    biggestPair = as.integer(names(theTable)[which.max(theTable)])
    matchI = intersect(which(pairsHessian$i == biggestPair), theNA)      
    matchJ = intersect(which(pairsHessian$j == biggestPair), theNA)      
    pairsHessian[matchI, 'other'] = pairsHessian[matchI, 'j']
    pairsHessian[matchJ, 'other'] = pairsHessian[matchJ, 'i']
    pairsHessian[c(matchI, matchJ),'pair'] = biggestPair

    Niter = Niter +1
    Nna = sum(is.na(pairsHessian$pair))
  }

  pairsHessian1 = as.data.frame(pairsHessian)[,c('pair','other')]
  names(pairsHessian1) = c('j','i')

  tableJ = as.data.frame(table(pairsHessian1$j))
  isOnes = pairsHessian1$j %in% as.integer(as.character(tableJ[tableJ$Freq == 1, 'Var1']))
  jOnes = pairsHessian1[isOnes,,drop=FALSE]
  pairsHessian2 = pairsHessian1[!isOnes, ,drop=FALSE]

  uniqueJ = unique(pairsHessian2[,'j'])
  pairsHessian2 = rbind(pairsHessian2, data.frame(j=uniqueJ, i=uniqueJ))
  pairsHessian2 = pairsHessian2[order(pairsHessian2$j, pairsHessian2$i), ]

  pairsHessian2 = as.list(pairsHessian2[order(pairsHessian2$j, pairsHessian2$i), c('i','j')])
  theTable = table(pairsHessian2$j)
  pairsHessian2$p = cumsum(c(0, theTable))
  pairsHessian2$j = as.integer(names(theTable))
  pairsHessian2$jLong = rep(pairsHessian2$j , diff(pairsHessian2$p))

  pairsHessian2 <- lapply(pairsHessian2, as.integer)

  pairsHessian2$onesJ = jOnes$j
  pairsHessian2$onesI = jOnes$i

  pairsHessian2$SjNotPair = setdiff(seq(0, len=Ngamma),unlist(pairsHessian2[c('i','j')]))
  pairsHessian2$SjNotPairOrDiag = setdiff(pairsHessian2$SjNotPair, unlist(pairsHessian2[c('onesI', 'onesJ')]))
  pairsHessian2$SiNotJ = sort(setdiff(pairsHessian2$i, c(pairsHessian2$j, pairsHessian2$onesJ)))



  sparsity = list(
    random = pairsHessian2, 
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

# which strata contribute to third
Nstrata = ncol(data$cc_matrixTp)
configForPairs = configForDiag
configForPairs$dense = FALSE
configForPairs$verbose = FALSE
configForPairs$num_threads = 1
configForPairs$sparsity = sparsity
configForPairs$sparsity$third$pairs$Nstrata = 1
configForPairs$sparsity$third$pairs$pStrata = seq(0, len=nrow(configForPairs$sparsity$third$pairs))
configForPairs$sparsity$third$pairs$pStrataEnd = configForPairs$sparsity$third$pairs$pStrata +1



pairsStrataList = parallel::mcmapply(
  pairsStrataFun,
Dstrata = seq(0, len=Nstrata),
MoreArgs = list(config = configForPairs, data=data, x = x),
mc.cores = config$num_threads, 
SIMPLIFY=FALSE )

pairsStrata = do.call(rbind, pairsStrataList)
pairsStrataN = aggregate(data.frame(N=rep(1, nrow(pairsStrata))), pairsStrata[,c('i','j')], sum)
# if over 80%, set to zero (so all strata are used)
pairsStrataN[pairsStrataN$N > ceiling(0.8*Nstrata), 'N'] = 0

pairsStrataN$pair = paste(pairsStrataN$i, pairsStrataN$j, sep='_')
pairsStrata$pair = paste(pairsStrata$i, pairsStrata$j, sep='_')
pairsStrataIndex = pairsStrata[pairsStrata$pair %in% pairsStrataN[pairsStrataN$N >0, 'pair'], , drop=FALSE]

#pairsStrataN = pairsStrataN[order(pairsStrataN$j, pairsStrataN$i), ]
pairsStrataIndex = pairsStrataIndex[order(pairsStrataIndex$j, pairsStrataIndex$i), ]


pairsStrataN$p = match(pairsStrataN$pair, pairsStrataIndex$pair)-1
pairsStrataN$pEnd = mapply(function(x, y) min(c(Inf, y[which(y > x)]), na.rm=TRUE), x=pairsStrataN$p, MoreArgs = list(y=pairsStrataN$p))
pairsStrataN[is.na(pairsStrataN$p), 'p'] = -1
pairsStrataN[pairsStrataN$pEnd==Inf, 'pEnd'] = -1


sparsity$third$pairs$pair = paste(sparsity$third$pairs$i, sparsity$third$pairs$j, sep='_')

theMatch = match(sparsity$third$pairs$pair, pairsStrataN$pair)
sparsity$third$pairs$pStrata = as.integer(pairsStrataN[theMatch, 'p'])
sparsity$third$pairs$pStrataEnd = as.integer(pairsStrataN[theMatch, 'pEnd'])
sparsity$third$pairs$Nstrata = as.integer(pairsStrataN[theMatch, 'N'])

sparsity$third$strata = as.integer(pairsStrata$strata)


  sparsity

}





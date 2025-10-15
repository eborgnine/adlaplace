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

  # sometimes zero parameter values lead to zeros in the hessian.
  x[x==0] = 1e-2

  denseHessian = objectiveFunctionC(
    x, data, configForDiag
  )$denseHessian

  if(FALSE) {
    resThirdDiag = thirdDiagonals(
      x, data, configForDiag
    ) 
    denseHessian = resThirdDiag$second
  }


  
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
  ijk = as.data.frame(ijk1[Nunique == 3, ])
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

  pairsHessian2 = as.data.frame(pairsHessian)[,c('pair','other')]
  names(pairsHessian2) = c('j','i')

  tableJ = as.data.frame(table(pairsHessian2$j))
  isOnes = pairsHessian2$j %in% as.integer(as.character(tableJ[tableJ$Freq == 1, 'Var1']))
  jOnes = pairsHessian2[isOnes,,drop=FALSE]
  pairsHessian2 = pairsHessian2[!isOnes, ,drop=FALSE]

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



  sparsity

}





#' @export
thirdDeriv = function(x, data, config) {

   # computing T_kii and H_ki, columnns are i, rows are k 
  resThirdDiag = thirdDiagonals(
    x, data, config
  ) 

  resThirdOffDiag = thirdOffDiagonals(
    x, data, config
  ) 

  if(identical(config$dense, TRUE)) {
          # T_iik, doubles are columns of resThirdDiag$diag
    fullHessian = Matrix::Matrix(resThirdDiag$second)
    thirdDiag = data.frame(
      i = rep(seq(0, len=nrow(resThirdDiag$diag)), each=ncol(resThirdDiag$diag)),
      k = rep(seq(0, len=ncol(resThirdDiag$diag)), nrow(resThirdDiag$diag)),
      x = 2*as.vector(resThirdDiag$diag)
    )
    thirdDiag$j = thirdDiag$i

    thirdNonDiag = config$sparsity$third$pairs[
    rep(1:nrow(config$sparsity$third$pairs), each=nrow(resThirdOffDiag)),c('i','j')]
    thirdNonDiag$k = rep(seq(0, len=nrow(resThirdOffDiag)), ncol(resThirdOffDiag))
    thirdNonDiag$taylor3 = as.vector(resThirdOffDiag)
    thirdNonDiag = thirdNonDiag[apply(thirdNonDiag[,c('i','j','k')], 1, lengthUnique)==3, ]
  } else {
    thirdDiag = data.frame(
      i = config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$i,
      k = config$sparsity$second$nonSymmetric$j,
      x = 2*(resThirdDiag$diag)
    )
    fullHessian = Matrix::forceSymmetric(
      Matrix::sparseMatrix(
        j = config$sparsity$second$nonSymmetric$j,
        i = config$sparsity$second$nonSymmetric$i,
        x = drop(resThirdDiag$second),
        dims = rep(length(resThirdDiag$first), 2), index1=FALSE
      ))

    thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
    thirdNonDiag$taylor3 = drop(resThirdOffDiag)
  }

  # 3rd taylor is   T_iik/2 + T_jjk/2 + T_ijk
  # pairs are i and j
  thirdDiag$ik = apply(thirdDiag[,c('i','k')],1,paste,collapse='_')
  thirdNonDiag$ik = apply(thirdNonDiag[,c('i','k')],1,paste,collapse='_')
  thirdNonDiag$jk = apply(thirdNonDiag[,c('j','k')],1,paste,collapse='_')
  matchIik = match(thirdNonDiag$ik, thirdDiag$ik)
  matchJjk = match(thirdNonDiag$jk, thirdDiag$ik)
  thirdNonDiag$Tiik = thirdDiag[matchIik, 'x']
  thirdNonDiag$Tjjk = thirdDiag[matchJjk, 'x']

  thirdNonDiag$x = thirdNonDiag$taylor3 - 0.5*(thirdNonDiag$Tiik  + thirdNonDiag$Tjjk)
  
# might be some duplicated
  nonDiagUnique = t(apply(thirdNonDiag[,c('i','j','k')], 1, sort))  
  thirdNonDiag = thirdNonDiag[!duplicated(nonDiagUnique), ]
# thirdNonDiag[apply(thirdNonDiag[,c('i','j','k')], 1, function(xx) all(xx %in% (c(2,3,4)-1))),]       


  theCols = c('i','j','k', 'x')
  third = rbind(
    thirdDiag[,theCols],
    thirdNonDiag[,theCols]
  )
#third[apply(third[,c('i','j','k')], 1, function(xx) all(xx %in% (c(2,3,4)-1))),]       

  thirdList = mapply(
    thirdTensor,
    k = seq(from=0, len=length(x)), 
    MoreArgs = list(third=third, N=length(x))
  )

list(fullHessian=fullHessian, thirdList =thirdList, first = resThirdDiag$first)

}
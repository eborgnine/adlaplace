#' @export
formatHpolData = function(data) {
	for(Ddata in c('X', 'A', 'cc_matrix')) {
      data[[paste0(Ddata, 'Tp')]] = as(Matrix::t(data[[Ddata]]), 'CsparseMatrix')
  	}
	if(!'Qdiag' %in% names(data)) {
		data$Qdiag = Matrix::diag(data$Q)
	}
	if(!'QsansDiag' %in% names(data)) {
		QsansDiag = as(Matrix::forceSymmetric(data$Q), 'TsparseMatrix')
		diag(QsansDiag) = 0
		data$QsansDiag = QsansDiag
	} else {
		if(!any(class(data$QsansDiag == 'TsparseMatrix'))) {
			data$QsansDiag = as(Matrix::forceSymmetric(data$QsansDiag), 'TsparseMatrix')
		}
	}
	
	ccMatrixToCheck = data$cc_matrix[seq(1, min(c(100, length(data$cc_matrix))))]

	data$y = as.integer(data$y)
	if(length(data$y) != ncol(data$ATp)) warning("data wrong length")
	 	
	return(data[setdiff(names(data), c('X','A','cc_matrix'))])
}

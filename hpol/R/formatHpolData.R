#' @export
formatHpolData = function(data) {
	for(Ddata in c('X', 'A', 'cc_matrix')) {
      data[[paste0(Ddata, 'Tp')]] = as(Matrix::t(data[[Ddata]]), 'dgCMatrix')
  	}
	if(!'Qdiag' %in% names(data)) {
		data$Qdiag = Matrix::diag(data$Q)
	}
	if(!'QsansDiag' %in% names(data)) {
		QsansDiag = as(Matrix::forceSymmetric(data$Q), 'dsTMatrix')
		diag(QsansDiag) = 0
		data$QsansDiag = QsansDiag
	} else {
		if(!any(class(data$QsansDiag == 'dsTMatrix'))) {
			data$QsansDiag = as(Matrix::forceSymmetric(data$QsansDiag), 'dsTMatrix')
		}
	}
	
	ccMatrixToCheck = cc_matrix[seq(1, min(c(100, length(cc_matrix))))]

	 	
	return(data)
}

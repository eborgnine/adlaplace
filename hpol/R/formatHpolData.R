#' @export
formatHpolData = function(data) {
	requireNamespace('Matrix')
	newMatrixClasses = c(X='dgRMatrix', A='dgRMatrix', cc_matrix = 'dgRMatrix')
	for(Ddata in names(newMatrixClasses)) {
		if(!any(class(data[[Ddata]] == newMatrixClasses[Ddata]))) {
			data[[Ddata]] = as(data[[Ddata]], newMatrixClasses[Ddata])
		}
	}
	if(!'Qdiag' %in% names(data)) {
		data$Qdiag = diag(data$Q)
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
	return(data)
}

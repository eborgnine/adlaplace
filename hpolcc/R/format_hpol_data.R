#' @export
formatHpolData <- function(data) {

  result <- NULL
  for (Ddata in c("X", "A")) {
    data[[paste0(Ddata, "Tp")]] <- 
      methods::as(
        methods::as(Matrix::t(data[[Ddata]]), "CsparseMatrix"),
        "dMatrix"
      )
    if(any(is.na(data[[paste0(Ddata, "Tp")]]@x))) {
      warning("missing data in ", Ddata)
    }
  }

  data$elgm_matrix <- methods::as(data$elgm_matrix, "CsparseMatrix")
  if (!"Qdiag" %in% names(data)) {
    data$Qdiag <- Matrix::diag(data$Q)
  }
  if (!"QsansDiag" %in% names(data)) {
    QsansDiag <- methods::as(Matrix::forceSymmetric(data$Q), "TsparseMatrix")
    diag(QsansDiag) <- 0
    data$QsansDiag <- QsansDiag
  } else {
    if (!any(class(data$QsansDiag == "TsparseMatrix"))) {
      data$QsansDiag <- methods::as(
        Matrix::forceSymmetric(data$QsansDiag),
        "TsparseMatrix"
      )
    }
  }
  data$y <- as.integer(data$y)
  if (length(data$y) != ncol(data$ATp)) warning("data wrong length")

  result = data[setdiff(names(data), c("X", "A"))]

  result
}

#' @export
formatHpolData <- function(data) {

  result = NULL
  for (Ddata in c("X", "A")) {
    data[[paste0(Ddata, "Tp")]] <- as(Matrix::t(data[[Ddata]]), "CsparseMatrix")
  }

  data$elgm_matrix <- as(data$elgm_matrix, "CsparseMatrix")
  if (!"Qdiag" %in% names(data)) {
    data$Qdiag <- Matrix::diag(data$Q)
  }
  if (!"QsansDiag" %in% names(data)) {
    QsansDiag <- as(Matrix::forceSymmetric(data$Q), "TsparseMatrix")
    diag(QsansDiag) <- 0
    data$QsansDiag <- QsansDiag
  } else {
    if (!any(class(data$QsansDiag == "TsparseMatrix"))) {
      data$QsansDiag <- as(
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

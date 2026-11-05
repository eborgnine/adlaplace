# Verbatim copies of Ziang Zhang / BayesGP seasonal-GP reference helpers
# for regression tests only. Source:
#   https://github.com/Bayes-GP/BayesGP (master)
#   R/01_utility.R: Compute_Q_sB, Compute_B_sB, compute_d_step_sgpsd
# Not part of the adlaplaceFem runtime API. FEM tests skip without Suggests: fda.

Compute_Q_sB <- function(a, k, region, accuracy = 5000, boundary = TRUE) {
  ss <- function(M) {
    Matrix::forceSymmetric(M + Matrix::t(M))
  }
  if ((accuracy %% 1) == 0) {
    x <- seq(min(region), max(region), length.out = accuracy)
  } else {
    x <- seq(min(region), max(region), by = accuracy)
  }
  if (boundary == TRUE) {
    B_basis <- suppressWarnings(fda::create.bspline.basis(
      rangeval = c(min(region), max(region)),
      nbasis = k,
      norder = 4,
      dropind = c(1, 2)
    ))
  } else {
    B_basis <- suppressWarnings(fda::create.bspline.basis(
      rangeval = c(min(region), max(region)),
      nbasis = k,
      norder = 4
    ))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj = 0, returnMatrix = TRUE)
  B1matrix <- fda::eval.basis(x, B_basis, Lfdobj = 1, returnMatrix = TRUE)
  B2matrix <- fda::eval.basis(x, B_basis, Lfdobj = 2, returnMatrix = TRUE)
  cos_matrix <- cos(a * x)
  sin_matrix <- sin(a * x)
  Bcos <- as(apply(Bmatrix, 2, function(x) x * cos_matrix), "dgCMatrix")
  B1cos <- as(apply(B1matrix, 2, function(x) x * cos_matrix), "dgCMatrix")
  B2cos <- as(apply(B2matrix, 2, function(x) x * cos_matrix), "dgCMatrix")
  Bsin <- as(apply(Bmatrix, 2, function(x) x * sin_matrix), "dgCMatrix")
  B1sin <- as(apply(B1matrix, 2, function(x) x * sin_matrix), "dgCMatrix")
  B2sin <- as(apply(B2matrix, 2, function(x) x * sin_matrix), "dgCMatrix")

  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0, x)))), "dgCMatrix")

  ### T
  T00 <- Matrix::t(Bcos) %*% Numerical_I %*% Bcos
  T10 <- Matrix::t(B1cos) %*% Numerical_I %*% Bcos
  T11 <- Matrix::t(B1cos) %*% Numerical_I %*% B1cos
  T20 <- Matrix::t(B2cos) %*% Numerical_I %*% Bcos
  T21 <- Matrix::t(B2cos) %*% Numerical_I %*% B1cos
  T22 <- Matrix::t(B2cos) %*% Numerical_I %*% B2cos

  ### L
  L00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bsin
  L10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bsin
  L11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1sin
  L20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bsin
  L21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1sin
  L22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2sin

  ### I
  I00 <- Matrix::t(Bsin) %*% Numerical_I %*% Bcos
  I10 <- Matrix::t(B1sin) %*% Numerical_I %*% Bcos
  I11 <- Matrix::t(B1sin) %*% Numerical_I %*% B1cos
  I20 <- Matrix::t(B2sin) %*% Numerical_I %*% Bcos
  I21 <- Matrix::t(B2sin) %*% Numerical_I %*% B1cos
  I22 <- Matrix::t(B2sin) %*% Numerical_I %*% B2cos

  ### Inner product involving B spline:
  Bmatrix <- as(Bmatrix, "dgCMatrix")
  B1matrix <- as(B1matrix, "dgCMatrix")
  B2matrix <- as(B2matrix, "dgCMatrix")

  BB <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bmatrix
  B2B2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  BB2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2matrix
  BS <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bsin
  BC <- Matrix::t(Bmatrix) %*% Numerical_I %*% Bcos
  BS1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1sin
  BC1 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B1cos
  BS2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2sin
  BC2 <- Matrix::t(Bmatrix) %*% Numerical_I %*% B2cos
  B2S <- Matrix::t(B2matrix) %*% Numerical_I %*% Bsin
  B2C <- Matrix::t(B2matrix) %*% Numerical_I %*% Bcos
  B2S1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1sin
  B2C1 <- Matrix::t(B2matrix) %*% Numerical_I %*% B1cos
  B2S2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2sin
  B2C2 <- Matrix::t(B2matrix) %*% Numerical_I %*% B2cos

  ## G = <phi,phj>
  G <- rbind(
    cbind(T00, Matrix::t(I00), Matrix::t(BC)),
    cbind(I00, L00, Matrix::t(BS)),
    cbind(BC, BS, BB)
  )

  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2 * a * ss(I21) - (a^2) * ss(T20) + 2 * (a^3) * ss(I10) +
    4 * (a^2) * L11 + (a^4) * T00
  C22 <- L22 + 2 * a * ss(I21) - (a^2) * ss(L20) - 2 * (a^3) * ss(I10) +
    4 * (a^2) * T11 + (a^4) * L00
  C12 <- I22 + 2 * a * T21 - (a^2) * ss(I20) - 2 * a * Matrix::t(L21) -
    4 * (a^2) * I11 + 2 * (a^3) * L10 - 2 * (a^3) * Matrix::t(T10) + (a^4) * I00

  C13 <- Matrix::t(B2C2) - 2 * a * Matrix::t(B2S1) - (a^2) * Matrix::t(B2C)
  C23 <- Matrix::t(B2S2) + 2 * a * Matrix::t(B2C1) - (a^2) * Matrix::t(B2S)
  C33 <- B2B2
  C <- rbind(
    cbind(C11, C12, C13),
    cbind(Matrix::t(C12), C22, C23),
    cbind(Matrix::t(C13), Matrix::t(C23), C33)
  )

  ## M = <phi,D^2phj>
  M11 <- Matrix::t(T20) - (2 * a) * Matrix::t(I10) - (a^2) * T00
  M12 <- Matrix::t(I20) + (2 * a) * Matrix::t(T10) - (a^2) * I00
  M21 <- Matrix::t(I20) - (2 * a) * Matrix::t(L10) - (a^2) * I00
  M22 <- Matrix::t(L20) + (2 * a) * Matrix::t(I10) - (a^2) * L00

  M13 <- Matrix::t(B2C)
  M23 <- Matrix::t(B2S)
  M31 <- BC2 - (2 * a) * BS1 - (a^2) * BC
  M32 <- BS2 + (2 * a) * BC1 - (a^2) * BS
  M33 <- BB2

  M <- rbind(cbind(M11, M12, M13), cbind(M21, M22, M23), cbind(M31, M32, M33))

  ### Compute the final precision matrix: Q
  Q <- (a^4) * G + C + (a^2) * ss(M)
  Matrix::forceSymmetric(Q)
}

Compute_B_sB <- function(x, a, k, region, boundary = TRUE) {
  if (boundary) {
    B_basis <- suppressWarnings(fda::create.bspline.basis(
      rangeval = c(min(region), max(region)),
      nbasis = k,
      norder = 4,
      dropind = c(1, 2)
    ))
  } else {
    B_basis <- suppressWarnings(fda::create.bspline.basis(
      rangeval = c(min(region), max(region)),
      nbasis = k,
      norder = 4
    ))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj = 0, returnMatrix = TRUE)
  cos_matrix <- cos(a * x)
  sin_matrix <- sin(a * x)
  Bcos <- apply(Bmatrix, 2, function(x) x * cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x * sin_matrix)
  cbind(Bcos, Bsin, Bmatrix)
}

compute_d_step_sgpsd <- function(d, a) {
  sqrt((1 / (a^2)) * ((d / 2) - (sin(2 * a * d) / (4 * a))))
}

#' Build grouped Hessian sparsity structures
#'
#' Constructs Hessian sparsity templates (global/outer and inner-gamma)
#' and per-group maps from shard sparsity patterns.
#'
#' @param sparsity_list List of per-group sparsity shards from \code{getAdFun()}.
#' @param Nbeta Number of fixed-effect parameters.
#' @param Ngamma Number of random-effect parameters.
#' @param Ntheta Number of variance parameters.
#'
#' @return A list with \code{outer}, \code{inner}, optional \code{chol_inner},
#'   \code{map_outer}, \code{map_inner}, and \code{sizes} (named
#'   \code{beta}/\code{gamma}/\code{theta}; suitable for
#'   \code{adlaplace_attach_hessian()}).
#'
#' @export
hessian_map <- function(sparsity_list, Nbeta, Ngamma, Ntheta) {
  if (Nbeta < 0L || Ngamma < 0L || Ntheta < 0L) {
    stop("Nbeta, Ngamma, and Ntheta must be non-negative")
  }

  Ngroups <- length(sparsity_list)
  Nparams <- Nbeta + Ngamma + Ntheta
  Sgamma0 <- seq.int(Nbeta, length.out = Ngamma)

  make_sparse_df <- function(shard_id, row, col) {
    n <- length(col)
    data.frame(
      shard = rep.int(as.integer(shard_id), n),
      row = as.integer(row),
      col = as.integer(col),
      local = seq.int(0L, length.out = n),
      stringsAsFactors = FALSE
    )
  }

  bind_sparse_df <- function(lst) {
    lst <- lst[lengths(lst) > 0L]
    if (length(lst) == 0L) {
      return(data.frame(
        shard = integer(0),
        row = integer(0),
        col = integer(0),
        local = integer(0),
        stringsAsFactors = FALSE
      ))
    }
    do.call(rbind, lst)
  }

  list_outer <- vector("list", Ngroups)
  list_inner <- vector("list", Ngroups)

  for (D in seq_len(Ngroups)) {
    shard <- sparsity_list[[D]]
    list_outer[[D]] <- make_sparse_df(
      shard_id = D - 1L,
      row = shard$row_hess,
      col = shard$col_hess
    )
    list_inner[[D]] <- make_sparse_df(
      shard_id = D - 1L,
      row = shard$row_hess_inner,
      col = shard$col_hess_inner
    )
  }

  sparsityOuter <- bind_sparse_df(list_outer)
  sparsityInner <- bind_sparse_df(list_inner)

  sparsityOuter$cell <- sparsityOuter$row + sparsityOuter$col * Nparams
  sparsityOuter$cellSparse <- as.integer(factor(sparsityOuter$cell)) - 1L

  sparsityInner$rowInner <- match(sparsityInner$row, Sgamma0) - 1L
  sparsityInner$colInner <- match(sparsityInner$col, Sgamma0) - 1L
  if (anyNA(sparsityInner$rowInner) || anyNA(sparsityInner$colInner)) {
    stop("inner hessian index out of gamma range")
  }
  sparsityInner$cell <- sparsityInner$rowInner + Ngamma * sparsityInner$colInner
  sparsityInner$cellSparse <- as.integer(factor(sparsityInner$cell)) - 1L

  outerMat <- sparsityOuter[!duplicated(sparsityOuter$cell), , drop = FALSE]
  innerMat <- sparsityInner[!duplicated(sparsityInner$cell), , drop = FALSE]

  hessian_outer <- Matrix::sparseMatrix(
    i = outerMat$row,
    j = outerMat$col,
    x = outerMat$cellSparse,
    repr = "C",
    symmetric = TRUE,
    index1 = FALSE,
    dims = c(Nparams, Nparams)
  )

  hessian_inner <- Matrix::sparseMatrix(
    i = innerMat$rowInner,
    j = innerMat$colInner,
    x = innerMat$cellSparse,
    repr = "C",
    symmetric = TRUE,
    index1 = FALSE,
    dims = c(Ngamma, Ngamma)
  )

  hessian_outer_map <- Matrix::sparseMatrix(
    i = sparsityOuter$cellSparse,
    j = sparsityOuter$shard,
    x = sparsityOuter$local,
    index1 = FALSE,
    repr = "C",
    dims = c(length(unique(sparsityOuter$cellSparse)), Ngroups)
  )

  hessian_inner_map <- Matrix::sparseMatrix(
    i = sparsityInner$cellSparse,
    j = sparsityInner$shard,
    x = sparsityInner$local,
    index1 = FALSE,
    repr = "C",
    dims = c(length(unique(sparsityInner$cellSparse)), Ngroups)
  )

  result_map <- list(
    outer = list(
      p = as.integer(hessian_outer_map@p),
      local = as.integer(hessian_outer_map@x),
      global = as.integer(hessian_outer_map@i),
      dims = dim(hessian_outer_map)
    ),
    inner = list(
      p = as.integer(hessian_inner_map@p),
      local = as.integer(hessian_inner_map@x),
      global = match(hessian_inner_map@i, as.integer(hessian_inner@x)) - 1L,
      dims = c(length(hessian_inner@x), ncol(hessian_inner_map))
    )
  )

  chol_inner <- NULL
  if (Ngamma > 0L && nrow(innerMat) > 0L) {
    chol_inner <- tryCatch(
      {
        dummy_x <- ifelse(innerMat$rowInner == innerMat$colInner, 10, 1)
        hi_dummy <- Matrix::sparseMatrix(
          i = innerMat$rowInner,
          j = innerMat$colInner,
          x = dummy_x,
          symmetric = TRUE,
          index1 = FALSE,
          dims = c(Ngamma, Ngamma)
        )
        Matrix::Cholesky(hi_dummy, perm = TRUE, LDL = TRUE)
      },
      error = function(e) NULL
    )
  }

  L <- Matrix::expand2(chol_inner)$L1
  Linv <- as(Matrix::solve(L), "nMatrix")

  list(
    outer = hessian_outer,
    inner = hessian_inner,
    chol_inner = chol_inner,
    Linv = Linv,
    map_outer = result_map$outer,
    map_inner = result_map$inner,
    sizes = c(beta = Nbeta, gamma = Ngamma, theta = Ntheta)
  )
}

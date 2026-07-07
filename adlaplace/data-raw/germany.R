# Maintainer script: requires INLA to build package data (not a runtime dependency).
if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("Install INLA to rebuild the germany dataset.", call. = FALSE)
}

data(Germany, package = "INLA")
g_file <- system.file("demodata/germany.graph", package = "INLA")
graph <- INLA::inla.read.graph(g_file)
adj <- INLA::inla.graph2matrix(graph)
adj <- methods::as(adj, "CsparseMatrix")

ref_id <- which.max(apply(adj, 1, sum))

deg <- Matrix::rowSums(adj)
Q_icar <- Matrix::Diagonal(nrow(adj), deg) - adj
Q_icar <- methods::as(methods::as(Q_icar, "generalMatrix"), "CsparseMatrix")

Q_scaled <- INLA::inla.scale.model(
  Q_icar,
  constr = list(A = matrix(1, 1, nrow(adj)), e = 0)
)
Q_scaled <- methods::as(methods::as(Q_scaled, "generalMatrix"), "CsparseMatrix")

ev <- eigen(as.matrix(Q_scaled), symmetric = TRUE, only.values = TRUE)$values
ev <- sort(ev, decreasing = TRUE)
log_det_gen <- sum(log(ev[seq_len(nrow(adj) - 1L)]))

germany <- list(
  Y = as.numeric(Germany$Y),
  E = as.numeric(Germany$E),
  region = as.integer(Germany$region),
  adj = adj,
  Q_scaled = Q_scaled,
  prec = list(
    Q = Q_scaled,
    log_det = log_det_gen,
    rank = nrow(adj) - 1L
  )
)

usethis::use_data(germany, overwrite = TRUE, internal = FALSE)

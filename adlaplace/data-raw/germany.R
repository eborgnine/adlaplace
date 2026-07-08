# Maintainer script: requires INLA to build package data (not a runtime dependency).
if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("Install INLA to rebuild the germany dataset.", call. = FALSE)
}

library(Matrix)
library(terra)
source(system.file("demodata/Bym-map.R", package = "INLA"))

id_seq <- unlist(lapply(germany, function(xx) xx[1, 1]))

xy <- as.data.frame(do.call(rbind, germany))
names(xy) <- c("id", "x", "y")

xy$newregion <- c(1, diff(xy$id)) != 0
xy$new_subregion <- xy$newregion + is.na(xy$x)
xy$subregion <- cumsum(xy$new_subregion)
xy <- xy[!is.na(xy$x) & !is.na(xy$y), ]

xy_split <- split(xy, xy$subregion)
vect_list <- lapply(xy_split, function(x) {
  terra::vect(cbind(c(x$x, x$x[1]), c(x$y, x$y[1])), type = "polygon")
})
gmap_1 <- terra::vect(vect_list)

gmap_1$id <- unlist(lapply(xy_split, function(x) x[1, "id"]))
gmap_1$subregion <- unlist(lapply(xy_split, function(x) x[1, "subregion"]))
bob <- Matrix::Matrix(terra::relate(gmap_1, gmap_1, "contains"))
Matrix::diag(bob) <- 0

duplicated_polys <- Matrix::Matrix(Matrix::as.matrix(bob) == 1 & t(Matrix::as.matrix(bob)) == 1)

gmap_1$duplicated <- rowSums(duplicated_polys) > 0
gmap_1$has_hole <- rowSums(bob) > 0
gmap_1$is_hole <- colSums(bob) > 0
gmap_1$neither <- !gmap_1$has_hole & !gmap_1$is_hole
gmap_1$is_hole_but_
gmap_1$area <- terra::expanse(gmap_1)

id_with_holes <- unique(gmap_1$id[gmap_1$has_hole])
id_are_holes <- unique(gmap_1$id[gmap_1$is_hole])

gmap_1$is_hole_and_has_hole <- gmap_1$is_hole & (gmap_1$id %in% id_with_holes)

values(gmap_1)[gmap_1$id %in% gmap_1$id[c(40, 43, 44)], ]

# gmap_1 = gmap_1[!gmap_1$is_hole_and_has_hole,]


plot(gmap_1[40:45, ])
text(centroids(gmap_1[40:45, ]), labels = 40:45)
plot(gmap_1[gmap_1$neither, ], col = "yellow", add = TRUE)
plot(gmap_1[gmap_1$has_hole, ], add = TRUE, col = "#00FF0040")
plot(gmap_1[gmap_1$is_hole, ], add = TRUE, col = "#FF000040")

gmap_2 <- terra::erase(gmap_1[gmap_1$has_hole, ], gmap_1[gmap_1$is_hole, ])
gmap_no_holes <- gmap_1[!gmap_1$has_hole, ]
gmap_3 <- terra::vect(c(gmap_2, gmap_no_holes))

gmap_1$in_merged <- gmap_1$id %in% gmap_3$id
gmap_4 <- terra::vect(c(gmap_1[!gmap_1$in_merged, ], gmap_3))



gmap <- terra::aggregate(gmap_4, "id")
gmap <- gmap[match(unique(xy$id), gmap$id), ]

adj <- Matrix::Matrix(terra::adjacent(gmap, type = "intersects", symmetric = TRUE, pairs = FALSE))
Matrix::diag(adj) <- 1

g_file <- system.file("demodata/germany.graph", package = "INLA")
graph <- INLA::inla.read.graph(g_file)
adj2 <- INLA::inla.graph2matrix(graph)
adj2 <- methods::as(adj2, "CsparseMatrix")

adj_diff_1 <- adj - adj2

table(adj_diff_1@x)

adj_diff <- Matrix::Matrix(as.matrix(adj_diff_1))
the_max <- min(which(Matrix::rowSums((adj_diff)) > 0))
which_touch <- setdiff(which(adj[, the_max] != 0), the_max)
which_touch_inla <- setdiff(which(adj2[, the_max] != 0), the_max)

plot(gmap[c(the_max, which_touch, which_touch_inla), ])
plot(gmap, add = TRUE)
plot(gmap[the_max, ], add = TRUE, col = "blue")
plot(gmap[which_touch, ], add = TRUE, col = "yellow")
plot(gmap[which_touch_inla, ], add = TRUE, col = "#FF000030")

c(the_max, setdiff(which_touch, which_touch_inla))

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

data("Germany", package = "INLA")
germany <- list(
  Y = as.numeric(Germany$Y),
  E = as.numeric(Germany$E),
  region = as.integer(Germany$region),
  adj = adj,
  prec = list(
    Q = Q_scaled,
    log_det = log_det_gen,
    rank = nrow(adj) - 1L
  ),
  map = terra::wrap(gmap)
)

save(germany, file = "adlaplace/data/germany.rda", compress = "xz")

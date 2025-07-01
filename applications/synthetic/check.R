
res <- readRDS("applications/synthetic/coverage.rds")
l <- sapply(res, length) |> max()
res <- sapply(res, \(x) c(x, rep(NA, length(x)-l))) |> rowMeans()

plot(res)
abline(h=.8, lty=2)


M <- matrix(res, nrow=3, byrow=T)
plot(M[1,])
abline(h=.8, lty=2)

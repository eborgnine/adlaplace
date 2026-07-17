test_that("ad_shards partitions a modest design matrix", {
  set.seed(11)
  A <- matrix(rnorm(120), nrow = 30, ncol = 4)
  G <- ad_shards(A, num_shards = 3L)
  expect_s4_class(G, "dgCMatrix")
  expect_equal(nrow(G), 30L)
  expect_gte(ncol(G), 1L)
  expect_lte(ncol(G), 3L)
})

test_that("ad_shards source uses max(dim(ATp)) for large-matrix warning", {
  src_path <- normalizePath(
    file.path(testthat::test_path(".."), "..", "R", "ad_shards.R"),
    mustWork = FALSE
  )
  skip_if_not(file.exists(src_path), "package source not available")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  expect_match(src, "max\\(dim\\(ATp\\)\\)\\s*>\\s*1e5")
  expect_match(src, "max\\(dim\\(ATp\\)\\)\\s*>\\s*100")
})


test_that("vignettes/test.Rmd renders", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("adlaplace")
  vignette <- testthat::test_path("../../vignettes/test.Rmd")
  skip_if_not(file.exists(vignette), "vignette file missing")
  output <- tempfile(fileext = ".html")
  expect_error(
    rmarkdown::render(
      vignette,
      output_file = output,
      quiet = TRUE,
      envir = new.env(parent = globalenv())
    ),
    NA
  )
})

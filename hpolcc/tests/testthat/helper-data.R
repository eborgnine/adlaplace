#' Build minimal synthetic data for hpolcc integration tests
make_hpolcc_test_data <- function(n_strata = 6L, n_per = 4L) {
  set.seed(0)
  n <- n_strata * n_per
  data <- data.frame(
    count = rpois(n, lambda = 2),
    hum = rnorm(n),
    region = rep(1:2, length.out = n),
    date = rep(seq_len(n_strata), each = n_per)
  )
  list(
    data = data,
    formula = count ~ hum + f(date, model = "iid"),
    cc_design = hpolcc:::ccDesign(time_var = "date", strat_vars = "region")
  )
}

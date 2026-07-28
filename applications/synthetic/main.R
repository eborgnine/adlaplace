# Packages (load from repo when developing) --------------------------------
if (file.exists("adlaplace/DESCRIPTION")) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all("adlaplace", quiet = TRUE)
    devtools::load_all("adlaplaceHgp", quiet = TRUE)
    devtools::load_all("hpolcc", quiet = TRUE)
  }
}
library(adlaplace)
library(adlaplaceHgp)
library(hpolcc)
library(data.table)
library(ggplot2)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(5)
}

# Data --------------------------------------------------------------------
genPM <- function(date) {
  l <- length(date) + 1
  date_int <- as.integer(date)
  date_int <- c(date_int[1] - 1, date_int)
  pm <- 10
  pm <- pm + 5 * sin(2 * pi * (date_int + rnorm(l, 0, 1)) / 365.25)
  pm <- pm + rnorm(l, 0, 1)
  pm <- (2 * pm[-l] + pm[-1]) / 3
  pmax(0, pm)
}

genHum <- function(date) {
  l <- length(date) + 1
  date_int <- as.integer(date)
  date_int <- c(date_int[1] - 1, date_int)
  hum <- 50
  hum <- hum + 5 * cos(2 * pi * (date_int + rnorm(l, 0, 1)) / 365.25)
  hum <- hum + rnorm(l, 0, 1)
  hum <- (2 * hum[-l] + hum[-1]) / 3
  hum <- pmax(0, hum)
  hum / max(hum)
}

region_effect1 <- rnorm(10, 1, 0.1)
region_effect2 <- rnorm(10, 1, 0.1)
genPmEffect <- function(pm, r1, r2) {
  r1 * pm / 10 + pm^(r2 / 3)
}
genCount <- function(hum, pm, r1 = 1, r2 = 1) {
  l <- length(pm)
  rpois(l, exp(hum + genPmEffect(pm, r1, r2)))
}

data <- lapply(1:2, \(i) {
  data <- data.table(date = as.Date(1:5000), region = i)
  data$hum <- genHum(data$date)
  data$pm <- genPM(data$date)
  data$count <- genCount(data$hum, data$pm, region_effect1[i], region_effect2[i])
  data
}) |> rbindlist()

data$monthDow <- format(data$date, "%Y-%m-%a")

# Model -------------------------------------------------------------------
knots_pm <- seq(0, ceiling(max(data$pm) / 5) * 5, 2.5)
ref_values <- list(pm = 10)

# Global PM smooth (iwp). For regional hierarchical PM, use hiwp(..., by = "region")
# once the hiwp AD stack is stable for your sample size.
cc_formula <- dirichlet_multinom(
  count,
  by = c("region", "monthDow")
) ~ hum + iwp(
  pm,
  p = 2,
  ref_value = 10,
  knots = knots_pm
)

fit <- hnlm(
  formula = cc_formula,
  data = data,
  config = list(
    transform_theta = TRUE,
    num_threads = 4L,
    num_shards = 100L,
    num_sim = 200L,
    verbose = FALSE
  ),
  control = list(maxit = 200, trace = 1),
  control_inner = list(maxit = 100, report.level = 0)
)

print(fit)
summary(fit)

# PM effect curves (conditional simulation) -------------------------------
xseq <- seq(knots_pm[1], rev(knots_pm)[1], by = 0.5)
newx <- list(pm = data.frame(pm = xseq, hum = 0))

sim <- adlaplaceHgp::cond_sim_iwp(
  fit = fit$extra,
  model_data = list(
    terms = fit$call$terms,
    term_data = list(info = fit$info)
  ),
  newx = newx,
  n = 500,
  probs = c(0.1, 0.5, 0.9)
)

pm_sim <- sim$sim$pm
pm_x <- sim$x$pm
pm_q <- sim$quantiles$common$pm

plot_df <- data.frame(
  var_value = rep(pm_x, times = ncol(pm_sim)),
  effect_value = as.vector(pm_sim),
  draw = rep(seq_len(ncol(pm_sim)), each = length(pm_x))
)

ggplot(plot_df, aes(x = var_value, y = effect_value, group = draw)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1, col = "gray75") +
  geom_line(alpha = 0.15) +
  geom_line(
    data = data.frame(
      var_value = pm_x,
      effect_value = pm_q[, "0.5"]
    ),
    linewidth = 1
  ) +
  labs(x = "pm", y = "RR", title = "Population-average PM effect (median)")

# True generating effects by region -------------------------------------
true_df <- do.call(
  rbind,
  lapply(unique(data$region), function(i) {
    x <- genPmEffect(xseq, region_effect1[i], region_effect2[i])
    k <- which(xseq == ref_values$pm)
    data.frame(
      variable = "pm",
      var_value = xseq,
      group = i,
      true = x - x[k]
    )
  })
)

ggplot(true_df, aes(x = var_value, y = exp(true), colour = as.factor(group))) +
  theme_bw() +
  geom_line() +
  labs(x = "pm", y = "RR", title = "Simulated regional PM effects (generating model)")

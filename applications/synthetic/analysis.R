# Posterior effect curves (run after main.R has created `fit`, `data`, `knots_pm`, `ref_values`)

library(adlaplaceHgp)

xseq <- seq(knots_pm[1], rev(knots_pm)[1], by = 0.01)
newx <- list(pm = expand.grid(
  pm = xseq,
  hum = 0,
  region = unique(data$region)[1]
))

sim <- cond_sim_iwp(
  fit = fit$laplace,
  model_data = fit$model_data,
  newx = newx,
  n = 1000,
  probs = c(0.01, 0.5, 0.99)
)

pm_q <- sim$quantiles$common$pm
res1 <- data.frame(
  variable = "pm",
  var_value = sim$x$pm,
  q_0.01 = pm_q[, "0.01"],
  q_0.5 = pm_q[, "0.5"],
  q_0.99 = pm_q[, "0.99"]
)

if (!is.null(sim$quantiles$group) && !is.null(sim$quantiles$group$pm)) {
  res1 <- do.call(
    rbind,
    lapply(names(sim$quantiles$group$pm), function(g) {
      gq <- sim$quantiles$group$pm[[g]]
      data.frame(
        variable = "pm",
        var_value = sim$x$pm,
        group = g,
        q_0.01 = gq[, "0.01"],
        q_0.5 = gq[, "0.5"],
        q_0.99 = gq[, "0.99"]
      )
    })
  )
}

for (i in unique(res1$group)) {
  x <- genPmEffect(res1$var_value[res1$group == i], region_effect1[i], region_effect2[i])
  k <- which(res1$var_value[res1$group == i] == ref_values$pm)
  res1$true[res1$group == i] <- x - x[k]
}

library(ggplot2)
ggplot(res1[res1$variable == "pm", ], aes(x = var_value, y = q_0.5, colour = as.factor(group))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col = "gray75") +
  geom_line() +
  geom_ribbon(aes(ymin = q_0.01, ymax = q_0.99, fill = as.factor(group)), col = NA, alpha = 0.3) +
  geom_line(aes(y = true), linetype = 3, colour = "black") +
  facet_grid(~variable, scales = "free_x")

res1 <- constructEffect(fit, exposure_var = "pm", 
                        group_var = "region", group = unique(data$region), 
                        values = seq(knots_pm[1],rev(knots_pm)[1],.01),
                        ref_values = ref_values,
                        pars = samp, probs = c(.01,.5,.99))
res2 <- constructEffect(fit, exposure_var = "pm", 
                        group_var = "region", group = NA, 
                        values = knots_pm[1]:rev(knots_pm)[1],
                        ref_values = ref_values,
                        pars = samp, probs = c(.01,.5,.99))

for(i in unique(res1$group)){
  x <- genPmEffect(res1$var_value[res1$group == i], region_effect1[i], region_effect2[i])
  k <- which(res1$var_value[res1$group == i] == ref_values$pm)
  res1$true[res1$group == i] <- x - x[k]
}


ggplot(res1[res1$variable == "pm",], aes(x=var_value, y=q_0.5, colour=as.factor(group))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, col="gray75") +
  geom_line() +
  geom_ribbon(aes(ymin = q_0.01, ymax = q_0.99, fill = as.factor(group)), col=NA, alpha=.3) + 
  geom_line(data=res1[res2$variable == "pm",][,-4], aes(y=true, group=group), linetype=3, colour="black") +
  geom_line(data=res2[res2$variable == "pm",][,-4], linetype=2, colour="black") +
  facet_grid(~variable, scales = "free_x")

fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "beta"]
fit$theta_info
fit$obj$env$last.par.best[names(fit$obj$env$last.par.best) == "theta"]



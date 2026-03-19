res = readRDS("~/research/healthcanada/pkg/cancc/vignettes/res.rds")
devtools::load_all("~/research/adlaplace/hpolcc", recompile=FALSE)
sim = cond_sim_iwp(
  res
)

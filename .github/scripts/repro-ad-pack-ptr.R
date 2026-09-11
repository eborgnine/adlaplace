#!/usr/bin/env Rscript
# Isolate the Windows hard-abort path seen after five gamSim prints in
# test-format-parameters.R ("ad_pack from ptr has empty info").
# Breadcrumbs flush so the last printed line is the failing call.

say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%OS3"), " ", paste(..., collapse = " "))
  cat(msg, "\n", sep = "")
  flush.console()
  invisible(NULL)
}

quit_ok <- function() {
  say("DONE ok")
  quit(save = "no", status = 0L)
}

say("start", R.version.string)
say("platform", R.version$platform)

if (!requireNamespace("mgcv", quietly = TRUE)) {
  say("SKIP: mgcv not installed")
  quit_ok()
}
if (!requireNamespace("adlaplace", quietly = TRUE)) {
  stop("adlaplace is not installed", call. = FALSE)
}

say("library(adlaplace)", as.character(utils::packageVersion("adlaplace")))
library(adlaplace)

say("gamSim(6)")
dat <- NULL
utils::capture.output(dat <- mgcv::gamSim(6, n = 80, scale = 0.2, dist = "poisson"))
stopifnot(!is.null(dat))

say("model_data")
md <- adlaplace::model_data(
  adlaplace::nbinom(y, lower = 1e-9) ~ x1 + adlaplace::iid(fac, init = 0.25),
  data = dat,
  verbose = FALSE
)

say("ad_pack(md) warmup")
af_md <- adlaplace::ad_pack(
  md,
  config = list(
    transform_theta = TRUE,
    num_shards = 4L,
    num_threads = 1L,
    verbose = FALSE
  )
)
say("ad_pack(md) ok", length(methods::slotNames(af_md)))
rm(af_md)
gc(FALSE)

config <- list(
  beta = md$term_data$info$beta$init,
  theta = md$term_data$info$theta$init,
  gamma = rep(0, nrow(md$term_data$info$gamma)),
  transform_theta = TRUE,
  verbose = FALSE
)
config$theta <- adlaplace::apply_theta_log(md$term_data$info$theta, cols = "init")$init
shards <- unname(c(md$observations, md$random, md$parameters))
say("n_shards", length(shards))

ptrs <- vector("list", length(shards))
for (i in seq_along(shards)) {
  sh <- shards[[i]]
  dens <- tryCatch(sh@density, error = function(e) class(sh)[1L])
  say(sprintf("ad_pack_ptr[%d/%d] density=%s", i, length(shards), dens))
  ptrs[[i]] <- adlaplace::ad_pack_ptr(sh, config = config)
  say(sprintf("ad_pack_ptr[%d] ok", i))
}

say("combine via do.call(c, ptrs)")
ptr <- do.call(c, ptrs)
say("combine ok")

say("ad_pack(ptr)")
af <- adlaplace::ad_pack(ptr)
say("ad_pack(ptr) ok", length(af@info))
stopifnot(identical(af@info, list()))

quit_ok()

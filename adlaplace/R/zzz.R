.onLoad <- function(libname, pkgname) {
  .Call(`_adlaplace_register_callables`)
}

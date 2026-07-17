# Compile smoke test: verify RCppAD headers work with R CMD SHLIB.
# Run during R CMD check (all .R files in tests/ are sourced).

has_cpp_compiler <- function() {
  cxx <- Sys.getenv("CXX", unset = NA_character_)
  if (is.na(cxx) || !nzchar(cxx)) {
    cxx <- Sys.which("g++")
  }
  if (!nzchar(cxx)) {
    cxx <- Sys.which("clang++")
  }
  nzchar(cxx)
}

rcppad_include_dir <- function() {
  # Source tree during development
  if (file.exists("DESCRIPTION") &&
      file.exists(file.path("inst", "include", "cppad", "cppad.hpp"))) {
    return(normalizePath(file.path("inst", "include"), mustWork = TRUE))
  }
  # Installed package: headers live under include/, not inst/include/
  installed <- system.file("include", package = "RCppAD", mustWork = FALSE)
  if (nzchar(installed) &&
      file.exists(file.path(installed, "cppad", "cppad.hpp"))) {
    return(normalizePath(installed, mustWork = TRUE))
  }
  stop("Could not locate RCppAD include directory", call. = FALSE)
}

run_cppad_compile_smoke <- function() {
  include <- rcppad_include_dir()

  tmp <- tempfile("rcppad_smoke_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  writeLines(
    c(
      "#include <cppad/cppad.hpp>",
      "extern \"C\" void rcppad_smoke(void) {",
      "  CppAD::AD<double> x = 1.0;",
      "  CppAD::AD<double> y = x * x;",
      "  (void) y;",
      "}"
    ),
    file.path(tmp, "smoke.cpp")
  )

  # PKG_*FLAGS via Makevars in the temp dir (reliable across platforms).
  writeLines(
    c(
      paste0("PKG_CPPFLAGS=-I", include),
      "PKG_CXXFLAGS=-std=c++17"
    ),
    file.path(tmp, "Makevars")
  )

  old_wd <- setwd(tmp)
  on.exit(setwd(old_wd), add = TRUE)

  # Clear R_TESTS: during R CMD check it is set to startup.Rs (relative),
  # and the child R would fail to source it from this tempdir.
  old_r_tests <- Sys.getenv("R_TESTS", unset = NA_character_)
  Sys.unsetenv("R_TESTS")
  on.exit({
    if (is.na(old_r_tests)) {
      Sys.unsetenv("R_TESTS")
    } else {
      Sys.setenv(R_TESTS = old_r_tests)
    }
  }, add = TRUE)

  status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "smoke.cpp"),
    stdout = TRUE,
    stderr = TRUE
  )

  so_candidates <- list.files(
    tmp,
    pattern = "\\.(so|dll|dylib)$",
    full.names = TRUE
  )
  if (length(so_candidates) == 0L) {
    stop(
      "RCppAD compile smoke failed (no shared library produced).\n",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

if (!has_cpp_compiler()) {
  message("Skipping RCppAD compile smoke: no C++ compiler found")
} else {
  run_cppad_compile_smoke()
  message("RCppAD compile smoke: OK")
}

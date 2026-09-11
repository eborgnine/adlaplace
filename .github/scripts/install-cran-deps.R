#!/usr/bin/env Rscript
# Install the union of Depends/Imports/LinkingTo/Suggests from local packages,
# excluding the local packages themselves (those are built/installed later in CI).
# Transitive Suggests are not installed (dependencies = NA) so unavailable CRAN
# Suggests of our Suggests do not break the solve.
# geostatsp is installed from R-universe (newer than CRAN); see below.

# Optional CLI args restrict which local package DESCRIPTIONs are scanned
# (e.g. Rscript install-cran-deps.R RCppAD adlaplace).
cli_locals <- commandArgs(trailingOnly = TRUE)
locals <- if (length(cli_locals)) {
  cli_locals
} else {
  c(
    "RCppAD",
    "adlaplace",
    "adlaplaceExample",
    "adlaplaceHgp",
    "adlaplaceFem"
  )
}

dep_types <- c("Depends", "Imports", "LinkingTo", "Suggests")

parse_dep_field <- function(raw) {
  if (is.na(raw) || !nzchar(raw)) {
    return(character())
  }
  parts <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  names <- sub("\\s*\\(.*\\)$", "", parts)
  names <- trimws(names)
  names[nzchar(names) & names != "R"]
}

pkgs <- character()
for (dir in locals) {
  desc <- file.path(dir, "DESCRIPTION")
  if (!file.exists(desc)) {
    stop("DESCRIPTION not found: ", desc, call. = FALSE)
  }
  dcf <- read.dcf(desc)
  for (ty in dep_types) {
    if (ty %in% colnames(dcf)) {
      pkgs <- c(pkgs, parse_dep_field(dcf[1, ty]))
    }
  }
}

base_pkgs <- rownames(installed.packages(priority = "base"))
pkgs <- sort(unique(setdiff(pkgs, c(locals, base_pkgs))))
message(
  "Installing ", length(pkgs), " CRAN deps (local + base packages excluded):\n  ",
  paste(pkgs, collapse = ", ")
)

pak_lib <- Sys.getenv("R_LIB_FOR_PAK", unset = "")
if (nzchar(pak_lib)) {
  library(pak, lib.loc = pak_lib)
} else {
  library(pak)
}

# macOS only: binaries, no source. New terra (and similar) releases often
# reach CRAN before macOS binaries exist; compiling them needs Homebrew GDAL.
# Windows already uses binaries; Linux uses P3M binaries + can compile.
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  plat <- tryCatch(pak::system_r_platform(), error = function(e) "macos")
  options(pkg.platforms = plat)
  options(pkgType = "binary")
  options(install.packages.compile.from.source = "never")
  message("macOS: binaries only (pkg.platforms = ", plat, ")")
}

# adlaplaceFem needs geostatsp newer than CRAN. Use the package-specific
# R-universe repo so other CRAN packages (terra, ...) stay on CRAN binaries.
runiverse_geostatsp <- "https://eborgnine.r-universe.dev/geostatsp"
if ("geostatsp" %in% pkgs) {
  message("Installing geostatsp from R-universe: ", runiverse_geostatsp)
  old_repos <- getOption("repos")
  options(repos = c(geostatsp = runiverse_geostatsp, old_repos))
  pak::pkg_install("geostatsp", dependencies = NA)
  options(repos = old_repos)
  pkgs <- setdiff(pkgs, "geostatsp")
}

# Hard deps only for transitive packages. Our Suggests remain direct targets.
if (length(pkgs)) {
  pak::pkg_install(pkgs, dependencies = NA)
}
message("CRAN dependency install complete.")

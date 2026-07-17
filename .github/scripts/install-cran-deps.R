#!/usr/bin/env Rscript
# Install the union of Depends/Imports/LinkingTo/Suggests from local packages,
# excluding the local packages themselves (those are built/installed later in CI).
# Transitive Suggests are not installed (dependencies = NA) so unavailable CRAN
# Suggests of our Suggests (e.g. geostatsp -> RandomFields) do not break the solve.

locals <- c(
  "RCppAD",
  "adlaplace",
  "adlaplaceExample",
  "adlaplaceHgp",
  "adlaplaceGrf",
  "hpolcc"
)

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

# Hard deps only for transitive packages: keeps geostatsp installable without
# its archived Suggests (RandomFields). Our Suggests remain direct targets.
pak::pkg_install(pkgs, dependencies = NA)
message("CRAN dependency install complete.")

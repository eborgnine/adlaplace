#!/usr/bin/env Rscript
# Install the union of Depends/Imports/LinkingTo/Suggests from local packages,
# excluding the local packages themselves (those are built/installed later in CI).

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

pkgs <- sort(unique(setdiff(pkgs, locals)))
message(
  "Installing ", length(pkgs), " CRAN deps (local packages excluded):\n  ",
  paste(pkgs, collapse = ", ")
)

pak_lib <- Sys.getenv("R_LIB_FOR_PAK", unset = "")
if (nzchar(pak_lib)) {
  library(pak, lib.loc = pak_lib)
} else {
  library(pak)
}

pak::pkg_install(pkgs, dependencies = TRUE)
message("CRAN dependency install complete.")

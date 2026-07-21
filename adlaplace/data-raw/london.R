#' Rebuild the bundled London mortality / PM10 dataset
#'
#' Downloads the Tobias et al. case-crossover example CSV and writes
#' \code{adlaplace/data/london.rda}.
#'
#' Run from the repository root:
#' \code{source("adlaplace/data-raw/london.R")}
#'
#' @keywords internal
NULL

url <- "https://raw.githubusercontent.com/aureliotobias/casecrossover/main/london.csv"
london <- as.data.frame(data.table::fread(url))
london$date <- as.Date(ISOdate(london$year, london$month, london$day))
london$dow <- factor(format(london$date, "%a"),
  levels = format(ISOdate(2000, 1, seq(from = 3, length.out = 7)), "%a")
)
london$month <- factor(format(london$date, "%b"),
  levels = format(ISOdate(2000, 1:12, 1), "%b")
)

dir.create("adlaplace/data", recursive = TRUE, showWarnings = FALSE)
save(london, file = "adlaplace/data/london.rda", compress = "xz")

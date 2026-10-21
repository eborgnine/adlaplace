#' @importFrom KFAS KFS SSMcustom SSModel
NULL

#' Initial state variance used by CAPTAIN's Kalman filter.
#'
#' `kalmanfis` documents a default initial covariance of `1e5 * I`.
#' @noRd
CAPTAIN_P0 <- 1e5

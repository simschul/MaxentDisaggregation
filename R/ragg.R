#' Title
#'
#' @param n
#' @param mean
#' @param sd
#' @param min
#' @param max
#'
#' @return
#' @export
#'
#' @examples
ragg <- function(n, mean, sd = NULL, min = 0, max = Inf, meanlog = NULL, sdlog = NULL) {
  if (!is.null(mean) & !is.null(sd) & min == -Inf & max == Inf) {
    # Normal distribution
    rnorm(n, mean, sd)
  } else if (!is.null(mean) & !is.null(sd)) {
    # Truncated normal
    rtruncnorm(n, a = min, b = max, mean = mean, sd = sd)
  } else if (!is.null(mean) & is.null(sd) & min == 0 & max == Inf) {
    # Exponential
    rexp(n, rate = 1 / mean)
  } else if (is.null(mean) & is.null(sd) & is.finite(min) & is.finite(max)) {
    # uniform
    runif(n, min = min, max = max)
  } else if (is.null(mean) & is.null(sd) & !is.null(meanlog) & !is.null(sdlog)) {
    rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  } else {
    stop('Case not implemented atm.')
  }
}

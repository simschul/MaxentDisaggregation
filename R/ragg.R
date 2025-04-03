#'
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
ragg <- function(n, mean, sd = NULL, min = 0, max = Inf, meanlog = NULL, sdlog = NULL,
                 log = FALSE) {

  #TODO: check mean, sd other args
  if (!is.null(mean) & !is.null(sd) & min == -Inf & max == Inf & isFALSE(log)) {
    # Normal distribution
    rnorm(n, mean, sd)
  } else if (!is.null(mean) & !is.null(sd) & min <= 0 & max == Inf & isTRUE(log)) {
    # Lognormal distribution
    rlnorm2(n, mean, sd)
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


#' Random number generation for the lognormal distribution parameterised to the `mean` and `sd` on the **linear** scale.
#' Taken from: https://stackoverflow.com/questions/56821688/sample-a-lognormal-distribution-to-an-exact-mean-and-sd
#'
#' @param n
#' @param mean the mean on the linear scale
#' @param sd the standard deviation on the linear scale
#'
#' @return
#' @export
#'
#' @examples
rlnorm2 <- function(n, mean, sd) {
  m2 <-  log(mean^2 / sqrt(sd^2 + mean^2))
  sd2 <- sqrt(log(1 + (sd^2 / mean^2)))
  r <- c(scale(rnorm(n)))
  return(exp(sd2*r+m2))
}

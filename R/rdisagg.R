#' Generate Random Disaggregates
#'
#' Creates a random sample of disaggregates based on user-specified information
#' about the total aggregate and the shares. Internally, it first samples
#' the aggregate (via [ragg()]) and then samples the shares (via [rshares()]).
#' The final output is the element-wise product of these two samples.
#'
#' @param n Integer. Sample size (number of random draws).
#' @param mean_0 Numeric. The best guess (mean) of the aggregate on the *linear* scale
#'   unless `log = TRUE`.
#' @param sd_0 Numeric. The standard deviation of the aggregate on the linear scale
#'   unless `log = TRUE`. Defaults to `NULL`.
#' @param min Numeric. Lower boundary for the aggregate distribution. Defaults to `0`.
#' @param max Numeric. Upper boundary for the aggregate distribution. Defaults to `Inf`.
#' @param shares Numeric vector of best-guess shares. Must sum to 1 if no \code{NA}
#'   values are present. See [rshares()] for details.
#' @param sds Numeric vector of standard deviations for the shares, or `NULL` if none
#'   are available (default).
#' @param log Logical. If `TRUE`, the aggregate distribution will be sampled in log-space
#'   (lognormal or truncated lognormal, etc.) depending on the parameters. Defaults to `FALSE`.
#'
#' @return A numeric matrix of size \code{n} times \code{length(shares)}, where each row
#'   sums to the corresponding entry in the aggregate sample. Specifically, each row
#'   is the element-wise product of the sampled aggregate (scalar in each row) and
#'   the sampled shares (vector in each row).
#'
#' @seealso
#' - [ragg()] for details on how the aggregate sample is generated.
#' - [rshares()] for details on how the shares sample is generated.
#'
#' @export
#'
#' @examples
#' # Example 1: Simple exponential aggregate, uniform shares
#' set.seed(123)
#' agg_mean <- 100
#' n_samples <- 5
#' shares_vec <- c(0.3, 0.4, 0.3)
#'
#' rdisagg(n = n_samples, mean_0 = agg_mean, shares = shares_vec)
#'
#' # Example 2: Normal aggregate, Dirichlet shares w/ MaxEnt
#' rdisagg(n = 5, mean_0 = 100, sd_0 = 10, shares = c(0.1, 0.2, 0.7))
rdisagg <- function(n, mean_0, sd_0 = NULL, min = 0, max = Inf,
                    shares, sds = NULL, log = FALSE) {

  sample_agg <- ragg(
    n    = n,
    mean = mean_0,
    sd   = sd_0,
    min  = min,
    max  = max,
    log  = log
  )

  sample_shares <- rshares(n = n, shares = shares, sds = sds)
  sample_disagg <- sample_shares * sample_agg
  return(sample_disagg)
}

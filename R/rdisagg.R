#' Generate random disaggregates
#'
#' Creates a random sample of disaggregates based on the information provided.
#' The aggregate and the shares are sampled independently. The distribution
#' from which to sample is determined internally based on the information
#' provided by the user.
#'
#' @param n sample size
#' @param mean_0 the best guess of the aggregate
#' @param sd_0 the standard deviation of the aggregate
#' @param min the lower boundary
#' @param max the upper boundary
#' @param shares best guesses for the shares (need to sum to one unless they contain NA's). See `?rshares`
#' @param sds standard deviations of the shares (can contain NA). Set to NULL (default) if none available.
#'
#' @return
#' @export
#'
#' @examples
rdisagg <- function(n, mean_0, sd_0 = NULL, min = 0, max = Inf,
                    shares, sds = NULL, log = FALSE) {
  sample_agg <- ragg(n = n, mean = mean_0, sd = sd_0, min = min, max = max, log = log)
  sample_shares <- rshares(n = n, shares = shares, sds = sds)
  sample_disagg <- sample_shares * sample_agg

  #names(sample_disagg) <- names(alpha)
  return(sample_disagg)
}

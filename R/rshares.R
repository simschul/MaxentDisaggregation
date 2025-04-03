#' Generate random numbers from Dirichlet distribution. The concenctration parameter
#' $\gamma$ (or sometimes also noted $\alpha_0$) is estimated based on the
#' the maximum entropy principle.
#'
#' See: `?find_gamma_maxent2`
#'
#' @param n number of samples
#' @param shares Vector containing a best-guess (mean) values for the shares. Must sum to 1!
#' @param ... other parameters passed to `find_gamma_maxent`
#'
#' @return a matrix with n rows and length(shares) cols, each containing a single Dirichlet random deviate
#' @export
#'
#' @examples
rdir_maxent <- function(n, shares, ...) {
  out <- find_gamma_maxent(shares, eval_f = dirichlet_entropy, ...)
  sample <- rdirichlet(n, shares = shares, gamma = out$solution)
  colnames(sample) <- names(shares)
  return(sample)
}



#' Generate uniform random number from a Dirichlet distribution.
#' Uniform means each alpha is 1:https://en.wikipedia.org/wiki/Dirichlet_distribution#When_each_alpha_is_1
#'
#'
#' @param n number of samples
#' @param length the number of variables
#' @param names optional: a character vector of length `length` used as column names for the output sample
#'
#' @return a matrix with n rows and `length` cols, each containing a single Dirichlet random deviate
#' @export
#'
#' @examples
rdir_uniform <- function(n, length, names = NULL) {
  sample <- rdirichlet_standard(n, rep(1, length))
  colnames(sample) <- names
  return(sample)
}



#' Generate uniform random number from a Generalised Dirichlet distribution.
#'
#' @param n sample size
#' @param shares Vector containing a best-guess (mean) values for the shares. Must sum to 1!
#' @param sds Vector of same length as `shares` containing the SDs of the shares.
#'
#' @references Plessis, Sylvain, Nathalie Carrasco, and Pascal Pernot. “Knowledge-Based Probabilistic Representations of Branching Ratios in Chemical Networks: The Case of Dissociative Recombinations.” The Journal of Chemical Physics 133, no. 13 (October 7, 2010): 134110. https://doi.org/10.1063/1.3479907.
#'
#' @return
#' @export
#'
#' @examples
rdir_generalised <- function(n, shares, sds) {

  if (!isTRUE(all.equal(names(shares), names(sds)))) {
    stop('shares and sds need to have the same column names. can also be both NULL')
  }

  alpha2 <- (shares / sds) ^ 2
  beta2 <- shares / (sds) ^ 2
  k <- length(alpha2)
  x <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    x[, i] <- rgamma(n, shape = alpha2[i], rate = beta2[i])
  }
  sample <- x / rowSums(x)
  colnames(sample) <- names(shares)
  return(sample)
}



#' Create Random Dirichlet distributed numbers.
#' Adds the gamma parameter compared to the standard Dir variant.
#' For all variables `i` with shares below `threshold`, the `rate` argument of
#' the `rgamma` functions is set to `1 / shares[i]`, and `shape` is set `1`.
#' This ensures less extreme values
#' for those variables which presumably originate from the know issues of
#' generating random gamma distributed number with a very small `shape` parameter.
#'
#' The reason for setting rate and shape parameters for small alphas like this is the following:
#' The gamma distribution (similar to the beta distribution), generates well behaved
#' random numbers only if the variables SD is not greater than the MEAN.
#' The gamma distribution can also be parametrised by a mean and sd using the
#' following relation between mean/sd and the shape and rate parameters:
#' $$
#' rate = \frac{\mu}{\sigma^2} \\
#' shape = (\frac{\mu}{\sigma})^2
#' $$
#' (derived from: https://en.wikipedia.org/wiki/Gamma_distribution#Mean_and_variance)
#'
#'
#' See `?rgamma`: "Note that for smallish values of shape (and moderate scale) a
#' large parts of the mass of the Gamma distribution is on values of x so
#' near zero that they will be represented as zero in computer arithmetic.
#' So rgamma may well return values which will be represented as zero."
#'
#' Note: the current implementation is rather a pragmatic workaround for this
#' issue. Both the threshold parameter and the formula to set the `rate` parameter
#' are somehow arbitrarily chosen. However, since only variables with a very small
#' mean are affected this shouldn't affect overall results too much.
#'
#' @param n
#' @param shares vector containing a best-guess (mean) values for the shares. Must sum to 1!
#' @param gamma
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
rdirichlet <- function(n, shares, gamma, threshold = 1E-2) {
  alpha <- gamma * shares
  l <- length(alpha)
  rate <- rep(1, l)
  #rate[alpha < threshold] <- sqrt(1/alpha[alpha < threshold])
  rate[alpha < threshold] <- 1 / alpha[alpha < threshold]
  alpha[alpha < threshold] <- 1
  x <- matrix(rgamma(l * n, alpha, rate), ncol = l, byrow = TRUE)
  # x <- matrix(rgamma2(n, alpha, rate), ncol = l, byrow = FALSE)
  sm <- rowSums(x)
  sample <- x/sm
  colnames(sample) <- names(alpha)
  return(sample)
}

#' Functions to generate random deviates from the Dirichlet distribution.
#'
#' Like `gtools::rdirichlet`, but with an adjustment when `alpha` contains small
#' values below a given `threshold`.
#'
#' Basically the same as `rdirichlet`, but keeping the parametrisation of `gtools::rdirichlet`
#' to allow easy switching between the two.
#'
#' Details see `?rdirichlet`.
#'
#' @param n
#' @param alpha
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
rdirichlet_standard <- function(n, alpha, threshold = 1E-2) {
  l <- length(alpha)
  rate <- rep(1, l)
  rate[alpha < threshold] <- 1 / alpha[alpha < threshold]
  alpha[alpha < threshold] <- 1
  x <- matrix(rgamma(l * n, alpha, rate), ncol = l, byrow = TRUE)
  sm <- rowSums(x)
  sample <- x/sm
  colnames(sample) <- names(alpha)
  return(sample)
}




#' Generate random shares based on the information provided.
#' Samples alyways sum to 1.
#'
#' Compared to generalised Dirichlet the function can handle mixed information,
#' e.g. SDs/means only for some of the shares.
#'
#' The function uses a nested sampling approach:
#' 1. sample those shares for which both, mean and SD is available (either from generalised Dirichlet `dirg()` if mean and sd are given for ALL shares, or from beta distiributions `rbeta3()` when mean and sd are only given for some shares.)
#' 2. sample reminder from standard dirichlet with gamma fitted to MaxEnt (`rdir_maxent`)
#'
#' @param n sample size
#' @param shares vector containing a best-guess (mean) values for the shares. Can contain NA's if no information is available.
#' @param sds Vector of same length as `shares` containing the SDs of the shares. Can contain NA's if no information is available.
#' @param na_action what to do with NA's in the `shares` argument. "fill" (default) fills them proportional with the difference `1-sum(shares)`. "remove" removes them from the vector.
#' @param max_iter see `?rbeta3`
#' @param ... other argeuments passed to `rbeta3`
#'
#' @return
#' @export
#'
#' @examples
rshares <- function(n, shares, sds = NULL,
                    na_action = 'fill', max_iter = 1E3, ...) {
  if (is.null(sds)) {
    sds <- rep(NA, length(shares))
    names(sds) <- names(shares)
  }
  if (na_action == 'remove') {
    sds <- sds[!is.na(shares)]
    shares <- shares[!is.na(shares)]
  } else if (na_action == 'fill') {
    shares[is.na(shares)] <- (1 - sum(shares, na.rm = TRUE)) / length(shares[is.na(shares)])
  } else {
    stop('na_action must be either "remove" or "fill"!')
  }

  if (!isTRUE(all.equal(sum(shares), 1))) {
    stop('shares must sum to one! If you have NAs in your shares consider setting "na_action" to "fill". ')
  }

  K <- length(shares)
  have_mean_only <- is.finite(shares) & !is.finite(sds)
  have_sd_only <- is.finite(sds) & !is.finite(shares)
  have_both <- is.finite(shares) & is.finite(sds)

  if (all(have_both)) {
    # Generalised dirichlet
    sample <- rdir_generalised(n, shares, sds)
  } else if (all(have_mean_only)) {
    # Dirichlet with maxent fitted gamma
    sample <- rdir_maxent(n, shares)
  } else {
    # partial info
    stop("Please either provide positive numeric values for all elements in `sds` for for which you also provided positive numeric values for `shares`. Or provide no `sds` at all (can be also a vector of NA's). ")
  }
  return(sample)
}


# rshares(100, c(3,3,3)) %>% boxplot
# rshares(100, c(0.1, 0.3, 0.6)) %>% boxplot
# rshares(100, c(0.1, 0.3, 0.6), c(0.1, 0.01, 0.2)) %>% boxplot
# rshares(100, shares = c(0.1, 0.3, 0.6), sds = c(0.1, 0.3, 0.003)) %>% boxplot
# rshares(100, shares = c(NA, NA, 0.6), sds = c(NA, NA,0.1)) %>% boxplot


#' Generate random numbers from the beta distribution. Similar to `base::rbeta`, but:
#' - vectorized
#' - parametrised with mean and sd (derived from: https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_variance)
#' - Ensures that the sum of the samples do not exceed 1.
#'
#' It is similar to a generalised Dirichlet, however the main functional difference
#' is that the samples do NOT sum to 1 (but are <= 1).
#'
#' The constraint that the sample sum is <= 1 is currently achieved by sampling
#' as long until `n` samples exist that satisfy the constraint. Depending on the
#' parameter combination this can be quite computational intensive.
#'
#' Note that the beta distribution is only defined if (mean * (1-mean)) < sd^2
#'
#' @param n number of samples
#' @param shares a vector of shares (between 0 and 1, sum(means) <= 1)
#' @param sds vector of standard devaitions of the shares.
#' @param fix if set to TRUE (default) the SD of all variables for which the
#' parameter combination of shares and sd is undedefined ((shares * (1-shares)) < sd^2).
#' If set to FALSE the function will throw an error message if there are undefined
#' parameter combinations.
#' @param max_iter the maximum number of iterations to sample.
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
rbeta3 <- function(n, shares, sds, fix = TRUE, max_iter = 1E3) {
  var <- sds^2
  undef_comb <- (shares * (1-shares)) <= var # beta distribution is not defined
  if (!all(!undef_comb)) {
    if (isTRUE(fix)) {
      var[undef_comb] <- (shares[undef_comb] * (1 - shares[undef_comb])) - 1E-1
      #shares[undef_comb]^2
    } else stop('The beta distribution is not defined for the
                                    parameter combination you provided!
                                    sd must be smaller or equal sqrt(shares*(1-shares))')

  }
  # from: https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_variance
  alpha <- shares * (((shares * (1-shares)) / var) - 1)
  beta <- (1-shares)*(((shares * (1-shares)) / var) - 1)
  k<- length(shares)
  x <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    x[, i] <- rbeta(n, alpha[i], beta[i])
  }
  larger_one <- rowSums(x) > 1
  count <- 0
  while(sum(larger_one) > 0) {
    for (i in 1:k) {
      x[larger_one, i] <- rbeta(sum(larger_one), alpha[i], beta[i])
    }
    larger_one <- rowSums(x) > 1
    count <- count + 1
    if (count > max_iter) stop('max_iter is reached. the combinations of shares
    and sds you provided does allow to generate `n` random samples that are
                               not larger than 1. Either increase max_iter, or
                               change parameter combination. ')
  }
  colnames(x) <- names(shares)
  return(x)
}


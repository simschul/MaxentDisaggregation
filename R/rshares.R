#' Generate random numbers from a Dirichlet distribution via Maximum Entropy
#'
#' This function estimates the concentration parameter \eqn{\gamma} (sometimes noted
#' \eqn{\alpha_0}) by applying the maximum entropy principle via
#' [find_gamma_maxent()]. The resulting parameter is then used to generate random
#' vectors from a Dirichlet distribution.
#'
#' For details on how \eqn{\gamma} is fitted, see [find_gamma_maxent()] or the
#' references therein. The entropy is evaluated using [dirichlet_entropy()].
#'
#' @param n Integer. Number of random vectors to generate.
#' @param shares Numeric vector. A best-guess (mean) of the shares. Must sum to 1.
#' @param ... Other parameters passed on to [find_gamma_maxent()].
#'
#' @return A numeric matrix of dimension \code{n} times \code{length(shares)}, where
#' each row is a random deviate from the Dirichlet distribution.
#'
#' @seealso [find_gamma_maxent()], [dirichlet_entropy()], [rdirichlet_standard()]
#' @export
#'
#' @examples
#' rdirichlet_maxent(5, c(0.2, 0.3, 0.5))
rdirichlet_maxent <- function(n, shares, ...) {
  out <- find_gamma_maxent(shares, eval_f = dirichlet_entropy, ...)
  sample <- rdirichlet(n, shares = shares, gamma = out$solution)
  colnames(sample) <- names(shares)
  return(sample)
}


#' Generate uniform random numbers from a Dirichlet distribution
#'
#' Generates random vectors from a Dirichlet distribution with each \eqn{\alpha = 1}.
#' This implies a uniform distribution over the simplex.
#'
#' For details, see
#' [Wikipedia on the Dirichlet distribution](https://en.wikipedia.org/wiki/Dirichlet_distribution#When_each_alpha_is_1).
#'
#' @param n Integer. Number of random vectors to generate.
#' @param length Integer. The number of variables (columns) in the output.
#' @param names Optional character vector of length \code{length} giving column names.
#'
#' @return A numeric matrix of dimension \code{n} times \code{length}, with each row
#' a Dirichlet random deviate drawn from a uniform (each \eqn{\alpha=1}) distribution.
#'
#' @seealso [rdirichlet_standard()]
#' @export
#'
#' @examples
#' rdirichlet_uniform(5, 3)
rdirichlet_uniform <- function(n, length, names = NULL) {
  sample <- rdirichlet_standard(n, rep(1, length))
  colnames(sample) <- names
  return(sample)
}


#' Generate random numbers from a Generalised Dirichlet distribution
#'
#' Generates random vectors from a Generalised Dirichlet distribution using the mean
#' (`shares`) and standard deviations (`sds`) of each component. The shape and rate
#' parameters (\eqn{\alpha_2} and \eqn{\beta_2}) are derived from the provided means
#' and SDs, following the approach of Plessis et al. (2010).
#'
#' Please note that unrealistic combinations of `shares` and `sds`, in particular
#' components with `sds/shares > 0.5`, might result in biased samples.
#' Against this background, after the draw, the function checks the sample mean of every
#' component against the target share.  If the *relative bias*
#' \deqn{\text{bias}_i = \frac{\bar x_i - \mu_i}{\mu_i}}
#' exceeds `max_rel_bias`, a warning is issued.
#' If you really think that the uncertainty of your `shares` is very high
#' (in the order of a coefficient of variation of 0.5 or more), consider using
#' the Maximum Entropy Dirichlet distribution instead by setting `sds=NULL` (the default).
#' This way the uncertainties are fitted such that they are maximally high (for all components at the same time).
#'
#' @param n Integer. Number of random vectors to generate.
#' @param shares Numeric vector of best-guess (mean) values for the shares. Must sum to 1.
#' @param sds Numeric vector (same length as \code{shares}) containing the standard
#'   deviations of the shares.
#' @param max_rel_bias Numeric scalar (default `0.1`).  Threshold for the
#'   largest acceptable relative bias (10 % by default) between the
#'   simulated mean and the target share.  Components whose simulated bias
#'   exceeds this value are listed in a warning
#'
#' @references
#' Plessis, Sylvain, Nathalie Carrasco, and Pascal Pernot.
#' "Knowledge-Based Probabilistic Representations of Branching Ratios in Chemical Networks:
#' The Case of Dissociative Recombinations." *The Journal of Chemical Physics* 133, no. 13 (2010): 134110.
#' doi:{10.1063/1.3479907}
#'
#' @return A numeric matrix of dimension \code{n} times \code{length(shares)}, with each
#' row a single Generalised Dirichlet random deviate.
#'
#' @seealso [rgamma()], [rdirichlet()]
#' @export
#'
#' @examples
#' rdirichlet_generalised(5, c(0.2, 0.3, 0.5), c(0.05, 0.07, 0.06))
rdirichlet_generalised <- function(n, shares, sds, max_rel_bias = 0.1) {
  if (!isTRUE(all.equal(names(shares), names(sds)))) {
    stop("`shares` and `sds` must have the same names. They can also both be NULL.")
  }

  alpha2 <- (shares / sds) ^ 2
  beta2  <- shares / (sds) ^ 2
  k <- length(alpha2)
  x <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    x[, i] <- rgamma(n, shape = alpha2[i], rate = beta2[i])
  }
  sample <- x / rowSums(x)
  colnames(sample) <- names(shares)

  # bias check
  rel_bias <- (colMeans(sample) - shares) / shares
  idx_bad  <- which(rel_bias > max_rel_bias)
  if (length(idx_bad) > 0) {
    raw_names <- names(rel_bias)[idx_bad]                    # may be NULL
    labels    <- ifelse(is.null(raw_names) || raw_names == "",
                        sprintf("shares[%d]", idx_bad),      # fallback
                        raw_names)                           # proper name
    detail_lines <- sprintf("%s: bias of %.3f", labels, rel_bias[idx_bad])

    warning(paste0(
      "The samples of the following components show a relative bias of > ",
      max_rel_bias, ":\n",
      paste(detail_lines, collapse = "\n"),
      "\nMost likely this is because the standard deviations you provided are too high for those components. You might decrease the SD of those components, or if you really think that the uncertainty on that components is that high consider using the Maximum Entropy Dirichlet distribution instead by setting `sds = NULL`. "
    ), call. = FALSE)
  }


  return(sample)
}


#' Generate random numbers from a Dirichlet distribution (adjusted for small alpha)
#'
#' Extends a standard Dirichlet random generator by adding a concentration parameter
#' (\eqn{\gamma}).
#' For each variable \eqn{i} whose mean value (\eqn{\alpha_i = \gamma \cdot share_i})
#' is below a specified `threshold`, a fallback parametrization of the Gamma distribution
#' is applied to avoid zero or near-zero sampling. This is especially useful for very
#' small shape parameters, which can cause numerical issues in [rgamma()].
#'
#' This pragmatic workaround sets:
#'
#' - \eqn{\alpha_i = 1} (shape) for shares below `threshold`
#' - \eqn{rate = 1 / \alpha_i}
#'
#' ensuring less extreme values.
#'
#' For more details, see the discussion in [rgamma()] under "small shape values" and
#' the references there. This approach helps mitigate issues where numeric precision
#' can push small Gamma-distributed values to zero.
#'
#' @param n Integer. Sample size.
#' @param shares Numeric vector of best-guess (mean) values for the shares. Must sum to 1.
#' @param gamma Numeric. Concentration parameter (\eqn{\gamma}).
#' @param threshold Numeric. All values in \code{shares} multiplied by \code{gamma}
#'   that are below this threshold use the fallback parametrization (default: \code{1e-2}).
#'
#' @return A numeric matrix of dimension \code{n} times \code{length(shares)}, with
#'   each row a single Dirichlet random deviate.
#'
#' @seealso [rdirichlet_standard()], [rgamma()]
#' @export
#'
#' @examples
#' rdirichlet(5, c(0.2, 0.3, 0.5), gamma = 10, threshold = 1e-2)
rdirichlet <- function(n, shares, gamma, threshold = 1e-2) {
  alpha <- gamma * shares
  l <- length(alpha)
  rate <- rep(1, l)
  # If alpha[i] < threshold, set alpha[i] = 1 and rate[i] = 1 / alpha[i]
  # to avoid extremely small shape parameters.
  rate[alpha < threshold] <- 1 / alpha[alpha < threshold]
  alpha[alpha < threshold] <- 1

  x <- matrix(rgamma(l * n, alpha, rate), ncol = l, byrow = TRUE)
  sm <- rowSums(x)
  sample <- x / sm
  colnames(sample) <- names(alpha)
  return(sample)
}


#' Generate random numbers from a Dirichlet distribution (standard parametrization)
#'
#' Similar to [gtools::rdirichlet()], but when \eqn{\alpha_i} is below `threshold`,
#' a fallback parametrization is used to avoid zero or near-zero outcomes.
#'
#' This function keeps the standard Dirichlet parametrization (\eqn{\alpha_i}) but
#' adjusts small \eqn{\alpha_i} to 1 and uses \eqn{rate = 1/\alpha_i} for the
#' underlying Gamma draws. For details, see [rdirichlet()].
#'
#' @param n Integer. Number of random vectors to generate.
#' @param alpha Numeric vector of shape parameters.
#' @param threshold Numeric. Threshold below which fallback parameters are used.
#'   Defaults to \code{1e-2}.
#'
#' @return A numeric matrix of dimension \code{n} times \code{length(alpha)}, with
#'   each row a Dirichlet random deviate.
#'
#' @seealso [gtools::rdirichlet()], [rdirichlet()]
#' @export
#'
#' @examples
#' rdirichlet_standard(5, c(0.1, 0.1, 0.8))
rdirichlet_standard <- function(n, alpha, threshold = 1e-2) {
  l <- length(alpha)
  rate <- rep(1, l)
  rate[alpha < threshold] <- 1 / alpha[alpha < threshold]
  alpha[alpha < threshold] <- 1

  x <- matrix(rgamma(l * n, alpha, rate), ncol = l, byrow = TRUE)
  sm <- rowSums(x)
  sample <- x / sm
  colnames(sample) <- names(alpha)
  return(sample)
}


#' Generate random shares that always sum to 1
#'
#' This function combines different approaches for generating shares (Dirichlet,
#' Generalised Dirichlet, Beta) in a single workflow, depending on the information
#' provided by `shares` and `sds`. All generated vectors sum to 1.
#'
#' The logic is:
#'
#' 1. If both mean and SD are available for **all** shares, use [rdirichlet_generalised()].
#' 2. If only mean values (no SDs) are available for **all** shares, use [rdirichlet_maxent()].
#' 3. If partial information in the mean values is present (i.e. there are \code{NA}'s in `shares`)
#' the output depends on `na_action`: By default, if `na_action` is "fill", any \code{NA} in `shares` is replaced
#' proportionally so that the new `shares` still sum to 1. If `na_action` is "remove",
#' elements with \code{NA} are dropped.
#' 4. Otherwise, throw an error if partial information is present (e.g., some but not all
#'    SDs are defined).
#'
#'
#'
#' @param n Integer. Sample size.
#' @param shares Numeric vector containing a best-guess (mean) values for the shares.
#'   May contain \code{NA}, which is handled depending on `na_action`. Must sum to 1
#'   after adjustment.
#' @param sds Numeric vector (same length as `shares`) for the standard deviations of
#'   the shares. Can also contain \code{NA}.
#' @param na_action Character. Either `"fill"` (default) to fill \code{NA} shares
#'   proportionally, or `"remove"` to drop them entirely.
#' @param max_iter Numeric. Used by [rbeta2()] if partial Beta sampling is performed
#'   (not used here unless you modify the source).
#' @param max_rel_bias Numeric scalar in (0, 1).  Maximum *relative* deviation
#'   that the internal moment–matching routine is allowed to have between the
#'   requested means in `shares` and the theoretical means implied by the
#'   fitted (generalised) Dirichlet parameters.  It is computed as
#'   \code{max(abs(mu_target - mu_fitted) / mu_target)} across all components. Defaults to 0.1.
#' @param ... Additional arguments passed on to [rbeta2()] or other internal functions.
#'
#' @return A numeric matrix of dimension \code{n} times \code{length(shares)}, with
#'   each row summing to 1.
#'
#' @seealso [rdirichlet_generalised()], [rdirichlet_maxent()], [rbeta2()]
#' @export
#'
#' @examples
#' # Example with no SDs -> uses rdirichlet_maxent
#' rshares(5, c(0.1, 0.3, 0.6))
#'
#' # Example with means & SDs for all -> uses rdirichlet_generalised
#' rshares(5, c(0.1, 0.3, 0.6), c(0.05, 0.07, 0.02))
#'
#' # 'na_action' example: fill missing shares to sum to 1
#' rshares(5, c(NA, 0.2, 0.3) , na_action = "fill")
rshares <- function(n, shares, sds = NULL,
                    na_action = "fill", max_iter = 1e3, max_rel_bias = 0.1, ...) {

  if (is.null(sds)) {
    sds <- rep(NA, length(shares))
    names(sds) <- names(shares)
  }

  if (na_action == "remove") {
    sds    <- sds[!is.na(shares)]
    shares <- shares[!is.na(shares)]
  } else if (na_action == "fill") {
    shares[is.na(shares)] <-
      (1 - sum(shares, na.rm = TRUE)) / length(shares[is.na(shares)])
  } else {
    stop('na_action must be either "remove" or "fill"!')
  }

  if (!isTRUE(all.equal(sum(shares), 1))) {
    stop("shares must sum to 1! If you have NAs, consider setting na_action = 'fill'.")
  }
  K <- length(shares)

  have_mean_only <- is.finite(shares) & !is.finite(sds)
  have_both      <- is.finite(shares) & is.finite(sds)

  if (all(have_both)) {
    # All shares have corresponding SDs: use generalised Dirichlet
    sample <- rdirichlet_generalised(n, shares, sds, max_rel_bias = max_rel_bias)
  } else if (all(have_mean_only)) {
    # All shares have means, no SDs: use Dirichlet with maxent for gamma
    sample <- rdirichlet_maxent(n, shares)
  } else {
    # nested Dirichlet approach
    sample <- .rshares_nested(
      n, shares, sds,
      max_rel_bias   = max_rel_bias,
      max_iter_bias_fix = 20,        # <- tweak if needed
      max_iter_rbeta2   = max_iter
    )
  }
  return(sample)
}

# Helper for the nested Dirichlet logic with iterative bias correction
.rshares_nested <- function(n, shares, sds, max_rel_bias = 0.10,
                            max_iter_bias_fix = 20, max_iter_rbeta2 = 1e3) {

  iter <- 0
  repeat {
    iter <- iter + 1
    if (iter > max_iter_bias_fix) {
      stop("Bias could not be reduced below ", max_rel_bias,
           " after ", max_iter_bias_fix, " iterations.\n",
           "Consider providing smaller SDs or more NA means.")
    }

    have_both      <- is.finite(shares) & is.finite(sds)
    have_mean_only <- is.finite(shares) & !is.finite(sds)

    K      <- length(shares)
    sample <- matrix(0, nrow = n, ncol = K)
    colnames(sample) <- names(shares)

    ## 1) Components with mean *and* SD  -->  Beta-truncated sampling
    if (any(have_both)) {
      sample[, have_both] <- rbeta2(
        n,
        shares = shares[have_both],
        sds    = sds[have_both],
        max_iter = max_iter_rbeta2
      )
    }

    ## 2) Components with mean only --> MaxEnt Dirichlet (rescaled afterwards)
    if (any(have_mean_only)) {
      alpha2        <- shares[have_mean_only] / sum(shares[have_mean_only])
      sample_dm     <- rdirichlet_maxent(n, alpha2)
      sample[, have_mean_only] <- sample_dm * (1 - rowSums(sample))
    }

    ## 3) Bias check on the “have_both” block
    if (any(have_both)) {
      emp_mean  <- colMeans(sample[, have_both, drop = FALSE])
      rel_bias  <- abs(emp_mean - shares[have_both]) / shares[have_both]

      if (all(rel_bias <= max_rel_bias)) {
        return(sample)                     # Success – exit the loop
      }

      # → At least one component is too far off: mark its SD as NA
      sds[have_both][rel_bias > max_rel_bias] <- NA
      warning(
        "Relative bias exceeded ", max_rel_bias,
        " for component(s): ",
        paste(names(rel_bias)[rel_bias > max_rel_bias], collapse = ", "),
        ".  Their SDs have been set to NA and will be handled with ",
        "a Maximum-Entropy Dirichlet on the next iteration.",
        call. = FALSE
      )
    } else {
      # No components left with SDs – nothing more to fix
      return(sample)
    }
  }
}

# .rdirichlet_nested <- function(n, shares, sds) {
#   K <- length(shares)
#   have_both      <- is.finite(shares) & is.finite(sds)
#
#   # nested Dirichlet approach
#   sample <- matrix(0, nrow = n, ncol = K)
#   colnames(sample) <- names(shares)
#
#   repeat{
#     # sample all shares with both mean + sd
#     sample[, have_both] <- rbeta2(n, shares = shares[have_both],
#                                   sds = sds[have_both], max_iter = max_iter)
#     # check possible bias due to rejection
#     emp_mean <- colMeans(sample[,have_both, drop = FALSE])
#     diff_mean  <- abs(emp_mean - shares[have_both]) / shares[have_both]
#     bias <- any(diff_mean > max_rel_bias)
#     if (isFALSE(bias)) {
#       break
#     } else {
#       # set problematic SDs to NA (so that Maxent dirichlet is used for those)
#       sds[have_both][diff_mean > max_rel_bias] <- NA
#       sample <- .rdirichlet_nested(n, shares, sds)
#     }
#   }
#
#     return(sample)
#
# }



#' Generate random numbers from the Beta distribution via mean and SD
#'
#' This function is similar to [stats::rbeta()], but vectorized and parametrized
#' directly in terms of mean and standard deviation. It also ensures that the sum
#' of sampled values across variables does not exceed 1 by discarding and
#' re-sampling any draws whose row-sum is > 1. Depending on the parameter
#' configuration, this can be computationally expensive.
#'
#' The Beta distribution is only well-defined if \eqn{\mathrm{Var}(X) \leq \mu(1-\mu)},
#' i.e., if \eqn{\mathrm{sd}(X)^2 \leq \mu(1-\mu)}.
#'
#' @param n Integer. Number of samples to generate.
#' @param shares Numeric vector of means (between 0 and 1). Must satisfy `sum(shares) <= 1`.
#' @param sds Numeric vector of standard deviations of the same length as `shares`.
#' @param fix Logical. If `TRUE` (default), any undefined parameter combinations
#'   (where \eqn{\mathrm{sd}^2 > \mu(1-\mu)}) are reset to a slightly smaller
#'   variance to make them valid. If `FALSE`, an error is thrown.
#' @param max_iter Integer. Maximum number of re-sampling iterations allowed when
#'   discarding samples whose row-sum is > 1.
#'
#' @return A numeric matrix of size \code{n} times \code{length(shares)}. Each row
#'   has \eqn{\sum_i x_i \leq 1}, but they do not necessarily sum exactly to 1.
#'
#' @seealso [stats::rbeta()]
#' @export
#'
#' @examples
#' # Means that sum to less than 1
#' rbeta2(5, c(0.2, 0.3, 0.4), c(0.05, 0.04, 0.03))
rbeta2 <- function(n, shares, sds, fix = TRUE, max_iter = 1e2) {
  var <- sds^2
  undef_comb <- (shares * (1 - shares)) <= var  # Beta undefined if var > mu*(1-mu)

  if (!all(!undef_comb)) {
    if (isTRUE(fix)) {
      var[undef_comb] <- (shares[undef_comb] * (1 - shares[undef_comb])) - 1e-1
    } else {
      stop(
        "The beta distribution is not defined for the parameter combination you provided!\n",
        "sd must satisfy: sd^2 <= shares * (1 - shares)."
      )
    }
  }

  # Derive alpha and beta from mean and variance
  alpha <- shares * (((shares * (1-shares)) / var) - 1)
  beta  <- (1 - shares) * (((shares * (1-shares)) / var) - 1)

  k <- length(shares)
  x <- matrix(0, nrow = n, ncol = k)
  for (i in seq_len(k)) {
    x[, i] <- rbeta(n, alpha[i], beta[i])
  }

  # Re-sample any rows whose sum exceeds 1
  larger_one <- rowSums(x) > 1
  count <- 0
  while (sum(larger_one) > 0) {
    for (i in seq_len(k)) {
      x[larger_one, i] <- rbeta(sum(larger_one), alpha[i], beta[i])
    }
    larger_one <- rowSums(x) > 1
    count <- count + 1
    if (count > max_iter) {
      stop(
        "`max_iter` reached. The combination of shares and sds does not allow ",
        "generating `n` valid samples that satisfy sum(x_i) <= 1."
      )
    }
  }

  colnames(x) <- names(shares)
  return(x)
}


# -------------------------------------------------------------------- #
# Helper: warn if empirical mean / sd differ too much from targets     #
# -------------------------------------------------------------------- #
.check_sample <- function(sample, target_means, target_sds,
                          thr_mean = 0.10, thr_sd = 0.20) {

  emp_mean <- colMeans(sample)
  emp_sd   <- apply(sample, 2, sd)

  ## mean check (only where target is finite & > 0)
  valid_mean <- is.finite(target_means) & target_means > 0
  diff_mean  <- abs(emp_mean[valid_mean] - target_means[valid_mean]) /
    target_means[valid_mean]

  if (any(diff_mean > thr_mean)) {
    idx <- which(diff_mean > thr_mean)
    warning(
      "Sample mean deviates by more than ",
      thr_mean * 100, "% from target for components: ",
      paste(names(target_means[valid_mean])[idx], collapse = ", ")
    )
  }

  ## sd check
  valid_sd  <- is.finite(target_sds) & target_sds > 0
  diff_sd   <- abs(emp_sd[valid_sd] - target_sds[valid_sd]) /
    target_sds[valid_sd]

  if (any(diff_sd > thr_sd)) {
    idx <- which(diff_sd > thr_sd)
    warning(
      "Sample SD deviates by more than ",
      thr_sd * 100, "% from target for components: ",
      paste(names(target_sds[valid_sd])[idx], collapse = ", ")
    )
  }
}

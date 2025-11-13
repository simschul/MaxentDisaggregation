#' Aggregate Random Number Generator
#'
#' Generates random numbers from a variety of distributions, chosen by the
#' combination of provided parameters:
#'
#' - **Normal** distribution (if `mean` and `sd` are given, `min = -Inf`, `max = Inf`, and `log = FALSE`)
#' - **Lognormal** distribution (if `mean` and `sd` are given, `min <= 0`, `max = Inf`, and `log = TRUE`)
#' - **Truncated Normal** distribution (if `mean` and `sd` are given, and finite `min`, `max` not both infinite)
#' - **Exponential** distribution (if `mean` is given, `sd` is `NULL`, and `min = 0`, `max = Inf`)
#' - **Uniform** distribution (if both `mean` and `sd` are `NULL`, and finite `min` and `max`)
#' - **Lognormal** via meanlog/sdlog** (if `meanlog` and `sdlog` are not `NULL` and `mean`/`sd` are `NULL`)
#'
#' If none of these distribution signatures is met, an error is raised.
#'
#' @param n Integer. Number of random values to generate.
#' @param mean Numeric. Mean of the distribution on the *linear* scale (for Normal,
#'   Truncated Normal, Lognormal parameterization, or Exponential).
#'   Set to `NULL` if not used.
#' @param sd Numeric. Standard deviation of the distribution on the *linear* scale
#'   (for Normal, Truncated Normal, or Lognormal parameterization).
#'   Set to `NULL` if not used.
#' @param min Numeric. Lower bound (for Truncated Normal or Uniform). Defaults to \code{0}.
#' @param max Numeric. Upper bound (for Truncated Normal or Uniform). Defaults to \code{Inf}.
#' @param meanlog Numeric. Mean on the log scale for a standard lognormal draw (when
#'   `mean`/`sd` are `NULL`). Defaults to `NULL`.
#' @param sdlog Numeric. Standard deviation on the log scale for a standard lognormal
#'   draw (when `mean`/`sd` are `NULL`). Defaults to `NULL`.
#' @param log Logical. If `TRUE`, interpret the distribution as lognormal
#'   (given `mean` and `sd` on the linear scale). Defaults to `FALSE`.
#'
#' @return
#' A numeric vector of length \code{n}, drawn from the specified distribution.
#' The distribution choice is determined by the argument combination.
#'
#' @seealso [rnorm()], [rlnorm()], [rlnorm2()], [rtruncnorm()], [rexp()], [runif()]
#'
#' @export
#'
#' @examples
#' # 1) Normal distribution
#' x1 <- ragg(10, mean = 5, sd = 1)
#'
#' # 2) Lognormal distribution (via linear-scale mean and sd)
#' x2 <- ragg(10, mean = 1, sd = 0.3, log = TRUE)
#'
#' # 3) Truncated Normal distribution
#' x3 <- ragg(10, mean = 5, sd = 1, min = 3, max = 7)
#'
#' # 4) Exponential distribution
#' x4 <- ragg(10, mean = 2, min = 0)
#'
#' # 5) Uniform distribution
#' x5 <- ragg(10, min = 0, max = 10)
#'
#' # 6) Standard lognormal (using meanlog, sdlog)
#' x6 <- ragg(10, meanlog = 0, sdlog = 1)
ragg <- function(n, mean = NULL, sd = NULL, min = -Inf, max = Inf, meanlog = NULL, sdlog = NULL,
                 log = FALSE) {

  # Case 1: Normal distribution
  if (!is.null(mean) & !is.null(sd) & min == -Inf & max == Inf & isFALSE(log)) {
    return(rnorm(n, mean, sd))

    # Case 2: Lognormal distribution (via linear-scale mean/sd)
  } else if (!is.null(mean) & !is.null(sd) & min <= 0 & max == Inf & isTRUE(log)) {
    return(rlnorm2(n, mean, sd))

    # Case : Truncated normoal with moment matching (if lower bound = 0 and no upper bound)
  } else if (!is.null(mean) & !is.null(sd) & min <= 0 & max == Inf & isFALSE(log)) {
    out <- tnorm_params_from_moments(mean, sd)
    mean2 <- out$mu
    sd2 <- out$sigma
    return(rtruncnorm(n, a = min, b = max, mean = mean2, sd = sd2))

    # Case 3: Truncated Normal (without moment matching)
  } else if (!is.null(mean) & !is.null(sd)) {
    
    return(rtruncnorm(n, a = min, b = max, mean = mean, sd = sd))


    # Case 4: Exponential
  } else if (!is.null(mean) & is.null(sd) & min == 0 & max == Inf) {
    return(rexp(n, rate = 1 / mean))

    # Case 5: Uniform
  } else if (is.null(mean) & is.null(sd) & is.finite(min) & is.finite(max)) {
    return(runif(n, min = min, max = max))

    # Case 6: Standard lognormal (meanlog, sdlog given)
  } else if (is.null(mean) & is.null(sd) & !is.null(meanlog) & !is.null(sdlog)) {
    return(rlnorm(n, meanlog = meanlog, sdlog = sdlog))

  } else {
    stop("No matching distribution case found. Check your arguments or see ?ragg.")
  }
}


#' Lognormal Random Number Generation with Mean and SD on the Linear Scale
#'
#' Generates lognormal random values with a specified mean and standard deviation
#' given on the *linear* scale. Internally, this is done by converting the linear-scale
#' mean and SD to the appropriate log-scale parameters for [rlnorm()].
#'
#' This approach is based on the relationship:
#'
#' \deqn{
#'   m_2 = \log\!\bigl(\frac{\mu^2}{\sqrt{\sigma^2 + \mu^2}}\bigr), \quad
#'   s_2 = \sqrt{\log\!\bigl(1 + \frac{\sigma^2}{\mu^2}\bigr)},
#' }
#'
#' where \eqn{\mu} is the mean (linear scale) and \eqn{\sigma} the SD (linear scale),
#' while \eqn{m_2} and \eqn{s_2} are the mean and SD on the *log scale*.
#'
#' @param n Integer. Number of random values to generate.
#' @param mean Numeric. The mean of the desired distribution on the linear scale (must be > 0).
#' @param sd Numeric. The standard deviation of the desired distribution on the linear scale.
#'
#' @return
#' A numeric vector of length \code{n}, containing lognormal random deviates.
#'
#' @seealso [rlnorm()], [ragg()]
#'
#' @references
#' This parameterization is discussed in [this StackOverflow post](https://stackoverflow.com/questions/56821688).
#'
#' @export
#'
#' @examples
#' # Generate lognormal samples with mean=1, sd=0.3 (on linear scale)
#' set.seed(123)
#' x <- rlnorm2(10, mean = 1, sd = 0.3)
#' summary(x)
rlnorm2 <- function(n, mean, sd) {
  m2  <- log(mean^2 / sqrt(sd^2 + mean^2))
  sd2 <- sqrt(log(1 + (sd^2 / mean^2)))
  r   <- rnorm(n)
  return(exp(sd2 * r + m2))
}

#' Truncated normal (lower=0) parameters from mean m and sd s
#'
#' Given target moments (m, s) of X ~ Normal(mu, sigma^2) truncated at [0, +inf),
#' return (mu, sigma). Robust and efficient.
#'
#' @param m target mean (m > 0)
#' @param s target standard deviation (s > 0)
#' @param tol absolute tolerance for root finding
#' @param maxit maximum iterations for Newton fallback
#' @return list(mu=..., sigma=..., alpha=..., method="uniroot"|"newton")
tnorm_params_from_moments <- function(m, s, tol = 1e-12, maxit = 50) {
  if (!is.finite(m) || !is.finite(s) || m <= 0 || s <= 0) {
    stop("m and s must be finite and strictly positive.")
  }

  # Stable inverse Mills ratio: lambda(a) = phi(a) / (1 - Phi(a))
  lambda <- function(a) {
    # log phi and log upper tail to avoid under/overflow
    log_phi <- dnorm(a, log = TRUE)
    log_tail <- pnorm(a, lower.tail = FALSE, log.p = TRUE)
    exp(log_phi - log_tail)
  }

  # Derivative: lambda'(a) = lambda(a) * (lambda(a) - a)
  lambda_prime <- function(a, lam = NULL) {
    if (is.null(lam)) lam <- lambda(a)
    lam * (lam - a)
  }

  # Root function in alpha (a): g(a) = s^2 (lam - a)^2 - m^2 (1 + a*lam - lam^2)
  gfun <- function(a) {
    lam <- lambda(a)
    s^2 * (lam - a)^2 - m^2 * (1 + a * lam - lam^2)
  }

  # Derivative g'(a)
  gprime <- function(a) {
    lam <- lambda(a)
    l1  <- lam * (lam - a) # lambda'(a)
    2 * s^2 * (lam - a) * (l1 - 1) - m^2 * (lam + a * l1 - 2 * lam * l1)
  }

  # Try to bracket a root for uniroot (safe & fast). Expand symmetrically.
  # We cap |alpha| to ~37 to stay well within numeric stability of pnorm's log tail.
  cap <- 37
  L <- -8; U <- 8
  gL <- gfun(L); gU <- gfun(U)

  expand_attempts <- 0
  while (sign(gL) == sign(gU) && expand_attempts < 10) {
    L <- max(-cap, 2 * L)
    U <- min( cap, 2 * U)
    gL <- gfun(L); gU <- gfun(U)
    if ((L <= -cap && U >= cap)) break
    expand_attempts <- expand_attempts + 1
  }

  # If we got a sign change, use uniroot (preferred).
  if (is.finite(gL) && is.finite(gU) && sign(gL) != sign(gU)) {
    r <- uniroot(function(a) gfun(a), lower = L, upper = U, tol = tol)
    a  <- r$root
    lam <- lambda(a)
    sigma <- m / (lam - a)
    mu    <- -a * sigma
    return(list(mu = as.numeric(mu),
                sigma = as.numeric(sigma),
                alpha = as.numeric(a),
                method = "uniroot",
                iters = r$iter,
                bracket = c(L, U)))
  }

  # Fallback: damped Newton starting from a0 = 0 (half-normal pivot).
  a <- 0
  ga <- gfun(a)
  for (k in 1:maxit) {
    gpa <- gprime(a)
    if (!is.finite(ga) || !is.finite(gpa)) break
    if (gpa == 0) break
    step <- -ga / gpa

    # Damping if the step doesn't reduce |g|
    new_a <- a + step
    new_ga <- gfun(new_a)
    damp <- 0
    while ((!is.finite(new_ga) || abs(new_ga) > abs(ga)) && damp < 20) {
      step <- step / 2
      new_a <- a + step
      new_ga <- gfun(new_a)
      damp <- damp + 1
    }

    a <- new_a
    ga <- new_ga
    if (abs(ga) <= tol) break
    # Keep alpha within safe numeric range
    if (a >  cap) a <-  cap
    if (a < -cap) a <- -cap
  }

  lam <- lambda(a)
  sigma <- m / (lam - a)
  mu    <- -a * sigma

  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    stop("Failed to solve for (mu, sigma). Check inputs (m, s).")
  }

  list(mu = as.numeric(mu),
       sigma = as.numeric(sigma),
       alpha = as.numeric(a),
       method = "newton",
       iters = k)
}


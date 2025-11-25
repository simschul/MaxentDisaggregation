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
  } else if (!is.null(mean) & !is.null(sd) & min == 0 & max == Inf & isTRUE(log)) {
    return(rlnorm2(n, mean, sd))

    # Case 3: Truncated normoal with moment matching (if lower bound = 0 and no upper bound)
  } else if (!is.null(mean) & !is.null(sd)) {
    out <- tnorm_ab_params_from_moments(mean, sd, a = min, b = max)
    mean2 <- out$mu
    sd2 <- out$sigma
    return(rtruncnorm(n, a = min, b = max, mean = mean2, sd = sd2))

    
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

#' Maximum-Entropy Truncated Normal Random Generation
#'
#' Draws random values from the **maximum-entropy** distribution on the interval
#' \eqn{[a, b]} subject to a fixed mean and standard deviation on the *truncated*
#' scale.
#'
#' Under the constraints:
#' \itemize{
#'   \item support restricted to \eqn{[a, b]} (with \eqn{a \le b}),
#'   \item target mean \eqn{E[X] = \code{mean}},
#'   \item target standard deviation \eqn{\mathrm{sd}(X) = \code{sd}},
#' }
#' the maximum-entropy distribution is a Normal law, possibly truncated:
#' \itemize{
#'   \item if \code{a = -Inf} and \code{b = Inf}, the solution is the ordinary
#'     Normal \eqn{X \sim \mathcal{N}(\code{mean}, \code{sd}^2)};
#'   \item otherwise, the solution is a **truncated normal**, defined by some
#'     parent Normal \eqn{Y \sim \mathcal{N}(\mu, \sigma^2)} truncated to
#'     \eqn{[a, b]}, where \eqn{(\mu, \sigma)} are chosen such that the truncated
#'     distribution has mean \code{mean} and sd \code{sd}.
#' }
#'
#' This function:
#' \enumerate{
#'   \item uses \code{\link{tnorm_ab_params_from_moments}} to solve for the
#'     underlying Gaussian parameters \eqn{(\mu, \sigma)} compatible with
#'     \code{mean}, \code{sd}, \code{a}, and \code{b};
#'   \item generates random values using \code{truncnorm::rtruncnorm()} on
#'     \eqn{[a, b]} with mean \eqn{\mu} and sd \eqn{\sigma};
#'   \item optionally checks the **empirical** mean and sd of the sample against
#'     the requested \code{mean} and \code{sd} via \code{.check_sample()}, and
#'     throws an error if the deviations exceed user-specified thresholds.
#' }
#'
#' In particular, when the requested pair \code{(mean, sd)} is **not feasible**
#' for any truncated normal on \eqn{[a, b]} (for example, demanding a very large
#' sd near a boundary), the parameter solver may not be able to match the target
#' moments. With \code{check = TRUE} (the default), such cases are detected by
#' comparing the empirical moments of the generated sample to the targets, and
#' a descriptive error is raised.
#'
#' @param n Integer. Number of random values to generate.
#' @param a Numeric scalar. Lower truncation bound (can be \code{-Inf} for no
#'   lower truncation).
#' @param b Numeric scalar. Upper truncation bound (can be \code{Inf} for no
#'   upper truncation). Must satisfy \code{a <= b}.
#' @param mean Numeric scalar. Target mean of the distribution on the truncated
#'   domain \eqn{[a, b]}.
#' @param sd Numeric scalar. Target standard deviation of the distribution on
#'   the truncated domain \eqn{[a, b]}. Must be strictly positive.
#' @param check Logical. If \code{TRUE} (default), the function validates the
#'   generated sample using \code{.check_sample()}, comparing the empirical mean
#'   and sd to the requested \code{mean} and \code{sd}. If the relative
#'   deviations exceed \code{thr_mean} and/or \code{thr_sd}, an error is thrown.
#'   If \code{FALSE}, no post-hoc validation is performed.
#' @param thr_mean Numeric scalar. Relative tolerance for the sample mean. A
#'   value of \code{0.10} means that empirical means deviating by more than
#'   10\% from the target are considered unacceptable (when \code{check = TRUE}).
#' @param thr_sd Numeric scalar. Relative tolerance for the sample standard
#'   deviation. A value of \code{0.20} means that empirical standard deviations
#'   deviating by more than 20\% from the target are considered unacceptable
#'   (when \code{check = TRUE}).
#'
#' @return
#' A numeric vector of length \code{n}, containing random draws from the
#' maximum-entropy distribution consistent with the supplied support
#' \code{[a, b]} and target moments \code{(mean, sd)}, **provided** that those
#' moments are attainable by a truncated normal. If \code{check = TRUE} and the
#' empirical moments deviate too much from the targets, an error is raised.
#'
#' @details
#' When \code{a = -Inf} and \code{b = Inf}, no truncation is applied and the
#' function simply calls \code{stats::rnorm(n, mean, sd)}.
#'
#' For truncated cases (\code{a > -Inf} or \code{b < Inf}), the internal solver
#' \code{\link{tnorm_ab_params_from_moments}} solves the nonlinear moment-matching
#' problem
#' \deqn{
#'   X = Y \mid a \le Y \le b, \quad Y \sim \mathcal{N}(\mu,\sigma^2)
#' }
#' such that \eqn{E[X] = \code{mean}} and \eqn{\mathrm{sd}(X) = \code{sd}}.
#' Not all combinations of \code{(mean, sd, a, b)} are feasible for a truncated
#' normal distribution (for instance, the mean must lie in \code{[a, b]} when
#' both bounds are finite, and the sd cannot exceed the maximum possible sd on
#' the interval). In such infeasible cases, either the parameter solver will
#' error directly, or, if parameters are found but the resulting sample moments
#' are too far from the targets, \code{rtruncnorm_maxent()} will error when
#' \code{check = TRUE}.
#'
#' The post-hoc validation uses \code{.check_sample()} on a single-column matrix
#' of draws. Internally, the sample is considered acceptable if the relative
#' deviation of the sample mean is at most \code{thr_mean} and the relative
#' deviation of the sample sd is at most \code{thr_sd}. Otherwise a descriptive
#' error is thrown indicating that the target moments cannot be achieved by a
#' truncated normal on the specified interval.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{tnorm_ab_params_from_moments}} for the underlying
#'     parameter solver;
#'   \item \code{truncnorm::rtruncnorm()} for the truncated normal RNG;
#'   \item \code{\link{ragg}} for a higher-level random-number generator that
#'     chooses distributions based on argument combinations;
#'   \item \code{.check_sample()} for the sample-based validation routine.
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'
#'   ## 1) Untruncated Normal: mean = 5, sd = 2
#'   x1 <- rtruncnorm_maxent(1e4, a = -Inf, b = Inf, mean = 5, sd = 2)
#'   mean(x1); sd(x1)
#'
#'   ## 2) Lower-truncated at 0: mean = 1, sd = 0.7
#'   x2 <- rtruncnorm_maxent(1e4, a = 0, b = Inf, mean = 1, sd = 0.7)
#'   range(x2); mean(x2); sd(x2)
#'
#'   ## 3) Two-sided truncation [0, 3]: mean = 1.4, sd = 0.6
#'   x3 <- rtruncnorm_maxent(1e4, a = 0, b = 3, mean = 1.4, sd = 0.6)
#'   range(x3); mean(x3); sd(x3)
#'
#'   ## 4) Infeasible target: will error (mean too high relative to sd and bounds)
#'   x4 <- rtruncnorm_maxent(1e4, a = 0, b = 15, mean = 10, sd = 5)
#' }
rtruncnorm_maxent <- function(n, a = -Inf, b = Inf, mean = 0, sd = 1,
                              check = TRUE,
                              thr_mean = 0.10,
                              thr_sd   = 0.20) {
  # Basic input checks
  if (!is.numeric(n) || length(n) != 1L || n <= 0 || !is.finite(n)) {
    stop("n must be a positive finite scalar.")
  }
  n <- as.integer(n)

  if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean)) {
    stop("mean must be a finite numeric scalar.")
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("sd must be a finite numeric scalar > 0.")
  }
  if (!is.numeric(a) || length(a) != 1L ||
      !is.numeric(b) || length(b) != 1L) {
    stop("a and b must be numeric scalars.")
  }
  if (a > b) {
    stop("Require a <= b.")
  }

  # Case 1: no truncation
  if (is.infinite(a) && is.infinite(b)) {
    samp <- stats::rnorm(n, mean = mean, sd = sd)
  } else {
    # Case 2: truncated case -> solve for underlying (mu, sigma)
    out <- tnorm_ab_params_from_moments(m = mean, s = sd, a = a, b = b)
    samp <- truncnorm::rtruncnorm(n, a = a, b = b, mean = out$mu, sd = out$sigma)
  }

  if (check) {
    # capture any warning emitted by .check_sample
    w <- NULL
    withCallingHandlers(
      {
        .check_sample(
          sample       = matrix(samp, ncol = 1,
                                dimnames = list(NULL, "target")),
          target_means = c(target = mean),
          target_sds   = c(target = sd),
          thr_mean     = thr_mean,
          thr_sd       = thr_sd
        )
      },
      warning = function(cond) {
        # store the warning and silence it
        w <<- cond
        invokeRestart("muffleWarning")
      }
    )

    # if .check_sample produced a warning, treat it as infeasible and error
    if (!is.null(w)) {
      stop(
        sprintf(
          paste0(
            "Requested (mean = %.6g, sd = %.6g) could not be matched ",
            "on [a, b] = [%.6g, %.6g] by a truncated normal.\n%s"
          ),
          mean, sd, a, b, conditionMessage(w)
        ),
        call. = FALSE
      )
    }
  }

  samp
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


#' Truncated Normal parameters from target mean m and sd s
#'
#' Given target moments (m, s) for X = Y | a <= Y <= b, Y ~ Normal(mu, sigma^2),
#' return (mu, sigma) such that the truncated distribution matches (m, s).
#' Handles all cases: untruncated, one-sided, and two-sided truncation.
#'
#' @param m numeric scalar, target mean of truncated variable X (must be finite)
#' @param s numeric scalar, target sd of truncated variable X (s > 0, finite)
#' @param a numeric scalar, lower bound (can be -Inf)
#' @param b numeric scalar, upper bound (can be  Inf); must have a < b
#' @param tol absolute tolerance for the residual norm in the 2D solve
#' @param maxit maximum iterations for the 2D Newton solver
#' @return list(mu=..., sigma=..., alpha=..., beta=..., method="closed|shift|upper|2d-newton", iters=...)
#' @export
tnorm_ab_params_from_moments <- function(m, s, a, b, tol = 1e-10, maxit = 100) {
  if (!is.finite(m)) stop("m must be finite.")
  if (!is.finite(s) || s <= 0) stop("s must be finite and strictly positive.")
  if (!(a < b)) stop("Require a < b.")

  # Untruncated: closed-form solution
  if (is.infinite(a) && is.infinite(b)) {
    return(list(mu = as.numeric(m),
                sigma = as.numeric(s),
                alpha = -Inf,
                beta  = Inf,
                method = "closed",
                iters = 0L))
  }

  # Finite two-sided case: existing checks here (if you added them)
  if (is.finite(a) && is.finite(b)) {
    # mean in [a,b]
    if (m < a || m > b) {
      stop("Infeasible: target mean m is outside [a, b].")
    }
    # Popoviciu bound: sd <= (b - a)/2
    max_sd <- (b - a) / 2
    if (s > max_sd + sqrt(.Machine$double.eps)) {
      stop(sprintf("Infeasible: target sd s (%.6g) exceeds max possible on [a,b] (%.6g).",
                   s, max_sd))
    }
  }

  # ---- One-sided lower truncation [a, +Inf) ----
  if (is.finite(a) && is.infinite(b)) {
    if (m < a) {
      stop("Infeasible: for [a, Inf) truncation, require mean m >= a.")
    }
    # Nonnegative support bound: sd <= m - a
    max_sd <- m - a
    if (s > max_sd + sqrt(.Machine$double.eps)) {
      stop(sprintf("Infeasible: for [a, Inf) with mean m, sd must satisfy sd <= m - a. Got sd=%.6g, m-a=%.6g.",
                   s, max_sd))
    }
    # Solve lower-truncation-at-0 for shifted mean (m - a)
    out <- tnorm_params_from_moments(m - a, s)
    return(list(mu = out$mu + a,
                sigma = out$sigma,
                alpha = (a - (out$mu + a)) / out$sigma,
                beta  = Inf,
                method = paste0("shift+", out$method),
                iters = out$iters))
  }

  # ---- One-sided upper truncation (-Inf, b] ----
  if (is.infinite(a) && is.finite(b)) {
    if (m > b) {
      stop("Infeasible: for (-Inf, b] truncation, require mean m <= b.")
    }
    # Reflect to [0, Inf): W = b - X, so mean(b - m), same sd.
    max_sd <- b - m
    if (s > max_sd + sqrt(.Machine$double.eps)) {
      stop(sprintf("Infeasible: for (-Inf, b] with mean m, sd must satisfy sd <= b - m. Got sd=%.6g, b-m=%.6g.",
                   s, max_sd))
    }
    out <- tnorm_params_from_moments(b - m, s)
    mu_x <- b - out$mu
    sig  <- out$sigma
    return(list(mu = mu_x,
                sigma = sig,
                alpha = -Inf,
                beta  = (b - mu_x) / sig,
                method = paste0("upper+", out$method),
                iters = out$iters))
  }

  # ---- Two-sided truncation [a,b]: 2D solve in (alpha, beta) ----
  # Helper: numerically stable pieces
  logspace_sub <- function(logx, logy) { # assumes x>y>0
    if (logy > logx) stop("logspace_sub: require logx >= logy")
    logx + log1p(-exp(logy - logx))
  }

  # Standardized truncated-N(0,1) moments on [alpha, beta]
  tnorm01_moments <- function(alpha, beta) {
    # guard: alpha < beta (strict)
    if (!(alpha < beta)) return(list(ok = FALSE))
    lPhi_a <- pnorm(alpha, log.p = TRUE)
    lPhi_b <- pnorm(beta,  log.p = TRUE)
    # den = Phi(beta) - Phi(alpha), computed in log-space
    if (lPhi_b <= lPhi_a) return(list(ok = FALSE))  # guard numeric weirdness
    lden <- logspace_sub(lPhi_b, lPhi_a)
    den  <- exp(lden)

    la <- dnorm(alpha, log = TRUE);  lb <- dnorm(beta, log = TRUE)
    phi_a <- exp(la);                phi_b <- exp(lb)

    # eta = E[Z | alpha<=Z<=beta] for Z~N(0,1)
    eta <- (phi_a - phi_b) / den

    # v = Var[Z | alpha<=Z<=beta]
    v <- 1 + (alpha * phi_a - beta * phi_b) / den - eta^2

    list(ok = is.finite(eta) && is.finite(v) && v > 0, eta = eta, v = v, den = den,
         phi_a = phi_a, phi_b = phi_b)
  }

  # Residuals:
  # r1 = alpha - (eta + ((a - m) * sqrt(v))/s)
  # r2 = beta  - (eta + ((b - m) * sqrt(v))/s)
  residuals <- function(alpha, beta) {
    mom <- tnorm01_moments(alpha, beta)
    if (!isTRUE(mom$ok)) return(c(Inf, Inf))
    rt <- sqrt(mom$v)
    r1 <- alpha - (mom$eta + ((a - m) * rt) / s)
    r2 <- beta  - (mom$eta + ((b - m) * rt) / s)
    c(r1, r2)
  }

  # Numerical Jacobian (central diff)
  jacobian <- function(alpha, beta, h = 1e-5) {
    r0 <- residuals(alpha, beta)
    # alpha direction
    rpa <- residuals(alpha + h, beta)
    rma <- residuals(alpha - h, beta)
    # beta direction
    rpb <- residuals(alpha, beta + h)
    rmb <- residuals(alpha, beta - h)
    cbind((rpa - rma) / (2 * h), (rpb - rmb) / (2 * h))
  }

  # Initial guess: "weak truncation" heuristic
  # Start near standardized (a-m)/s, (b-m)/s
  alpha <- (a - m) / s
  beta  <- (b - m) / s
  # Make sure alpha < beta (strict); nudge if equal
  if (!(alpha < beta)) {
    eps <- 1e-3
    alpha <- alpha - eps
    beta  <- beta  + eps
  }

  # Clamp to a numerically safe range for pnorm/dnorm log-space math
  clamp <- function(x, lo = -37, hi = 37) pmax(lo, pmin(hi, x))

  # Damped Newton
  iters <- 0
  for (k in 1:maxit) {
    iters <- k
    alpha <- clamp(alpha); beta <- clamp(beta)
    r  <- residuals(alpha, beta)
    nr <- sqrt(sum(r * r))
    if (!is.finite(nr)) break
    if (nr <= tol) break

    J <- jacobian(alpha, beta)
    if (!all(is.finite(J))) break

    # Solve J * step = -r
    step <- tryCatch(solve(J, -r), error = function(e) rep(0, 2))
    if (!all(is.finite(step))) step <- rep(0, 2)

    # Backtracking line search
    t <- 1
    best_alpha <- alpha; best_beta <- beta; best_nr <- nr
    for (ls in 1:20) {
      a_try <- alpha + t * step[1]
      b_try <- beta  + t * step[2]
      if (a_try >= b_try) { t <- t / 2; next }
      a_try <- clamp(a_try); b_try <- clamp(b_try)
      r_try <- residuals(a_try, b_try)
      nr_try <- sqrt(sum(r_try * r_try))
      if (is.finite(nr_try) && nr_try < best_nr) {
        best_alpha <- a_try; best_beta <- b_try; best_nr <- nr_try
        break
      }
      t <- t / 2
    }
    # If no improvement, stop
    if (best_nr >= nr) break

    alpha <- best_alpha; beta <- best_beta
  }

  # Final check / build parameters
  mom <- tnorm01_moments(alpha, beta)
  if (!isTRUE(mom$ok)) {
    stop("Failed to solve for (mu, sigma) on [a,b]; try different (m,s) or adjust tol/maxit.")
  }
  sigma <- s / sqrt(mom$v)
  mu    <- m - sigma * mom$eta

  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    stop("Solution produced non-finite or non-positive sigma.")
  }

  list(mu = as.numeric(mu),
       sigma = as.numeric(sigma),
       alpha = as.numeric(alpha),
       beta  = as.numeric(beta),
       method = "2d-newton",
       iters = iters)
}

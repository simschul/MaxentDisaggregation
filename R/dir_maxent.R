#' Derivative of the Dirichlet Entropy Function (Experimental)
#'
#' Computes an analytical derivative of the [dirichlet_entropy()] function with
#' respect to the concentration parameter \eqn{x}, given a vector of shares.
#'
#' \strong{Note:} This function was derived using the `Deriv` and `autodiffr`
#' packages. However, the correctness of this expression is \strong{not guaranteed},
#' and it is currently \emph{not} used in the optimization workflow (see also
#' the note in the source code).
#'
#' @param x Numeric value (scalar). Concentration parameter (sometimes called \eqn{\gamma}).
#' @param shares Numeric vector of shares that sum to 1.
#'
#' @return A numeric value corresponding to the (tentative) derivative of the
#'   Dirichlet entropy function evaluated at \eqn{x} with given \code{shares}.
#'
#' @seealso [dirichlet_entropy()], [beta2()]
#'
#' @examples
#' # Use with caution; for demonstration only:
#' x_val <- 5
#' sh_vec <- c(0.2, 0.3, 0.5)
#' dirichlet_entropy_grad(x_val, sh_vec)
dirichlet_entropy_grad <- function(x, shares)
{
  .e1 <- shares * x
  .e2 <- sum(.e1)
  -(((1 - prod(gamma(.e1)) / (beta2(.e1) * gamma(.e2))) * digamma(.e2) +
       (.e2 - length(.e1)) * trigamma(.e2)) * sum(shares) -
      sum(shares * ((.e1 - 1) * trigamma(.e1) + digamma(.e1))))
}


#' Generalized Beta Function for Multiple Arguments
#'
#' Extends the usual two-argument Beta function to multiple arguments, i.e.,
#' \eqn{B(\alpha_1, \ldots, \alpha_k) = \frac{\prod_i \Gamma(\alpha_i)}{\Gamma(\sum_i \alpha_i)}}.
#'
#' @param alpha A non-negative numeric vector of parameters.
#'
#' @return A single numeric value, the generalized Beta function of the \code{alpha} vector.
#'
#' @seealso [base::beta()], [gamma()]
#' @export
#'
#' @examples
#' # Example with a 3-argument vector:
#' beta2(c(1, 2, 3))
beta2 <- function(alpha) {
  prod(gamma(alpha)) / gamma(sum(alpha))
}


#' Dirichlet Entropy (Negative) for Given Shares and Concentration
#'
#' Calculates the (negative) entropy of a Dirichlet distribution parameterized by
#' `shares * x`. The formula is from
#' [Wikipedia](https://en.wikipedia.org/wiki/Dirichlet_distribution#Entropy):
#'
#' \deqn{
#'   -\Bigl[\ln B(\alpha) + (\alpha_0 - K)\psi(\alpha_0) - \sum_{i=1}^{K}(\alpha_i - 1)\psi(\alpha_i)\Bigr]
#' }
#'
#' where \eqn{\alpha_i = x \cdot \mathrm{shares}[i]},
#' \eqn{\alpha_0 = \sum_i \alpha_i}, and \eqn{B(\alpha)} is the
#' generalized Beta function ([beta2()]).
#'
#' @param x Numeric. The concentration parameter (\eqn{\gamma}).
#' @param shares Numeric vector of shares (must sum to 1).
#'
#' @return A single numeric value: the negative of the Dirichlet entropy.
#'
#' @seealso [beta2()], [digamma()], [dirichlet_entropy_grad()]
#' @export
#'
#' @examples
#' dirichlet_entropy(x = 10, shares = c(0.2, 0.3, 0.5))
dirichlet_entropy <- function(x, shares) {
  alpha <- x * shares
  K <- length(alpha)
  psi <- digamma(alpha)
  alpha0 <- sum(alpha)
  -( log(beta2(alpha)) +
       (alpha0 - K) * digamma(alpha0) -
       sum((alpha - 1) * psi) )
}


#' Find the \eqn{\gamma} that Maximizes Dirichlet Entropy
#'
#' Uses numerical optimization via the **nloptr** package to find the concentration
#' parameter \eqn{\gamma} that maximizes the entropy of a Dirichlet distribution
#' given a vector of shares. Internally, it calls [dirichlet_entropy()] or a
#' user-specified function for the objective.
#'
#' If the chosen `x0` (initial guess) leads to non-finite objective values, random
#' attempts are made to find a valid starting value. If no valid starting point is
#' found after `x0_n_tries`, shares below `shares_lb` are excluded and re-scaled
#' until a valid configuration is found or the function errors out.
#'
#' @param shares Numeric vector of positive reals that sum to 1.
#' @param eval_f Function that returns the value of the objective (by default, the
#'   negative Dirichlet entropy). Must accept `x` (the \eqn{\gamma} parameter) and
#'   `shares`.
#' @param x0 Numeric. Initial guess for \eqn{\gamma}. Defaults to `1`.
#' @param x0_n_tries Integer. Number of attempts to find a valid initial `x0` if
#'   the objective is non-finite at the current guess. Defaults to `100`.
#' @param bounds Numeric vector of length 2 giving lower and upper bounds on \eqn{\gamma}.
#'   Defaults to `c(0.001, 172)`.
#' @param shares_lb Numeric. Lower bound for shares. Any share below this threshold
#'   is removed and the remaining shares are re-scaled to sum to 1. Defaults to `0`.
#' @param local_opts List of local optimization options passed to [nloptr()].
#'   Defaults to `list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1.0e-4)`.
#' @param opts List of global optimization options for [nloptr()]. Defaults to
#'   `list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1.0e-4,"maxeval"=1000,"local_opts"=local_opts,"print_level"=0)`.
#'
#' @return An object of class `nloptr`, containing (among other components):
#' \itemize{
#'   \item `solution`: The best-fit \eqn{\gamma} that maximizes the Dirichlet entropy.
#'   \item `objective`: The value of the objective function at that solution.
#' }
#'
#' @seealso [nloptr::nloptr()], [dirichlet_entropy()]
#' @export
#'
#' @examples
#' # Simple usage with default objective = dirichlet_entropy
#' sol <- find_gamma_maxent(c(0.2, 0.3, 0.5))
#' sol$solution     # Optimal gamma
#' sol$objective    # The maximum negative entropy
find_gamma_maxent <- function(
    shares,
    eval_f = dirichlet_entropy,
    # eval_grad_f = dirichlet_entropy_grad, # optionally specify a gradient function
    x0 = 1,
    x0_n_tries = 100,
    bounds = c(0.001, 172),
    shares_lb = 0,
    local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
    opts = list("algorithm" = "NLOPT_GN_DIRECT",
                "xtol_rel" = 1.0e-4,
                "maxeval" = 1e3,
                "local_opts" = local_opts,
                "print_level" = 0)
) {
  if (!isTRUE(all.equal(sum(shares), 1))) {
    stop("`shares` must sum to 1, but sum(shares) = ", sum(shares), ".")
  }

  # Exclude shares below shares_lb and re-scale
  shares <- shares[shares > shares_lb]
  shares <- shares / sum(shares)

  lb <- bounds[1]
  ub <- bounds[2]

  count  <- 0
  count2 <- 0
  while(!is.finite(eval_f(x = x0, shares = shares))) {
    if (count > x0_n_tries) {
      if (count2 == 1) {
        stop(
          "Could not find an initial `x0` for `eval_f`.\n",
          "Try increasing `x0_n_tries`, adjusting `bounds`, or raising `shares_lb`.\n"
        )
      }
      # Second attempt: remove small shares again and re-scale
      shares <- shares[shares > shares_lb]
      shares <- shares / sum(shares)
      count  <- 0
      count2 <- 1
    } else {
      x0 <- runif(1, min = lb, max = ub)
      count <- count + 1
    }
  }

  # Run optimizer
  res <- nloptr(
    x0          = x0,
    eval_f      = eval_f,
    # eval_grad_f = eval_grad_f, # uncomment if derivative is correct/desired
    lb          = lb,
    ub          = ub,
    opts        = opts,
    shares      = shares
  )

  return(res)
}

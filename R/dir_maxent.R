#' Derivative of the Dirichlet entropy function.
#' Derivated using the `Deriv` package and `autodiffr`.
#' For derivation see script:
#'
#' NOTE: Derivation does not seem to be correct. Is not used anymore in the optimisation.
#' @param x
#' @param shares
#'
#' @return
#'
#' @examples
dirichlet_entropy_grad <- function (x, shares)
{
  #
  .e1 <- shares * x
  .e2 <- sum(.e1)
  -(((1 - prod(gamma(.e1))/(beta2(.e1) * gamma(.e2))) * digamma(.e2) +
       (.e2 - length(.e1)) * trigamma(.e2)) * sum(shares) -
      sum(shares * ((.e1 - 1) * trigamma(.e1) + digamma(.e1))))
}

#' Generalized beta function
#'
#' A generalization of the base `beta` function with any number of variables
#'
#' @param alpha a non-negative numeric vector
#'
#' @return
#' @export
#'
#' @examples
beta2 <- function(alpha) {
  prod(gamma(alpha)) / gamma(sum(alpha))
}




#' Dirichlet entropy function
#'
#' calculates the negative entropy of a dirichlet distribution given `shares` (sum to 1) and a concentration parameter `x`, so that `alpha = shares * x`
#' from:  https://en.wikipedia.org/wiki/Dirichlet_distribution#Entropy
#' @param x value of the concentration parameter (gamma)
#' @param shares the shares (must sum to 1)
#'
#' @return
#' @export
#'
#' @examples
dirichlet_entropy <- function(x, shares) {
  alpha <- x * shares
  K <- length(alpha)
  psi <- digamma(alpha)
  alpha0 <- sum(alpha)
  (-(log(beta2(alpha)) + (alpha0 - K) * digamma(alpha0) - sum((alpha - 1) * psi)))
}


#' Finds the Gamma value which maximises the entropy of a Dirichlet distribution
#' with given alphas (=shares), using the nloptr optimazation package.
#'
#'
#' @param shares a vector of positive reals which sum to 1
#' @param x0 initial value. If not defined by dirichlet_entropy and/or eval_grad_f there is a automatical procedure to find a valid initial value within the `bounds`. If none is found, the `shares_lb` arguments is used to remove all shares below that threshold (default is 0)
#' @param bounds a vector of size 2 giving lower and upper bounds, default: `c(0.001, 172)`
#' @param shares_lb lower bound of shares. Only is of relevance if no starting value x0 could found be which is defined by dirichlet_entropy and/or eval_grad_f. set this to e.g. 1E-4) to increase chance of convergence
#' @param local_opts see `?nloptr`
#' @param opts passed to `nloptr`, check out `?nloptr`
#' @param eval_f function that returns the value of the objective function (the entropy function of the dirichlet distribution)
#' @param x0_n_tries number of tries to find an initial value of gamma, default: 100
#'
#' @return
#' @export
#'
#' @examples
find_gamma_maxent <- function(shares,
                               eval_f = eval_f,
                               #eval_grad_f = dirichlet_entropy_grad,
                               x0 = 1, # initial value of gamma
                               x0_n_tries = 100,
                               bounds = c(0.001, 172),
                               shares_lb = 0,
                               local_opts = list( "algorithm" = "NLOPT_LD_MMA", # optim options
                                                  "xtol_rel" = 1.0e-4 ),
                               opts = list( "algorithm"= "NLOPT_GN_DIRECT",
                                            "xtol_rel"= 1.0e-4,
                                            "maxeval"= 1E3,
                                            "local_opts" = local_opts,
                                            "print_level" = 0 )

) {




  if (!isTRUE(all.equal(sum(shares), 1))) {
    stop('Shares must sum to 1. But `sum(shares)` gives ', sum(shares))
  }
  shares <- shares[which(shares > shares_lb)]
  shares <- shares/sum(shares)

  # lower and upper bounds
  lb <- bounds[1]
  ub <- bounds[2]

  count <- 0
  count2 <- 0
  while(!is.finite(eval_f(x = x0, shares = shares))) {
    if (count > x0_n_tries) {

      if (count2 == 1) {
        stop('Error: Could not find an initial value x0 which is defined by eval_f and/or eval_grad_f. \n
                You have different options:\n
                1. increase the lower bound for the shares (under which all shares are excluded). E.g. shares_lb = 1E-3 \n
                2. increase x0_n_tries (default: 100)\n
                3. increase the parameter space with the bounds argument (e.g. bounds = c(1E-5, 1E4). Note: might take longer.\n
                If none works, rise an issue on Github. ')
        break
      }

      #warning('Warning: Could not find an initial value x0 which is defined by eval_f and/or eval_grad_f. Either increase x0_n_tries (defaul: 100), or increase the parameter space with the bounds argument. Q&D solution: x0 is set to 1')
      #x0 <- 1

      # remove all very small shares
      shares <- shares[shares > shares_lb]
      # rescale to 1
      shares <- shares/sum(shares)

      count <- 0
      count2 <- 1

    } else {
      x0 <- runif(1, min = lb, max = ub)
      count <- count + 1
    }
  }

  # Run the optimizer
  res <- nloptr ( x0 = x0,
                  eval_f = eval_f,
                  #eval_grad_f = eval_grad_f, # uncomment if wanna use the derivative
                  lb = lb,
                  ub = ub,
                  opts = opts,
                  shares = shares
  )

  return(res)
}


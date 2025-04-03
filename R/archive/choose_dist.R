
#' Choose an Appropriate Distribution Based on Aggregated Statistics
#'
#' This function selects an appropriate distribution based on the input mean, standard deviation,
#' and boundary parameters 'min' and 'max'. It handles various cases and constraints for the parameters,
#' providing warnings and defaulting values as necessary. It aims to guide the user towards a suitable
#' distribution choice based on the provided statistics.
#'
#' @param mean The mean of the distribution, can be NULL or NA.
#' @param sd The standard deviation of the distribution, can be NULL or NA.
#' @param min The lower bound of the distribution, can be NULL or NA.
#' @param max The upper bound of the distribution, can be NULL or NA.
#'
#' @return A character string indicating the type of distribution selected based on the inputs,
#' or NA if the input values do not allow for a proper distribution choice.
#' Possible return values include 'truncnorm', 'norm', 'exp', 'unif', or NA.
#' @export
#'
#' @examples
#' # To choose a distribution with known mean and standard deviation:
#' dist_type <- choose_dist_agg(mean = 5, sd = 2, min = NULL, max = NULL)


choose_dist_agg <- function(mean, sd, min, max) {
  if (is.null(mean)) mean <- NA
  if (is.null(sd)) sd <- NA
  if (is.null(min)) min <- NA
  if (is.null(max)) max <- NA

  if (is.finite(min) & is.finite(max) & min >= max) {
    warning('min must be < max')
    return(NA)
  }
  if (is.finite(sd) & sd <=0) {
    warning('sd must be positive')
    return(NA)
  }
  if (is.finite(mean) & is.finite(min) & mean <= min) {
    warning('mean must be > min')
    return(NA)
  }
  if (is.finite(mean) & is.finite(max) & mean >= max) {
    warning('mean must be < max')
    return(NA)
  }

  if (is.finite(mean)) {
    if (is.finite(sd)) {
      if (is.finite(min) | is.finite(max)) {
        # return(generate_distribution('rtruncnorm', mean = mean, sd = sd,
        #                              min = min, max = max))

        return('truncnorm')
      } else {
        # return(generate_distribution('rnorm', mean = mean, sd = sd))
        return('norm')
      }
    } else if (isTRUE(all.equal(min, 0)) & !is.finite(max)) {
      return('exp')
    } else {
      return(NA)
    }
  } else if (is.finite(min) & is.finite(max) & !is.finite(sd)) {
    return('unif')
  } else {
    NA
  }
}

#' Choose an Appropriate Disaggregated Distribution
#'
#' This function determines the appropriate disaggregated distribution type based on shares and sds parameters.
#' It considers various scenarios and constraints related to the input parameters to select from specific distribution types.
#'
#' @param shares A numeric vector representing the shares parameter(s) for the distribution.
#' @param sds A numeric vector representing the sds parameter(s) for the distribution, can be NULL or NA.
#'
#' @return A character string indicating the type of disaggregated distribution selected based on the inputs,
#' Possible return values include 'dir1' for Dirichlet 1, 'dir_maxent' for Dirichlet Maximum Entropy, 'gdir' for Generalized Dirichlet,
#' or NA if the inputs do not match any predefined distribution type.
#' @export
#'
#' @examples
#' # To determine a distribution with known shares values:
#' dist_type <- choose_dist_disagg(shares = c(1,2,3), sds = NULL)
choose_dist_disagg <- function(shares, sds) {
  if (all(is.na(sds))) sds <- NULL

  if (isTRUE(all.equal(var(shares), 0)) & is.null(sds)) {
    # Dirichlet 1
    return('dir1')
  } else if (is.null(sds)) {
    # Dirichlet MaxEnt
    return('dir_maxent')
  } else if (!is.null(sds)) {
    # Gen. Dirichlet
    return('gdir')
  } else {
    return(NA)
  }
}


#' Generate Distribution Function
#'
#' Generates a distribution function based on specified parameters and distribution type.
#' The function allows dynamic creation of distribution functions from a preset list of supported distributions.
#' It handles numeric, character, and named vector parameters.
#'
#' @param fun The distribution function to be generated (runif, rnorm, rlnorm, rweibull, rdir, etc.).
#' @param ... Additional parameters for the distribution function such as mean, sd for normal distribution,
#'        or shares for dirichlet distribution. Supports named vectors for parameter values.
#' @return A new function that takes a single argument 'n' and returns a sample of size 'n' from the specified distribution.
#'         The function is self-contained and includes all necessary parameters as part of its definition.
#' @examples
#' # Generate a normal distribution function with mean=10 and sd=2
#' my_normal <- generate_distribution(rnorm, mean = 10, sd = 2)
#' sample <- my_normal(100)  # Generate 100 samples
#'
#' # Generate a dirichlet distribution function with named shares vector
#' my_dirichlet <- generate_distribution(rdir, shares = c('3' = 1))
#' sample <- my_dirichlet(50)  # Generate 50 samples
#'
#' @export

generate_distribution <- function(fun, ...) {
  options_allowed <- c('runif', 'rnorm', 'rlnorm', 'rexp', 'rtruncnorm',
                       'rdir', 'rdirg')
  fun_name <- as.character(substitute(fun))

  if (!(fun_name %in% options_allowed)) {
    stop(paste0("'fun' should be one of ", toString(options_allowed)))
  }

  # Capture the additional arguments as a list
  args <- list(...)

  # Create a string representation of the arguments with values
  # args_with_values <- sapply(names(args), function(arg) {
  #   value <- args[[arg]]
  #   if (is.numeric(value) && length(value) == 1) {
  #     return(paste0(arg, " = ", value))
  #   } else if (is.character(value) && length(value) == 1) {
  #     return(paste0(arg, " = '", value, "'"))
  #   } else if (is.numeric(value) || is.character(value)) {
  #     # Properly format vectors
  #     vector_string <- toString(value)
  #     return(paste0(arg, " = c(", vector_string, ")"))
  #   } else {
  #     # Add more conditions as necessary for other types of arguments
  #     return(paste0(arg, " = ", value))
  #   }
  # })

  # args_with_values <- sapply(names(args), function(arg) {
  #   value <- args[[arg]]
  #   if (is.numeric(value) && length(value) == 1) {
  #     return(paste0(arg, " = ", value))
  #   } else if (is.character(value) && length(value) == 1) {
  #     return(paste0(arg, " = '", value, "'"))
  #   } else if (is.numeric(value) || is.character(value)) {
  #     # Format named vectors properly
  #     elements <- if (is.null(names(value))) {
  #       # Unnamed vector
  #       paste(value, collapse = ", ")
  #     } else {
  #       # Named vector, ensuring names and values are quoted appropriately
  #       paste(sprintf("'%s' = %s", names(value), value), collapse = ", ")
  #     }
  #     return(paste0(arg, " = c(", elements, ")"))
  #   } else {
  #     # Add more conditions as necessary for other types of arguments
  #     return(paste0(arg, " = ", value))
  #   }
  # })

  args_with_values <- sapply(names(args), function(arg) {
    value <- args[[arg]]
    # Checking if it's a named vector
    if (is.numeric(value) || is.character(value)) {
      if (!is.null(names(value))) {
        # Named vector: Create a string representation preserving the names
        named_elements <- paste(sprintf("'%s' = %s", names(value), value), collapse = ", ")
        return(paste0(arg, " = c(", named_elements, ")"))
      } else {
        # Unnamed vector or single value: Convert directly to string
        return(paste0(arg, " = c(", toString(value), ")"))
      }
    } else {
      # Handle other types of arguments if necessary
      return(paste0(arg, " = ", value))
    }
  })

  # Create the function body as a string
  fun_body_string <- sprintf("function(n) %s(n, %s)",
                             fun_name,
                             paste(args_with_values, collapse = ", "))

  # Create a new function from the string
  eval(parse(text = fun_body_string), envir = parent.frame())
}


#' Generate Aggregated Distribution Based on Parameters
#'
#' This function generates a random sample from an aggregated distribution determined by specified mean, standard deviation,
#' and boundary parameters. It uses the 'choose_dist_agg' function to select an appropriate distribution based on the input parameters
#' and then generates a random sample from the chosen distribution.
#'
#' @param mean The mean for the distribution, default is NA.
#' @param sd The standard deviation for the distribution, default is NA.
#' @param min The lower boundary for the distribution, default is -Inf.
#' @param max The upper boundary for the distribution, default is Inf.
#'
#' @return A random sample from the chosen aggregated distribution or a stop message if no suitable distribution is found.
#' @export
#'
#' @examples
#' # Generate a distribution with specified mean and sd:
#' sample <- generate_distribution_agg(mean = 5, sd = 2)
generate_distribution_agg <- function(mean = NA, sd = NA, min = -Inf, max = Inf) {
  dist <- choose_dist_agg(mean = mean, sd = sd, min = min, max= max)

  if (dist == 'truncnorm') {
    if (is.null(min) || is.na(min)) min <- -Inf
    if (is.null(max) || is.na(max)) max <- Inf
    return(generate_distribution(rtruncnorm, mean = mean, sd = sd, a= min, b=max))
  } else if (dist == 'norm') {
    return(generate_distribution(rnorm, mean = mean, sd = sd))
  } else if (dist == 'exp') {
    return(generate_distribution(rexp, rate = 1/mean))
  } else if (dist == 'unif') {
    return(generate_distribution(runif, min = min, max = max))
  } else {
    stop('no distribution found for parameter combintaion')
  }
}

#' Title
#'
#' @param shares
#' @param sds
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
generate_distribution_shares <- function(shares, sds, ...) {
  dist <- choose_dist_disagg(shares = shares, sds = sds)

  if (dist == 'dir1') {
    return(generate_distribution(rdir, shares = shares, gamma = length(shares)))
  } else if (dist == 'dir_maxent') {
    out <- find_gamma_maxent2(shares, eval_f = dirichlet_entropy, ...)
    return(generate_distribution(rdir, shares = shares, gamma = out$solution))
  } else if (dist == 'gdir') {
    return(generate_distribution(rdirg, shares = shares, sds = sds))
  } else {
    stop('no distribution found for parameter combintaion')
  }
}

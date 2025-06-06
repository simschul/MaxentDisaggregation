% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dir_maxent.R
\name{dirichlet_entropy_grad}
\alias{dirichlet_entropy_grad}
\title{Derivative of the Dirichlet Entropy Function (Experimental)}
\usage{
dirichlet_entropy_grad(x, shares)
}
\arguments{
\item{x}{Numeric value (scalar). Concentration parameter (sometimes called \eqn{\gamma}).}

\item{shares}{Numeric vector of shares that sum to 1.}
}
\value{
A numeric value corresponding to the (tentative) derivative of the
Dirichlet entropy function evaluated at \eqn{x} with given \code{shares}.
}
\description{
Computes an analytical derivative of the \code{\link[=dirichlet_entropy]{dirichlet_entropy()}} function with
respect to the concentration parameter \eqn{x}, given a vector of shares.
}
\details{
\strong{Note:} This function was derived using the \code{Deriv} and \code{autodiffr}
packages. However, the correctness of this expression is \strong{not guaranteed},
and it is currently \emph{not} used in the optimization workflow (see also
the note in the source code).
}
\examples{
# Use with caution; for demonstration only:
x_val <- 5
sh_vec <- c(0.2, 0.3, 0.5)
dirichlet_entropy_grad(x_val, sh_vec)
}
\seealso{
\code{\link[=dirichlet_entropy]{dirichlet_entropy()}}, \code{\link[=beta2]{beta2()}}
}

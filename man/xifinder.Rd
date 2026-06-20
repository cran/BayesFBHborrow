\name{xifinder}
\alias{xifinder}
\title{Find the root of the $p_0$ mixture probability equation}
\description{
  Finds the value of the root (in the interval \code{[0, 2]}) for the function \code{p0tp} involving mixture probability \code{p_0} and hyperparameters \code{b_tau} and \code{d_tau}.
}
\usage{
 xifinder(b_tau, d_tau, p_0)
}
\arguments{
  \item{b_tau}{Shape (or scale) parameter for the first mixture component.}
  \item{d_tau}{Shape (or scale) parameter for the second mixture component.}
  \item{p_0}{Mixture proportion; probability for the first component, in (0, 1).}
}
\details{
  Uses \code{uniroot()} to numerically solve for the root of the function \code{p0tp} on the interval \code{[0, 2]}. The function \code{p0tp} should accept arguments \code{x}, \code{b_tau}, \code{d_tau}, and \code{p_0}. This is typically used for calibrating the mixture prior or solving for a relevant quantile in hierarchical Bayesian models.

  \strong{Note:} The auxiliary function \code{p0tp} must be defined in the environment.
}
\value{
  Numeric value; the root of the \code{p0tp} equation within the interval.
}
\references{
  Scott DAV and Lewin A (2026). Borrowing from historical
  control data in a Bayesian time-to-event model with flexible baseline 
  hazard function. Biostatistics 27(1)

  Scott DAV, Axillus S, Lewin A and Izmirlian G. 
  BayesFBHborrow: An R Package for Bayesian borrowing for time-to-event
  data from a flexible baseline hazard function (2026)
  <doi:10.48550/arXiv.2408.04327> 
}
\seealso{\code{plot}, \code{summary}, \code{update} methods for class
  \code{"BayesFBHborrow"} objects
}
\examples{
## Suppose p0tp is defined elsewhere
## Find a root for specific parameters:
  \dontrun{
    xifinder(b_tau = 0.1, d_tau = 1, p_0 = 0.75)
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

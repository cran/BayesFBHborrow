\name{summary.BayesFBHborrow}
\alias{summary.BayesFBHborrow}
\title{
 Summary Method for the BayesFBHborrow Class
}
\description{
 Summary method for the BayesFBHborrow class
}
\usage{
  \method{summary}{BayesFBHborrow}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{BayesFBHborrow}.}
  \item{...}{Additional arguments to generic summary}
}
\details{
  Every call to the function runs the sampler and then produces the
  conditional estimand and marginal contrast for the
  treatment/pretreatment variable effects. The argument \code{G_compute}
  in the \code{BayesFBHborrow} call tells the print, plot, summary and
  coef methods which version the user wants displayed.
}
\value{
  The result of the summary method is an object of class
  \code{BayesFBHborrow} identical to the input with the addition of two
  new components -- the ones shown by the print method.
    \item{coefficients}{If the original model involved p pre-treatment
      covariates and a treatment indicator, this will be a 2p + 1 by 4
      matrix showing median of the logged hazard ratio, its
      exponentiated value and its 95 percent credible interval.}
    \item{surv_summary}{Shows arm specific conditional or A 8 by 6 matrix showing arm label C/I,
      information fraction Description of 'comp2'}
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
\author{
  Darren A. V. Scott, Sophia Axillus and Grant Izmirlian  
}
\seealso{
  \code{\link{BayesFBHborrow}},
}
\examples{
    \dontrun{
        my_coefs_cndl <- coef(my_bfhb)
        my_coefs_mgnl <- coef(update(my_bfhb, G_compute=TRUE))
    }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

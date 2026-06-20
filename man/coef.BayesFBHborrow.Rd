\name{coef.BayesFBHborrow}
\alias{coef.BayesFBHborrow}
\alias{coefficients}
\title{
The coefficients extractor method for the \code{BayesFBHborrow} class
}
\description{
The coefficients extractor method for the \code{BayesFBHborrow} class
with functionality dependent upon whether you select the conditional
estimand or marginal contrast.
}
\usage{
coef.BayesFBHborrow(object, ...)
}
\arguments{
  \item{object}{
    An object of class \code{BayesFBHborrow}
}
  \item{\dots}{
    Not used    
}
}
\details{
  Every call to the function runs the sampler and then produces the
  conditional estimand and marginal contrast for the
  treatment/pretreatment variable effects. The argument \code{G_compute}
  in the \code{BayesFBHborrow} call tells the print, plot, summary and
  coef methods which version the user wants displayed.
}
\value{
  \code{G_compute=FALSE}, then \code{coef()} will extract the
  conditional logged hazard ratio coefficients.
  
  If \code{G_compute=TRUE} then a table of showing the marginal
  treatment contrast e.g. the difference between the log-logged marginal
  survival fucntion at landmark information times, 0.25, 0.50, 0.75,
  1.00, with each line showing median (MTE: marginal treatment effect)
  its exponentiated value and its 95 percent credicble interval.
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

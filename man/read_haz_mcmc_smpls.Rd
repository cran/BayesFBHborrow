\name{read_haz_mcmc_smpls}
\alias{read_haz_mcmc_smpls}
\title{
  Reads in Fine Mesh Grid Hazards MCMC Samples
}
\description{
  Reads in fine mesh grid conditional baseline hazard and arm specific
  marginal hazards of MCMC samples from the posterior that were written
  out sample at a time during iteration.  
}
\usage{
  read_haz_mcmc_smpls(obj)
}
\arguments{
  \item{obj}{
    An object of class \code{BayesFBHborrow}. 
  }
}
\details{
  MCMC samples of the baseline hazard are natively represented as number
  of split-points, location of split-points and piecewise
  hazards. However the gibbs/MH updating done internally invovles
  integrals over the time to event CDF is done with respect to a fine
  mesh grid. Thus all \code{iter} samples must be snapped to this common
  fine mesh grid of length \code{max_grid}. Row by row updates of this
  matrix within iterations of the sampler is impractical as multiple
  copies and persistence within a loop will quickly consume
  resources. For this reason the\code{max_grid} by \code{iter}
  matrix representations of the conditional baseline hazard, marginal
  control arm hazard and marginal treated arm hazard are written out to
  three corresponding binary files: \code{cnd_blh_*.bin},
  \code{mgnl_haz0_*.bin}, and \code{mgnl_haz1_*.bin}. The files are
  written out sample by sample during iterations and then read back in
  at the end for summarizing. The files are given a unique random hex
  suffix for each call to the function and are left in the working
  directory so that the user may use them to calculate any desired
  functional of the posterior distribution. If you delete the *.bin
  files it will, of course, not work. The names of the binary files
  associated with a given result are returned with the object so if
  desired you can clean house.  
}
\value{
  A list with components
    \item{cnd_haz_0}{samples of the conditional baseline hazard
      posterior as \code{iter} by \code{max_grid} matrix}
    \item{mgnl_haz_0}{samples of the marginal control arm hazard
      posterior as \code{iter} by \code{max_grid} matrix} 
    \item{mgnl_haz_1}{samples of the marginal treated arm hazard
      posterior as \code{iter} by \code{max_grid} matrix} 
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
  \code{\link{BayesFBHborrow}}
}
\examples{
  \dontrun{
      read_haz_mcmc_smpls(fit_my_bfbhb)
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

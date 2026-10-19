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
  of split-points, location of split-points and corresponding values of
  the piecewise hazards. However the gibbs/MH updating done internally
  invovles integrals over the time to event CDF is done with respect to
  a fine mesh grid. Thus all \code{iter} samples must be snapped to a
  common fine mesh grid of length \code{max_grid}. Row by row updates of
  this matrix within iterations of the sampler is impractical based upon
  the way R handles updating arrays involving multiple copies of it and
  persistence of this memory usage within a loop until completion. This
  will quickly consume free memory resulting in diskcaching and
  ultimately machine freeze. For this reason binary file writes, per
  iteration of the sampling loop, are done for the snapped to grid
  versions of the following fucntions: the control arm conditional hazard,
  conditional cumulative hazard, arm specific logged negative logged
  marginal survival functions for which the difference is the marginal
  treatment effect as a function of time, and the marginal hazards,
  computed as the finite difference quotient of the their exponentiated
  values. These files are then read back in after the sampling loop
  for the summarization step. The corresponding six file names in
  order given above are: \code{cnd_haz_0_*.bin}, \code{cnd_cumHz0_*.bin},
  \code{lml_mgnl_surv_0_*.bin}, \code{lml_mgnl_surv_1_*.bin},
  \code{mgnl_haz_0_*.bin}, and \code{_mgnl_haz_1_*.bin}. The default
  behavior is to delete them when the main routine,
  \code{BayesFBHborrow}, exits. However, if the user wishes to use them
  in downstream calculation of functionals of the posterior distribution
  of these quantities, the argument \code{keep_binary_files} should be
  set to \code{TRUE}. Then one calls this function with the single
  argument \code{object}, the result of \code{BayesFBHborrow}. 
  The names of the binary files associated with a given result are
  returned with the object.
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

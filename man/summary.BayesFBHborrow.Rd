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
  \item{\dots}{Additional arguments to generic summary}
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
  \donttest{
    n_cc_1 <- 200
    n_cc_0 <- 100
    n_hst <- 100
    int_cc <- int_hst <- -log(3)
    B_x_cc <- B_x_hst <- c(-0.3,0.5,0.25,-0.5)
    B_trt <- log(0.55)

    dat_lst <- 
        genBFBHBdat(n_cc_1=200, n_cc_0=100, n_hst=100,
                    B_trt=B_trt, B_x_cc=B_x_cc, B_x_hst=B_x_hst, 
                    int_cc=int_cc, int_hst=int_hst, shape=2,
                    t_er=0.5, t_fin=1.5, X_fact_levs=3)

    DAT_cc <- trial_data(dat_lst)
    DAT_hst <- historical_data(dat_lst)

    my_bfbhb <- BayesFBHborrow(Surv(tte, event) ~ X_trt + X_01 + X_02 + X_03,
                               data = DAT_cc, data_hist = DAT_hst,
                               model_choice = 'mix', iter = 10, warmup_iter = 5,
                               verbose = FALSE)

    my_bfbhb.s <- summary(my_bfbhb)
    my_bfbhb.s
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

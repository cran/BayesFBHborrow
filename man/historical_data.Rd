\name{historical_data}
\alias{historical_data}
\title{
  Historical controls data extractor for a BFBHBdat object
}
\description{
  Historical controls data extractor for object of class BFBHBdat as returned
  by the simulation routine genBFBHBdat.
}
\usage{
historical_data(x)
}
\arguments{
  \item{x}{
    An object of class BFBHdat as returned by the  simulation routine genBFBHBdat
  }
}
\value{
   The historical controls data.frame with variables 
    \item{tte}{Time to event}
    \item{evnet}{Event indicator}
    \item{X_01 ... X_0n}{Pre-treatment covariates as directed in the
      call to genBFBHBdat}
}
\author{
  Darren A. V. Scott, Sophia Axillus and Grant Izmirlian
}
\seealso{
  \code{\link{trial_data}}, \code{\link{genBFBHBdat}}
}
\examples{
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
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

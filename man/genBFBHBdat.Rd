\name{genBFBHBdat}
\alias{genBFBHBdat}
\title{
Generate Current Trial/Historical Controls via Weibull Model
}
\description{
  Data genereration to simultaneously generate from a proportional
  hazards Weibull model, two datasets, one to take the role as the
  current trial and the other to take the role as the historical
  controls.
}
\usage{
genBFBHBdat(n_cc_1, n_cc_0, n_hst, B_trt, B_x_cc, int_cc, B_x_hst,
            int_hst, X_fact_levs = NULL, shape = 1, t_er, t_fin,
            dbg = FALSE) 
}
\arguments{
  \item{n_cc_1}{
    Size of treated arm in current trial.
}
  \item{n_cc_0}{
    Size of control arm in current trial.
}
  \item{n_hst}{
    Size of historical controls dataset.
}
  \item{B_trt}{
    Treatment effect (logged hazard ratio).
}
  \item{B_x_cc}{
    Coefficients (logged hazard ratio per unit) for effect of
    pretreatment covariates in current trial. The leading end of the
    vector specifies continuous covariates and the tail end of the
    vector specifies coefficients for levels of factor variables. See
    \code{X_fact_levs} below. 
}
  \item{int_cc}{
    Intercept for current trial model.
}
  \item{B_x_hst}{
    Coefficients (logged hazard ratio per unit) for effect of
    pretreatment covariates in historical coontrols. See X_fact_levs.  
}
  \item{int_hst}{
    Intercept for historical control model.
}
  \item{X_fact_levs}{
    A numeric vector of length 'k' which specifies both the number, 'k',
    of factor variables among the list of pre-treatment covariates, and
    the number of levels of each factor variable in their components.
    For example if we desire 3 continous variables and two factor
    variables having 3 levels and 4 levels, then the coefficient vectors
    described above \code{B_x_} must be of length 8 (3 + 2 + 3), and
    \code{X_fact_levs=c(3,4)}. Any combination from 0/0, n/0, 0/m of
    continuous and factor variables is possible. 
}
  \item{shape}{
    The shape parameter for the Weibull distribution--currently the same
    for the current trial and historical controls. 
}
  \item{t_er}{
    The time at which randomization ends for use in administrative
    censoring for both the current trial and historical controls. 
}
  \item{t_fin}{
    The time of maximum follow-up for use in administrative censoring for
    both the current trial and historical controls. 
}
  \item{dbg}{
    If run at the command line, set \code{dbg=TRUE} to run the command
    line debugger to step through the function line at a time. 
}
}
\value{
  A list of two component data.frames, \code{DAT_cc} and \code{DAT_hst},
  the current trial and the historical controls. Both current trial and
  historical controls set have time to event and event indicator and any
  pretreatment covariates specified via the arguments \code{B_x_} and
  \code{X_fact_levs}. In addition the current trial, \code{DAT_cc}, has
  a treatment indicator. Variables in the data.frames are named as
  follows.
  \item{tte}{Time to event, both data.frames.}
  \item{event}{Event indicator, both data.frames.}
  \item{X_trt}{Treatment indicator, current trial only.}
  \item{X_01...X_m}{Pre-treatment covariates both data.frames.}
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

    DAT_cc <- dat_lst$DAT_cc
    DAT_hst <- dat_lst$DAT_hst
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

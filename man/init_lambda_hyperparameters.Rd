\name{init_lambda_hyperparameters}
\alias{init_lambda_hyperparameters}
\title{Initialize lambda hyperparameters}
\description{
  Propose lambda hyperparameters for the choice of 
  initial values for lambda
}
\usage{
  init_lambda_hyperparameters(group_data, s, w = 0.5)
}
\arguments{
  \item{group_data}{group level data}
  \item{s}{split points}
  \item{w}{weight}
}
\value{
  shape and rate for the estimated lambda distribution
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
    set.seed(111)

    ## simulate an example data and set your initial values and hyper parameters 
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
    
    ## Say we want to know the group level data for the following split points
    s <- with(DAT_cc, quantile(tte, c(0, 0.45, 0.65, 1), names = FALSE))

    grpd_data <- 
      with(DAT_cc, group_summary(tte, event, NULL, s))

    init_lambda_hyperparameters(grpd_data, s)
}

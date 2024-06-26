% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs.R
\name{plot.BayesFBHborrow}
\alias{plot.BayesFBHborrow}
\title{Plot the MCMC results}
\usage{
\method{plot}{BayesFBHborrow}(x, x_lim, x_pred = NULL, ...)
}
\arguments{
\item{x}{object of class "BayesFBHborrow" to be visualized}

\item{x_lim}{x-axis to be used for plot, set to NULL to use default from MCMC sampling}

\item{x_pred}{vector of chosen predictors}

\item{...}{other plotting arguments, see .plot_matrix()
for more information}
}
\value{
nested list of 'plots' (posterior predictive hazard, survival, 
and hazard ratio) as well as their samples.
}
\description{
S3 object which produces predictive probabilities of the survival,
hazard, and hazard ratio for a given set of predictors
}
\examples{
data(weibull_cc, package = "BayesFBHborrow")

# Set your tuning parameters
tuning_parameters <- list("Jmax" = 5,
                          "pi_b" = 0.5,
                          "cprop_beta" = 0.5)
                          
# run the MCMC sampler
out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters = tuning_parameters,
                      iter = 3, warmup_iter = 1)

# for the treatment group
plots <- plot(out$out, out$out$time_grid, x_pred = c(1))
}

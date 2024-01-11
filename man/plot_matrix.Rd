% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs.R
\name{plot_matrix}
\alias{plot_matrix}
\title{Plot smoothed baseline hazards}
\usage{
plot_matrix(
  x_lim,
  y,
  percentiles = c(0.05, 0.95),
  title = "",
  xlab = "",
  ylab = "",
  color = "blue",
  fill = "blue",
  linewidth = 1,
  alpha = 0.2
)
}
\arguments{
\item{x_lim}{time grid}

\item{y}{samples}

\item{percentiles}{percentiles to include in plot, default is c(0.025, 0.975)}

\item{title}{optional, add title to plot}

\item{xlab}{optional, add xlabel}

\item{ylab}{optional, add ylabel}

\item{color}{color of the mid line, default is blue}

\item{fill}{color of the percentiles, default is blue}

\item{linewidth}{thickness of the plotted line, default is 1}

\item{alpha}{opacity of the percentiles, default is 0.2}
}
\value{
a ggplot2 object
}
\description{
Plot mean and given quantiles of a matrix. Can also be used to 
plot derivatives of the baseline hazard, such as estimated cumulative hazard 
and survival function.
}
\examples{
data(weibull_cc, package = "BayesFBHborrow")

# Set your initial values and hyper parameters
s <- c(0, quantile(weibull_cc$tte, c(0.5, 1.0), names = FALSE))
initial_param <- list("J" = 1, 
                      "s_r" = s[2], # split points only (length J) 
                      "mu" = 0, 
                      "sigma2" = 2,
                      "lambda" = c(1, 1), 
                       "beta" = 0)
               
# run the MCMC sampler CHANGE THIS BEFORE FINAL CHECK
out <- BayesFBHborrow(weibull_cc, NULL, initial_param,
                iter = 10, warmup_iter = 1)

# Visualize the smoothed baseline hazard
time_grid <- seq(0, max(s), length.out = 2000)
gg <- plot_matrix(time_grid, out$out_slam, 
                  title = "Example plot of smoothed baseline hazard",
                  xlab = "time", ylab = "baseline hazard")
}
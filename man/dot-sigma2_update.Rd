% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_updates.R
\name{.sigma2_update}
\alias{.sigma2_update}
\title{Calculate sigma2 posterior update}
\usage{
.sigma2_update(mu, lambda_0, Sigma_s, J, a_sigma, b_sigma)
}
\arguments{
\item{mu}{mean.}

\item{lambda_0}{Baseline hazard.}

\item{Sigma_s}{VCV matrix (j + 1) x (j + 1).}

\item{J}{Number of split point.}

\item{a_sigma}{Hyperparameter a.}

\item{b_sigma}{Hyperparameter b.}
}
\value{
sigma2 draw from IG
}
\description{
Calculate sigma2 posterior update
}

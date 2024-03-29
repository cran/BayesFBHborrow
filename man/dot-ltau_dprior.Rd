% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RJMCMC.R
\name{.ltau_dprior}
\alias{.ltau_dprior}
\title{Calculate log density tau prior}
\usage{
.ltau_dprior(tau, a_tau, b_tau, c_tau = NULL, d_tau = NULL, p_0 = NULL, type)
}
\arguments{
\item{tau}{current value(s) of tau}

\item{a_tau}{tau hyperparameter}

\item{b_tau}{tau hyperparameter}

\item{c_tau}{tau hyperparameter}

\item{d_tau}{tau hyperparameter}

\item{p_0}{mixture ratio}

\item{type}{choice of borrowing, "mix", "uni", or any other string for 
borrowing on every baseline hazard without mixture}
}
\value{
log density of tau
}
\description{
Calculate log density tau prior
}

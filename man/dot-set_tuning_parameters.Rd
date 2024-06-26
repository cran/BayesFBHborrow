% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.set_tuning_parameters}
\alias{.set_tuning_parameters}
\title{Set tuning parameters}
\usage{
.set_tuning_parameters(tuning_parameters = NULL, borrow, X, X_0 = NULL)
}
\arguments{
\item{tuning_parameters}{list of tuning_parameters, could contain any combination
of the listed tuning parameters}

\item{borrow}{choice of borrow, could be TRUE or FALSE}

\item{X}{design matrix for concurrent trial}

\item{X_0}{design matrix for historical trial}
}
\value{
filled list of tuning_parameters
}
\description{
Set tuning parameters
}

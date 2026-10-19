\name{rel_eff}
\alias{rel_eff}
\title{
  Ad-hoc Relative Efficiency for Model With and Without Borrowing
}
\description{
  Calculates an ad-hoc relative efficiency for a model with and
  without borrowing, both being fitted BayesFBHborrow models.
  This is done by computing the squared ratio of the widths
  of their 95\% credible intervals about the estimated
  treatment effect.
}
\usage{
  rel_eff(wb, nb, G_compute = FALSE)
}
\arguments{
  \item{wb}{
    A fitted BayesFBHborrow model, with borrowing
  }
  \item{nb}{
    The same fitted BayesFBHborrow model, without borrowing
  }
  \item{G_compute}{
    If \code{G_compute} is \code{FALSE} (default), conditional
    estimates are compared. If it is set to \code{TRUE}, marginal
    estimates are compared.
  }
}
\value{
  A double precision numeric.
}
\author{
    Darren A. V. Scott, Sophia Axillus and Grant Izmirlian
}

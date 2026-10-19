\name{MTE}
\alias{MTE}
\title{
 Marginal Treatment Effect Extractor
}
\description{
 Marginal Treatment Effect Extractor. Works on any BayesFBHborrow object
 having one or more pretreatment covariates.
}
\usage{
MTE(object, ...)
}
\arguments{
  \item{object}{
     an object of class BayesFBHborrow e.g. result of the main method. 
  }
  \item{\dots}{
    other arguments for the update method
  }
}
\value{
  A table (data.frame) showing the marginal treatment at 25\%, 50\%, 75\% and 100\%
  \item{InfFrac}{Information fraction: portion of total events observed
    0.25, 0.50, 0.75 and 1.00 defining landmark times at which the
    marginal treament effect is calculated.}
  \item{Time}{The landmark times at which the indicated information
    fractions are attained.}
  \item{MTE}{The calculated marginal treatment effect, e.g. the
    difference of arm specific logged negative logged marginal survival
    functions; on the scale of logged relative risk}
  \item{exp(MTE)}{Exponentiated marginal treatment effect; on the scale
    of relative risk.}
  \item{MTE_L}{The lower bound of the 95\% credible interval for MTE.}
  \item{MTE_U}{The upper bound of the 95\% credible interval for MTE.}
}
\author{
  Darrent Scott; Sophia Axillus and Grant Izmirlian
}

\name{Combine}
\alias{Combine}
\title{
An Method for Easily Combining Two Similar BayesFBHborrow Plots
}
\description{
  Combines two survival/hazard plots side by side or two TrtEff density
  plots atop one another.
}
\usage{
  Combine(p1, p2, lgnd = NULL)
}
\arguments{
  \item{p1}{
    First plot
  }
  \item{p2}{
    Second plot
  }
  \item{lgnd}{
    For the side by side hazard/survival plots the default is
    that one is conditional and one is marginal. If your use case is
    different, then put a legend consisting of two character strings
    to label the plots from left to right. See examples. 
  }
}
\details{
  Experimental at this stage. 
}
\value{
  A ggplot2 object
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
  \code{\link{plot.BayesFBHborrow}}, \code{\link{BayesFBHborrow}}
}
\examples{
  \dontrun{
      p_haz_nb_wb <- Combine(p_haz_nb, p_haz_wb)
      p_gbcs_TrtEff <- Combine(p_gbcs_nb_TrtEff, p_gbcs_wb_TrtEff, lgnd=c("NB","WB"))
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

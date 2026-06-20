\name{plot.BayesFBHborrow}
\alias{plot.BayesFBHborrow}
\title{
Plot Method for the BayesFBHborrow Class
}
\description{
The plot method for the BayesFBHborrow class which allows quick
convenient plots of important model summaries. 
}
\usage{
\method{plot}{BayesFBHborrow}(x, type = c("survival", "hazard", "trace", "TrtEff"),
                              col = NULL, InfFrac = NULL, ylim = NULL, ...)
}
\arguments{
  \item{x}{
    An object of class \code{BayesFBHborrow}.
  }
  \item{type}{
    The type of plot desired. Can be one of the following "survival", "hazard",
    "trace", "TrtEff".  
  }
  \item{col}{
    When \code{type="trace"} is specified you must name a column of the
    \code{samples} data.frame, which is a component of the output. 
  }
  \item{InfFrac}{
    When \code{type="TrtEff"} is specified and the fitted object has
    been cast to its marginal estimand form, e.g. \code{update(x,
      G_compute=TRUE)}, then you must specify 0.25, 0.50, 0.75 or 1
    being the information fraction time at which you are deriving the
    marginal treatment effect. Defaults to 1.    
  }
  \item{ylim}{
    Sometimes in \code{type="hazard"} plots the credible interval
    ribbons become obnoxious. Hence, this argument. Set to c(0, Mx). 
  }
  \item{...}{
    Additional arguments to generic plot function
  } 
}
\value{
  A ggplot2 object.
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
  \code{\link{BayesFBHborrow}}, \code{\link{Combine}}
}
\examples{
  \dontrun{ 
    p_haz_wb <- plot(fit_sim_wb, type="hazard")
    p_gbcs_wb_trace <- plot(fit_gbcs_wb, type="trace", col="tamoxifen")
    p_gbcs_wb_TrtEff <- plot(fit_gbcs_wb, type="TrtEff")
    p_gbcs_wb_TrtEff_mgnl <- plot(update(fit_gbcs_wb, G_compute=TRUE), type="TrtEff")
    p_gbcs_wb_TrtEff_mgnl_75 <- plot(update(fit_gbcs_wb, G_compute=TRUE),
                                     type="TrtEff", InfFrac=0.75)
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}
\concept{plot methods}

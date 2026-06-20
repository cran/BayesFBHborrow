\name{update.BayesFBHborrow}
\alias{update.BayesFBHborrow}
\title{
 The Update Method for the \code{BayesFBHborrow} Class
}
\description{
 The Update Method for the \code{BayesFBHborrow} Class
}
\usage{
update.BayesFBHborrow(object, formula = NULL, data, data_hist = NULL,
                      subset, na.action, CntlOnly = FALSE, model_choice = "mix",
                      tuning_parameters = NULL, hyperparameters = NULL, iter = 1000,
                      warmup_iter = 100, refresh = 0, verbose = FALSE, max_grid = 2000,
                      standardise = TRUE, G_compute = FALSE, preprocess = TRUE,
                      fnm_cnd_blh = "cnd_blh", fnm_mgnl_haz0 =  "mgnl_haz0",
                      fnm_mgnl_haz1 = "mgnl_haz1", dbg = FALSE, ...) 
}
\arguments{
  \item{object}{
    An object of class \code{BayesFBHborrow}
  }
  \item{formula}{
    Optional but most likely preferred. When using the formula/data
    interface, specify the formula in the form typically used in other
    survival models e.g.     Surv(TOS, Event) ~ Trt + Oth + Covs + etc
    The first listed covariate must be the treatment variable. The
    rest can be continuous or factor.
  }
  \item{data}{
    Required. The current trial/study data.frame with a treated and a
    control arm. For the formula/data interface, there is no naming
    convention for variables as the user specifies these with the
    formula. For the data only interface, you must name the time to
    event variable 'tte', and the event variable 'event'. All
    covariates including the treatment indicator must have the 'X_'
    prefix.
  }
  \item{data_hist}{
    When fitting the full model with borrowing, the historical controls
    data.frame with variables named to agree with the current trial
    data.frame specified in the 'data' argument. Comments above
    regarding the formula/data versus data only interface apply here as
    well.
  }
  \item{subset}{
    When using the formula/data interface, you may also use
    \code{subset} and \code{na.action} as you do in other statistical
    routines in R.
  }
  \item{na.action}{
    When using the formula/data interface, you may also use
    \code{subset} and \code{na.action} as you do in other statistical
    routines in R.
  }
  \item{CntlOnly}{Set to \code{TRUE} when running the model with on
    the historical controls dataset only. The default value is
    \code{FALSE}.
  }
  \item{model_choice}{
    Borrowing method. Determine the borrowing per split point by
    specifying the option "mix" (default) for a mixture of
    inverse gamma densities or "uni" for a single inverse gamma
    density. Otherwise, determine the borrowing once for all split
    points by specifying the option "all".
  }
  \item{tuning_parameters}{list with components \code{cprop_beta}
    (\code{cprop_beta_0} for historical data), \code{a_lambda},
    \code{b_lambda}, \code{pi_b} and \code{alpha}. Default is
    \code{list(cprop_beta=1.1, cprop_beta_0=1.1, a_lambda=0.01,
    b_lambda=0.01, pi_b=0.5, alpha=0.4)}. The most important for the user
    are \code{cprop_beta} and \code{cprop_beta_0}. Starting at the
    default values, the user must adjust these values up or down in a
    series of test runs (\code{iter=200}) so that \code{beta_acc_ratio}
    and \code{beta_0_acc_ratio} (in the output list) are between 0.3 and
    0.4. The relationship between \code{cprop_beta} parameters and the 
    acceptance ratios is inverse.
  }
  \item{hyperparameters}{list with components \code{beta_prior},
    \code{beta_0_prior}, \code{a_tau}, \code{b_tau}, \code{c_tau},
    \code{d_tau}, \code{type}, \code{p_0}, \code{a_sigma},
    \code{b_sigma}, \code{phi}, \code{clam_smooth},
    \code{Jmax}. Default is \code{list(beta_prior=10^2,
    beta_0_prior=10^2, a_tau=1, b_tau=0.001, c_tau=1, d_tau=5,
    type="mix", p_0=0.8, a_sigma=1, b_sigma=1, phi=3, clam_smooth=0.7,
    Jmax=5)}. 
  }
  \item{iter}{
    Number of iterations for MCMC sampler.
  }
  \item{warmup_iter}{
    Number of warmup iterations (burn-in) for MCMC sampler.
  }
  \item{refresh}{
    Number of iterations between printed screen updates.
  }
  \item{verbose}{
    Whether or not to print the progress to the console
    every \code{refresh} iterations. Defaults to \code{FALSE}. 
  }
  \item{max_grid}{
    Grid size for the smoothed baseline hazard. Default value is 2000.
  }
  \item{standardise}{
    Whether or not to normalize the current trial baseline hazard to
    unity and apply the same normalization constant to the historical
    data. Default is \code{TRUE}. 
  }
  \item{G_compute}{
    If set to \code{TRUE} when treatment effect and pretreatment
    covariates are entered into the model (current trial only or current
    trial with borrowing), results are displayed as a marginal
    estimand computed via g-computation. This setting affects the
    behavior of the `coef` and `summary` methods and thus the `print`
    method. The default is \code{FALSE} which results in display of
    a conditional estimand.
  }
  \item{preprocess}{
    Controls the interpretation of the 'baseline' hazard consequent to
    recoding of the input dataset(s). When set to \code{TRUE}, all
    continuous variables will be centered at their mean over the dataset
    in the \code{data} argument, and factor variables will be centered
    to their mean under a hypothetical design having balanced levels for
    each factor. This defaults to \code{TRUE} and works with either the
    formula/data interface or the data only interface. If complete
    control over recoding is desired, one must use the data only
    interface together with \code{preprocess=FALSE}. This is the only
    way to include interactions.
  }
  \item{fnm_cnd_blh}{
    Binary files generated during run time which are given a default name
    prefix with a random hex suffix. The collection of
    "snapped to common grid" samples from the posterior distribution of
    the baseline hazards is a matrix of dimension \code{max_grid} by
    \code{iter}. Adding updated columns to this very large matrix would
    result in memory caching and for this reason, these are written out
    to binary files during the run of the MCMC, and then read back in
    during the final sumarization. This file is for the conditional
    baseline hazard. These files are available to the user for computing
    any functional of the posterior distribution. See the documentation
    for the \code{read_haz_mcmc_smpls} function.
  }
  \item{fnm_mgnl_haz0}{
    Binary files generated during run time which are given a default name
    prefix with a random hex suffix. The collection of
    "snapped to common grid" samples from the posterior distribution of
    the baseline hazards is a matrix of dimension \code{max_grid} by
    \code{iter}. Adding updated columns to this very large matrix would
    result in memory caching and for this reason, these are written out
    to binary files during the run of the MCMC, and then read back in
    during the final sumarization. This file is for the control arm
    log-log marginal hazard. These files are available to the user for
    computing any functional of the posterior distribution. See the
    documentation for the \code{read_haz_mcmc_smpls} function.
  }
  \item{fnm_mgnl_haz1}{
    Binary files generated during run time which are given a default name
    prefix with a random hex suffix. The collection of
    "snapped to common grid" samples from the posterior distribution of
    the baseline hazards is a matrix of dimension \code{max_grid} by
    \code{iter}. Adding updated columns to this very large matrix would
    result in memory caching and for this reason, these are written out
    to binary files during the run of the MCMC, and then read back in
    during the final sumarization. This file is for the treated arm
    log-log marginal hazard. These files are available to the user for
    computing any functional of the posterior distribution. See the
    documentation for the \code{read_haz_mcmc_smpls} function.
  }
  \item{dbg}{
    When running in command line mode, set to \code{TRUE} in order to
    enable line by line debugging.
  }
  \item{\dots}{
    Other arguments as needed.
  }
}
\details{
  Note that the marginal contrast is always computed in addition to the
  conditional estimand. If \code{G_compute=FALSE}, (the default) in the
  original call, the print, coef, summary and plot methods will
  produce results reflecting the conditional estimand, whereas these
  will reflect the marginal contrast otherwise. NOTE: calling 
  \code{update(my_bfbhb, G_compute=TRUE)} doesn't require rerunning the
  sampler as the results are already there. 
}
\value{
  An object of class \code{BayesFBHborrow}. See the documentation.
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
\seealso{\code{plot}, \code{summary}, \code{update} methods for class
  \code{"BayesFBHborrow"} objects
}
\examples{
    \dontrun{
       update(my_bfbhb, G_compute=TRUE)
    }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

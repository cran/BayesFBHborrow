\name{BayesFBHborrow}
\alias{BayesFBHborrow}
\title{
Fit the Bayesian Flexible Baseline Hazard Model with Borrowing via MCMC/MH
}
\description{
  Fits the Bayesian Flexible Baseline Hazard Model with Borrowing via
  MCMC/MH. Each sample in the Gibbs sampler is a piecewise exponential
  model, with prior over number and location of split
  points. Accommodates borrowing over a historical control arm as well
  as G-computation. 
}
\usage{
BayesFBHborrow(formula = NULL, data, data_hist = NULL, subset, na.action,
               CntlOnly = FALSE, model_choice = "mix", tuning_parameters = NULL,
               hyperparameters = NULL, iter = 1000, warmup_iter = 100,
               verbose = FALSE, max_grid = 2000, standardise = TRUE,
               G_compute = FALSE, preprocess = TRUE, fnm_cnd_haz_0 = "cnd_haz_0",
               fnm_cnd_cumHz_0 = "cnd_cumHz_0", fnm_lml_mgnl_surv_0 = "lml_mgnl_surv_0",
               fnm_lml_mgnl_surv_1 = "lml_mgnl_surv_1", fnm_mgnl_haz_0 = "mgnl_haz_0",
               fnm_mgnl_haz_1 = "mgnl_haz_1", sfx = NULL,
               keep_bin_files = FALSE, dbg = FALSE, ...)
}
\arguments{
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
  \item{verbose}{
    Whether or not you want to see the status messages.
    Defaults to \code{FALSE}. 
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
  \item{fnm_cnd_haz_0}{
    name prefix with a random hex suffix. The collection of
    "snapped to common grid" samples from the posterior distribution of
    the baseline hazards is a matrix of dimension \code{max_grid} by
    \code{iter}. Adding updated columns to this very large matrix would
    result in memory caching and for this reason, these are written out
    to binary files during the run of the MCMC, and then read back in
    during the final sumarization. These files are written to a hidden
    subdirectory named \code{.xxxxxx} (dot randome hex string) in the
    current working directory. Their complete paths are stored in the
    resulting object and are available to the user for computing any
    downstream functional of the posterior distribution. See the
    documentation for the \code{read_haz_mcmc_smpls} function. This file
    contains samples of the conditional baseline hazard. These are given
    default names so there's no need to specify them. Default behavior
    is to delete them unless the user specifies otherwise (see below).
  }
  \item{fnm_cnd_cumHz_0}{
    (see above) This file contains samples of the conditional baseline
    cumulative hazard.
  }
  \item{fnm_lml_mgnl_surv_0}{
    (see above) This file contains samples of the logged negative logged
    control arm marginal survival function.
  }
  \item{fnm_lml_mgnl_surv_1}{
    (see above) This file contains samples of the logged negative logged
    treated arm marginal survival function.
  }
  \item{fnm_mgnl_haz_0}{
    (see above) This file contains samples of the control arm marginal
    hazard which have been derived as the finite difference quotient of
    the exponentiated values in the lml_mgnl_surv_0 file.
  }
  \item{fnm_mgnl_haz_1}{
    (see above) This file contains samples of the treated arm marginal
    hazard which have been derived as the finite difference quotient of
    the exponentiated values in the lml_mgnl_surv_0 file. 
  }
  \item{sfx}{
    The name of the subdirectory created in the current working
    directory to contain these binary files. Also used in the suffix of
    the file names. The default behavior is to give this a random 6
    digit hex string so there's no need to specify this argument either.
  }
  \item{keep_bin_files}{
    Defaults to \code{FALSE} as the default behavior is to delete the
    entire temp binary file directory when the run concludes. If the
    user wishes to perform downstream calculation using the full set
    of MCMC/MH sampled functions snapped to the finely spaced grid
    then set this to \code{TRUE} and then use
    \code{read_haz_mcmc_samples} (see below). 
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
  BayesFBHborrow implements our semiparametric Bayesian borrowing model
  allowing inferences on a current trial with commensurate Bayesian
  borrowing from a set of historical controls. The model supports
  covariate adjusted borrowing, which can reduce prior-data conflict and
  improve power when differences in outcomes are attributable to changes
  in covariate distribution. As the treatment effect estimator is
  non-collapsible, the marginal hazard ratio can be estimated via
  Bayesian G-computation, while still permitting an adjusted analysis to
  account for control group drift. The package now fully implements
  parallel computation using parallelly, future and future.apply.
}
\value{
  \item{Y}{Times on study in the 'data' argument data.frame, scaled to
    unit event rate.}
  \item{I}{Event indicators in the 'data' argument data.frame}  
  \item{X}{Treatment indicator and pretreatment covariates in the 'data'
    argument data.frame, recoded if 'prepocess' is TRUE.}  
  \item{Y_0}{Times on study in the 'data_hist' argument data.frame,
    scaled to unit event rate.} 
  \item{I_0}{Event indicators in the 'data_hist' argument data.frame}
  \item{X_0}{Pretreatment covariates in the 'data_hist' argument
    data.frame, recoded if 'prepocess' is TRUE.} 
  \item{samples}{An 'iter' by ncol(X) + 3*J + 4 matrix containing all
    post warmup samples from the MCMC sampler: the model coefficients,
    the number of split points, hyperparameters mu and sigma2, the
    locations of all split points and the values of the piecewise
    exponential hazards corresponding to the current trial and to the
    historical controls data distributions.}
  \item{surv_dat}{Median together with 95 percent credible interval of
    posterior samples of the per arm conditional hazards/survival
    functions, per arm marginal hazard/survival functions, 
    and marginal treatment effect at each time point on the fine meshed
    grid of length 'max_grid'.}
  \item{MTE_dens_dfs}{kernel density estimates for samples of the
    posterior distribution of the marginal treatment effect at
    information fraction 25, 50, 75 and 100 percent.}
  \item{T.if}{Times on study at 25, 50, 75 and 100 percent information.} 
  \item{fnms}{File names of binary files of snapped to fine meshed grid
    hazard function samples from the posterior. Provided for calculating
    custom functionals of the hazard function posterior distribution.}
  \item{beta_acc_ratio}{Metropolis-Hastings algorithm acceptance ratio
    for current data model coefficients. Find value of cprop_beta in
    tuning_parameters which results in this being between 0.30 and 0.40.}
  \item{beta_0_acc_ratio}{Metropolis-Hastings algorithm acceptance ratio
    for historical data model coefficients. Find value of cprop_beta in 
    tuning_parameters which results in this being between 0.30 and 0.40.} 
  \item{call}{The call which produced the object.}
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
\note{
  This runs in parallel using the future and future.apply. You need do
  nothing special if you're on a mac, linux PC, WSL or Windows. If
  you're working on an HPC cluster at the command line, then you must
  request a specially enabled ineractive node via the following:
  \preformatted{
    srun --partition=interactive \
         --nodes=1 \
         --ntasks=1 \
         --cpus-per-task=16 \
         --mem=16G \
         --time=02:00:00 \
         --pty /bin/bash
  }
  Once in the interactive node, you must load R:
  \preformatted{
    module load R
  }
  or if you want a tricked out version of R, you need to lookup
  the name of the custom module created at your institution.
  You can look these up via the command:
  \preformatted{
    module spider R
  }
  Then once you've identified an appropriate module, then
  execute the appropriate \code{module load} command, e.g.
  at Astrazeneca:
  \preformatted{
    module load R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1
  }
  Otherwise if you have in browser workbench access to the HPC
  cluster, you will use the menu guided options to select the multiple
  cpus (or nodes) per task. 
}
\seealso{\code{plot}, \code{summary}, \code{update} methods for class
  \code{"BayesFBHborrow"} objects
}
\examples{
n_cc_1 <- 200
n_cc_0 <- 100
n_hst <- 100
int_cc <- int_hst <- -log(3)
B_x_cc <- B_x_hst <- c(-0.3,0.5,0.25,-0.5)
B_trt <- log(0.55)

dat_lst <- 
    genBFBHBdat(n_cc_1=200, n_cc_0=100, n_hst=100,
                B_trt=B_trt, B_x_cc=B_x_cc, B_x_hst=B_x_hst, 
                int_cc=int_cc, int_hst=int_hst, shape=2,
                t_er=0.5, t_fin=1.5, X_fact_levs=3)

DAT_cc <- trial_data(dat_lst)
DAT_hst <- historical_data(dat_lst)

## using default settings for hyperparameters and tuning parameters
## formula interface
## cntls only; no cov's
  \donttest{
    test_f0 <- BayesFBHborrow(Surv(tte, event)~1, data = DAT_hst,
                              CntlOnly=TRUE, iter = 200, warmup_iter = 50,
                              verbose = FALSE)

    ## cntls only; with cov's
    test_f1 <- BayesFBHborrow(Surv(tte, event) ~ X_01 + X_02 + X_03,
                              data = DAT_hst, CntlOnly=TRUE,
                              iter = 200, warmup_iter = 50, verbose = FALSE)
    
    ##  trtd+cntls; w/o borrow; no cov's: 
    test_f2 <- BayesFBHborrow(Surv(tte, event) ~ X_trt, data = DAT_cc,
                              iter = 200, warmup_iter = 50, verbose = FALSE)
    
    ## trtd+cntls; w borrow, no cov's:
    test_f3 <- BayesFBHborrow(Surv(tte, event) ~ X_trt, data = DAT_cc,
                              data_hist = DAT_hst, model_choice = 'mix',
                              iter = 200, warmup_iter = 50, verbose = FALSE)

    ##  trtd+cntls; w/o borrow; w cov's: 
    test_f4 <- BayesFBHborrow(Surv(tte, event) ~ X_trt + X_01 + X_02 + X_03,
                              data = DAT_cc, model_choice = 'mix',
                              iter = 200, warmup_iter = 50, verbose = FALSE)
    
    ## trtd+cntls; w borrow, w cov's:
    test_f5 <- BayesFBHborrow(Surv(tte, event) ~ X_trt + X_01 + X_02 + X_03,
                              data = DAT_cc, data_hist = DAT_hst,
                              model_choice = 'mix', iter = 10, warmup_iter = 5,
                              verbose = FALSE)
  }
}
\concept{Bayesian borrowing}
\concept{flexible baseline hazard}
\concept{time-to-event}
\concept{covariate adjustment}
\concept{G-computation}

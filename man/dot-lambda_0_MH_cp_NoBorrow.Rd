% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda_updates.R
\name{.lambda_0_MH_cp_NoBorrow}
\alias{.lambda_0_MH_cp_NoBorrow}
\title{Lambda_0 MH step, proposal from conditional conjugate posterior}
\usage{
.lambda_0_MH_cp_NoBorrow(
  df_hist,
  Y_0,
  I_0,
  X_0 = NULL,
  s,
  beta_0 = NULL,
  mu,
  sigma2,
  lambda_0,
  bp_0 = 0,
  J,
  clam,
  a_lam = 0.01,
  b_lam = 0.01,
  lambda_0_count = 0,
  lambda_0_move = 0
)
}
\arguments{
\item{df_hist}{data.frame from dataframe_fun()}

\item{Y_0}{historical trial data}

\item{I_0}{historical trial censoring indicator}

\item{X_0}{historical trial design matrix}

\item{s}{split point locations, (J+2)}

\item{beta_0}{parameter value for historical covariates}

\item{mu}{prior mean for baseline hazard}

\item{sigma2}{prior variance hyperparameter for baseline hazard}

\item{lambda_0}{baseline hazard}

\item{bp_0}{number of covariates, length(beta_0)}

\item{J}{number of split points}

\item{clam}{controls neighbor interactions, in range (0, 1)}

\item{a_lam}{lambda hyperparameter, default is 0.01}

\item{b_lam}{lambda hyperparameter, default is 0.01}

\item{lambda_0_count}{number of total moves for lambda_0}

\item{lambda_0_move}{number of accepted moves for lambda_0}
}
\value{
list of updated (if accepted) lambda_0 and data.frames, as well as the 
number of accepted moves
}
\description{
Lambda_0 MH step, proposal from conditional conjugate posterior
}

#' @title Loglikelihood ratio calculation for beta parameters
#' @description Compute log likelihood for beta update
#'
#' @param df data.frame from dataframe_fun()
#' @param beta beta values
#' @param beta_new proposed beta values
#'
#' @return likelihood ratio
.llikelihood_ratio_beta <- function(df, beta, beta_new) {
  
  X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
  xdpb <- X %*% beta
  xdpb_new <- X %*% beta_new

  llikelihood_ratio <- sum((xdpb_new  - xdpb) * df$I - 
                           ((df$Y - df$tstart) * df$lambda) * 
                           (exp(xdpb_new) - exp(xdpb)))
  
  return(llikelihood_ratio)
}


#' @title Beta Metropolis-Hastings Random walk move
#' @description Update beta via a Metropolis-Hastings Random Walk move
#'
#' @param df data.frame from dataframe_fun()
#' @param beta beta values
#' @param bp number of covariates
#' @param cprop_beta hyperparameter for beta proposal standard deviation
#' @param beta_count number of moves done for beta
#'
#' @return beta, either old or new move
.beta_MH_RW <- function(df, beta, bp, cprop_beta, beta_count) {
  for (k in 1:bp) {
    beta_new <- beta
    beta_prop <- stats::rnorm(1, beta[k], cprop_beta)
    beta_new[k] <- beta_prop

    logacc <- .llikelihood_ratio_beta(df, beta, beta_new)

    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      beta_count[k] <- beta_count[k] + 1
    }

  }

  return(list("beta" = beta, "beta_count" = beta_count))

}

#' @title Mean for MALA using derivative for beta proposal
#'
#' @param df Data frame with indicators
#' @param k index for beta
#' @param beta vector of parameters
#' @param bp number of covariates
#' @param cprop_beta proposal standard dev
#'
#' @return proposal mean
.beta_mom <- function(df, k, beta, bp, cprop_beta) {
  X <- as.matrix(df[, paste0("X", 1:bp)])
  xdpb <- X %*% beta
  x <- X[, k]
  
  D1 <- sum(df$I * x - (df$Y - df$tstart) * df$lambda * x * exp(xdpb))  
  mu_prop <- beta[k] + (cprop_beta[k]**2) / 2 * D1
  
  return(mu_prop) 
}


#' @title Log density of proposal for MALA
#'
#' @param beta_prop proposal beta
#' @param mu mean of proposal distribution
#' @param cprop_beta proposal standard dev
#'
#' @return log density
.lprop_density_beta <- function(beta_prop, mu, cprop_beta) {
  ldens <- (-1 / (2 * cprop_beta**2)) * (beta_prop - mu)**2
}



#' @title Proposal beta with a Metropolis Adjusted Langevin (MALA)
#'
#' @param df Data frame with indicators
#' @param beta vector of parameters
#' @param bp number of covariates
#' @param cprop_beta proposal variance standard deviation
#' @param beta_count count number of accepts
#'
#' @return updated beta vector
.beta_MH_MALA <- function(df, beta, bp, cprop_beta, beta_count) {
  for(k in 1:bp){
    beta_new <- beta
    mu_prop <- .beta_mom(df, k, beta, bp, cprop_beta)
    beta_prop <- stats::rnorm(n = 1, mean = mu_prop, sd = cprop_beta[k])
    beta_new[k] <- beta_prop
    
    mu_old <- .beta_mom(df, k, beta_new, bp, cprop_beta)
    
    log_prop_ratio <- .lprop_density_beta(mu_prop, beta[k], cprop_beta[k]) - 
      .lprop_density_beta(beta_prop, mu_old, cprop_beta[k])
    target_ratio <- .llikelihood_ratio_beta(df, beta, beta_new)
    
    logacc <- target_ratio - log_prop_ratio  
    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      beta_count[k] <- beta_count[k] + 1
    }
  }
  
  return(list("beta" = beta, "beta_count" = beta_count))
  
}



#' @title Fit frequentist piecewise exponential model for MLE and information matrix of beta 
#' @description Compute MLE for PEM 
#'
#' @param df Data frame with time-to-event, censoring indicator and covariates
#'
#' @return beta MLE and inverse of information matrix
.glmFit <- function(df){
 
  lenp <- length(df[1, grepl("X", names(df))])
  lab <- paste0("X", 1:lenp)
  splits <- length(unique(df$tstart))
  
  if(splits == 1){
    alllab <- paste(paste(lab, collapse= "+"), "+ offset(log(Y))")
    fmla <- stats::as.formula(paste("I ~", paste(paste(lab, collapse= "+"), "+ offset(log(Y))")))
  }else{
    df$tstart <- as.factor(df$tstart)
    alllab <- paste(paste(lab, collapse= "+"), "+ offset(log(Y))")
    fmla <- stats::as.formula(paste("I ~ tstart +", paste(paste(lab, collapse= "+"), "+ offset(log(Y))")))
  }
  
  #fit PWE model
  fit <- stats::glm(fmla, data = df)
  ss.x <- grepl("X", names(fit$coefficients))
  
  beta.mu <- fit$coefficients[ss.x]
  beta.vcov <- stats::vcov(fit)[ss.x, ss.x]

  return(list("beta.mu" = beta.mu, "beta.vcov" = beta.vcov))
  
}



#' @title Beta MH RW sampler from freq PEM fit
#' @description Sample beta from RW sampler  
#'
#' @param df Data frame with indicators
#' @param beta vector of parameters
#' @param beta_count count number of accepted proposals
#' @param cprop_beta proposal scalar
#' 
#' @return beta, either old or new move
.beta.MH.RW.glm <- function(df, beta, beta_count, cprop_beta){
 
  glm.mom <- .glmFit(df)
  beta.new <- beta
  cd2 <- cprop_beta**2 / length(beta)
  
  beta.prop <- as.vector(mvtnorm::rmvnorm(1, mean = beta, sigma = cd2 * glm.mom$beta.vcov))
  
  for (k in 1:length(beta.new)) {
    
    beta.new[k] <- beta.prop[k]
    logacc <- .llikelihood_ratio_beta(df, beta, beta.new)
    
    if (logacc > log(stats::runif(1))) {
      beta <- beta.new
      beta_count[k] <- beta_count[k] + 1
    }
    
  }    
  
  return(list("beta" = beta, "beta_count" = beta_count))
  
}


#' @title log Gaussian proposal density for Newton Raphson proposal
#'
#' @param beta.prop beta proposal
#' @param mu_old density mean
#' @param var_old density variance
#' 
#' @return log Gaussian density
.lprop.dens.beta.NR <- function(beta.prop, mu_old, var_old){
  ldens <- (-1 / (2 * var_old)) * (beta.prop - mu_old)**2
  return(ldens)
}




#' @title First and second derivative of target for mode and variance of proposal
#'
#' @param df Data frame with indicators
#' @param k index 
#' @param beta vector of parameters
#' @param bp number of covariates
#' @param cprop_beta proposal variance standard deviation
#'
#' @return First and second derivative mode and variance 
.beta_mom.NR.fun <- function(df, k, beta, bp, cprop_beta) {

  bp <- length(beta)
  X <- as.matrix(df[, paste0("X", 1:bp)])
  xdpb <- X %*% beta
  x <- X[, k]
  
  D1 <- sum(df$I * x - (df$Y - df$tstart) * df$lambda * x * exp(xdpb))  
  D2 <- - sum(df$lambda * (df$Y - df$tstart) * x**2  * exp(xdpb))
  mu <- beta[k] - D1 / D2
  var <- -cprop_beta**2 / D2
  
  return(list("D1" = D1, "D2" = D2, "mu" = mu, "var" = var))  
  
}



#' @title Newton Raphson MH move 
#' @description Sample beta from RW sampler
#'
#' @param df Data frame with indicators
#' @param beta vector of parameters
#' @param bp number of covariates
#' @param cprop_beta proposal scalar
#' @param beta_count count number of accepts
#'
#' @return updated beta
.beta_MH_NR <- function(df, beta, bp, cprop_beta, beta_count){
  
  for(k in 1:bp){
    
    beta.new <- beta
    
    mom.prop <- .beta_mom.NR.fun(df, k, beta, bp, cprop_beta)
    beta.prop <- stats::rnorm(n = 1, mean = mom.prop$mu, sd = sqrt(mom.prop$var))
    beta.new[k] <- beta.prop
    
    mom.old <- .beta_mom.NR.fun(df, k, beta.new, bp, cprop_beta)
    
    log_prop_ratio <- .lprop.dens.beta.NR(beta[k], mom.prop$mu, mom.prop$var) - .lprop.dens.beta.NR(beta.prop,  mom.old$mu,  mom.old$var)
    target_ratio <- .llikelihood_ratio_beta(df, beta, beta.new)
    
    logacc <- target_ratio - log_prop_ratio  
    
    if(logacc > log(stats::runif(1))){
      beta[k] <- beta.prop
      beta_count[k] <- beta_count[k] + 1
    }
    
  }
  
  return(list("beta" = beta, "beta_count" = beta_count))
  
}









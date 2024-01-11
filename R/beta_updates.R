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
#' @param count_beta number of moves done for beta
#'
#' @return beta, either old or new move
.beta_MH_RW <- function(df, beta, bp, cprop_beta, count_beta) {
  for (k in 1:bp) {
    beta_new <- beta
    beta_prop <- stats::rnorm(1, beta[k], cprop_beta)
    beta_new[k] <- beta_prop

    logacc <- .llikelihood_ratio_beta(df, beta, beta_new)

    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      count_beta[k] <- count_beta[k] + 1
    }

  }

  return(list("beta" = beta, "count_beta" = count_beta))

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



#' @title Proposal beta with a Metrolopis Adjusted Langevin (MALA)
#'
#' @param df Data frame with indicators
#' @param beta vector of parameters
#' @param bp number of covariates
#' @param cprop_beta proposal variance standard deviation
#' @param count_beta count number of accepts
#'
#' @return updated beta vector
.beta_MH_MALA <- function(df, beta, bp, cprop_beta, count_beta) {
  for(k in 1:bp){
    
    beta_new <- beta
    
    mu_prop <- .beta_mom(df, k, beta, bp, cprop_beta)
    beta_prop <- stats::rnorm(n = 1, mean = mu_prop, sd = cprop_beta[k])
    beta_new[k] <- beta_prop
    
    mu_old <- .beta_mom(df, k, beta_new, bp, cprop_beta)
    
    log_prop_ratio <- .lprop_density_beta(beta[k], mu_prop, cprop_beta[k]) - 
      .lprop_density_beta(beta_prop, mu_old, cprop_beta[k])
    target_ratio <- .llikelihood_ratio_beta(df, beta, beta_new)
    
    logacc <- target_ratio - log_prop_ratio  
    
    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      count_beta[k] <- count_beta[k] + 1
    }
    
  }
  
  return(list("beta" = beta, "count_beta" = count_beta))
  
}
#' @title Calculate sigma2 posterior update
#'
#' @param mu mean.
#' @param lambda_0 Baseline hazard.
#' @param Sigma_s VCV matrix (j + 1) x (j + 1).
#' @param J Number of split point.
#' @param a_sigma Hyperparameter a.
#' @param b_sigma Hyperparameter b.
#'
#' @return sigma2 draw from IG
.sigma2_update <- function(mu, lambda_0, Sigma_s, J, a_sigma, b_sigma) {
  a_post <- a_sigma + (J + 1) / 2
  
  one <- rep(1, J + 1)
  cp <- t(mu * one - log(lambda_0)) %*% solve(Sigma_s) %*% (mu * one - log(lambda_0))
  b_post <- b_sigma + cp / 2   
  
  sigma2 <- invgamma::rinvgamma(1, shape = a_post, rate = b_post)
}

#' @title Calculate mu posterior update
#'
#' @param Sigma_s VCV matrix (j + 1) x (j + 1).
#' @param lambda_0 Baseline hazard.
#' @param sigma2 Scale variance.
#' @param J  Number of split point.
#'
#' @return mu update from Normal.
.mu_update <- function(Sigma_s, lambda_0, sigma2, J) {
  one <- rep(1, J + 1)
  mu_num <- (t(one) %*% solve(Sigma_s) %*% log(lambda_0)) 
  mu_den <- t(one) %*% solve(Sigma_s) %*% one
  mu_mu <- mu_num / mu_den
  
  mu_var <- sigma2 / mu_den
  
  mu <- stats::rnorm(1, mean = mu_mu, sd = sqrt(mu_var))
}

#' @title Calculate covariance matrix in the MVN-ICAR
#'
#' @param s split points, J + 2
#' @param J number of split points
#' @param clam controls neighbor interactions, in range (0, 1)
#'
#' @return Sigma_s = (I - W)^(-1) * Q, W, Q
.ICAR_calc <- function(s, J, clam) {
  
  W <- matrix(rep(0,(J + 1) * (J + 1)), nrow = J + 1)
  Q <- matrix(rep(0,(J + 1) * (J + 1)), nrow = J + 1)
  
  interval_length <- diff(s[!(is.na(s))])
  
  if (J < 2) {
    if (J == 1) {
      
      W[1, 2] <- clam * (interval_length[1] + interval_length[2]) / (2 * interval_length[1] + interval_length[2])
      W[2, 1] <- clam * (interval_length[2]+ interval_length[1]) / (interval_length[1] + 2 * interval_length[2])
      Q[1, 1] <- 2 / (2 * interval_length[1] + interval_length[2])
      Q[2, 2] <- 2 / (interval_length[1] + 2 * interval_length[2])
      Sigma_s <- solve(diag(J + 1) - W) %*% Q
      
    } else {
      
      Sigma_s <- as.matrix(1)
      
    }  
  } else {
    
    for (j in 2:J) {
      
      W[j, j + 1] <- clam * (interval_length[j] + interval_length[j + 1]) / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])
      W[j, j-1] <- clam * (interval_length[j] + interval_length[j-1]) / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])
      Q[j, j] <- 2 / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])
      
    }
    
    Q[j + 1, j + 1] <- 2 / (interval_length[J] + 2 * interval_length[j + 1])
    Q[1, 1] <- 2 / (2 * interval_length[1] + interval_length[2])
    W[1, 2] <- clam * (interval_length[1] + interval_length[2]) / (2 * interval_length[1] + interval_length[2])
    W[j + 1, J] <- clam * (interval_length[j + 1] + interval_length[J]) / (interval_length[J]+ 2 * interval_length[j + 1])
    
    Sigma_s <- solve(diag(j + 1) - W) %*% Q
    
  }
  
  return(list("Sigma_s" = Sigma_s, "W" = W, "Q" = Q))
  
}

#' @title Sample tau from posterior distribution
#'
#' @param lambda_0 historical baseline hazard
#' @param lambda baseline hazard
#' @param J number of split points
#' @param s split point locations, J + 2
#' @param a_tau Inverse Gamma hyperparameter
#' @param b_tau Inverse Gamma hyperparameter
#' @param c_tau Inverse Gamma hyperparameter
#' @param d_tau Inverse Gamma hyperparameter
#' @param p_0 mixture ratio
#' @param type choice of borrowing, "mix", "uni", or any other string for 
#' borrowing on every baseline hazard without mixture
#'
#' @return list containing tau and new mixture ratio
.tau_update <- function(lambda_0, lambda, J, s, 
                    a_tau, b_tau, c_tau = NULL, d_tau = NULL, 
                    p_0 = NULL, type) {
  
  sq_diff <- (log(lambda_0) - log(lambda))**2 
  
  if (type == "mix") {
    
    #compute prob on the log scale
    lw_0_num <- a_tau * log(b_tau) + lgamma(a_tau + 0.5) 
    lw_0_den <- (0.5 + a_tau) * log((sq_diff / 2) + b_tau) + lgamma(a_tau) 
    lw_1_num <- c_tau * log(d_tau) + lgamma(c_tau + 0.5) 
    lw_1_den <- (0.5 + c_tau) * log((sq_diff / 2) + d_tau) + lgamma(c_tau)
    
    p_0_new <- log(p_0) +  lw_0_num - lw_0_den
    p_1_new <- log(1 - p_0) +  lw_1_num - lw_1_den
    
    probability_mat <- cbind(p_0_new, p_1_new)
    
    # normalize with log sum exp trick - avoid overflow
    p_new <- apply(probability_mat, 1, .normalize_prob)
    
    # sample mixture 
    comp <- apply(p_new, 2, sample, x = 1:2, size = 1, replace = F)
    
    # hyperparameters
    ac  <- matrix(rep(c(0.5 + a_tau, 0.5 + c_tau), J + 1), nrow = J + 1, byrow = T)
    bd <- cbind(sq_diff / 2 + b_tau, sq_diff / 2 + d_tau)  
    
    call <- cbind(1:(J + 1), comp)
    tau <- invgamma::rinvgamma(n = J + 1, shape = ac[call], rate = bd[call])
    
  }else if (type == "uni") {
    
    shape_tau <- 0.5 + a_tau
    rate_tau <- sq_diff / 2 + b_tau
    
    # placeholder
    p_new <- 1
    
    tau <- invgamma::rinvgamma(n = J + 1, shape = shape_tau, rate = rate_tau)
    
    
  } else {
    
    sq_diff_all <- sum(sq_diff)  
    
    # compute prob on the log scale
    lw_0_num <- a_tau * log(b_tau) + lgamma(a_tau + (J + 1) / 2) 
    lw_0_den <- ((J + 1) / 2 + a_tau) * log((sq_diff_all / 2) + b_tau) + lgamma(a_tau) 
    lw_1_num <- c_tau * log(d_tau) + lgamma(c_tau + (J + 1) / 2)
    lw_1_den <- ((J + 1) / 2 + c_tau) * log((sq_diff_all / 2) + d_tau) + lgamma(c_tau) 
    
    p_0_new <- log(p_0) +  lw_0_num - lw_0_den
    p_1_new <- log(1 - p_0) +  lw_1_num - lw_1_den
    
    probability_mat <- cbind(p_0_new, p_1_new)
    
    # normalize with log sum exp trick - avoid overflow
    p_new <- apply(probability_mat, 1, .normalize_prob)
    
    # sample mixture 
    mix <- sample(x = 1:2, size = 1, replace = F, prob = p_new)
    
    # hyperparameters
    ac  <- c((J + 1) / 2 + a_tau, (J + 1) / 2 + c_tau)
    bd <- c(sq_diff_all / 2 + b_tau, sq_diff_all / 2 + d_tau)  
    
    tau <- invgamma::rinvgamma(n = 1, shape = ac[mix], rate =  bd[mix])

  }
  
  return(list("tau" = tau, "p_new" = p_new))
  
}

#' @title Metropolis Hastings step: shuffle the split point locations (with 
#' Bayesian borrowing)
#'
#' @param df_hist dataframe containing historical trial data and parmaeters
#' @param df_curr data.frame containing current trial data and parameters
#' @param Y_0 historical trial data
#' @param I_0 historical trial censoring indicator
#' @param X_0 historical trial design matrix
#' @param lambda_0 historical baseline hazard
#' @param beta_0 historical parameter vector
#' @param Y data
#' @param I censoring indicator
#' @param X design matrix
#' @param lambda baseline hazard
#' @param beta parameter vector
#' @param s split point locations, J + 2
#' @param J number of split points
#' @param bp number of covariates in current trial
#' @param bp_0 number of covariates in historical trial
#' @param clam_smooth neighbor interactions, in range (0, 1), for ICAR update
#' @param maxSj the smallest of the maximal time points, min(max(Y), max(Y_0))
#'
#' @return list containing new split points, updated Sigma_s and data.frames
#' for historic and current trial data 
.shuffle_split_point_location <- function(df_hist, df_curr, Y_0, I_0, X_0, lambda_0, 
                                    beta_0, Y, I, X, lambda, beta, s, J, bp_0, 
                                    bp, clam_smooth, maxSj) {
  Sigma_s <- .ICAR_calc(s, J, clam_smooth)$Sigma_s
  
  #star individual proposals
  #prop vector of proposals and current
  for (j in 1:J) {
    if (j == J) {
      s_star <- stats::runif(1, min = s[j], max = maxSj)
    } else {
      s_star <- stats::runif(1, min = s[j], max = s[j + 2])
    }
    s_prop <- s 
    s_prop[j + 1] <- s_star
    
    # Create new data.frames for s_prop
    df_hist_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda_0, bp = bp_0, J = J)
    df_curr_prop <- .dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda, bp = bp, J = J)
    
    # Update ICAR
    Sigma_s_prop <- .ICAR_calc(s_prop, J, clam_smooth)$Sigma_s 
    
    # Probability of accepting
    if (j == J) {
      lprior_num <- log(maxSj - s_star) + log(s_star - s[j]) 
      lprior_den <- log(maxSj - s[j + 1]) + log(s[j + 1] - s[j]) 
    } else {
      lprior_num <- log(s[j + 2] - s_star) + log(s_star - s[j]) 
      lprior_den <- log(s[j + 2] - s[j + 1]) + log(s[j + 1] - s[j])  
    }  
    
    llike_num <- .log_likelihood(df_hist_prop, beta_0) + .log_likelihood(df_curr_prop, beta)
    llike_den <- .log_likelihood(df_hist, beta_0) + .log_likelihood(df_curr, beta)
    
    # Acceptance ratio
    logacc <- llike_num - llike_den + lprior_num - lprior_den 

    if (logacc > log(stats::runif(1))) {
      Sigma_s <- Sigma_s_prop 
      df_hist <- df_hist_prop
      df_curr <- df_curr_prop
      s <- s_prop
    }
    
  }
  
  return(list("s" = s, "Sigma_s" = Sigma_s, "df_hist" = df_hist, "df_curr" = df_curr))
  
}


#' @title Metropolis Hastings step: shuffle the split point locations (without 
#' Bayesian borrowing)
#'
#' @param df dataframe containing trial data and parameters
#' @param Y_0 data
#' @param I_0 censoring indicator
#' @param X_0 design matrix
#' @param lambda_0 baseline hazard
#' @param beta_0 parameter vector
#' @param s split point locations, J + 2
#' @param J number of split points
#' @param bp_0 number of covariates in historical trial
#' @param clam_smooth neighbor interactions, in range (0, 1), for ICAR update
#'
#' @return list containing new split points, updated Sigma_s and data.frames
#' for historic and current trial data 
.shuffle_split_point_location_NoBorrow <- function(df, Y_0, I_0, X_0, 
                     lambda_0, beta_0, s, J, 
                     bp_0, clam_smooth) {
  
  Sigma_s <- .ICAR_calc(s, J, clam_smooth)$Sigma_s
  
  #star individual proposals
  #prop vector of proposals and current
  
  for (j in 1:J) {
    
    s_star <- stats::runif(1, min = s[j], max = s[j + 2])
    s_prop <- s 
    s_prop[j + 1] <- s_star
    
    ##like
    df_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda_0, bp = bp_0, J = J)
    
    #ICAR
    Sigma_s_prop <- .ICAR_calc(s_prop, J, clam_smooth)$Sigma_s 
    
    #Prob of accepting
    lprior_num <- log(s[j + 2] - s_star) + log(s_star - s[j]) 
    lprior_denom <- log(s[j + 2] - s[j + 1]) + log(s[j + 1] - s[j]) 
    
    llike_num <- .log_likelihood(df_prop, beta_0)
    llike_den <- .log_likelihood(df, beta_0) 
    
    #Prob 
    logacc <- llike_num - llike_den + lprior_num - lprior_denom 
    
    if (logacc > log(stats::runif(1))) {
      Sigma_s <- Sigma_s_prop 
      df <- df_prop
      s <- s_prop
    }
    
  }
  
  return(list("s" = s, "Sigma_s" = Sigma_s, "df_all" = df))
}
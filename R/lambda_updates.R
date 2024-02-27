#' @title Propose lambda from a gamma conditional conjugate posterior proposal
#'
#' @param df data.frame from dataframe_fun()
#' @param beta parameter value for beta
#' @param j current split point
#' @param bp number of covariates
#' @param alam lambda hyperparameter, default set to 0.01
#' @param blam lambda hyperparameter, default set to 0.01
#'
#' @return list containing proposed lambda, shape and rate parameters
.lambda_conj_prop <- function(df, beta, j, bp,  alam = 0.01, blam = 0.01) {
  
  indx <- unique(df$tstart)
  df_ss <- df[df$tstart == indx[j],  ]
  
  if(!is.null(beta)) {
    X <- as.matrix(df_ss[, paste0("X", 1:bp)])
    xdpb <- X %*% beta
    rate_prop <- blam + sum((df_ss$Y - df_ss$tstart) * exp(xdpb))
  }else{
    rate_prop <- blam + sum((df_ss$Y - df_ss$tstart))
  } 
  shape_prop <- alam + sum(df_ss$I) 
  
  lambda_prop <- 0
  # if sum(df_ss) = 0, this causes lambda_prop --> 0, which will cause NA in logacc
  while (lambda_prop == 0) {
  lambda_prop <- stats:: rgamma(1, shape = shape_prop, rate = rate_prop)
  }

  return(list("lambda_prop" = lambda_prop, "shape_prop" = shape_prop, "rate_prop" = rate_prop))
  
}

#' @title Log likelihood for lambda / lambda_0 update
#'
#' @param df data.frame from dataframe_fun()
#' @param df_prop proposal data.frame
#' @param beta parameter value for beta
#'
#' @return log likelihood ratio for lambda
.llikelihood_ratio_lambda <- function(df, df_prop, beta) {
  if(!is.null(beta)) {
    X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
    xdpb <- X %*% beta
    
    llikelihood_ratio <- sum((log(df_prop$lambda) - log(df$lambda)) * df$I - 
                         (df$Y - df$tstart) * (df_prop$lambda - df$lambda) * exp(xdpb))  
    
  }else{
    llikelihood_ratio <- sum((log(df_prop$lambda) - log(df$lambda)) * df$I  - 
                         (df$Y - df$tstart) * (df_prop$lambda - df$lambda)) 
  }
  return(llikelihood_ratio)
}

#' @title Calculates nu and sigma2 for the Gaussian Markov random field prior, 
#' for a given split point j
#'
#' @param j current split point
#' @param lambda_0 historical baseline hazard
#' @param mu prior mean for baseline hazard
#' @param sigma2 prior variance hyperparameter for baseline hazard
#' @param W influence from right and left neighbors
#' @param Q individual effect of neighborhood
#' @param J number of split points
#'
#' @return nu and sigma2
.nu_sigma_update <- function(j, lambda_0, mu, sigma2, W, Q, J) {
  
  nu <- mu + (lambda_0 - rep(mu, nrow(W))) %*% W[j,]
  
  if(J > 0) {
    sigma2j <- sigma2 * Q[j,j]
  }else{
    sigma2j <- sigma2 * 1
  }
  
  return(list("nu" = nu, "sigma2j" = sigma2j))
  
}

#' @title Calculate log gamma ratio for two different parameter values
#'
#' @param x1 old parameter value
#' @param x2 proposed parameter value
#' @param shape shape parameter
#' @param rate rate parameter
#'
#' @return log gamma ratio
.lgamma_ratio <- function(x1, x2, shape, rate) {
  (shape - 1) * log(x1) - rate * x1 - (shape - 1) * log(x2) + rate * x2
}

#' @title Lambda_0 MH step, proposal from conditional conjugate posterior
#'
#' @param df_hist data.frame from dataframe_fun()
#' @param Y_0 historical trial data
#' @param I_0 historical trial censoring indicator
#' @param X_0 historical trial design matrix
#' @param s split point locations, (J+2)
#' @param beta_0 parameter value for historical covariates
#' @param mu prior mean for baseline hazard
#' @param sigma2 prior variance hyperparameter for baseline hazard
#' @param lambda baseline hazard
#' @param lambda_0 historical baseline hazard
#' @param tau borrowing parameter
#' @param bp_0 number of covariates, length(beta_0)
#' @param J number of split points
#' @param clam controls neighbor interactions, in range (0, 1)
#' @param a_lam lambda hyperparameter, default is 0.01
#' @param b_lam lambda hyperparameter, default is 0.01
#' @param lambda_0_count number of total moves for lambda_0
#' @param lambda_0_move number of accepted moves for lambda_0
#'
#' @return list of updated (if accepted) lambda_0 and data.frames, as well as the 
#' number of accepted moves 
.lambda_0_MH_cp <- function(df_hist, Y_0, I_0, X_0 = NULL, s, beta_0 = NULL,
                           mu, sigma2,  lambda, lambda_0, tau, bp_0 = 0, J,
                           clam, a_lam = 0.01, b_lam = 0.01, lambda_0_count = 0,
                           lambda_0_move = 0) {
  ICAR <- .ICAR_calc(s, J, clam) 
  Sigma_s <- ICAR$Sigma_s
  Q <- ICAR$Q
  W <- ICAR$W
  
  for (j in 1:(J + 1)) {
    
    lambda_0_new <- lambda_0   
    lambda_0_prop_all <- .lambda_conj_prop(df_hist, beta = beta_0, j, bp = bp_0, alam = a_lam,
                                           blam = b_lam)  
    lambda_0_prop <- lambda_0_prop_all$lambda_prop
    
    lambda_0_new[j] <- lambda_0_prop 
    shape_prop <- lambda_0_prop_all$shape_prop
    rate_prop <- lambda_0_prop_all$rate_prop
    
    df_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0_new, bp = bp_0, J = J)  

    nu_sigma <-.nu_sigma_update(j, lambda_0, mu, sigma2, W, Q, J)
    
    llikelihood_ratio <- .llikelihood_ratio_lambda(df_hist, df_prop, beta_0)

    log_prop_ratio <- .lgamma_ratio(x1 = lambda_0[j], x2 = lambda_0_prop, shape = shape_prop, rate = rate_prop)

    target_num <- stats::dnorm(log(lambda_0_prop), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                  log(lambda_0_prop)                                                                                                                           
    
    target_den <- stats::dnorm(log(lambda_0[j]), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                  log(lambda_0[j])

    if(length(tau) > 1) {
      target_num <- target_num + 
        stats::dnorm(log(lambda[j]), log(lambda_0_prop), sqrt(tau[j]), log = T)
        
      target_den <- target_den +
        stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau[j]), log = T)
    }else{
      target_num <- target_num +
        stats::dnorm(log(lambda[j]), log(lambda_0_prop), sqrt(tau), log = T)
        
      target_den <- target_den +
        stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau), log = T)
    }

    logacc <- llikelihood_ratio + target_num - target_den + log_prop_ratio
    
    if(logacc > log(stats::runif(1))) {
      lambda_0 <- lambda_0_new
      df_hist <- df_prop 
      lambda_0_move <- lambda_0_move + 1
    }
    
    lambda_0_count <- lambda_0_count + 1
  }
  
  return(list("lambda_0" = lambda_0, "df_hist" = df_hist, "lambda_0_count" = lambda_0_count, "lambda_0_move" = lambda_0_move))
  
}

#' @title Lambda_0 MH step, proposal from conditional conjugate posterior
#'
#' @param df_hist data.frame from dataframe_fun()
#' @param Y_0 historical trial data
#' @param I_0 historical trial censoring indicator
#' @param X_0 historical trial design matrix
#' @param s split point locations, (J+2)
#' @param beta_0 parameter value for historical covariates
#' @param mu prior mean for baseline hazard
#' @param sigma2 prior variance hyperparameter for baseline hazard
#' @param lambda_0 baseline hazard
#' @param bp_0 number of covariates, length(beta_0)
#' @param J number of split points
#' @param clam controls neighbor interactions, in range (0, 1)
#' @param a_lam lambda hyperparameter, default is 0.01
#' @param b_lam lambda hyperparameter, default is 0.01
#' @param lambda_0_count number of total moves for lambda_0
#' @param lambda_0_move number of accepted moves for lambda_0
#'
#' @return list of updated (if accepted) lambda_0 and data.frames, as well as the 
#' number of accepted moves 
.lambda_0_MH_cp_NoBorrow <- function(df_hist, Y_0, I_0, X_0 = NULL, s,
                                    beta_0 = NULL, mu, sigma2,  lambda_0, 
                                    bp_0 = 0, J, clam, a_lam = 0.01,
                                    b_lam = 0.01, lambda_0_count = 0,
                                    lambda_0_move = 0) {
  ICAR <- .ICAR_calc(s, J, clam) 
  Sigma_s <- ICAR$Sigma_s
  Q <- ICAR$Q
  W <- ICAR$W
  
  for (j in 1:(J + 1)) {
    
    lambda_0_new <- lambda_0   
    lambda_0_prop_all <- .lambda_conj_prop(df_hist, beta = beta_0, j, bp = bp_0, alam = a_lam,
                                           blam = b_lam)  
    lambda_0_prop <- lambda_0_prop_all$lambda_prop
    
    lambda_0_new[j] <- lambda_0_prop 
    shape_prop <- lambda_0_prop_all$shape_prop
    rate_prop <- lambda_0_prop_all$rate_prop
    
    df_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0_new, bp = bp_0, J = J)  
    
    nu_sigma <-.nu_sigma_update(j, lambda_0, mu, sigma2, W, Q, J)
    
    llikelihood_ratio <- .llikelihood_ratio_lambda(df_hist, df_prop, beta_0)
      log_prop_ratio <- stats::dgamma(lambda_0[j], shape = shape_prop, rate = rate_prop, log = T) -
                        stats::dgamma(lambda_0_prop, shape = shape_prop, rate = rate_prop, log = T)  
      
      target_num <- .log_likelihood(df_prop, beta_0) +
                      stats::dnorm(log(lambda_0_prop), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                      log(lambda_0_prop)
      
      target_den <- .log_likelihood(df_hist, beta_0) +
                      stats::dnorm(log(lambda_0[j]), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) - 
                      log(lambda_0[j])
      
      logacc <- target_num - target_den + log_prop_ratio
    
    if(logacc > log(stats::runif(1))) {
      lambda_0 <- lambda_0_new
      df_hist <- df_prop 
      lambda_0_move <- lambda_0_move + 1
    }
    
      lambda_0_count <- lambda_0_count + 1
  }
  
  return(list("lambda_0" = lambda_0, "df_hist" = df_hist, "lambda_0_count" = lambda_0_count, "lambda_0_move" = lambda_0_move))
  
}

#' @title Lambda MH step, proposal from conditional conjugate posterior
#'
#' @param df_hist data.frame from dataframe_fun()
#' @param df_curr data.frame from dataframe_fun()
#' @param Y  data
#' @param I  censoring indicator
#' @param X  design matrix
#' @param s split point locations, J + 2
#' @param beta parameter value for covariates
#' @param beta_0 parameter value for historical covariates
#' @param mu 
#' @param mu prior mean for baseline hazard
#' @param sigma2 prior variance hyperparameter for baseline hazard
#' @param lambda baseline hazard
#' @param lambda_0 historical baseline hazard
#' @param tau borrowing parameter
#' @param bp number of covariates, length(beta)
#' @param bp_0 number of covariates, length(beta_0)
#' @param J number of split points
#' @param a_lam lambda hyperparameter
#' @param b_lam lambda hyperparameter
#' @param lambda_move number of accepted lambda moves
#' @param lambda_count total number of lambda moves
#' @param alpha power parameter
#'
#' @return list of updated (if accepted) lambda and data.frames, as well as the 
#' number of accepted moves 
.lambda_MH_cp <- function(df_hist, df_curr, Y, I, X, s, beta, beta_0 = NULL, mu, sigma2, lambda, lambda_0, tau, 
                         bp, bp_0 = 0, J, a_lam = 0.01, b_lam = 0.01, lambda_move = 0,
                         lambda_count = 0, alpha = 0.3) {

  for (j in 1:(J + 1)) {
    
    lambda_new <- lambda
    lambda_prop_cc <- .lambda_conj_prop(df_curr, beta, j, bp = bp, alam = a_lam, 
                                        blam = b_lam) 
    cc_shape_prop <- lambda_prop_cc$shape_prop
    cc_rate_prop <- lambda_prop_cc$rate_prop

    lambda_prop_hist <- .lambda_conj_prop(df_hist, beta_0, j, bp = bp_0, alam = a_lam, 
                                          blam = b_lam) 
    hist_shape_prop <- lambda_prop_hist$shape_prop
    hist_rate_prop <- lambda_prop_hist$rate_prop

    shape_prop <- a_lam + cc_shape_prop + alpha * hist_shape_prop
    rate_prop <- b_lam + cc_rate_prop + alpha * hist_rate_prop
    
    lambda_prop <- stats:: rgamma(1, shape = shape_prop, rate = rate_prop)
    lambda_new[j] <- lambda_prop
    
    df_prop <- .dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda_new, bp = bp, J = J)        
    
    llikelihood_ratio <- .llikelihood_ratio_lambda(df_curr, df_prop, beta)
    log_prop <- .lgamma_ratio(x1 = lambda[j], x2 = lambda_prop, shape = shape_prop, rate = rate_prop)
    
    target_num <- (- log(lambda_prop))  
    target_den <- (- log(lambda[j]))
    
    # Adjust for non piecewise tau
    if(length(tau) > 1) {
      target_num <- target_num + stats::dnorm(log(lambda_prop), log(lambda_0[j]), sqrt(tau[j]), log = T)
      target_den <- target_den + stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau[j]), log = T)
    }else{
      target_num <- target_num + stats::dnorm(log(lambda_prop), log(lambda_0[j]), sqrt(tau), log = T)
      target_den <- target_den + stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau), log = T)
    }
    
    logacc <- llikelihood_ratio + target_num - target_den + log_prop  

    if(logacc > log(stats::runif(1))) {
      lambda <- lambda_new
      df_curr <- df_prop
      lambda_move <- lambda_move + 1  
    }
    
    lambda_count <- lambda_count + 1
    
  }
  
  return(list("lambda" = lambda, "df_curr" = df_curr, "lambda_count" = lambda_count, "lambda_move" = lambda_move))
  
}

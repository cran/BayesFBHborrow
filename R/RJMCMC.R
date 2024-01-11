#' @title Birth move in RJMCMC
#' @description Calculates new values of x when proposing another split point, based on a 
#' weighted mean, as x_new/x <- (1-U)/U
#'
#' @param U uniform random number
#' @param sj upcoming split point location, j
#' @param s_star new split point location, *
#' @param sjm1 previous split point location, j-1
#' @param x vector of parameter values, length J + 1
#' @param j split point
#'
#' @return vector with adjusted parameter values after additional split point,
#' length J + 2
.birth_move <- function(U, sj, s_star, sjm1, x, j) {
  lxj <- log(x[j]) - (sj - s_star) / (sj - sjm1) * log((1 - U) / U)

  lxjp1 <- log(x[j]) + (s_star - sjm1) / (sj - sjm1) * log((1 - U) / U)

  x_prop <- append(x, exp(c(lxj, lxjp1)) , after = j)[-j]

  return(x_prop)
}

#' @title Death move in RJMCMC
#' @description Calculates new values of x when proposing the death of a split point
#'
#' @param sjp1 upcoming split point location, J + 1
#' @param sj split point location to be removed, j
#' @param sjm1 previous split point location, j-1
#' @param x vector of parameter values, length J + 1
#' @param j split point
#'
#' @return vector with adjusted parameter values after removal of split point,
#' length J
.death_move <- function(sjp1, sj, sjm1, x, j) {
  lxj <- ((sj- sjm1) * log(x[j-1]) + (sjp1 - sj) * log(x[j])) / (sjp1 - sjm1)

  x_prop <-  append(x, exp(lxj), after = j)[-((j-1):j)]

  return(x_prop)
}

#' @title Calculate log density tau prior
#'
#' @param tau current value(s) of tau
#' @param a_tau tau hyperparameter
#' @param b_tau tau hyperparameter
#' @param c_tau tau hyperparameter
#' @param d_tau tau hyperparameter
#' @param p_0 mixture ratio
#' @param type choice of borrowing, "mix", "uni", or any other string for 
#' borrowing on every baseline hazard without mixture
#'
#' @return log density of tau
.ltau_dprior <- function(tau, a_tau, b_tau, c_tau = NULL, d_tau = NULL, p_0 = NULL, type) {
  if (type == "mix") {
    ldtau <- sum(log(p_0 * invgamma::dinvgamma(tau, shape = a_tau, rate = b_tau) + 
                       (1 - p_0) * invgamma::dinvgamma(tau, shape = c_tau, rate = d_tau))) 
  } else if (type == "uni") {
    ldtau <- sum(invgamma::dinvgamma(tau, shape = a_tau, rate = b_tau, log = T))
  }
}

#' @title RJMCMC (with Bayesian Borrowing)
#' @description Metropolis-Hastings Green Reversible Jump move, with Bayesian 
#' Borrowing
#'
#' @param df_hist data_frame containing historical data.
#' @param df_curr data_frame containing current trial data.
#' @param Y data.
#' @param Y_0 historical data.
#' @param I censoring indicator.
#' @param I_0 historical trial censoring indicator.
#' @param X design matrix.
#' @param X_0 historical trial design matrix.
#' @param lambda baseline hazard.
#' @param lambda_0 historical trial baseline hazard.
#' @param beta current trial parameters.
#' @param beta_0 historical trial parameters.
#' @param mu prior mean for baseline hazard.
#' @param sigma2 prior variance hyperparameter for baseline hazard.
#' @param tau borrowing parameter.
#' @param s split point locations, J + 2.
#' @param J number of split points.
#' @param Jmax maximum number of split points.
#' @param bp number of covariates in current trial.
#' @param bp_0 number of covariates in historical trial.
#' @param clam_smooth neighbor interactions, in range (0, 1), for ICAR update.
#' @param a_tau tau hyperparameter.
#' @param b_tau tau hyperparameter.
#' @param c_tau tau hyperparameter.
#' @param d_tau tau hyperparameter.
#' @param type choice of borrowing, "mix", "uni", or any other string for 
#' borrowing on every baseline hazard without mixture.
#' @param p_0 mixture ratio.
#' @param phi J hyperparameter.
#' @param pi_b probability of birth move.
#' @param maxSj maximal time point, either current or historic.
#'
#' @return list of proposed J and s, with adjusted values of lambda, lambda_0, 
#' tau, Sigma_s, and data_frames for historical and current trial data.
#'
.J_RJMCMC <- function(df_hist, df_curr, Y, Y_0, I, I_0, 
                     X, X_0, lambda, lambda_0, 
                     beta, beta_0, 
                     mu, sigma2,
                     tau,  
                     s, J, Jmax, bp, bp_0, 
                     clam_smooth,
                     a_tau = NULL, b_tau = NULL, c_tau = NULL, d_tau = NULL, type,
                     p_0 = NULL, phi, pi_b,
                     maxSj) {
  sindx <- 1:(J + 2)
  Sigma_s <- .ICAR_calc(s, J, clam_smooth)$Sigma_s

  #Birth or death move
  if (J==0) {
    move <- 0 
    pi_b <- 1 
    pi_d <- 1
  } else if (J == Jmax) {
    move <- 2
    pi_b <- 1
    pi_d <- 1
  } else {
    move <- stats::runif(1)
    pi_d <- 1 - pi_b
  }

  if (move < pi_b) {
    # Birth move, update split point locations
    s_star <-stats::runif(1, s[1], maxSj)
    s_max <- max(s)
    jlow <- max(sindx[s < s_star])
    jup <- min(sindx[s > s_star])
    slow <- s[jlow]
    sup <- s[jup]
    s_prop <- append(s, s_star , after = jlow)
    J_prop <- J + 1

    U2 <- stats::runif(1)
    U3 <- stats::runif(1)
    U4 <- stats::runif(1)

    # lambda proposal
    lambda0_prop <- .birth_move(U = U3, sj = sup , s_star = s_star, sjm1 = slow, x = lambda_0, j = jlow)
    lambda_prop <-  .birth_move(U = U2, sj = sup , s_star = s_star, sjm1 = slow, x = lambda, j = jlow)

    # Update data.frames and calculate llikelihood ratio
    df_hist_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    df_curr_prop <- .dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda_prop, bp = bp, J = J_prop)  

    llike_num <- .log_likelihood(df_hist_prop, beta_0) + .log_likelihood(df_curr_prop, beta)
    llike_den <- .log_likelihood(df_hist, beta_0) + .log_likelihood(df_curr, beta)

    # Calculate lpriors
    Sigma_s_prop <- .ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s

    lprior_num <- mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J + 2), Sigma_s_prop * sigma2, log = T) + 
      log(s_star - slow) + log(sup - s_star) + log(2 * J + 3) + log(2 * J + 2) + 
      stats::dpois(J_prop, phi, log = T)

    lprior_den <- mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) + 
      stats::dpois(J, phi, log = T) + log(sup - slow) + 2 * log(s_max)

    # Adjust for scalar tau if J == 0 and non-piecewise tau
    if (type %in% c("uni", "mix")) {
      tau_prop <- .birth_move(U = U4, sj = sup , s_star = s_star, sjm1 = slow, x = tau, j = jlow)
      tau_star <- tau_prop[c(jlow, jup)]
      tau_curr <- tau[jlow]

      lprior_num <- lprior_num + mvtnorm::dmvnorm(log(lambda_prop), log(lambda0_prop), diag(tau_prop), log = T) +
        .ltau_dprior(tau_star, a_tau, b_tau, c_tau, d_tau, p_0, type) 
      lprior_den <- lprior_den + .ltau_dprior(tau_curr, a_tau, b_tau, c_tau, d_tau, p_0, type) 

      if (J == 0) {
        lprior_den <- lprior_den + stats::dnorm(log(lambda), log(lambda_0), sqrt(tau), log = T)
      } else {
        lprior_den <- lprior_den + mvtnorm::dmvnorm(log(lambda), log(lambda_0), diag(tau), log = T)
      }
    } else {
      lprior_num <- lprior_num + mvtnorm::dmvnorm(log(lambda_prop), log(lambda0_prop), diag(tau, J + 2, J + 2), log = T)
      if (J == 0) {
        lprior_den <- lprior_den + stats::dnorm(log(lambda), log(lambda_0), sqrt(tau), log = T)
      } else {
        lprior_den <- lprior_den + mvtnorm::dmvnorm(log(lambda), log(lambda_0), diag(tau, J + 1, J + 1), log = T)
      }  
    }

    # Proposal
    lprop <- log(pi_d) - log(J + 1) - log(pi_b) + log(maxSj)

    # Jacobian
    ljac <- -log(U2) - log(1 - U2) - log(U3) - log(1 - U3)
    if (type %in% c("uni", "mix")) {
      ljac <-  ljac - log(U4) - log(1 - U4)
    }

    # Prob 
    logacc <- llike_num - llike_den + lprior_num - lprior_den + lprop + ljac
    
    if (logacc > log(stats::runif(1))) {
      lambda_0 <- lambda0_prop
      lambda <- lambda_prop
      s <- s_prop
      J <- J_prop
      Sigma_s <- Sigma_s_prop 
      if (type %in% c("uni", "mix")) {
        tau <- tau_prop
      }
      df_hist <- df_hist_prop
      df_curr <- df_curr_prop
    } 

    
  } else {
    #Death move
    if (J >= 2) {
      j <- sample(sindx[-c(1, J + 2)], 1)
    } else {
      j <- 2
    }
    s_max <- max(s)
    sj <- s[j]
    slow <- s[j - 1]
    sup <- s[j + 1]

    s_prop <- s[-j]

    J_prop <- J - 1

    U2 <- stats::runif(1)
    U3 <- stats::runif(1)
    U4 <- stats::runif(1)
    
    # lambda proposal
    lambda0_prop <- .death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda_0, j = j)
    lambda_prop <- .death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda, j = j)
    
    # likelihood
    df_hist_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    df_curr_prop <- .dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda_prop, bp = bp, J = J_prop)  
    
    llike_num <- .log_likelihood(df_hist_prop, beta_0) + .log_likelihood(df_curr_prop, beta)
    llike_den <- .log_likelihood(df_hist, beta_0) + .log_likelihood(df_curr, beta)
    
    # lprior calculations
    Sigma_s_prop <- .ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s
    
    lprior_num <- stats::dpois(J_prop, phi, log = T)  +
      mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J), Sigma_s_prop * sigma2, log = T) + 
      2 * log(s_max) + log(sup - slow)
    
    lprior_den <- stats::dpois(J, phi, log = T)  +  
      mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) + 
      log(sj - slow) + log(sup - sj) + log(2 * J + 1) + log(2 * J) 
    
    # Account for scalar / piecewise tau
    if (type %in% c("uni", "mix")) {
      tau_prop <- .death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = tau, j = j)
      tau_star <- tau_prop[j-1]
      tau_curr <- tau[(j-1):j]
      
      lprior_num <- lprior_num + .ltau_dprior(tau_star, a_tau, b_tau, c_tau, d_tau, p_0, type) 
      lprior_den <- lprior_den + mvtnorm::dmvnorm(log(lambda), log(lambda_0), diag(tau), log = T) +
      .ltau_dprior(tau_curr, a_tau, b_tau, c_tau, d_tau, p_0, type) 
      if (J == 1) {
        lprior_num<- lprior_num + stats::dnorm(log(lambda_prop), log(lambda0_prop), sqrt(tau_prop), log = T)
      } else {
        lprior_num <- lprior_num + mvtnorm::dmvnorm(log(lambda_prop), log(lambda0_prop), diag(tau_prop), log = T)
      }  
      
    } else { #Non piecewise tau (no proposal for tau)
      lprior_den <- lprior_den + mvtnorm::dmvnorm(log(lambda), log(lambda_0), diag(tau, J + 1, J + 1), log = T)
      if (J == 1) {
        lprior_num <- lprior_num + stats::dnorm(log(lambda_prop), log(lambda0_prop), sqrt(tau), log = T)
      } else {
        lprior_num <- lprior_num + mvtnorm::dmvnorm(log(lambda_prop), log(lambda0_prop), diag(tau, J, J), log = T)
      }  
    }

    # Proposal
    lprop <- log(pi_b) -  log(maxSj) - log(pi_d) + log(J)

    # Jacobian
    ljac <- log(U2) + log(1 - U2) + log(U3) + log(1 - U3)
    if (type %in% c("uni", "mix")) {
      ljac <-  ljac + log(U4) + log(1 - U4)
    }

    # Acceptance ratio
    logacc <- llike_num - llike_den + lprior_num - lprior_den + lprop + ljac

    if (logacc > log(stats::runif(1))) {
      lambda_0 <- lambda0_prop
      lambda <- lambda_prop
      s <- s_prop
      J <- J_prop
      Sigma_s <- Sigma_s_prop 
      if (type %in% c("uni", "mix")) {
        tau <- tau_prop
      }
      df_hist <- df_hist_prop
      df_curr <- df_curr_prop
    } 

  }

  return(list("J" = J, "s" = s, "lambda" = lambda, "lambda_0" = lambda_0, 
              "tau" = tau, "Sigma_s" = Sigma_s, "df_hist"= df_hist, "df_curr" = df_curr))
} 

#' @title RJMCMC (without Bayesian Borrowing)
#' @description Metropolis-Hastings Green Reversible Jump move, without Bayesian 
#' Borrowing
#'
#' @param df data_frame 
#' @param Y_0 data
#' @param I_0 censoring indicator
#' @param X_0 design matrix
#' @param lambda_0 baseline hazard
#' @param beta_0 historical trial parameters
#' @param mu prior mean for baseline hazard
#' @param sigma2 prior variance hyperparameter for baseline hazard
#' @param s split point locations, J + 2
#' @param J number of split points
#' @param Jmax maximum number of split points
#' @param bp_0 number of covariates in historical trial
#' @param clam_smooth neighbor interactions, in range (0, 1), for ICAR update
#' @param phi J hyperparameter
#' @param pi_b probability of birth move
#'
#' @return list of proposed J and s, with adjusted values of lambda, lambda_0, 
#' tau, Sigma_s, and data_frames for historical and current trial data
.J_RJMCMC_NoBorrow <- function(df, Y_0, I_0, X_0, lambda_0, beta_0, mu, sigma2, 
                              s, J, Jmax, bp_0, clam_smooth, phi, pi_b) {
  sindx <- 1:(J + 2)
  Sigma_s <- .ICAR_calc(s, J, clam_smooth)$Sigma_s 
  
  #Birth or death
  if (J==0) {
    move <- 0 
    pi_b <- 1 
    pi_d <- 1
  } else if (J == Jmax) {
    move <- 2
    pi_b <- 1
    pi_d <- 1
  } else {
    move <- stats::runif(1)
    pi_d <- 1 - pi_b 
  }
  
  if (move < pi_b) {
    s_star <-stats::runif(1, s[1], s[J + 2])
    s_max <- max(s)
    jlow <- max(sindx[s < s_star])
    jup <- min(sindx[s > s_star])
    slow <- s[jlow]
    sup <- s[jup]
    s_prop <- append(s, s_star , after = jlow)
    
    J_prop <- J + 1
    U2 <- stats::runif(1)
    
    #lambda proposal
    lambda0_prop <- .birth_move(U = U2, sj = sup , s_star = s_star, sjm1 = slow, x = lambda_0, j = jlow)
    
    ##Likelihood
    df_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    
    llike_num <- .log_likelihood(df_prop, beta_0) 
    llike_den <- .log_likelihood(df, beta_0) 
    
    ##Prior
    Sigma_s_prop <- .ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s
    
    lprior_num <- mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J + 2), Sigma_s_prop * sigma2, log = T) +
      log(s_star - slow) + log(sup - s_star) + log(2 * J + 3) + log(2 * J + 2) + stats::dpois(J + 1, phi, log = T)
    
    lprior_denom <- mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) + 
      stats::dpois(J, phi, log = T) +  log(sup - slow) + 2 * log(s_max) 
    
    
    ##Proposal
    lprop <- log(pi_d)  - log(pi_b) - log(J + 1) + log(s_max)
    
    ##Jacobian
    ljac <- -log(U2) - log(1 - U2) 
    
    #Prob 
    logacc <- llike_num - llike_den + lprior_num - lprior_denom + lprop + ljac
    
    if (logacc > log(stats::runif(1))) {
      lambda_0 <- lambda0_prop
      s <- s_prop
      J <- J_prop
      Sigma_s <- Sigma_s_prop
      df <- df_prop
    } 
    
  } else {
    
    if (J >= 2) {
      j <- sample(sindx[-c(1, J + 2)], 1)
    } else {
      j <- 2
    }
    
    s_max <- max(s)
    sj <- s[j]
    slow <- s[j - 1]
    sup <- s[j + 1]
    
    s_prop <- s[-j]
    
    J_prop <- J - 1
    U2 <- stats::runif(1)
    
    #lambda proposal
    lambda0_prop <- .death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda_0, j = j)
    
    ##Likelihood
    df_prop <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    
    llike_num <- .log_likelihood(df_prop, beta_0) 
    llike_den <- .log_likelihood(df, beta_0) 
    
    ##Prior
    Sigma_s_prop <- .ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s
    
    lprior_num <- stats::dpois(J_prop, phi, log = T) + mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J), Sigma_s_prop * sigma2, log = T) + 
      2 * log(s_max) + log(sup- slow) 
    
    lprior_denom <- stats::dpois(J, phi, log = T) + mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) + 
      log(sj - slow) + log(sup - sj) + log(2 * J + 1) + log(2 * J) 
    
    
    ##Proposal
    lprop <- log(pi_b) - log(pi_d) + log(J) - log(s_max)
    
    ##Jacobian
    ljac <- log(U2) + log(1 - U2)
    
    #Prob 
    logacc <- llike_num - llike_den + lprior_num - lprior_denom + lprop + ljac
    
    if (logacc > log(stats::runif(1))) {
      lambda_0 <- lambda0_prop
      s <- s_prop
      J <- J_prop 
      Sigma_s <- Sigma_s_prop
      df <- df_prop
    } 
  }
  
  return(list("J" = J, "s" = s, "lambda_0" = lambda_0, "Sigma_s" = Sigma_s, "df_all" = df))
} 
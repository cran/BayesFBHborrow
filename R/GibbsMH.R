#' @title S3 generic, calls the correct GibbsMH sampler
#' 
#' @description An MCMC sampler for Bayesian borrowing with time-to-event data. 
#' We obtain a flexible baseline hazard function by making the split points 
#' random within a piecewise exponential model and using a Gaussian Markov 
#' random field prior to smooth the baseline hazards. Only calls the sampler and
#' does not run any input checks. Best practice is to call BayesFBHborrow(), if the 
#' user is not familiar with the model at hand.
#'
#' @param Y data
#' @param I event indicator
#' @param X design matrix
#' @param Y_0 historical data, default is NULL
#' @param I_0 historical event indicator, default is NULL
#' @param X_0 historical design matrix, default is NULL
#' @param tuning_parameters list of "cprop_beta", "cprop_beta_0", "alpha", "Jmax",
#' and "pi_b"
#' @param initial_values list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "tau", "lambda_0", "lambda", "beta_0", "beta") (optional)
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#'  "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma", "Jmax", "clam_smooth",
#'  "cprop_beta", "phi", "pi_b"). Default is list("a_tau" = 1,"b_tau" = 1,"c_tau" = 1,
#' "d_tau" = 0.001, "type" = "mix", "p_0" = 0.5, "a_sigma" = 2, "b_sigma" = 2,
#' "Jmax" = 20, "clam_smooth" = 0.8, "cprop_beta" = 0.3, "phi" = 3, "pi_b" = 0.5)
#' @param lambda_hyperparameters contains two hyperparameters (a_lambda and b_lambda) used 
#' for the update of lambda and lambda_0
#' @param iter number of iterations for MCMC sampler, excluding warmup, 
#' default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler, 
#' default is 2000
#' @param refresh number of iterations between printed screen updates, 
#' default is 500
#' @param max_grid grid size for the smoothed baseline hazard, default is 2000
#'
#' @return depending on if the user wishes to borrow; returns a list with values
#'  after each iteration for parameters: out_fixed (J, mu, sigma2, beta), lambda,
#'  lambda_0, tau, s, as well as tuning values of the total number of accepts:
#'  lambda_move, lambda_0_move and beta_move. Also included is the out_slam which
#'  contains the shrunk estimate of the baseline hazard.
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load example data and set your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' # The datasets consists of 3 (2) columns named "tte", "event" and "X" 
#' # (only for concurrent). To explicitly run the sampler, extract the samples as
#' # following
#' Y <- weibull_cc$tte
#' I <- weibull_cc$event
#' X <- matrix(weibull_cc$X_trt)
#' 
#' Y_0 <- weibull_hist$tte
#' I_0 <- weibull_hist$event
#' X_0 <- NULL
#' 
#' # Specify hyperparameters and tuning parameters
#' hyper <-  list("a_tau" = 1, 
#'                "b_tau" = 0.001,
#'                "c_tau" = 1,
#'                "d_tau" = 1, 
#'                "type" = "all",
#'                "p_0" = 0.5, 
#'                "a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "clam_smooth" = 0.5,
#'                "phi" = 3)
#' 
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5,
#'                           "alpha" = 0.4)
#'                           
#' # Set initial values to 'NULL' for default settings
#' output <- GibbsMH(Y, I, X, Y_0, I_0, X_0,
#'                   tuning_parameters, initial_values = NULL, hyper, 
#'                   iter = 5, warmup_iter = 1)
#'                   
GibbsMH <- function(Y, I, X, Y_0 = NULL, I_0 = NULL, X_0 = NULL, tuning_parameters,
                    initial_values, hyperparameters,
                    lambda_hyperparameters, iter, warmup_iter, refresh,
                    max_grid) {
  checkmate::assert_numeric(Y)
  Y_0 <- Y_0
  if (is.null(Y_0)) {
    class(Y) <- c("numeric", "NoBorrow")
  } else {
    class(Y) <- c("numeric", "WBorrow")
  }
  UseMethod("GibbsMH", Y)
}


#' @title GibbsMH sampler, with Bayesian Borrowing
#' @description An MCMC sampler for Bayesian borrowing with time-to-event data. 
#' We obtain a flexible baseline hazard function by making the split points 
#' random within a piecewise exponential model and using a Gaussian Markov 
#' random field prior to smooth the baseline hazards. Only calls the sampler and
#' does not run any input checks. Best practice is to call BayesFBHborrow(), if the 
#' user is not familiar with the model at hand.
#' 
#' @param Y data
#' @param I event indicator
#' @param X design matrix
#' @param Y_0 historical data
#' @param I_0 historical event indicator
#' @param X_0 historical design matrix
#' @param tuning_parameters list of "cprop_beta", "cprop_beta_0", "alpha", "Jmax",
#' and "pi_b"
#' @param initial_values list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "tau", "lambda_0", "lambda", "beta_0", "beta") (optional)
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#'  "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma", "Jmax", "clam_smooth",
#'  "cprop_beta", "phi", "pi_b"). Default is list("a_tau" = 1,"b_tau" = 1,"c_tau" = 1,
#' "d_tau" = 0.001, "type" = "mix", "p_0" = 0.5, "a_sigma" = 2, "b_sigma" = 2,
#' "Jmax" = 20, "clam_smooth" = 0.8, "cprop_beta" = 0.3, "phi" = 3, "pi_b" = 0.5)
#' @param lambda_hyperparameters contains two hyperparameters (a_lambda and b_lambda)
#' used for the update of lambda and lambda_0. Default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler, excluding warmup, 
#' default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler, 
#' default is 2000
#' @param refresh number of iterations between printed screen updates, 
#' default is 500
#' @param max_grid grid size for the smoothed baseline hazard, default is 2000
#'
#' @return list with values after each iteration for parameters: out_fixed (J, 
#'  mu, sigma2, beta), lambda, lambda_0, tau, s, as well as tuning values of the total number 
#'  of accepts: lambda_move, lambda_0_move and beta_move. Also included is the out_slam which 
#'  contains the shrunk estimate of the baseline hazard.
#'  
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load example data and set your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' # The datasets consists of 3 (2) columns named "tte", "event" and "X" 
#' # (only for concurrent). To explicitly run the sampler, extract the samples as
#' # following
#' Y <- weibull_cc$tte
#' I <- weibull_cc$event
#' X <- matrix(weibull_cc$X_trt)
#' 
#' Y_0 <- weibull_hist$tte
#' I_0 <- weibull_hist$event
#' X_0 <- NULL
#' 
#' # Specify hyperparameters and tuning parameters
#' hyper <-  list("a_tau" = 1, 
#'                "b_tau" = 0.001,
#'                "c_tau" = 1,
#'                "d_tau" = 1, 
#'                "type" = "all",
#'                "p_0" = 0.5, 
#'                "a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "clam_smooth" = 0.5,
#'                "phi" = 3)
#' 
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5,
#'                           "alpha" = 0.4)
#'                           
#' # Set initial values to 'NULL' for default settings
#' output <- GibbsMH(Y, I, X, Y_0, I_0, X_0,
#'                   tuning_parameters, initial_values = NULL, hyper, 
#'                   iter = 5, warmup_iter = 1)
GibbsMH.WBorrow <- function(Y, I, X, 
                            Y_0, I_0, X_0,
                            tuning_parameters,
                            initial_values = NULL,
                            hyperparameters = list(
                              "a_tau" = 1,
                              "b_tau" = 0.001,
                              "c_tau" = 1,
                              "d_tau" = 1,
                              "type" = "mix",
                              "p_0" = 0.8,
                              "a_sigma" = 1,
                              "b_sigma" = 1,
                              "phi" = 3, 
                              "clam_smooth" = 0.8),
                            lambda_hyperparameters = list(
                              "a_lambda" = 0.01,
                              "b_lambda" = 0.01
                            ),
                            iter = 150L,
                            warmup_iter = 10L,
                            refresh = 0,
                            max_grid = 2000L
) {
  ### Initialize parameters ###
  # count accept for MH beta
  lambda_0_count <- 0
  lambda_count <- 0
  lambda_move <- 0
  lambda_0_move <- 0
  
  # proposal prior 
  a_lambda <- lambda_hyperparameters$a_lambda
  b_lambda <- lambda_hyperparameters$b_lambda
  #alpha <- lambda_hyperparameters$alpha
  
  # hyperparameters
  a_sigma <- hyperparameters$a_sigma
  b_sigma <- hyperparameters$b_sigma
  clam <- hyperparameters$clam_smooth
  phi <- hyperparameters$phi
  
  if (is.null(initial_values)) {
    J <- phi
    sigma2 <- b_sigma / (a_sigma + 1)
    quantiles <- quantile(Y, probs = seq(0, 1, length.out = J + 2),
                          na.rm = TRUE, names = FALSE)
    s <- c(0, quantiles[1:J + 1], max(Y, Y_0))
    
    group_data_cc <- group_summary(Y[X[,1] == 0], I[X[,1] == 0], NULL, s)
    group_data_hist <- group_summary(Y_0, I_0, NULL, s)
    
    lambda_init <- init_lambda_hyperparameters(group_data_cc, s)
    lambda_init_hist <- init_lambda_hyperparameters(group_data_hist, s)
    
    lambda <- mapply(stats::rgamma, n = 1,
                     shape = lambda_init$shape,
                     rate = lambda_init$rate)
    
    lambda_0 <- mapply(stats::rgamma, n = 1,
                       shape = lambda_init_hist$shape,
                       rate = lambda_init_hist$rate)
    
    ### mu and sigma2
    # add the time exposed, etc.
    lambda_init_sum_hist <- init_lambda_hyperparameters(lapply(group_data_hist, sum), s[c(1, J + 2)])
    mu <- mean(log(mapply(stats::rgamma, n = 500,
                      shape = lambda_init_sum_hist$shape,
                      rate = lambda_init_sum_hist$rate)))
    
    sigma2 <- var(log(mapply(stats::rgamma, n = 500,
                          shape = lambda_init_sum_hist$shape,
                          rate = lambda_init_sum_hist$rate)))
    
    #Data and beta initial values
    bp <- ncol(X)
    df_curr <- .dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp = ncol(X), J = J)
    
    if(is.null(X_0)){
      bp_0 <- 0
      beta_0 <- NULL
    } else {
      bp_0 <- ncol(X_0)
      df_hist <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0, bp = bp_0, J = J)
      glm.mom_0 <- .glmFit(df_hist)
      beta_0 <- glm.mom_0$beta.mu
      beta_0_count <- numeric(bp_0)
    }  
    
    glm.mom <- .glmFit(df_curr)
    beta <- glm.mom$beta.mu
    bp <- ncol(X)
    beta_count <- numeric(bp)
    
    if(is.null(X_0)){
      df_hist <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0, bp = bp_0, J = J)
    }
    
    
  } else {
    J <- initial_values$J
    mu <- initial_values$mu
    sigma2 <- initial_values$sigma2
    s <- c(0, sort(initial_values$s_r), max(Y, Y_0))
    lambda <- initial_values$lambda
    lambda_0 <- initial_values$lambda_0
    beta <- initial_values$beta
    bp <- length(beta)
    beta_count <- numeric(bp)
    
    if(!is.null(X_0)) {
      beta_0 <- initial_values$beta_0
      bp_0 <- length(beta_0)
      beta_0_count <- rep(0, bp_0)
    } else {
      bp_0 <- 0
      beta_0 <- NULL
    }
    
    #Map lambda and introduce indicators. 
    df_hist <- .dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0, bp = bp_0, J = J)
    df_curr <- .dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp = bp, J = J)
    
  }
  
  a_tau <- hyperparameters$a_tau
  b_tau <- hyperparameters$b_tau
  type <- hyperparameters$type
  
  
  if (type == "mix") {
    c_tau <- hyperparameters$c_tau
    d_tau <- hyperparameters$d_tau
    p_0 <- hyperparameters$p_0
    tau <- p_0 * invgamma::rinvgamma(J + 1, shape = a_tau, rate = b_tau) + (1-p_0) *
           invgamma::rinvgamma(J + 1, shape = c_tau, rate = d_tau)
    
  } else if (type == "uni") {
    c_tau <- NULL
    d_tau <- NULL
    p_0 <- NULL
    tau <- invgamma::rinvgamma(J+1, shape = a_tau, rate = b_tau)
    
  } else {
    c_tau <- hyperparameters$c_tau
    d_tau <- hyperparameters$d_tau
    p_0 <- hyperparameters$p_0
    tau <- p_0 * invgamma::rinvgamma(1, shape = a_tau, rate = b_tau) + (1-p_0) *
      invgamma::rinvgamma(1, shape = c_tau, rate = d_tau)
  }
  
  # Tuning parameters
  Jmax <- tuning_parameters$Jmax
  cprop_beta <- tuning_parameters$cprop_beta
  pi_b <- tuning_parameters$pi_b
  alpha <- tuning_parameters$alpha
  if (!is.null(X_0)) {
    cprop_beta_0 <- tuning_parameters$cprop_beta_0
  }
  
  maxSj <- min(max(Y), max(Y_0))

  out_fixed <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = 3 + bp + bp_0))
  colnames(out_fixed)[1:3] <- c("J", "mu", "sigma2")
  colnames(out_fixed)[4:(3 + bp)] <- paste0("beta_", 1:bp)
  if(bp_0 > 0) {
    colnames(out_fixed)[(3 + bp + 1):(3 + bp + bp_0)] <- paste0("beta_0_", 1:bp_0)
  }
  
  out_lambda <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax + 1))
  out_lambda_0 <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax + 1))
  out_s <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax + 2))
  
  if(type %in% c("uni", "mix")) {
    out_tau <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax +1))
  }else{
    out_tau <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = 1))
    colnames(out_tau) <- NULL
  }
  
  #Max number of grid points
  time_grid <- seq(1e-8, max(Y), length.out =  max_grid)
  out_slam <-  data.frame(matrix(data = NA, nrow = 0, ncol =  max_grid))
  colnames(out_slam) <- time_grid
  
  
  ### MCMC START ###
  sample <- c(rep("(Warmup)", warmup_iter), rep("(Sampling)", iter))
  mess <- character(iter + warmup_iter)
  if (refresh != 0 && refresh < (iter + warmup_iter)) {
    mess[seq(refresh, iter + warmup_iter, refresh)] <- paste("Iteration:", 
                                                             seq(refresh, iter + warmup_iter, refresh), "/", (iter + warmup_iter), 
                                                             sample[seq(refresh, iter + warmup_iter, refresh)])
  } else if (refresh > (iter + warmup_iter)) {
    message("'refresh' is larger than number of iterations, using default value (0)")
    refresh <- 0
  }

  for(i in 1:(iter + warmup_iter)) {
    if(i%%refresh == 0 && refresh != 0){message(mess[i])}
    
    ICAR <- .ICAR_calc(s = s, J = J, clam = clam)
    Sigma_s <- ICAR$Sigma_s 
    
    # 1. Conjugate posterior updates [Gibbs]
    mu <- .mu_update(Sigma_s, lambda_0, sigma2, J)
    sigma2 <- .sigma2_update(mu, lambda_0, Sigma_s, J, a_sigma, b_sigma)
    tau_all <- .tau_update(lambda_0, lambda, J, s, a_tau, b_tau, c_tau, d_tau, p_0, type)
    tau <- tau_all$tau
    
    # 2. Update beta [MH NR] 
    beta_all <- .beta_MH_NR(df_curr, beta, bp, cprop_beta, beta_count)
    beta <- beta_all$beta
    beta_count <- beta_all$beta_count
    
    if(bp_0 > 0) {
      beta_0_all <- .beta_MH_NR(df = df_hist, beta = beta_0, bp = bp_0, 
                                cprop_beta = cprop_beta_0, beta_count = beta_0_count)
      beta_0 <- beta_0_all$beta
      beta_0_count <- beta_0_all$beta_count
    }
      
    # 3. Update lambda_0, propose new lambda_0 from conditional [MH]
    #   conjugate posterior
    df_lambda_0 <- .lambda_0_MH_cp(df_hist, Y_0, I_0, X_0, s,
                                   beta_0, mu, sigma2,  lambda, lambda_0, tau, 
                                   bp_0, J, clam, a_lam = a_lambda, b_lam = b_lambda,
                                   lambda_0_count, lambda_0_move)
    lambda_0 <- df_lambda_0$lambda_0
    lambda_0_move <- df_lambda_0$lambda_0_move
    lambda_0_count <- df_lambda_0$lambda_0_count
    df_hist <- df_lambda_0$df_hist
    
    # 4. Update lambda, propose new lambda from conditional 
    #   conjugate posterior [MH]
    df_lambda <- .lambda_MH_cp(df_hist, df_curr, Y, I, X, s,
                               beta, beta_0, mu, sigma2, lambda, lambda_0, tau, 
                               bp, bp_0, J, a_lam = a_lambda, b_lam = b_lambda, lambda_move,
                               lambda_count, alpha)
    lambda <- df_lambda$lambda
    lambda_move <- df_lambda$lambda_move
    lambda_count <- df_lambda$lambda_count
    df_curr <- df_lambda$df_curr
    
    # 5. Shuffle split point locations (accept/reject) [MH]
    if (J > 0) {
      swap_df <- .shuffle_split_point_location(df_hist, df_curr, Y_0, I_0, X_0, 
                                               lambda_0, beta_0, Y, I, X, lambda, beta, s, J, bp_0, bp, 
                                               clam, maxSj)
      s <- swap_df$s
      Sigma_s <- swap_df$Sigma_s
      df_hist <- swap_df$df_hist
      df_curr <- swap_df$df_curr
    }
    
    # 6. Propose a birth/death of a split point via a reversible jump step [MH-Green]
    # This will update lambda_0, lambda, tau, s and J (via weighted mean) if accepted
    rjmcmc_out <- .J_RJMCMC(df_hist, df_curr, Y, Y_0, I, I_0, X, X_0, 
                            lambda, lambda_0, beta, beta_0, 
                            mu, sigma2, tau, s, J, Jmax, 
                            bp, bp_0, 
                            clam,
                            a_tau, b_tau, c_tau, d_tau, type,
                            p_0, phi, pi_b, maxSj)
    
    J <- rjmcmc_out$J
    s <- rjmcmc_out$s
    lambda <- rjmcmc_out$lambda
    lambda_0 <- rjmcmc_out$lambda_0
    Sigma_s <- rjmcmc_out$Sigma_s
    if(type %in% c("uni", "mix")) {
      tau <- rjmcmc_out$tau
    }
    df_hist <- rjmcmc_out$df_hist
    df_curr <- rjmcmc_out$df_curr
    
    # 7. Save parameter values of iteration i
    out_fixed[i, 1] <- J
    out_fixed[i, 2] <- mu
    out_fixed[i, 3] <- sigma2
    out_fixed[i, 4:(3 + bp)] <- beta
    if(bp_0 > 0) {
      out_fixed[i, (3 + bp + 1):(3 + bp + bp_0)] <- beta_0
    }
    out_lambda[i, 1:(length(lambda))] <- lambda
    out_lambda_0[i, 1:(length(lambda_0))] <- lambda_0
    out_s[i, 1:(length(s))] <- s
    out_tau[i, 1:(length(tau))] <- tau
    
    #Grid of baseline hazards for shrunk estimate
    indx <- findInterval(time_grid, s, left.open = T)
    out_slam[i,] <- lambda[indx]
    
  }
  
  # Remove burn-in
  out_fixed <- out_fixed[-(1:warmup_iter),]
  out_lambda <- out_lambda[-(1:warmup_iter),]
  out_lambda_0 <- out_lambda_0[-(1:warmup_iter),]
  out_s <- out_s[-(1:warmup_iter),]
  out_tau <- out_tau[-(1:warmup_iter),]
  out_slam <- out_slam[-(1:warmup_iter),]
  
  if (type %in% c("uni", "mix")) {
    out_list <- list("out_fixed" = out_fixed, "lambda" = out_lambda, 
                     "lambda_0" = out_lambda_0, "s" = out_s, "tau" = out_tau,
                     "lambda_0_move" = lambda_0_move, 
                     "lambda_move" = lambda_move,
                     "beta_move" = beta_count,
                     "out_slam" = out_slam)
    class(out_list) <- c("BayesFBHborrow", "list")
  } else {
    out_fixed <- cbind(out_fixed, out_tau)
    out_fixed <- dplyr::rename_all(out_fixed, dplyr::recode, out_tau = "tau")
    out_list <- list("out_fixed" = out_fixed, "lambda" = out_lambda, 
                     "lambda_0" = out_lambda_0, "s" = out_s,
                     "lambda_0_move" = lambda_0_move, 
                     "lambda_move" = lambda_move,
                     "beta_move" = beta_count,
                     "out_slam" = out_slam)
    class(out_list) <- c("BayesFBHborrow", "list")
  }
  
  if(!is.null(beta_0)) {
    out_list$beta_0_move <- beta_0_count
  }
  
  return(out_list)
  
}


#' @title GibbsMH sampler, without Bayesian Borrowing
#' 
#' @description An MCMC sampler for time-to-event data, without Bayesian Borrowing.
#' We obtain a flexible baseline hazard function by making the split points 
#' random within a piecewise exponential model and using a Gaussian Markov 
#' random field prior to smooth the baseline hazards. Only calls the sampler and
#' does not run any input checks. Best practice is to call BayesFBHborrow(), if the 
#' user is not familiar with the model at hand.
#' 
#' @param Y data
#' @param I event indicator
#' @param X design matrix
#' @param Y_0 historical data, default is NULL
#' @param I_0 historical event indicator, default is NULL
#' @param X_0 historical design matrix, default is NULL
#' @param tuning_parameters list of "cprop_beta", "Jmax",
#' and "pi_b"
#' @param initial_values list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "lambda", beta") (optional)
#' @param hyperparameters list containing the hyperparameters c("a_sigma",
#' "b_sigma", "Jmax", "clam_smooth", "cprop_beta", "phi"). Default is 
#' list("a_sigma" = 2, "b_sigma" = 2, "Jmax" = 20, "clam_smooth" = 0.8, 
#' "cprop_beta" = 0.3, "phi" = 3)
#' @param lambda_hyperparameters contains two hyperparameters ("a" and "b") used for
#' the update of lambda, default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler, excluding warmup, 
#' default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler, 
#' default is 2000
#' @param refresh number of iterations between printed screen updates, 
#' default is 500
#' @param max_grid grid size for the smoothed baseline hazard, default is 2000
#'
#' @return list with values after each iteration for parameters: out_fixed (J, 
#'  mu, sigma2, beta), lambda, s, as well as tuning values of the total number 
#'  of accepts: lambda_move and beta_move. Also included is the out_slam which 
#'  contains the shrunk estimate of the baseline hazard.
#'  
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load example data and set your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' # The datasets consists of 3 (2) columns named "tte", "event" and "X".
#' # To explicitly run the sampler, extract the samples as following
#' Y <- weibull_cc$tte
#' I <- weibull_cc$event
#' X <- matrix(weibull_cc$X_trt)
#' 
#' # Specify hyperparameters and tuning parameters
#' hyper <-  list("a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "clam_smooth" = 0.5,
#'                "phi" = 3)
#' 
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5)
#'                           
#' # Set initial values to 'NULL' for default settings
#' output <- GibbsMH(Y, I, X, NULL, NULL, NULL,
#'                   tuning_parameters, initial_values = NULL, hyper, 
#'                   iter = 5, warmup_iter = 1)

                  
GibbsMH.NoBorrow <- function(Y, I, X = NULL, Y_0 = NULL, I_0 = NULL, X_0 = NULL,
                             tuning_parameters,
                             initial_values = NULL,
                             hyperparameters = list(
                               "a_sigma" = 1,
                               "b_sigma" = 1,
                               "phi" = 3, 
                               "clam_smooth" = 0.8),
                             lambda_hyperparameters = list(
                               "a_lambda" = 0.01,
                               "b_lambda" = 0.01
                             ),
                             iter = 1500L,
                             warmup_iter = 10L,
                             refresh = 0,
                             max_grid = 2000L
) {

  #Proposal prior 
  a_lambda <- lambda_hyperparameters$a_lambda
  b_lambda <- lambda_hyperparameters$b_lambda
  
  a_sigma <- hyperparameters$a_sigma
  b_sigma <- hyperparameters$b_sigma
  Jmax <- hyperparameters$Jmax
  clam_smooth <- hyperparameters$clam_smooth
  phi <- hyperparameters$phi
  pi_b <- hyperparameters$pi_b
  
  # Count accepts
  lambda_count <- 0
  lambda_move <- 0
  
  # hyperparameters
  a_sigma <- hyperparameters$a_sigma
  b_sigma <- hyperparameters$b_sigma
  clam <- hyperparameters$clam_smooth
  phi <- hyperparameters$phi
  
  
  if (is.null(initial_values)) {
    J <- phi
    sigma2 <- b_sigma / (a_sigma + 1)
    quantiles <- quantile(Y, probs = seq(0, 1, length.out = J + 2),
                          na.rm = TRUE, names = FALSE)
    s <- c(0, quantiles[1:J + 1], max(Y))
   
    group_data <- group_summary(Y, I, NULL, s)
    
    lambda_init <- init_lambda_hyperparameters(group_data, s)
    lambda <- mapply(stats::rgamma, n = 1,
                     shape = lambda_init$shape,
                     rate = lambda_init$rate)
    
    lambda_init_sum <- init_lambda_hyperparameters(lapply(group_data, sum), s[c(1, J + 2)])
    
    log_hazard_sample <- log(mapply(stats::rgamma, n = 500,
                                    shape = lambda_init_sum$shape,
                                    rate = lambda_init_sum$rate))
    mu <- mean(log_hazard_sample)
    sigma2 <- var(log_hazard_sample)
    
    if(is.null(X)) {
      beta <- NULL
      bp <- 0
    } else {
      bp <- ncol(X)
      df_all <- .dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp =  ncol(X), J = J)
      glm.mom <- .glmFit(df_all)
      beta <- glm.mom$beta.mu
      beta_count <- numeric(bp)
    }
    
    if (bp == 0) {
      df_all <- .dataframe_fun(Y = Y, I = I, X = NULL, s = s, lambda = lambda, bp =  0, J = J)
    }
    
    
    
  } else {
    
    J <- initial_values$J
    mu <- initial_values$mu
    sigma2 <- initial_values$sigma2
    s <- c(0, sort(initial_values$s_r), max(Y, Y_0))
    lambda <- initial_values$lambda
    
    if (!is.null(X)) {
      beta <- initial_values$beta
      bp <- length(beta)
      beta_count <- numeric(bp)
    } else {
      beta <- NULL
      bp <- 0
    }
    
    if (bp > 0) {
      df_all <- .dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp = bp, J = J)
    }else{
      df_all <- .dataframe_fun(Y = Y, I = I, X = NULL, s = s, lambda = lambda, bp =  0, J = J)
    }
  
  }
  
  # Tuning parameters
  Jmax <- tuning_parameters$Jmax
  cprop_beta <- tuning_parameters$cprop_beta
  pi_b <- tuning_parameters$pi_b
  
  #Output array 
  # J, mu, sigma2, beta
  out_fixed <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = 3 +  bp))
  colnames(out_fixed)[1:3] <- c("J", "mu", "sigma2")
  if (bp > 0) {
    colnames(out_fixed)[4:(3 +  bp)] <- paste0("beta_", 1: bp)
  }
  
  out_lambda <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax + 1))
  out_s <- data.frame(matrix(NA, nrow = iter + warmup_iter, ncol = Jmax + 2))
  
  #Max number of grid points
  t <- seq(1e-8, max(Y), length.out = max_grid)
  out_slam <-  data.frame(matrix(data = NA, nrow = 0, ncol = max_grid))
  colnames(out_slam) <- t
  
  
  ### MCMC START ###
  sample <- c(rep("(Warmup)", warmup_iter), rep("(Sampling)", iter))
  mess <- character(iter + warmup_iter)
  if (refresh != 0) {
    mess[seq(refresh, iter + warmup_iter, refresh)] <- paste("Iteration:", 
                                                             seq(refresh, iter + warmup_iter, refresh), "/", (iter + warmup_iter), 
                                                             sample[seq(refresh, iter + warmup_iter, refresh)])
  } else if (refresh > (iter + warmup_iter)) {
    message("Refresh is larger than number of iterations, using default value (0)")
    refresh <- 0
  }
  
  for (i in 1:(iter + warmup_iter)) {
    if(i%%refresh == 0 && refresh != 0){message(mess[i])}
    
    ICAR <- .ICAR_calc(s, J, clam_smooth)
    Sigma_s <- ICAR$Sigma_s 
    
    mu <- .mu_update(Sigma_s, lambda, sigma2, J)
    sigma2 <- .sigma2_update(mu, lambda, Sigma_s, J, a_sigma, b_sigma)
    
    #Map lambda and introduce indicators. 
    if (bp > 0) {
      beta_all <- .beta_MH_NR(df = df_all, beta = beta, bp =  bp, cprop_beta = cprop_beta, beta_count = beta_count)
      beta <- beta_all$beta
      beta_count <- beta_all$beta_count
    }
    
    #lambda, lambda and adjusted data frames
    dflambda <-.lambda_0_MH_cp_NoBorrow(df_all, Y, I, X, s, beta, mu, 
                                        sigma2, lambda,  bp, J, clam_smooth, 
                                        a_lam = a_lambda, b_lam = b_lambda, lambda_count,
                                        lambda_move)
    
    lambda <- dflambda$lambda_0
    df_all <- dflambda$df
    lambda_count <- dflambda$lambda_0_count 
    lambda_move <- dflambda$lambda_0_move
    
    #shuffle s
    if (J > 0) {
      swap_df <- .shuffle_split_point_location_NoBorrow(df_all, Y, I,  X, 
                                                        lambda, beta, s, J,  bp, clam_smooth)
      s <- swap_df$s
      Sigma_s <- swap_df$Sigma_s
      df_all <- swap_df$df_all
    }
    
    #J update
    rjmcmc_out <- .J_RJMCMC_NoBorrow(df_all, Y, I, X, lambda, beta, 
                                     mu, sigma2, s, J, Jmax,  bp, clam_smooth,
                                     phi, pi_b)
    J <- rjmcmc_out$J
    s <- rjmcmc_out$s
    lambda <- rjmcmc_out$lambda
    Sigma_s <- rjmcmc_out$Sigma_s
    df_all <- rjmcmc_out$df_all
    
    #J, mu, sigma2, beta, beta
    out_fixed[i, 1] <- J
    out_fixed[i, 2] <- mu
    out_fixed[i, 3] <- sigma2
    if (bp > 0) {
      out_fixed[i, 4:(3 + bp)] <- beta
    } 
    
    out_lambda[i, 1:(length(lambda))] <- lambda
    out_s[i, 1:(length(s))] <- s
    
    #Grid of baseline hazards for shrunk estimate
    indx <- findInterval(t, s, left.open = T)
    out_slam[i, ] <- lambda[indx]
  }
  
  # Remove burn-in
  out_fixed <- out_fixed[-(1:warmup_iter), ]
  out_lambda <- out_lambda[-(1:warmup_iter), ]
  out_s <- out_s[-(1:warmup_iter), ]
  out_slam <- out_slam[-(1:warmup_iter), ]
  
  out_list <- list("out_fixed" = out_fixed, "lambda" = out_lambda, "s" = out_s,
                   "out_slam" = out_slam,  "lambda_count" = lambda_count,
                   "lambda_move" = lambda_move)
  class(out_list) <- c("BayesFBHborrow", "list")
  
  if (!is.null(beta)) {
    out_list$beta_move <- beta_count
  }
  
  return(out_list)
}
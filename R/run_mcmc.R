#' @title BayesFBHborrow: Run MCMC for a piecewise exponential model
#'
#' @description Main function of the BayesFBHborrow package. This generic function 
#' calls the correct MCMC sampler for time-to-event Bayesian borrowing. 
#'
#' @param data data.frame containing atleast three vectors of "tte" (time-to-event)
#' and "event" (censoring), and covariates "X_i" (where i should be a number/
#' indicator of the covariate)
#' @param data_hist data.frame containing atleast three vectors of "tte" 
#'  (time-to-event) and "event" (censoring), with the option of adding covariates
#'  named "X_0_i" (where i should be a number/
#' indicator of the covariate), for historical
#'  data
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#' "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma", "Jmax", "clam_smooth",
#'  "cprop_beta", "phi", "pi_b")
#' @param initial_parameters list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "tau", "lambda_0", "lambda", "beta_0", "beta")
#' @param lambda_hyperparameters contains two (three) hyperparameters (a, b (,alpha))
#'  used for the update of lambda and lambda_0. alpha is the power parameter when 
#'  sampling for lambda (effects how much is borrowed)
#' @param iter number of iterations for MCMC sampler
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler.
#' @param refresh number of iterations between printed screen updates
#' @param verbose TRUE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grids size for the smoothed baseline hazard
#'
#' @export
#' 
#' @return list of samples for both fixed (can be found in $out_fixed) and 
#' multidimensional parameters (lambda, lambda_0, s, tau)
#'
#' @examples
#' set.seed(123)
#' # Load the example data and write your initial values and hyper parameters
#' data(piecewise_exp_cc, package = "BayesFBHborrow")
#' data(piecewise_exp_hist, package = "BayesFBHborrow")
#' 
#' # Set your initial values and hyper parameters
#' s <- c(0, quantile(piecewise_exp_cc$tte, c(0.25, 0.75, 1.0), names = FALSE))
#' group_data_cc <- group_summary(piecewise_exp_cc[piecewise_exp_cc$X_trt == 0,]$tte, 
#'                                piecewise_exp_cc[piecewise_exp_cc$X_trt == 0,]$event, 
#'                                NULL, s)
#' group_data_hist <- group_summary(piecewise_exp_hist$tte, piecewise_exp_hist$event, NULL, s)
#' lambda_init_cc <- init_lambda_hyperparameters(group_data_cc, s) 
#' lambda_init_hist <- init_lambda_hyperparameters(group_data_hist, s) 
#' initial_param <- list("J" = 2, 
#'                       "s_r" = s[2:3], # split points only (length J) 
#'                       "mu" = 0, 
#'                       "sigma2" = 2,
#'                       "tau" = c( 1, 1, 1),
#'                       "lambda_0" = mapply(stats::rgamma, n=1, 
#'                                           shape = lambda_init_hist$shape, 
#'                                           rate = lambda_init_hist$rate), 
#'                       "lambda" = mapply(stats::rgamma, n=1, 
#'                                         shape = lambda_init_cc$shape, 
#'                                         rate = lambda_init_cc$rate), 
#'                        "beta_0" = NULL,
#'                        "beta" = 0)
#'
#' hyper <-  list("a_tau" = 1, 
#'                "b_tau" = 0.001,
#'                "c_tau" = 1,
#'                "d_tau" = 10, 
#'                "type" = "mix",
#'                "p_0" = 0.5, 
#'                "a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "Jmax" = 5, 
#'                "clam_smooth" = 0.5, 
#'                "cprop_beta" = 0.3,
#'                "phi" = 2, 
#'                "pi_b" = 0.5)
#'
#' out <- BayesFBHborrow(piecewise_exp_cc, piecewise_exp_hist, initial_param, hyper,
#'                 iter = 5, warmup_iter = 1)
#' 
#' # Create a summary of the output
#' # summary(out, estimator = "out_fixed")
#' 
#' # Plot some of the estimates
#' # Do beta (trace), s (hist) and lambda (matrix)
#' trace <- plot(out, 1:5, estimator = "beta_1", type = "trace")
#' hist <- plot(out, estimator = "J", type = "hist")
#' smoothed_baseline_hazard <- plot(out, 1:2000, estimator = "out_slam", 
#'                                  type = "matrix")
BayesFBHborrow <- function(data, data_hist = NULL, initial_parameters, hyperparameters,
                     lambda_hyperparameters, iter, warmup_iter, refresh,
                     verbose, max_grid) {
  checkmate::assert_data_frame(data)
  data_hist = data_hist
  if (is.null(data_hist)) {
    class(data) <- c("data.frame", "NoBorrow")
  } else {
    class(data) <- c("data.frame", "WBorrow")
  }
  UseMethod("BayesFBHborrow", as.data.frame(data))
}


#' @title Run the MCMC sampler with Bayesian Borrowing
#' 
#' @description Main function of the BayesFBHborrow package. This generic function 
#' calls the correct MCMC sampler for time-to-event Bayesian borrowing. 
#' 
#' @param data data.frame containing atleast three vectors called "tte" 
#' (time-to-event), "event" (censoring), and covariates "X_i" (where i should be a number/
#' indicator of the covariate)
#' @param data_hist data.frame containing atleast two vectors called "tte" 
#'  (time-to-event) and "event" (censoring), with the option of adding covariates
#'  named "X_0_i" (where i should be a number/
#' indicator of the covariate), for historical data
#' @param initial_parameters list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "tau", "lambda_0", "lambda", "beta_0", "beta")
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#'  "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma", "Jmax", "clam_smooth",
#'  "cprop_beta", "phi", "pi_b"). Default is list("a_tau" = 1,"b_tau" = 1,"c_tau" = 1,
#' "d_tau" = 0.001, "type" = "mix", "p_0" = 0.5, "a_sigma" = 2, "b_sigma" = 2,
#' "Jmax" = 20, "clam_smooth" = 0.8, "cprop_beta" = 0.3, "phi" = 3, "pi_b" = 0.5)
#' @param lambda_hyperparameters contains three hyperparameters (a, b and alpha)
#' used for the update of lambda and lambda_0. Default is c(0.01, 0.01, 0.4)
#' @param iter number of iterations for MCMC sampler. Default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler. 
#' Default is 2000
#' @param refresh number of iterations between printed console updates. Default
#' is 0 
#' @param verbose TRUE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grid size for the smoothed baseline hazard. Default is 2000
#'
#' @return list of samples for both fixed (can be found in $out_fixed) and 
#' multidimensional parameters (lambda, lambda_0, s, tau)
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load the example data and write your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' # Set your initial values and hyperparameters
#' s <- c(0, quantile(weibull_cc$tte, c(0.25, 0.75, 1.0), names = FALSE))
#' group_data_cc <- group_summary(weibull_cc[weibull_cc$X_trt == 0,]$tte, 
#'                                weibull_cc[weibull_cc$X_trt == 0,]$event, 
#'                                NULL, s)
#' group_data_hist <- group_summary(weibull_hist$tte, weibull_hist$event, NULL, s)
#' lambda_init_cc <- init_lambda_hyperparameters(group_data_cc, s) 
#' lambda_init_hist <- init_lambda_hyperparameters(group_data_hist, s) 
#' initial_param <- list("J" = 2, 
#'                       "s_r" = s[2:3], # split points only (length J) 
#'                       "mu" = 0, 
#'                       "sigma2" = 2,
#'                       "tau" = c( 1, 1, 1),
#'                       "lambda_0" = mapply(stats::rgamma, n=1, 
#'                                           shape = lambda_init_hist$shape, 
#'                                           rate = lambda_init_hist$rate), 
#'                       "lambda" = mapply(stats::rgamma, n=1, 
#'                                         shape = lambda_init_cc$shape, 
#'                                         rate = lambda_init_cc$rate), 
#'                        "beta_0" = NULL,
#'                        "beta" = 0)
#'
#' hyper <-  list("a_tau" = 1, 
#'                "b_tau" = 0.001,
#'                "c_tau" = 1,
#'                "d_tau" = 10, 
#'                "type" = "mix",
#'                "p_0" = 0.5, 
#'                "a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "Jmax" = 5, 
#'                "clam_smooth" = 0.5, 
#'                "cprop_beta" = 0.3,
#'                "phi" = 2, 
#'                "pi_b" = 0.5)
#'
#' out <- BayesFBHborrow(weibull_cc, weibull_hist, initial_param, hyper,
#'                 iter = 5, warmup_iter = 1)
#' 
#' # Create a summary of the output
#' summary(out, estimator = "out_fixed")
#' 
#' # Plot some of the estimates
#' trace <- plot(out, 1:5, estimator = "beta_1", type = "trace")
#' hist <- plot(out, estimator = "J", type = "hist")
#' smoothed_baseline_hazard <- plot(out, 1:2000, estimator = "out_slam", 
#'                                  type = "matrix")
BayesFBHborrow.WBorrow <- function(data, data_hist, 
                 initial_parameters,
                 hyperparameters = list(
                   "a_tau" = 1,
                   "b_tau" = 0.001,
                   "c_tau" = 1,
                   "d_tau" = 1,
                   "type" = "mix",
                   "p_0" = 0.8,
                   "a_sigma" = 1,
                   "b_sigma" = 1,
                   "Jmax" = 5,
                   "clam_smooth" = 0.8,
                   "cprop_beta" = 0.3,
                   "phi" = 3, 
                   "pi_b" = 0.5), 
                 lambda_hyperparameters = list(
                   "a_lambda" = 0.01,
                   "b_lambda" = 0.01,
                   "alpha" = 0.4),
                 iter = 150, 
                 warmup_iter = 10,
                 refresh = 0,
                 verbose = FALSE,
                 max_grid = 2000) {
  checkmate::assert_flag(verbose) # user choice on output
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_names(
    names(data),
    must.include = c("tte", "event")
  )
  checkmate::assert_data_frame(data[grepl("X", names(data))], any.missing = FALSE, min.cols = 1)
  
  checkmate::assert_data_frame(data_hist, any.missing = FALSE)
  if (length(data_hist[grepl("X", names(data_hist))]) != 0) {
    checkmate::assert_data_frame(data_hist[grepl("X", names(data_hist))], any.missing = FALSE, min.cols = 1)
    X_0 <- as.matrix(data_hist[grepl("X", names(data_hist))]) 
  } else {
    X_0 <- NULL
    }
  checkmate::assert_names(
    names(data_hist),
    must.include = c("tte", "event")
  )
  
  ## Hyperparameters
  checkmate::assert_names(
    names(hyperparameters),
    must.include = c("a_tau", "b_tau", "c_tau", "d_tau","type", "p_0", "a_sigma",
                     "b_sigma", "Jmax", "clam_smooth", "cprop_beta", "phi", "pi_b")
  )
  
  ## Initial values
  checkmate::assert_names(
    names(initial_parameters),
    must.include = c("J", "s_r", "mu", "sigma2", "tau", "lambda_0", 
                     "lambda", "beta_0", "beta")
  )
  
  # Call the Gibbs_MH_WBorrow function with the input data
  X <- as.matrix(data[grepl("X", names(data))])
  
  # Call the .input_check function
  if (verbose) {
    s <- .input_check(data$tte, data_hist$tte, X, X_0, hyperparameters, initial_parameters)
    message((paste0(s, "\nInputs look ok\n")))
  } else {
    suppressMessages(
      .input_check(data$tte, data_hist$tte, X, X_0, hyperparameters, initial_parameters)
    )
  }
  if (verbose) {
    # print diagnostics
    message("Starting MCMC sampler")
  }
  
  out <- GibbsMH(Y = data$tte, I = data$event, X = X, 
                            Y_0 = data_hist$tte, I_0 = data_hist$event, X_0 = X_0,
                            initial_parameters = initial_parameters,
                            hyperparameters = hyperparameters, 
                            lambda_hyperparameters = lambda_hyperparameters,
                            iter = iter, 
                            warmup_iter = warmup_iter,
                            refresh = refresh,
                            max_grid = max_grid)
  if (verbose) {
    # print diagnostics
    beta_acc_ratio <- out$beta_move/(warmup_iter + iter)
    message("MCMC sampler complete")
    message(paste("Beta acceptance ratio: ", beta_acc_ratio))
  }

  class(out) <- c("BayesFBHborrow", "list")
  return(out)
}

#' @title Run the MCMC sampler without Bayesian Borrowing
#' 
#' @description Main function of the BayesFBHborrow package. This generic function 
#' calls the correct MCMC sampler for time-to-event without Bayesian borrowing. 
#'
#' @param data data.frame containing atleast three vectors of "tte" (time-to-event)
#' and "event" (event indicator), and covariates "X_i" (where i should be a number/
#' indicator of the covariate)
#' @param data_hist NULL (not used)
#' @param initial_parameters list containing the initial values of c("J", "s_r",
#' "mu", "sigma2", "lambda", beta")
#' @param hyperparameters list containing the hyperparameters c("a_sigma",
#' "b_sigma", "Jmax", "clam_smooth", "cprop_beta", "phi", "pi_b"). Default is 
#' list("a_sigma" = 2, "b_sigma" = 2, "Jmax" = 20, "clam_smooth" = 0.8, 
#' "cprop_beta" = 0.3, "phi" = 3, "pi_b" = 0.5)
#' @param lambda_hyperparameters contains two hyperparameters ("a" and "b") used for
#' the update of lambda, default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler. Default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler. 
#' Default is 2000
#' @param refresh number of iterations between printed console updates. Default
#' is 0 
#' @param verbose TRUE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grid size for the smoothed baseline hazard. Default is 2000
#'
#' @return list of samples for both fixed (can be found in $out_fixed) and 
#' multidimensional parameters (lambda, s, tau)
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load the example data and write your initial values and hyper parameters
#' # For this example we want to calculate the Effective Historical Sample Size
#' # (EHSS), and therefor we need to run the sampler for both with and without
#' # borrowing.
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' # Set your initial values and hyperparameters
#' s <- c(0, quantile(weibull_cc$tte, c(0.25, 0.75, 1.0), names = FALSE))
#' group_data_cc <- group_summary(weibull_cc[weibull_cc$X_trt == 0,]$tte, 
#'                                weibull_cc[weibull_cc$X_trt == 0,]$event, 
#'                                NULL, s)
#' group_data_hist <- group_summary(weibull_hist$tte, weibull_hist$event, NULL, s)
#' lambda_init_cc <- init_lambda_hyperparameters(group_data_cc, s) 
#' lambda_init_hist <- init_lambda_hyperparameters(group_data_hist, s) 
#' initial_param <- list("J" = 2, 
#'                       "s_r" = s[2:3], # split points only (length J) 
#'                       "mu" = 0, 
#'                       "sigma2" = 2,
#'                       "tau" = c( 1, 1, 1),
#'                       "lambda_0" = mapply(stats::rgamma, n=1, 
#'                                           shape = lambda_init_hist$shape, 
#'                                           rate = lambda_init_hist$rate), 
#'                       "lambda" = mapply(stats::rgamma, n=1, 
#'                                         shape = lambda_init_cc$shape, 
#'                                         rate = lambda_init_cc$rate), 
#'                        "beta_0" = NULL,
#'                        "beta" = 0)
#'
#' hyper <-  list("a_tau" = 1, 
#'                "b_tau" = 0.001,
#'                "c_tau" = 1,
#'                "d_tau" = 10, 
#'                "type" = "mix",
#'                "p_0" = 0.5, 
#'                "a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "Jmax" = 5, 
#'                "clam_smooth" = 0.5, 
#'                "cprop_beta" = 0.3,
#'                "phi" = 2, 
#'                "pi_b" = 0.5)
#'        
#' initial_param_NB <- list("J" = 2, 
#'                       "s_r" = s[2:3], # split points only (length J) 
#'                       "mu" = 0, 
#'                       "sigma2" = 2,
#'                       "lambda" = mapply(stats::rgamma, n=1, 
#'                                         shape = lambda_init_cc$shape, 
#'                                         rate = lambda_init_cc$rate), 
#'                        "beta" = 0)
#'
#' hyper_NB <-  list("a_sigma" = 2,
#'                "b_sigma" = 2,
#'                "Jmax" = 5, 
#'                "clam_smooth" = 0.5, 
#'                "cprop_beta" = 0.3,
#'                "phi" = 2, 
#'                "pi_b" = 0.5)         
#'                     
#' output_Noborrowing <- BayesFBHborrow(weibull_cc, NULL, initial_param_NB, hyper_NB,
#'                 iter = 10, warmup_iter = 1)
BayesFBHborrow.NoBorrow <- function(data, data_hist = NULL,
                             initial_parameters,
                             hyperparameters = list(
                               "a_sigma" = 1,
                               "b_sigma" = 1,
                               "Jmax" = 5,
                               "clam_smooth" = 0.8,
                               "cprop_beta" = 0.3,
                               "phi" = 3, 
                               "pi_b" = 0.5), 
                             lambda_hyperparameters = list(
                               "a_lambda" = 0.01,
                               "b_lambda" = 0.01),
                             iter = 150, 
                             warmup_iter = 10,
                             refresh = 0,
                             verbose = FALSE,
                             max_grid = 2000) {
  checkmate::assert_flag(verbose) # user choice on output
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_names(
    names(data),
    must.include = c("tte", "event")
  )
  checkmate::assert_data_frame(data[grepl("X", names(data))], any.missing = FALSE, min.cols = 1)

  ## Hyperparameters
  checkmate::assert_names(
    names(hyperparameters),
    must.include = c("a_sigma", "b_sigma", "Jmax", "clam_smooth", "cprop_beta",
                     "phi", "pi_b")
  )
  
  ## Initial values
  checkmate::assert_names(
    names(initial_parameters),
    must.include = c("J", "s_r", "mu", "sigma2", "lambda", "beta")
  )
  
  # Call the Gibbs_MH_WBorrow function with the input data
  X <- as.matrix(data[grepl("X", names(data))])
  
  # Call the .input_check function
  if (verbose) {
    s <- .input_check(data$tte, NULL, X, NULL, hyperparameters, initial_parameters)
    message((paste0(s, "\nInputs look ok\n")))
  } else {
    suppressMessages(
      .input_check(data$tte, NULL, X, NULL, hyperparameters, initial_parameters)
    )
  }

  out <- GibbsMH(Y = data$tte, I = data$event, X = X, 
                             initial_parameters = initial_parameters,
                             hyperparameters = hyperparameters,
                             lambda_hyperparameters = lambda_hyperparameters,
                             iter = iter,
                             warmup_iter = warmup_iter,
                             refresh = refresh,
                             max_grid = max_grid
                             )
  if (verbose) {
    # print diagnostics
    beta_acc_ratio <- out$beta_move/(warmup_iter + iter)
    message("MCMC sampler complete")
    message(paste("Beta acceptance ratio: ", beta_acc_ratio))
  }

  class(out) <- c("BayesFBHborrow", "list")
  return(out)
}
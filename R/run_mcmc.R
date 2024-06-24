#' @title BayesFBHborrow: Run MCMC for a piecewise exponential model
#'
#' @description Main function of the BayesFBHborrow package. This generic function 
#' calls the correct MCMC sampler for time-to-event Bayesian borrowing. 
#'
#' @param data data.frame containing atleast three vectors of "tte" (time-to-event)
#' and "event" (censoring), and covariates "X_i" (where i should be a number/
#' indicator of the covariate)
#' @param data_hist data.frame containing atleast two vectors of "tte" 
#'  (time-to-event) and "event" (censoring), with the option of adding covariates
#'  named "X_0_i" (where i should be a number/
#' indicator of the covariate), for historical
#'  data
#' @param borrow TRUE (default), will run the model with borrowing
#' @param model_choice choice of which borrowing model to use out of 'mix', 'uni'
#' or 'all'
#' @param tuning_parameters list of "cprop_beta" ("cprop_beta_0" for historical data), "alpha", "Jmax",
#' and "pi_b". Default is ("Jmax" = 5, "clam_smooth" = 0.8, "cprop_beta" = 0.5, 
#' "cprop_beta_0"  = 0.5, "pi_b" = 0.5)
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#'  "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma"). Default is list("a_tau" = 1,
#'  "b_tau" = 1,"c_tau" = 1, "d_tau" = 0.001, "type" = "mix", "p_0" = 0.5, 
#'  "a_sigma" = 2, "b_sigma" = 2, "phi" = 3)
#' @param lambda_hyperparameters contains two hyperparameters (a_lambda and b_lambda)
#' used for the update of lambda and lambda_0. Default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler.
#' @param refresh number of iterations between printed screen updates
#' @param verbose FALSE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grids size for the smoothed baseline hazard
#'
#' @export
#' 
#' @return a nested list of two items, 'out' and 'plots'. The list 'out' will
#' contain all the samples of the MCMC chain, as well as acceptance ratios. The
#' latter, 'plots', contains plots (and data) of the smoothed baseline hazard, 
#' smoothed survival, a histogram of the sampled number of split points, and the
#' trace plot of the treatment effect beta_1
#'
#' @examples
#' set.seed(123)
#' # Load the example data
#' data(piecewise_exp_cc, package = "BayesFBHborrow")
#' data(piecewise_exp_hist, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 3.25,
#'                           "alpha" = 0.4)
#'                           
#' # Set hyperparameters to default, with the borrowing model "mix"
#' out <- BayesFBHborrow(data = piecewise_exp_cc, data_hist = piecewise_exp_hist,
#'                       model_choice = 'mix', tuning_parameters = tuning_parameters,
#'                       iter = 2, warmup_iter = 0)
#' 
#' # Create a summary of the output
#' summary(out$out, estimator = "out_fixed")
#' 
#' # Plot the predictive curves for the treatment group
#' plots <- plot(out$out, out$out$time_grid, x_pred = c(1))
BayesFBHborrow <- function(data, data_hist = NULL, borrow = TRUE, model_choice,
                           tuning_parameters, hyperparameters, 
                           lambda_hyperparameters, iter, warmup_iter, refresh,
                           verbose, max_grid) {
  checkmate::assert_data_frame(data)
  data_hist = data_hist
  if (!borrow || is.null(data_hist)) {
    class(data) <- c("data.frame", "NoBorrow")
  } else if (borrow) {
    class(data) <- c("data.frame", "WBorrow")
  } else (stop("unrecognized 'borrow', should be ('borrow', 'no_borrow')"))
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
#' @param borrow TRUE (default), will run the model with borrowing
#' @param model_choice choice of which borrowing model to use out of 'mix', 'uni'
#' or 'all'
#' @param tuning_parameters list of "cprop_beta" ("cprop_beta_0" for historical data), "alpha", "Jmax",
#' and "pi_b". Default is ("Jmax" = 5, "clam_smooth" = 0.8, "cprop_beta" = 0.5, 
#' "cprop_beta_0"  = 0.5, "pi_b" = 0.5)
#' @param hyperparameters list containing the hyperparameters c("a_tau", "b_tau",
#'  "c_tau", "d_tau","type", "p_0", "a_sigma", "b_sigma"). Default is list("a_tau" = 1,
#'  "b_tau" = 1,"c_tau" = 1, "d_tau" = 0.001, "type" = "mix", "p_0" = 0.5, 
#'  "a_sigma" = 2, "b_sigma" = 2, "phi" = 3)
#' @param lambda_hyperparameters contains three hyperparameters (a_lambda, b_lambda)
#' used for the update of lambda and lambda_0. Default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler. Default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler. 
#' Default is 2000
#' @param refresh number of iterations between printed console updates. Default
#' is 0 
#' @param verbose FALSE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grid size for the smoothed baseline hazard. Default is 2000
#'
#' @return a nested list of two items, 'out' and 'plots'. The list 'out' will
#' contain all the samples of the MCMC chain, as well as acceptance ratios. The
#' latter, 'plots', contains plots (and data) of the smoothed baseline hazard, 
#' smoothed survival, a histogram of the sampled number of split points, and the
#' trace plot of the treatment effect beta_1
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load the example data
#' data(piecewise_exp_cc, package = "BayesFBHborrow")
#' data(piecewise_exp_hist, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 3.25,
#'                           "alpha" = 0.4)
#'                           
#' # Set hyperparameters to default, with the borrowing model "mix"
#' out <- BayesFBHborrow(data = piecewise_exp_cc, data_hist = piecewise_exp_hist,
#'                       model_choice = 'mix', tuning_parameters = tuning_parameters,
#'                       iter = 2, warmup_iter = 0)
#' 
#' # Create a summary of the output
#' summary(out$out, estimator = "out_fixed")
#' 
#' # Plot the predictive curves for the treatment group
#' plots <- plot(out$out, out$out$time_grid, x_pred = c(1))
BayesFBHborrow.WBorrow <- function(data, data_hist, borrow = TRUE,
                                   model_choice = "mix",
                                   tuning_parameters = NULL,
                                   hyperparameters = NULL,
                                   lambda_hyperparameters = list(
                                     "a_lambda" = 0.01,
                                     "b_lambda" = 0.01
                                   ),
                                   iter = 2000, 
                                   warmup_iter = 2000,
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
  # Call the Gibbs_MH_WBorrow function with the input data
  X <- as.matrix(data[grepl("X", names(data))])
  
  
  if (verbose) {
    # print diagnostics
    hyperparameters <- .set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice)
    tuning_parameters <- .set_tuning_parameters(tuning_parameters = tuning_parameters, 
                                                borrow, X, X_0)
    message("Starting MCMC sampler")
  } else {
    suppressMessages(
      hyperparameters <- .set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice))
    suppressMessages(
      tuning_parameters <- .set_tuning_parameters(tuning_parameters = tuning_parameters, 
                                                  borrow, X, X_0))}
  
  out <- GibbsMH(Y = data$tte, I = data$event, X = X, 
                            Y_0 = data_hist$tte, I_0 = data_hist$event, X_0 = X_0,
                            tuning_parameters = tuning_parameters,
                            hyperparameters = hyperparameters, 
                            lambda_hyperparameters = lambda_hyperparameters,
                            iter = iter, 
                            warmup_iter = warmup_iter,
                            refresh = refresh,
                            max_grid = max_grid)
  
  out$beta_acc_ratio <- out$beta_move/(warmup_iter + iter)
  out$beta_0_acc_ratio <- out$beta_0_move/(warmup_iter + iter)
  if (verbose) {
    # print diagnostics
    message("MCMC sampler complete")
    message(paste0("beta_",c(1:length(out$beta_acc_ratio)) ," acceptance ratio: ", out$beta_acc_ratio, collapse = "\n"))
  }
  
  ## Output plots
  bp <- dim(X)[2]
  beta_samples <- out$out_fixed[, 4:(3 + bp)] %>%
    as.matrix()
  
  # smoothed baseline hazard (treatment and control)
  plots <- list()
  plots$smooth_hazard_trt <- .smooth_hazard(out$out_slam, out$out_fixed$beta_1)
  plots$smooth_hazard_null <- .smooth_hazard(out$out_slam, NULL)
  plots$smooth_baseline_hazard <- .plot_matrix(out$time_grid, plots$smooth_hazard_trt, 
                                            title = "Smoothed posterior baseline hazard", 
                                            xlab = "Time", ylab = "Hazard function h(t|x)",
                                            y2 = plots$smooth_hazard_null)
  
  # smoothed baseline hazard (treatment and control)
  tm1 <- c(0,out$time_grid)[-(length(out$time_grid)+1)]
  grid_width <- out$time_grid-tm1
  plots$smooth_survival_mat <- .smooth_survival(grid_width, out$out_slam, beta_samples = beta_samples)
  plots$smooth_survival_null <- .smooth_survival(grid_width, out$out_slam, beta_samples = NULL)
  plots$smooth_survival  <- .plot_matrix(out$time_grid, plots$smooth_survival_mat,
                                      title = "Smoothed posterior survival", 
                                      xlab = "Time", ylab = "Survival function S(t|x)",
                                      y2 = plots$smooth_survival_null)
  
  # plots the density of beta_1
  plots$treatment_density <- ggplot(data.frame(x = out$out_fixed$beta_1), aes(x = .data$x)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = "Treatment density", x = "Treatment effect", y = "Density") +
    theme_minimal()
  
  # histogram of split points
  plots$split_points_hist <- .plot_hist(out$out_fixed$J, ylab = "Frequency", title = "Number of split points",
                                     xlab = "Number of split points", binwidth = 0.5)
  
  # trace of beta_1
  plots$treatment_trace <- .plot_trace(1:iter, out$out_fixed$beta_1, title = "Trace plot of treatment effect",
                                    xlab = "Iterations", ylab = "Treatment effect", linewidth = 0.3)

  class(out) <- c("BayesFBHborrow", "list")
  return(list("out" = out, "plots" = plots))
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
#' @param borrow FALSE (default), will run the model with borrowing
#' @param model_choice 'no_borrow' (default), for no borrowing
#' @param tuning_parameters list of "cprop_beta", "Jmax", and "pi_b". Default is 
#' ("Jmax" = 5, "cprop_beta" = 0.5, "pi_b" = 0.5)
#' @param hyperparameters list containing the hyperparameters c("a_sigma", "b_sigma", 
#' "phi", clam_smooth"). Default is list("a_sigma" = 2, "b_sigma" = 2, "phi" = 3 ,
#' "clam_smooth" = 0.8)
#' @param lambda_hyperparameters contains two hyperparameters ("a_lambda" and "b_lambda")
#' used for the update of lambda, default is c(0.01, 0.01)
#' @param iter number of iterations for MCMC sampler. Default is 2000
#' @param warmup_iter number of warmup iterations (burn-in) for MCMC sampler. 
#' Default is 2000
#' @param refresh number of iterations between printed console updates. Default
#' is 0 
#' @param verbose FALSE (default), choice of output, if TRUE will output 
#' intermittent results into console
#' @param max_grid grid size for the smoothed baseline hazard. Default is 2000
#'
#' @return a nested list of two items, 'out' and 'plots'. The list 'out' will
#' contain all the samples of the MCMC chain, as well as acceptance ratios. The
#' latter, 'plots', contains plots (and data) of the smoothed baseline hazard, 
#' smoothed survival, a histogram of the sampled number of split points, and the
#' trace plot of the treatment effect beta_1
#' @export
#'
#' @examples
#' set.seed(123)
#' # Load the example data
#' data(piecewise_exp_cc, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "cprop_beta" = 3.25)
#'                           
#' # Set initial values to default
#' out <- BayesFBHborrow(piecewise_exp_cc, NULL, borrow = FALSE, 
#'                       tuning_parameters = tuning_parameters,
#'                       iter = 2, warmup_iter = 0)
BayesFBHborrow.NoBorrow <- function(data, data_hist = NULL, borrow = FALSE,
                                    model_choice = "no_borrow",
                                    tuning_parameters = NULL,
                                    hyperparameters = NULL,
                                    lambda_hyperparameters = list(
                                      "a_lambda" = 0.01,
                                      "b_lambda" = 0.01
                                    ),
                                    iter = 2000, 
                                    warmup_iter = 2000,
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

  # Call the Gibbs_MH_WBorrow function with the input data
  X <- as.matrix(data[grepl("X", names(data))])
  
  if (verbose) {
    # print diagnostics
    hyperparameters <- .set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice)
    tuning_parameters <- .set_tuning_parameters(tuning_parameters = tuning_parameters, 
                                                borrow, X, NULL)
    message("Starting MCMC sampler")
  } else {
    suppressMessages(
      hyperparameters <- .set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice))
    suppressMessages(
      tuning_parameters <- .set_tuning_parameters(tuning_parameters = tuning_parameters, 
                                                  borrow, X, NULL))}
  
  out <- GibbsMH(Y = data$tte, I = data$event, X = X, 
                 tuning_parameters = tuning_parameters,
                 hyperparameters = hyperparameters,
                 lambda_hyperparameters = lambda_hyperparameters,
                 iter = iter,
                 warmup_iter = warmup_iter,
                 refresh = refresh,
                 max_grid = max_grid
                 )
  out$beta_acc_ratio <- out$beta_move/(warmup_iter + iter)
  if (verbose) {
    # print diagnostics
    message("MCMC sampler complete")
    message(paste0("beta_",c(1:length(out$beta_acc_ratio)) ," acceptance ratio: ", out$beta_acc_ratio, collapse = "\n"))
  }

  bp <- dim(X)[2]
  beta_samples <- out$out_fixed[, 4:(3 + bp)] %>%
    as.matrix()
  
  # Plots
  plots <- list()
  plots$smooth_hazard_mat <- .smooth_hazard(out$out_slam, out$out_fixed$beta_1)
  plots$smooth_baseline_hazard <- .plot_matrix(out$time_grid, plots$smooth_hazard_mat, 
                                            title = "Smoothed posterior baseline hazard", 
                                            xlab = "Time", ylab = "Hazard function h(t|x)")
  
  # smoothed baseline hazard (treatment and control)
  tm1 <- c(0,out$time_grid)[-(length(out$time_grid)+1)]
  grid_width <- out$time_grid-tm1
  plots$smooth_survival_mat <- .smooth_survival(grid_width, out$out_slam, beta_samples = beta_samples)
  plots$smooth_survival  <- .plot_matrix(out$time_grid, plots$smooth_survival_mat,
                                         title = "Smoothed posterior survival", 
                                         xlab = "Time", ylab = "Survival function S(t|x)")
  # plots the density of beta_1
  plots$treatment_density <- ggplot(data.frame(x = out$out_fixed$beta_1), aes(x = .data$x)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = "Treatment density", x = "Treatment effect", y = "Density") +
    theme_minimal()
  
  # histogram of split points
  plots$split_points_hist <- .plot_hist(out$out_fixed$J, ylab = "Frequency", title = "Number of split points",
                                     xlab = "Number of split points", binwidth = 0.5)
  
  # trace of beta_1
  plots$treatment_trace <- .plot_trace(1:iter, out$out_fixed$beta_1, title = "Trace plot of treatment effect",
                                       xlab = "Iterations", ylab = "Treatment effect", linewidth = 0.3)
  
  class(out) <- c("BayesFBHborrow", "list")
  return(list("out" = out, "plots" = plots))
}
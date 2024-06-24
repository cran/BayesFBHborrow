#' @title Summarize fixed MCMC results
#' 
#' @description S3 method for with borrowing. Returns summary of mean, median and given 
#' percentiles for the one dimensional parameters.
#'
#' @param object MCMC sample object from BayesFBHborrow()
#' @param estimator The type of estimator to summarize, could be "fixed", "lambda",
#' "lambda_0" or "s". The default is NULL and will print a summary of the output list.
#' @param percentiles Given percentiles to output, default is c(0.025, 0.25, 0.75, 0.975)
#' @param ... other arguments, see summary.default
#'
#' @importFrom magrittr %>%
#' @importFrom stats var
#' @importFrom stats quantile
#' @import dplyr
#'
#' @return summary of the given estimator
#' @export
#'
#' @examples
#' data(piecewise_exp_cc, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5)
#'                           
#' # run the MCMC sampler
#' out <- BayesFBHborrow(piecewise_exp_cc, NULL, tuning_parameters = tuning_parameters,
#'                       iter = 3, warmup_iter = 1)
#' 
#' # Create a summary of the output
#' summary(out$out, estimator = "out_fixed")
summary.BayesFBHborrow <- function(object, estimator = NULL, 
                             percentiles = c(0.025, 0.25, 0.75, 0.975), ...) {
  summary <- NULL
  if (is.null(estimator)) {
    summary <- summary.default(object)
  } else if (estimator == "out_fixed") {
    ## For out_fixed
    mean_val <- apply(object$out_fixed, 2, mean)
    variance_val <- apply(object$out_fixed, 2, var)
    
    # Given (or default) percentiles
    quantiles <- apply(object$out_fixed, 2, quantile, probs = percentiles)
    
    # Output structure
    summary <- tibble::tibble(
      id = names(object$out_fixed),
      Mean = mean_val,
      sd = sqrt(variance_val)
    )
    summary <-
      summary %>%
      dplyr::bind_cols(tibble::as_tibble(t(quantiles)))
    
  } else if (estimator == "lambda") {
    if (exists("lambda_0", object)) {
      lambda <- cbind(object$lambda, object$lambda_0)
      lam_names <- paste0("lam_X", 1:ncol(object$lambda))
      lam0_names <- paste0("lam0_X", 1:ncol(object$lambda))
      names(lambda) <- c(lam_names, lam0_names)
    } else {lambda <- object$lambda}
    
    mean_val <- apply(lambda, 2, mean, na.rm = TRUE)
    variance_val <- apply(lambda, 2, var, na.rm = TRUE)
    N_samples <- colSums(!is.na(lambda))
    quantiles <- apply(lambda, 2, stats::quantile, probs = percentiles, na.rm = TRUE)
    
    summary <- tibble::tibble(
      id = names(mean_val),
      samples = N_samples,
      mean = mean_val,
      sd = sqrt(variance_val)
    )
    summary <-
      summary %>%
      bind_cols(tibble::as_tibble(t(quantiles)))
    
  } else if (estimator == "s") {
    table_J <- table(object$out_fixed$J)
    most_frequent <- as.numeric(names(table_J[which.max(table_J)]))
    cat(paste0("Most frequent number of split points: ", most_frequent, "\n"))
    
    mean_val <- apply(object$s, 2, mean, na.rm = TRUE)
    variance_val <- apply(object$s, 2, stats::var, na.rm = TRUE)

    # Given (or default) percentiles
    N_samples <- colSums(!is.na(object$s))
    quantiles <- apply(object$s, 2, stats::quantile, probs = percentiles, na.rm = TRUE)
    
    summary <- tibble::tibble(
      id = names(mean_val),
      samples = N_samples,
      mean = mean_val,
      sd = sqrt(variance_val)
    )
    summary <-
      summary %>%
      dplyr::bind_cols(tibble::as_tibble(t(quantiles)))
    
  } else {stop("Type not recognized")}
  
  return(summary)
}

#' @title Extract mean posterior values
#' 
#' @description S3 method for class "BayesFBHborrow", returns the mean posterior values
#' for the fixed parameters
#'
#' @param object  MCMC sample object from BayesFBHborrow()
#' @param ... other arguments, see coef.default()
#'
#' @return mean values of given samples
#' @export
#'
#' @examples
#' data(weibull_cc, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5)
#'                           
#' # run the MCMC sampler
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters = tuning_parameters,
#'                       iter = 3, warmup_iter = 1)
#' 
#' # Plot the posterior mean values of the fixed parameters
#' coef(out$out)
coef.BayesFBHborrow <- function(object, ...) {
  return(apply(object$out_fixed, 2, mean))
}

#' @title Plot smoothed baseline hazards
#' 
#' @description Plot mean and given quantiles of a matrix. Can also be used to 
#' plot derivatives of the baseline hazard, such as estimated cumulative hazard 
#' and survival function.
#'
#' @param x_lim time grid
#' @param y samples
#' @param percentiles percentiles to include in plot, default is c(0.025, 0.975)
#' @param title optional, add title to plot
#' @param xlab optional, add xlabel
#' @param ylab optional, add ylabel
#' @param color color of the mid line, default is blue
#' @param fill color of the percentiles, default is blue
#' @param linewidth thickness of the plotted line, default is 1
#' @param alpha opacity of the percentiles, default is 0.2
#' @param y2 (optional) second set of samples for comparison
#' @param color2 (optional) color of the mid line, default is red
#' @param fill2 (optional) color of the percentiles, default is red
#' 
#' @import ggplot2
#'
#' @return a ggplot2 object
.plot_matrix <- function(x_lim, y, percentiles = c(0.05, 0.95), title = "", 
                        xlab = "", ylab = "", color = "blue", fill = "blue",
                        linewidth = 1, alpha = 0.2, y2 = NULL, color2 = "red", fill2 = "red") {
  mean_values <- apply(y, 2, mean)
  ql_values <- apply(y, 2, stats::quantile, probs = percentiles[1])
  qu_values <- apply(y, 2, stats::quantile, probs = percentiles[2])

  # data frame
  plot_data <- data.frame(
    time_grid = x_lim,
    mean = mean_values,
    ql = ql_values,
    qu = qu_values
  )
  
  if (is.null(y2)) {
  gg <- ggplot(plot_data, aes(x = .data$time_grid)) +
    geom_line(aes(y = .data$mean), color = color, linewidth = linewidth) +
    geom_ribbon(aes(ymin = .data$ql, ymax = .data$qu), fill = fill, alpha = alpha) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()
  
  } else {
    plot_data$Group = "Treatment"
    
    mean_values <- apply(y2, 2, mean)
    ql_values <- apply(y2, 2, stats::quantile, probs = percentiles[1])
    qu_values <- apply(y2, 2, stats::quantile, probs = percentiles[2])
    
    plot_data_null <- data.frame(
      time_grid = x_lim,
      mean = mean_values,
      ql = ql_values,
      qu = qu_values,
      Group = "Control"
    )
    
    gg <- ggplot() +
      geom_line(data = plot_data, aes(x = .data$time_grid, y = .data$mean, color = .data$Group), linewidth = linewidth) +
      geom_ribbon(data = plot_data, aes(x = .data$time_grid, ymin = .data$ql, ymax = .data$qu, fill = .data$Group), alpha = alpha) +
      geom_line(data = plot_data_null, aes(x = .data$time_grid, y = .data$mean, color = .data$Group), linewidth = linewidth) +
      geom_ribbon(data = plot_data_null, aes(x = .data$time_grid, ymin = .data$ql, ymax = .data$qu, fill = .data$Group), alpha = alpha) +
      labs(x = xlab, y = ylab, title = title) +
      scale_color_manual(values = c("Treatment" = color, "Control" = color2)) +
      scale_fill_manual(values = c("Treatment" = fill, "Control" = fill2)) +
      theme_minimal()
    
    
  }
  
  
  return(gg)
}

#' @title Plot MCMC trace
#' 
#' @description Creates a trace plot of given MCMC samples.
#'
#' @param x_lim x-axis of the plot
#' @param samples samples from MCMC
#' @param title optional, add title to plot
#' @param xlab optional, add xlabel
#' @param ylab optional, add ylabel
#' @param color color of the mid line, default is black
#' @param linewidth thickness of the plotted line, default is 1
#' 
#' @import ggplot2
#' 
#' @return a ggplot2 object
.plot_trace <- function(x_lim, samples, title = "", xlab = "", ylab = "", 
                       color = "black", linewidth = 1) {
  trace_data <- data.frame(iterations = x_lim, samples = samples)
  
  gg <- ggplot(trace_data, aes(x = .data$iterations, y = .data$samples)) +
    geom_line(color = color, linewidth = linewidth) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()
  
  
  return(gg)
}

#' @title Plot histogram from MCMC samples
#'
#' @description Plots a histogram of the given discrete MCMC samples
#'
#' @param samples data.frame containing the discrete MCMC samples
#' @param title title of the plot, default is none
#' @param xlab x-label of the plot, default is "Values"
#' @param ylab y-label of the plot, default is "Frequency"
#' @param color outline color for the bars, default is "black"
#' @param fill fill color, default is "blue"
#' @param binwidth width of the histogram bins, default is 0.5
#' @param scale_x option to scale the x-axis, suitable for discrete samples, 
#' default is FALSE
#' 
#' @import ggplot2
#'
#' @return a ggplot2 object
.plot_hist <- function(samples,  title = "", xlab = "Values", ylab = "Frequency", 
                          color = "black", fill = "blue", binwidth = 0.05,
                          scale_x = FALSE) {

  if (scale_x == TRUE) {
    gg <- ggplot(data.frame(values = samples), aes(x = .data$values)) +
      geom_histogram(binwidth = binwidth, fill = fill, color = color) +
      labs(x = xlab, y = ylab, title = title) +
      scale_x_continuous(breaks = seq(min(samples), max(samples), by = 1)) +
      theme_minimal()
  } else {
    gg <- ggplot(data.frame(values = samples), aes(x = .data$values)) +
      geom_histogram(binwidth = binwidth, fill = fill, color = color) +
      labs(x = xlab, y = ylab, title = title) +
      theme_minimal()
  }
  

  return(gg)
}

#' Predictive survival from BayesFBHborrow object
#'
#' @param grid_width size of time step
#' @param out_slam samples from the smoothed baseline hazard
#' @param x_pred set of predictors to be used for calculating the predictive survival
#' @param beta_samples samples of the covariates
#'
#' @return matrix of the predictive survival
.predictive_survival <- function(grid_width, out_slam, x_pred, beta_samples){
  ##predictive survival - returns a grid for plot
  if(length(x_pred) > 1){
    xbeta <- beta_samples %*% x_pred
  }else{
    xbeta <- beta_samples * x_pred
  }
  integrand <- as.matrix(-(out_slam * exp(xbeta))) %*% diag(grid_width)
  all_sp <- exp(t(apply(integrand, 1, cumsum)))
  return(all_sp)
}

#' Predictive hazard from BayesFBHborrow object
#'
#' @param out_slam samples from the smoothed baseline hazard
#' @param x_pred set of predictors to be used for calculating the predictive hazard
#' @param beta_samples samples of the covariates
#'
#' @return matrix of the predictive hazard
.predictive_hazard <- function(out_slam, x_pred, beta_samples){
  ##predictive hazard - returns a grid for plot
  if(length(x_pred) > 1){
    xbeta <- beta_samples %*% x_pred
  }else{
    xbeta <- beta_samples * x_pred
  }
  all_ph <- out_slam * exp(xbeta)
  return(all_ph)
}

#' Predictive hazard ratio (HR) from BayesFBHborrow object
#'
#' @param x_pred set of predictors to be used for calculating the predictive HR
#' @param beta_samples samples of the covariates
#'
#' @return posterior samples for expectation and credible intervals
.predictive_hazard_ratio <- function(x_pred, beta_samples){
  if(length(x_pred) > 1){
    xbeta <- beta_samples %*% x_pred
  }else{
    xbeta <- beta_samples * x_pred
  }
  phr <- exp(xbeta)
  return(phr)
}

#' Smoothed hazard function
#'
#' @param out_slam samples from GibbsMH of the baseline hazard
#' @param beta_samples samples from GibbsMH from the treatment effect
#'
#' @return smoothed function for the baseline hazard
.smooth_hazard <- function(out_slam, beta_samples = NULL){
  ##smooth hazard for trt / if no trt just average out_slam
  if (is.null(beta_samples)) {
    all_ph <- out_slam
  } else {
    all_ph <- out_slam * exp(beta_samples)
  }
  return(all_ph)
}


#' Smoothed survival curve
#'
#' @param grid_width step size
#' @param out_slam samples from GibbsMH of the baseline hazard
#' @param beta_samples samples from GibbsMH from the treatment effect
#'
#' @return smoothed survival function
.smooth_survival <- function(grid_width, out_slam, beta_samples = NULL){
  ##Smooth survival - beta_samples = NULL return control grid
  if(is.null(beta_samples)){
    integrand <- as.matrix(-out_slam) %*% diag(grid_width)
  }else{
    integrand <- as.matrix(-(out_slam * exp(beta_samples))) %*% diag(grid_width)
  }
  all_sp <- exp(t(apply(integrand, 1, cumsum)))
  return(all_sp)
}


#' @title Plot the MCMC results
#'
#' @description S3 object which produces predictive probabilities of the survival,
#' hazard, and hazard ratio for a given set of predictors
#'
#' @param x object of class "BayesFBHborrow" to be visualized
#' @param x_lim x-axis to be used for plot, set to NULL to use default from MCMC sampling
#' @param x_pred vector of chosen predictors
#' @param ... other plotting arguments, see .plot_matrix()
#' for more information
#'
#' @importFrom graphics plot.default
#' @import ggplot2
#' 
#' @return nested list of 'plots' (posterior predictive hazard, survival, 
#' and hazard ratio) as well as their samples.
#' 
#' @export
#'
#' @examples
#' data(weibull_cc, package = "BayesFBHborrow")
#' 
#' # Set your tuning parameters
#' tuning_parameters <- list("Jmax" = 5,
#'                           "pi_b" = 0.5,
#'                           "cprop_beta" = 0.5)
#'                           
#' # run the MCMC sampler
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters = tuning_parameters,
#'                       iter = 3, warmup_iter = 1)
#' 
#' # for the treatment group
#' plots <- plot(out$out, out$out$time_grid, x_pred = c(1))
plot.BayesFBHborrow <- function(x, x_lim, x_pred = NULL, ...) {
  if(is.null(x_pred)) {
    stop("Please specify set of predictors to use with 'x_pred', calling default")
  } else if (length(x_pred) != length(x$beta_move)) {
    stop("Length of 'x_pred' should be the same as the number of covariates")
  }
  
  bp <- length(x_pred)
  tm1 <- c(0,x_lim)[-(length(x_lim)+1)]
  grid_width <- x_lim-tm1
  # Need the beta posterior samples as a matrix
  beta_samples <- x$out_fixed[, 4:(3 + bp)] %>%
    as.matrix()
  
  
  p <- list()
  p$predictive_survival <- .predictive_survival(grid_width, x$out_slam, x_pred, beta_samples)
  p$predictive_hazard <- .predictive_hazard(x$out_slam, x_pred, beta_samples)
  p$predictive_hazard_ratio <- .predictive_hazard_ratio(x_pred, beta_samples)
  
  plots <- list()
  
  plots$predictive_survival <- .plot_matrix(x_lim = x_lim, y = p$predictive_survival, 
                            title = "Smoothed posterior predictive survival", xlab = "Time", 
                            ylab = "Survival function S(t|x)", ...)
  plots$predictive_hazard <- .plot_matrix(x_lim = x_lim, y = p$predictive_hazard, 
                                         title = "Smoothed posterior predictive hazard", xlab = "Time", 
                                         ylab = "Hazard function h(t|x)", ...)
  
  # Plot the kernel density for the HR
  plots$predictive_hazard_ratio <- ggplot(data.frame(x = p$predictive_hazard_ratio), aes(x = .data$x)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = "Posterior predictive hazard ratio density", x = "Hazard ratio", y = "Density") +
    theme_minimal()
  

  return(list("plots" = plots, "predictive_survival" = p$predictive_survival,
              "predictive_hazard" = p$predictive_hazard, 
              "predictive_hazard_ratio" = p$predictive_hazard_ratio))
}

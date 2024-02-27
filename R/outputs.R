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
#' out <- BayesFBHborrow(piecewise_exp_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Create a summary of the output
#' summary(out, estimator = "out_fixed")
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
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Plot the posterior mean values of the fixed parameters
#' coef(out)
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
#' 
#' @import ggplot2
#'
#' @return a ggplot2 object
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
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Visualize the smoothed baseline hazard
#' time_grid <- seq(0, max(weibull_cc$tte), length.out = 2000)
#' gg <- plot_matrix(time_grid, out$out_slam, 
#'                   title = "Example plot of smoothed baseline hazard",
#'                   xlab = "time", ylab = "baseline hazard")
plot_matrix <- function(x_lim, y, percentiles = c(0.05, 0.95), title = "", 
                        xlab = "", ylab = "", color = "blue", fill = "blue",
                        linewidth = 1, alpha = 0.2) {
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

  # Create the plot of survival
  gg <- ggplot(plot_data,  aes_string(x = "time_grid")) +
    geom_line(aes_string(y = "mean"), color = color, linewidth = linewidth) +
    geom_ribbon(aes_string(ymin = "ql", ymax = "qu"), fill = fill, alpha = alpha) +
    labs(x = xlab, y = ylab,
         title = title) +
    theme_minimal()
  
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
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Create a tarce plot of the treatment effect, beta_1
#' time_grid <- seq(0, max(weibull_cc$tte), length.out = 2000)
#' gg <- plot_trace(1:10, out$out_fixed$beta_1, 
#'                   title = "Example trace plot",
#'                   xlab = "iterations", ylab = "beta_1 (treatment effect)")
plot_trace <- function(x_lim, samples, title = "", xlab = "", ylab = "", 
                       color = "black", linewidth = 1) {
  trace_data <- data.frame(iterations = x_lim, samples = samples)

  gg <- ggplot(trace_data, aes_string(x = "iterations", y = "samples")) +
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
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Plot the frequency of the number of split points, J with a histogram
#' time_grid <- seq(0, max(weibull_cc$tte), length.out = 2000)
#' gg <- plot_hist(out$out_fixed$J, title = "Example histogram of J",
#'                 scale_x = TRUE)
plot_hist <- function(samples,  title = "", xlab = "Values", ylab = "Frequency", 
                          color = "black", fill = "blue", binwidth = 0.05,
                          scale_x = FALSE) {

  if (scale_x == TRUE) {
    gg <- ggplot(data.frame(values = samples), aes_string(x = "values")) +
      geom_histogram(binwidth = binwidth, fill = fill, color = color) +
      labs(x = xlab, y = ylab, title = title) +
      scale_x_continuous(breaks = seq(min(samples), max(samples), by = 1)) + 
      theme_minimal() 
  } else {
    gg <- ggplot(data.frame(values = samples), aes_string(x = "values")) +
      geom_histogram(binwidth = binwidth, fill = fill, color = color) +
      labs(x = xlab, y = ylab, title = title) +
      theme_minimal() 
  }

  return(gg)
}

#' @title Plot the MCMC results
#'
#' @description S3 object which produces different plots depending on the 
#' "type" variable
#'
#' @param x object of class "BayesFBHborrow" to be visualized
#' @param x_lim x-axis to be used for plot
#' @param estimator which estimate to be visualized
#' @param type The type of plot to be produced, 
#' "trace" will produce a trace plot of the "fixed" parameters, 
#' "hist" will give a histogram for the "fixed" parameters, 
#' and "matrix" will plot the mean and quantiles of a given sample.
#' @param ... other plotting arguments, see plot_trace(), plot_hist(), plot_matrix()
#' for more information
#'
#' @importFrom graphics plot.default
#' 
#' @return ggplot2 object
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
#' out <- BayesFBHborrow(weibull_cc, NULL, tuning_parameters, 
#'                       initial_values = NULL,
#'                       iter = 10, warmup_iter = 1)
#' 
#' # Now let's create a variety of plots
#' 
#' # Staring with a histogram of beta_1 (treatment effect)
#' gg_hist <- plot(out, NULL, estimator = "beta_1", type = "hist",
#'                 title = "Example histogram of beta_1")
#'
#' # And an accompanied trace plot of the same parameter                 
#' gg_trace <- plot(out, 1:10, estimator = "beta_1", type = "trace",
#'                   title = "Example trace plot", xlab = "iterations",
#'                   ylab = "beta_1 (treatment effect)")
#'                   
#' # Lastly. visualize the smoothed baseline hazard
#' time_grid <- seq(0, max(weibull_cc$tte), length.out = 2000)
#' gg_matrix <- plot(out, time_grid, estimator = "out_slam", type = "matrix",
#'                   title = "Example plot of smoothed baseline hazard",
#'                   xlab = "time", ylab = "baseline hazard")
plot.BayesFBHborrow <- function(x, x_lim, estimator = NULL, type = NULL, ...) {
  
  if (is.null(type)) {
    message("Please specify type of plot for BayesFBHborrow class object, calling default")
    gg_object <- plot.default(x_lim, x, ...)
    
  } else if (type == "trace") {
    samples = x[[1]]
    gg_object <- plot_trace(x_lim = x_lim, samples = samples[[estimator]], ...)
    
  } else if (type == "hist") {
    samples = x[[1]]
    gg_object <- plot_hist(samples = samples[[estimator]], ...)
      
  } else if (type == "matrix") {
    gg_object <- plot_matrix(x_lim = x_lim, y = x[[as.symbol(estimator)]], ...)
    
  } else {stop("'type' not recognized")}
  
  return(gg_object)
}
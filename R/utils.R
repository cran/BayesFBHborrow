#' @title Input checker
#' @description Checks inputs before Gibbs sampler is run
#'
#' @param Y current time-to-event data
#' @param Y_0 historical time-to-event data
#' @param X design Matrix
#' @param X_0 design Matrix for historical data
#' @param tuning_parameters list of tuning parameters
#' @param initial_values list of initial values (optional)
#' @param hyperparameters list of hyperparameters
#'
#' @return a print statement
.input_check <- function(Y, Y_0, X, X_0, tuning_parameters, initial_values = NULL, hyperparameters) {

  if (any(unlist(hyperparameters) < 0)) {
    stop((paste("in hyperparameters(s): ",
                names(which(unlist(hyperparameters) < 0)),
    "hyperparameters must be non negative")))

  } else if (max(hyperparameters$p_0, hyperparameters$clam_smooth,
                 tuning_parameters$pi_b) > 1) {
    listed <- c("p_0" = hyperparameters$p_0, "clam_smooth" = hyperparameters$clam_smooth,
              "pi_b" = tuning_parameters$pi_b)
    stop(paste("in hyperparameters:",
    names(which(unlist(listed) > 1)),
    ", should be in range [0, 1]"))

  }

  
  if (!is.null(Y_0)) {
    if (!is.null(initial_values)) {
      J <- initial_values$J
    if (hyperparameters$type != c("uni")) {
      borr_choice <- ifelse(hyperparameters$type == c("mix"), "mix", "all")
      s <- (paste("Choice of borrow:", borr_choice))
      if (any(hyperparameters[c("c_tau", "d_tau")] == "NULL")) {
        stop(paste("specify hyperparameters for",
                   names(which(hyperparameters[c("c_tau", "d_tau")] == "NULL"))))
      }
      if (borr_choice == "mix" && length(initial_values$tau) != J + 1) {
        stop("wrong dimension for tau, should be J+1")
      }
      else if (borr_choice == "all" && length(initial_values$tau) == J + 1) {
        stop("borrow is 'all' but several initial values for tau were provided")
      }
    } else {
      s <- "Choice of borrow: uni"
      if (any(hyperparameters[c("c_tau", "d_tau")] != "NULL")) {
        message("borrow is 'uni', choice of c_tau and d_tau will be ignored")
      }
      if (length(initial_values$tau) != J + 1) {
        stop("wrong dimension for tau, should be J+1")
      }
    }

      maxSj <- min(max(Y), max(Y_0))
      if (max(initial_values$s_r) > maxSj) {
        stop("all s_r must be < min(max(Y),max(Y_0))")
      } else if (any(initial_values$s_r < 0)) {
        stop("all s_r must be > 0")
      }
      
      if (any(sapply(initial_values[c("lambda", "lambda_0")], length) != J+1)) {
        stop(paste("dimension error in", names(which(sapply(
          initial_values[c("lambda", "lambda_0")], length) != J+1))))
      }
      
      if (any(c(initial_values$lambda, initial_values$lambda_0) < 0)) {
        stop("baseline hazard must be > 0")
      }
      
      if (is.null(X_0)) {
        if (length(initial_values$beta) != dim(X)[[2]] || length(initial_values$beta_0) != 0) {
          arg <- which(sapply(initial_values[c("beta", "beta_0")],
                              length) != c(ncol(X), ncol(X_0)))
          arg_dim <- dim(initial_values[names(arg)])
          trial <- switch(ncol(X) != length(initial_values$beta), dim(X), NULL)
          stop(paste0("dimension mismatch in ", names(arg), " with length ", arg_dim,
                      ", given design matrix has dimension ", trial))
        }
      } else {
        if (length(initial_values$beta) != dim(X)[[2]] || length(initial_values$beta_0) != dim(X_0)[[2]]) {
          arg <- which(sapply(initial_values[c("beta", "beta_0")],
                              length) != c(ncol(X), ncol(X_0)))
          arg_dim <- length(initial_values[names(arg)][[1]])
          trial <- dplyr::if_else(ncol(X) != length(initial_values$beta),
                                  dim(X)[2], dim(X_0)[2])
          stop(paste0("dimension mismatch in ", names(arg), " with length ",
                      arg_dim,", given design matrix has dimension ", trial))
        }
      }
      cprop_beta <- tuning_parameters$cprop_beta
      max_length <- max(length(initial_values$beta), length(initial_values$beta_0))
      if (max_length != length(cprop_beta)) {
        arg <- which(sapply(initial_values[c("beta", "beta_0")],
                                   length) == max_length)
        arg_length <- length(initial_values[names(arg)][[1]])
        stop(paste0("dimension mismatch in 'cprop_beta' with length ",
                    length(cprop_beta),", given the number of covariates in ",names(arg), " is ", arg_length, "\n"))
      }
      if (initial_values$sigma2 < 0) {
        stop("sigma2 must be > 0")
      }
      
      if (length(initial_values$s_r) != J) {
        stop("dimension error in s_r")
      }
    } else { # if no initial values are provided
      if (hyperparameters$type != c("uni")) {
      borr_choice <- ifelse(hyperparameters$type == c("mix"), "mix", "all")
      s <- (paste("Choice of borrow:", borr_choice))
      if (any(hyperparameters[c("c_tau", "d_tau")] == "NULL")) {
        stop(paste("specify hyperparameters for",
                   names(which(hyperparameters[c("c_tau", "d_tau")] == "NULL"))))
      }
    } else {
      s <- "Choice of borrow: uni"
      if (any(hyperparameters[c("c_tau", "d_tau")] != "NULL")) {
        message("borrow is 'uni', choice of c_tau and d_tau will be ignored")
        }
      }
    }
  } else if (!is.null(initial_values)) {
    s <- "No borrow"
    lambda <- initial_values[grepl("^lambda", names(initial_values))]
    beta <- initial_values[grepl("^beta", names(initial_values))]
    J <- initial_values$J
    
      if (max(initial_values$s_r) > max(Y)) {
        stop("all s_r must be < max(Y)")
      } else if (any(initial_values$s_r < 0)) {
        stop("all s_r must be > 0")
      }
      
      if (length(lambda[[1]]) != J+1) {
        stop(paste0("dimension error in ", names(lambda), " with length ", 
                    length(lambda),", should be of length ", J+1))
      }
      
      if (any(lambda[[1]] < 0)) {
        stop("baseline hazard must be > 0")
      }
    
      if (!is.null(X) && length(beta[[1]]) != ncol(X)) {
        stop(paste0("dimension error in beta with length ", length(beta), 
                    ", should be of length "), ncol(X))
      }
      if (initial_values$sigma2 < 0) {
        stop("sigma2 must be > 0")
      }
    
      if (length(initial_values$s_r) != J) {
        stop("dimension error in s_r")
      }

  } else {
    s <- "No borrow"
  }

  return(invisible(s))
}

#' @title Create data.frame for piecewise exponential models
#' 
#' @description Construct a split data.frame for updated split points
#'
#' @param Y time-to-event
#' @param I censor indicator
#' @param X design Matrix
#' @param s split point locations, including start and end (length J + 2)
#' @param lambda baseline Hazards (length J+1)
#' @param bp number of covariates
#' @param J number of split points
#'
#' @return data.frame with columns c(tstart, id, X1,..., Xp, Y, I, lambda)
#'
#' @import survival
.dataframe_fun <- function(Y, I, X, s, lambda, bp, J) {
  id <- seq_along(Y)
  if (bp > 0) {
    xcol <- 1:bp
    df_like <- as.data.frame(cbind(Y, I, id, X))
    colnames(df_like)[xcol + 3] <- paste0("X", xcol)
  } else {
    df_like <- as.data.frame(cbind(Y, I, id))
  }

  #Create indicators
  df_split <- survival::survSplit(formula = Surv(Y, I) ~ ., 
                                  data = df_like, cut = s)
  
  lam_mat <- as.data.frame(cbind(s[1:(J + 1)], lambda))
  colnames(lam_mat)[1] <- "tstart"

  df_all <- merge(df_split, lam_mat, by = "tstart", all.x = T)
}

#' Computes the logarithmic sum of an exponential
#'
#' @param x set of log probabilities
#'
#' @return the logarithmic sum of an exponential
.logsumexp <- function(x) {
  c <- max(x)
  p <- c + log(sum(exp(x - c)))
}

#' Normalize a set of probability to one, using the the log-sum-exp trick
#'
#' @param x set of log probabilities
#'
#' @return normalized set of log probabilities
.normalize_prob <- function(x) {
  exp(x - .logsumexp(x))
}

#' Log likelihood function
#'
#' @param df data.frame containing data, time split points, and lambda
#' @param beta coefficients for covariates
#'
#' @return log likelihood given lambdas and betas
.log_likelihood <- function(df, beta) {
  if (!is.null(beta)) {
    X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
    xdpb <- X %*% beta

    # df has the delta_ij indicator integrated within it. Each individual
    # is partitioned by interval with Y adjusted.
    # I is the nu_i indicator.
    llike <- sum(log(df$lambda) * df$I + xdpb * df$I - 
                   ((df$Y - df$tstart) * df$lambda) * exp(xdpb))
  } else {
    llike <- sum(log(df$lambda) * df$I + df$I - ((df$Y - df$tstart) * df$lambda))
  }

}

#' @title Create group level data
#' 
#' @description Aggregate individual level data into group level data
#'
#' @param Y data
#' @param I censoring indicator
#' @param X design matrix
#' @param s split points, J + 2
#'
#' @return list of group level data
#' @export
#'
#' @examples
#' set.seed(111)
#' # Load example data and set your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' Y <- weibull_cc$tte
#' I <- weibull_cc$event
#' X <- weibull_cc$X_trt
#' 
#' # Say we want to know the group level data for the following split points
#' s <- quantile(Y, c(0, 0.45, 0.65, 1), names = FALSE)
#' 
#' group_summary(Y, I, X, s)
group_summary <- function(Y, I, X, s) {
  data <- as.data.frame(cbind(seq_along(Y), Y, I))

  if (is.null(X)) {

    # Censoring count
    df_cens <- (data[data$I == 0, ])
    df_cens$C <- 1 - df_cens$I

    events_split <- survival::survSplit(formula = Surv(Y, I) ~ .,
                                   data = data, cut = s)
    cnsr_split <- survival::survSplit(formula = Surv(Y, C) ~ .,
                                   data = df_cens, cut = s)

    # Number of events in each time segment
    events <- sapply(1:(length(s) - 1), function(i) {
      sum(events_split[events_split$tstart == s[i], ]$I)})

    # Number of censored patients in each time segment
    cnsr <- sapply(1:(length(s) - 1), function(i) {
      sum(cnsr_split[cnsr_split$tstart == s[i], ]$C)})

    # Time exposed
    time_exposed <- sapply(1:(length(s) - 1), function(i) {
                  sum(events_split[events_split$tstart == s[i], ]$Y -
                        events_split[events_split$tstart == s[i], ]$tstart)})
    
    # Number of patients no events in each time segment
    num_at_risk <- sapply(1:(length(s) - 1), function(i) {
      nrow(events_split[events_split$tstart == s[i], ])})

    return(list("events" = events, "time_exposed" = time_exposed, 
                "num_at_risk" = num_at_risk, "num_cnsr" = cnsr))

  } else {
    # Concurrent to include covariate
    data <- as.data.frame(cbind(data, X))
    
    # Censoring count which includes covariate
    df_cens <- (data[data$I == 0, ])
    df_cens$C <- 1 - df_cens$I

    events_split <- survival::survSplit(formula = Surv(Y, I) ~ .,
                                   data = data, cut = s)
    csnr_split <- survival::survSplit(formula = Surv(Y, I) ~ .,
                                   data = df_cens, cut = s)

    # Number of events control
    events_c <- sapply(1:(length(s) - 1), function(i) {
      sum(events_split[events_split$tstart == s[i] & events_split$X == 0, ]$I)})

    num_at_risk_c <- sapply(1:(length(s) - 1), function(i) {
      nrow(events_split[events_split$tstart == s[i] & events_split$X == 0, ])})

    cens <- sapply(1:(length(s) - 1), function(i) {
      sum(csnr_split[csnr_split$tstart == s[i] & csnr_split$X == 0, ]$I)})

    # Number of events treatment
    events_trt <- sapply(1:(length(s) - 1), function(i) {
      sum(events_split[events_split$tstart == s[i] & events_split$X == 1, ]$I)})

    num_at_risk_trt <- sapply(1:(length(s) - 1), function(i) {
      nrow(events_split[events_split$tstart == s[i] & events_split$X == 1, ])})

    # Time exposed
    time_exposed_c <- sapply(1:(length(s) - 1), function(i) {
      sum(events_split[events_split$tstart == s[i] & events_split$X == 0, ]$Y -
            events_split[events_split$tstart == s[i] & events_split$X == 0, ]$tstart)})

    time_exposed_trt <- sapply(1:(length(s) - 1), function(i) {
      sum(events_split[events_split$tstart == s[i] & events_split$X == 1, ]$Y -
            events_split[events_split$tstart == s[i] & events_split$X == 1, ]$tstart)})

    return(list("events_c" = events_c, "events_trt" = events_trt,
                "time_c" = time_exposed_c, "time_trt" = time_exposed_trt,
                "num_at_risk_c" = num_at_risk_c, "num_at_risk_trt" = num_at_risk_trt))
  }

}

#' @title Initialize lambda hyperparameters
#' 
#' @description Propose lambda hyperparameters for the choice of 
#' initial values for lambda
#'
#' @param group_data group level data
#' @param s split points
#' @param w weight
#'
#' @return shape and rate for the estimated lambda distribution
#' @export
#'
#' @examples
#' set.seed(111)
#' # Load example data and set your initial values and hyper parameters
#' data(weibull_cc, package = "BayesFBHborrow")
#' data(weibull_hist, package = "BayesFBHborrow")
#' 
#' Y <- weibull_cc$tte
#' I <- weibull_cc$event
#' X <- weibull_cc$X_trt
#' 
#' # Say we want to know the group level data for the following split points
#' s <- quantile(Y, c(0, 0.45, 0.65, 1), names = FALSE)
#' 
#' group_data <- group_summary(Y, I, NULL, s)
#' init_lambda_hyperparameters(group_data, s)
init_lambda_hyperparameters <- function(group_data, s, w = 0.5) {
  h_star <- rep(0, length(s)- 1)
  
  # Set index
  idx <- 1:length(h_star)
  xi <- s[-1] - s[-length(s)]
  n_dash <- group_data$num_at_risk - group_data$num_cnsr / 2
  numerator <- ((n_dash - group_data$events/2) * xi)
  
  # Adjust index for 0
  idx <- idx[numerator!= 0]
  h_star[idx] <- (group_data$events / numerator)[idx]
  
  ingroup_data <- (1:length(xi))[group_data$events == 0]
  indnz <- (1:length(xi))[group_data$events != 0]
  
  rate <- NULL
  shape <-group_data$num_at_risk * w
  
  # Account for zero number of events
  rate[indnz] <- shape[indnz] / h_star[indnz]
  rate[ingroup_data] <- shape[ingroup_data] / rep(0.1, length(ingroup_data))
  
  # t2_gamma hyper prior for mu - using nonzero values only
  t2 <- (log(max(h_star) / min(h_star[h_star > 0]))) ** 2
  
  # Adjust for zero shape and rate
  if (any(shape == 0)) {
    shape[shape == 0] <- 0.001
  }
  
  if (any(rate == 0)) {
    rate[rate == 0] <- 0.001
  }
  
  return(list("shape" = shape, "rate" = rate, "t2" = t2))
}

#' Set tuning parameters
#'
#' @param tuning_parameters list of tuning_parameters, could contain any combination
#' of the listed tuning parameters
#' @param borrow choice of borrow, could be TRUE or FALSE
#' @param X design matrix for concurrent trial
#' @param X_0 design matrix for historical trial
#'
#' @return filled list of tuning_parameters
.set_tuning_parameters <- function(tuning_parameters = NULL, borrow, X, 
                                   X_0 = NULL) {
  tuning_parameters_out <- tuning_parameters
  
  if (borrow) {
    n_beta = ncol(X)
    defaults <- list("Jmax" = 5,
                     "pi_b" = 0.5,
                     "alpha" = 0.4,
                     "cprop_beta" = 0.5
    )
    
    if (!is.null(X_0)) {
      defaults$cprop_beta_0 <- 0.5
    }
    
    for (key in names(defaults)) {
      if (!key %in% names(tuning_parameters)) {
        tuning_parameters_out[[key]] <- defaults[[key]]
      }
    }
    
  } else {
    n_beta = ncol(X)
    defaults <- list("Jmax" = 5,
                     "pi_b" = 0.5,
                     "cprop_beta" = 0.5)
    
    for (key in names(defaults)) {
      if (!key %in% names(tuning_parameters)) {
        tuning_parameters_out[[key]] <- defaults[[key]]
      }
    }
  }
  
  set_to_default <- setdiff(names(tuning_parameters_out), names(tuning_parameters))
  if (length(set_to_default) > 0) {
    default_params_str <- paste(set_to_default, collapse = ", ")
    message("The following tuning_parameters were set to default: ", default_params_str)
  }
  return(tuning_parameters_out)
}

#' Set tuning parameters
#'
#' @param hyperparameters list of hyperparameters, could contain any combination
#' of the listed hyperparameters
#' @param model_choice choice of model, could be either of 'mix', 'uni' or 'all'
#'
#' @return filled list of tuning_parameters
.set_hyperparameters <- function(hyperparameters = NULL, model_choice) {
  hyperparameters_out <- hyperparameters
  if (model_choice == "mix") {
    defaults <- list("a_tau" = 1,
                    "b_tau" = 0.001,
                    "c_tau" = 1,
                    "d_tau" = 1,
                    "type" = "mix",
                    "p_0" = 0.8,
                    "a_sigma" = 1,
                    "b_sigma" = 1,
                    "phi" = 3, 
                    "clam_smooth" = 0.8)
    
    for (key in names(defaults)) {
      if (!key %in% names(hyperparameters)) {
        hyperparameters_out[[key]] <- defaults[[key]]
      }
    }
    
  } else if (model_choice == "all") {
    defaults <- list("a_tau" = 1,
                     "b_tau" = 0.001,
                     "c_tau" = 1,
                     "d_tau" = 1,
                     "type" = "all",
                     "p_0" = 0.8,
                     "a_sigma" = 1,
                     "b_sigma" = 1,
                     "phi" = 3, 
                     "clam_smooth" = 0.8)
    
    for (key in names(defaults)) {
      if (!key %in% names(hyperparameters)) {
        hyperparameters_out[[key]] <- defaults[[key]]
      }
    }
    
  } else if (model_choice == "uni") {
    defaults <- list("a_tau" = 1,
                     "b_tau" = 0.001,
                     "type" = "uni",
                     "a_sigma" = 1,
                     "b_sigma" = 1,
                     "phi" = 3, 
                     "clam_smooth" = 0.8)
  
  for (key in names(defaults)) {
    if (!key %in% names(hyperparameters)) {
      hyperparameters_out[[key]] <- defaults[[key]]
      }
  }
    
  } else if (model_choice == "no_borrow") {
    defaults <- list("a_sigma" = 1,
                     "b_sigma" = 1,
                     "phi" = 3, 
                     "clam_smooth" = 0.8)
    
    for (key in names(defaults)) {
      if (!key %in% names(hyperparameters)) {
        hyperparameters_out[[key]] <- defaults[[key]]
      }
    }
  }
  
  set_to_default <- setdiff(names(hyperparameters_out), names(hyperparameters))
  set_to_default <- set_to_default[set_to_default != "type"]
  if (length(set_to_default) > 0) {
    default_params_str <- paste(set_to_default, collapse = ", ")
    message("The following hyperparameters were set to default: ", default_params_str)
  }
  message("Choice of borrowing: ", model_choice)
  return(hyperparameters_out)
}
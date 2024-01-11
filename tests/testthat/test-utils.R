test_that(".dataframe_fun.R has correct formatting", {
  Y <-  c(9, 13, 13, 18, 23)
  I <- c(1 ,1 ,0, 1, 1)
  X <- NULL
  s_r <- 10
  s <- c(0, sort(s_r), max(Y))
  lambda <- c(0.1, 0.9)
  J <- 1

  order <- c("tstart","id","Y", "I", "lambda")
  
  # Order, no covariates
  df <- .dataframe_fun(Y, I, X, s, lambda, bp = 0, J)
  expect_equal(names(df), order)
  
  # Order, with covariates
  X <- matrix(c(stats::rbinom(length(Y),5,0.5), stats::rbinom(length(Y),5,0.5), 
         stats::rbinom(length(Y),5,0.5)), ncol = 3)
  df <- .dataframe_fun(Y, I, X, s, lambda, bp = 3, J)
  order <- c("tstart","id","X1", "X2", "X3","Y", "I", "lambda")
  expect_equal(names(df), order)
  
  # Data Consistency, number of events is equal before and after split
  expect_equal(1, sum(df[df$tstart == 0,]$I))
  expect_equal(3, sum(df[df$tstart == 10,]$I))
  
  # Events are the same even when data is not ordered
  Y_wrong <-  c(13, 18, 23, 13, 9)
  I_wrong <- c(1 ,1 ,1, 0, 1)
  expect_equal(1, sum(df[df$tstart == 0,]$I))
  expect_equal(3, sum(df[df$tstart == 10,]$I))
  
  # Expect error when s is defined wrong (needs to be >= 2)
  # Written specifically in code, should not occur
  
  # Missing Data Handling
  Y_missing <-  c(9, NA, 13, 18, NA)
  I_missing <- c(1 , NA, 0, 1, NA)
  df_missing <- .dataframe_fun(Y_missing, I_missing, X, s, lambda, bp = 0, J)
  expect_equal(sum(is.na(df_missing$Y)), sum(is.na(Y_missing)))
  
  # there is some weird behavior in survSplit where it will give different values of I
  # (NA or 0) every time you run it, not sure if this interferes with anything [Q]
  
  # Empty dataframe, should not be able to split data
  expect_error(data.frame.fun(c(), c(), c(), s, lambda, bp = 0, J))
  
  # Performance Test
  X_large <- matrix(nrow = 1000, ncol = 1000)
  Y_large <- 1:1000
  I_large <- rbinom(n=1000, size=1, prob=0.5)
  for (ii in 1:1000) {
    X_large[ii,] <- 1:1000
  }
  system.time(df_large <- .dataframe_fun(Y_large, I_large, X_large, s, lambda, bp = ncol(X_large), J))
  expect_true(ncol(df_large) == 5+1000)
})

test_that(".logsumexp operates as it should", {
  
  # Postive values
  input_pos <- c(0.5, 0.4, 0.1, 0.01)
  expected_result <- log(sum(exp(input_pos)))
  expect_equal(.logsumexp(input_pos), expected_result)
  
  # Negative values
  input_neg <- c(-0.5, -0.4, -0.1, -0.01)
  expected_result <- log(sum(exp(input_neg)))
  expect_equal(.logsumexp(input_neg), expected_result)
  
  # All zero
  input_zero <- c(0, 0, 0, 0, 0)
  expected_result <- log(sum(exp(input_zero)))
  expect_equal(.logsumexp(input_zero), expected_result)
  
  # NAs present, should they be ignored? Or checked before? [Q]
  input_NA <- c(0.5, 0.4, NA, 0.01)
  expected_result <- log(sum(exp(input_NA)))
  expect_equal(.logsumexp(input_NA), expected_result)
  
})

test_that(".normalize_prob operates as it should", {
  
  # Postive values
  input_pos <- c(0.5, 0.4, 0.1, 0.01)
  expected_result <- exp(input_pos - .logsumexp(input_pos))
  expect_equal(.normalize_prob(input_pos), expected_result)
  
  # Negative values
  input_neg <- c(-0.5, -0.4, -0.1, -0.01)
  expected_result <- exp(input_neg - .logsumexp(input_neg))
  expect_equal(.normalize_prob(input_neg), expected_result)
  
  # All zero
  input_zero <- c(0, 0, 0, 0, 0)
  expected_result <- exp(input_zero - .logsumexp(input_zero))
  expect_equal(.normalize_prob(input_zero), expected_result)
  
  # NAs present, should they be ignored? Or checked before? [Q]
  input_NA <- c(0.5, 0.4, NA, 0.01)
  expected_result <- exp(input_NA - .logsumexp(input_NA))
  expect_equal(.normalize_prob(input_NA), expected_result)
  
})

test_that(".log_likelihood operates as it should", {
  Y <-  c(9, 13, 13, 18, 23)
  I <- c(1 ,1 ,0, 1, 1)
  X <- NULL
  s_r <- 10
  s <- c(0, sort(s_r), max(Y))
  lambda <- c(0.1, 0.9)
  J <- 1
  df <- .dataframe_fun(Y, I, X, s, lambda, bp = 0, J)
  
  ## Without covariates
  # Empty data.frames
  expect_error(.log_likelihood(df = c(), beta = NULL))
  
  # NAs present, should it return NA? How should we handle NAs? [Q]
  Y_NA <-  c(NA, 13, 13, 18, 23)
  I_NA <- c(NA , 1, 0, 1, 1)
  df <- .dataframe_fun(Y_NA, I_NA, X, s, lambda, bp = 0, J)
  expect_equal(sum(is.na(.log_likelihood(df, beta = NULL))), 1)
  
  ## With covariates
  X <- matrix(c(stats::rbinom(length(Y),5,0.5), stats::rbinom(length(Y),5,0.5), 
                stats::rbinom(length(Y),5,0.5)), ncol = 3)
  beta = c(-2, 0.5, 1)
  df <- .dataframe_fun(Y, I, X, s, lambda, bp = 3, J)
  
  # NAs present, should it return NA? How should we handle NAs? [Q]
  Y_NA <-  c(NA, 13, 13, 18, 23)
  I_NA <- c(NA , 1, 0, 1, 1)
  df <- .dataframe_fun(Y_NA, I_NA, X, s, lambda, bp = 3, J)
  expect_equal(sum(is.na(.log_likelihood(df, beta = beta))), 1)
  
  # Message for zero/high/low?
  
})

test_that("Input check: valid inputs do not produce errors/warnings", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list("a_tau" = 1, 
                "b_tau" = 1, 
                "c_tau" = 1, 
                "d_tau" = 1, 
                "p_0" = 0.1, 
                "clam_smooth" = 0.7, 
                "pi_b" = 0.3, 
                "type" = "mix", 
                "a_sigma" = 1, 
                "b_sigma" = 2, 
                "cprop_beta" = 1, 
                "phi" = 2)
  initial_param <- list("J" = 2, 
                        "s_r" = c(5, 10), 
                        "mu" = 5, 
                        "sigma2" = 0.5, 
                        "tau" = c(0.1, 0.2, 0.3), 
                        "lambda_0" = c(0.1, 0.2, 0.3), 
                        "lambda" = c(0.1, 0.2, 0.3), 
                        "beta_0" = NULL, 
                        "beta" = -0.5) 
 expect_no_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param))
})

test_that("Input check: Negative hyperparameter values produce errors", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "uni", 
                a_sigma = -1, # error
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "negative")
})

test_that("Input check: Hyperparameters outside [0, 1] range produce errors", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 1.3, #error
                type = "abc", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "range")
})

test_that("Input check: borrowing type 'uni' with given c_tau and d_tau produces a warning", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "uni", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_message(.input_check(Y, Y_0, X, X_0, hyper, initial_param), "Borrowing is 'uni'")
})

test_that("Input check: s_r values greater than maxSj produce errors", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 50), #error
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "s_r must be <")
})

test_that("Input check: negative sigma2 value produces an error", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = -0.5, #error
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "sigma2 must be > 0")
})

test_that("Input check: incorrect dimensions for tau produce error", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3, 0.1), # error
                        lambda_0 = c(0.1, 0.2, 0.3, 0.1), # error
                        lambda = c(0.1, 0.2, 0.3, 0.1), # error
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "dimension")
})

test_that("Input check: incorrect dimensions for lambda produce error", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3, 0.1), # error
                        lambda = c(0.1, 0.2, 0.3),
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), regexp = "dimension")
})

test_that("Input check: Negative lambdas renders error", {
  Y <- c(2, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(-0.1, 0.2, 0.3), #error
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, 
                        beta = -0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), 
               regexp = "baseline hazard")
  initial_param$lambda <- c(-0.1, 0.2, 0.3)
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), 
               regexp = "baseline hazard")
})

test_that("Input check: throws error for wrong dimension of beta/beta_0 when adding covariates", {
  Y <- c(2, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(c(1,0, 1), nrow = 3, ncol = 1)
  X_0 <- matrix(c(1,0, 0), nrow = 3, ncol = 1)
  hyper <- list(a_tau = 1, 
                b_tau = 1, 
                c_tau = 1, 
                d_tau = 1, 
                p_0 = 0.1, 
                clam_smooth = 0.7, 
                pi_b = 0.3, 
                type = "mix", 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3), 
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = NULL, # ERROR
                        beta = 0.5) 
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), 
               regexp = "beta_0")
  initial_param$lambda <- c(-0.1, 0.2, 0.3)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        tau = c(0.1, 0.2, 0.3), 
                        lambda_0 = c(0.1, 0.2, 0.3),
                        lambda = c(0.1, 0.2, 0.3), 
                        beta_0 = 2, 
                        beta = c(1, 2)) # ERROR
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param), 
               regexp = "beta")
  
})

test_that("group_summary() is working as it should", {
  Y <- stats::rweibull(20, 10, 0.4)
  I <- stats::rbinom(20, 1, 0.2)
  
  #W Without covariates
  X <- NULL
  
  s <- c(0, quantile(Y, probs = c(0.25, 0.75)), max(Y))
  
  # Throws no error for normal input
  expect_no_error(group_no_covariates <- group_summary(Y, I, X, s))
  expect_true(any(!is.na(group_no_covariates)))
  
  # Has correct output
  expect_equal(sum(group_no_covariates$events), sum(I))
  expect_equal(sum(group_no_covariates$num_cnsr), sum(I == 0))
  
  #W With covariates
  X <- stats::rbinom(20, 1, 0.5)
  
  # Throws no error for normal input
  expect_no_error(group <- group_summary(Y, I, X, s))
  expect_true(any(!is.na(group)))
  
  # Has correct output
  expect_equal(sum(group$events_c), sum(I[X == 0])) # control
  expect_equal(sum(group$events_t), sum(I[X == 1])) # treatment
  expect_true(all(group$num_at_risk_c + group$num_at_risk_t == group_no_covariates$num_at_risk))
  
  
})

test_that("init_lambda_hyperparameters() is working as it should", {
  Y <- stats::rweibull(20, 10, 0.4)
  I <- stats::rbinom(20, 1, 0.2)
  
  # Without covariates
  X <- NULL
  
  s <- c(0, quantile(Y, probs = c(0.25, 0.75)), max(Y))
  
  # Throws no error for normal input
  group_data <- group_summary(Y, I, X, s)
  expect_no_error(lambdas <- init_lambda_hyperparameters(group_data, s, w = 0.5))
  
  # Shape is correct:
  expect_length(lambdas$shape, length(s) - 1)
  expect_length(lambdas$rate, length(s) - 1)

})

test_that("Input check runs for only historical/current data (NoBorrow)" ,{
  Y <- c(1, 20, 30)
  Y_0 <- NULL
  X <- matrix(1, nrow = 3, ncol = 1)
  X_0 <- NULL
  hyper <- list(clam_smooth = 0.7, 
                pi_b = 0.3, 
                a_sigma = 1, 
                b_sigma = 2, 
                cprop_beta = 1, 
                phi = 2)
  initial_param <- list(J = 2, 
                        s_r = c(5, 10), 
                        mu = 5, 
                        sigma2 = 0.5, 
                        lambda = c(0.1, 0.2, 0.3), # error
                        beta = -0.5) 
  expect_no_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param))
  
  s <- .input_check(Y, Y_0, X, X_0, hyper, initial_param)
  expect_equal(s, "No borrowing")
})

test_that("Input check for cprop_beta dimensions", {
  Y <- c(1, 20, 30)
  Y_0 <- c(1, 20, 31)
  X <- matrix(1, nrow = 3, ncol = 3)
  X_0 <- NULL
  hyper <- list("a_tau" = 1, 
                "b_tau" = 1, 
                "c_tau" = 1, 
                "d_tau" = 1, 
                "p_0" = 0.1, 
                "clam_smooth" = 0.7, 
                "pi_b" = 0.3, 
                "type" = "mix", 
                "a_sigma" = 1, 
                "b_sigma" = 2, 
                "cprop_beta" = c(1, 1, 1), 
                "phi" = 2)
  initial_param <- list("J" = 2, 
                        "s_r" = c(5, 10), 
                        "mu" = 5, 
                        "sigma2" = 0.5, 
                        "tau" = c(0.1, 0.2, 0.3), 
                        "lambda_0" = c(0.1, 0.2, 0.3), 
                        "lambda" = c(0.1, 0.2, 0.3), 
                        "beta_0" = NULL, 
                        "beta" = c(-1, 1, 1)) 
  expect_no_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param))
  
  hyper <- list("a_tau" = 1, 
                "b_tau" = 1, 
                "c_tau" = 1, 
                "d_tau" = 1, 
                "p_0" = 0.1, 
                "clam_smooth" = 0.7, 
                "pi_b" = 0.3, 
                "type" = "mix", 
                "a_sigma" = 1, 
                "b_sigma" = 2,
                "cprop_beta" = c(1, 1), #
                "phi" = 2)
  expect_error(.input_check(Y, Y_0, X, X_0, hyper, initial_param),
                 regexp = "dimension mismatch")
})
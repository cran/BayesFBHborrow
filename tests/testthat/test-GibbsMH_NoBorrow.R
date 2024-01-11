set.seed(2023)

# Simulate larger, historic data
n_0 <- 250
lambda <- c(1.0, 3.0)
id <- stats::runif(n_0)
segment1_times <- stats::rexp(sum(id <= 0.5), rate = lambda[1])
segment2_times <- stats::rexp(sum(id > 0.5), rate = lambda[2])

Y_0 <- c(segment1_times, segment2_times)
censoring <- stats::rbinom(n_0, size = 1, prob = 0.3)
I_0 <- ifelse(censoring == 1, 0, 1)
X_0 <- NULL
initial_parameters <- list("J" = 1, 
                           "s_r" = c(1.0),
                           "mu" = 1, 
                           "sigma2" = 2,
                           "lambda" = lambda,
                           "beta" = NULL)

test_that("Gibbs_NoBorrow: output parameters are correct", {
  hyperparameters <- list("a_sigma" = 2,
                           "b_sigma" = 2,
                           "Jmax" = 20, 
                           "clam_smooth" = 0.8, 
                           "cprop_beta" = 0.3, #controls proposal sd for beta, beta0, 
                           "phi" = 3, 
                           "pi_b" = 0.5)
  iter <- 100
  
  .input_check(Y_0, NULL, X_0, NULL, hyperparameters, initial_parameters)
  
  expect_no_error(GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL, 
                                          initial_parameters,
                                         hyperparameters,
                                         iter = iter))
  
  out <- GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL, 
                          initial_parameters,
                          hyperparameters,
                          iter = iter)
  # 10 outputs
  expect_equal(length(out), 6)
  
  # 4 "fixed" parameters: J, mu, sigma
  expect_length(out$out_fixed, 3)
  expect_named(out$out_fixed, c("J", "mu", "sigma2"))
  
  # Dimensions are correct for split-point dependent parameters
  expect_equal(dim(out$s), c(iter, hyperparameters$Jmax + 2))
  expect_equal(dim(out$lambda), c(iter, hyperparameters$Jmax + 1))
  
  # Sampling parameters are scalars and above zero
  expect_type(out$lambda_move, "double")
  expect_gt(out$lambda_move, 0)
})


test_that("Gibbs_MH_NoBorrow runs for default parameter values", {
  initial_parameters <- list("J" = 1, 
                             "s_r" = c(1.0),
                             "mu" = 1, 
                             "sigma2" = 2,
                             "lambda" = c(1.0, 3.0),
                             "beta" = NULL)
  
  expect_true(all(!is.na(GibbsMH(
                                         Y_0, I_0, X_0, NULL, NULL, NULL, 
                                         initial_parameters)$out_fixed)))
  
})

## split points
test_that("Gibbs_MH_NoBorrow runs for zero split points", {
  initial_parameters <- list("J" = 0, 
                             "s_r" = NULL,
                             "mu" = 1, 
                             "sigma2" = 2,
                             "lambda" = c(1.0),
                             "beta" = NULL)
  
  hyperparameters <- list("a_sigma" = 2,
                           "b_sigma" = 2,
                           "Jmax" = 20, 
                           "clam_smooth" = 0.8, 
                           "cprop_beta" = 0.3, #controls proposal sd for beta, beta0, 
                           "phi" = 3, 
                           "pi_b" = 0.5)
  
  expect_true(all(!is.na(GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL, 
                                          initial_parameters,
                                         hyperparameters, 
                                         iter = 10)$out_fixed)))
  
  out <- GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL, 
                         initial_parameters,
                         hyperparameters, 
                         iter = 100)
  
  expect_gt(out$lambda_move, 0)
})


test_that("Gibbs_MH_NoBorrow runs for covariates on historical data", {
  # rate proposal blows up for medium -> large values of beta
  X_0 <- matrix(c(stats::rnorm(length(Y_0),0,2), stats::rbinom(length(Y_0),1,0.5), 
                  stats::rnorm(length(Y_0),0,2)), ncol = 3)
  
  initial_parameters <- list("J" = 1, 
                             "s_r" = c(0.5),
                             "mu" = 1, 
                             "sigma2" = 2,
                             "lambda" = c(1.0, 3.0),
                             "beta" = c(0.1, 0.1, 0.1))
  
  hyperparameters <- list("a_sigma" = 2,
                           "b_sigma" = 2,
                           "Jmax" = 20, 
                           "clam_smooth" = 0.8, 
                           "cprop_beta" = c(0.1, 0.1, 0.1), #controls proposal sd for beta, beta0, 
                           "phi" = 3, 
                           "pi_b" = 0.5)
  
  expect_true(all(!is.na(GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL, 
                                         initial_parameters,
                                         hyperparameters,
                                         iter = 10)$out_fixed)))
  
  expect_no_error(out <- GibbsMH(Y_0, I_0, X_0, NULL, NULL, NULL,  initial_parameters,
                         hyperparameters, 
                         iter = 100))
  
  # three covariates in -> three out
  expect_gt(out$lambda_move, 0)
  #expect_true(all(out$beta_move > 0))
  expect_length(out$beta_move, ncol(X_0))
  expect_length(out$out_fixed[paste0("beta_",1:3)],3)
  
})

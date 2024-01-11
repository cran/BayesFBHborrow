cprop_beta =c(1, 2)

Y_0 <-  c(7, 25, 51, 13, 64, 11)
I_0 <- c(1 ,1 ,0, 1, 1, 1)
X_0 <- matrix(stats::rnorm(length(Y_0)*2, 1, 0.5), ncol = 2)
Y <-  c(9, 13, 13, 18, 23)
I <- c(1 ,1 ,0, 1, 1)
X <- matrix(stats::rbinom(length(Y), 1, 0.5))
#s_r <- c(10, 15, 35) # one split point larger than Y_max, could this occur?
s_r <- c(10, 15)
s <- c(0, s_r, max(c(Y_0,Y)))
J <- 2

# without covariates, lambda_0
bp_0 <- ncol(X_0)
lambda_0 <- c(1.5, 1.4, 1.0)
beta_0 <- c(2, 5)

bp <- ncol(X)
lambda <- c(1.5, 1.4, 1.0)
beta <- 5

beta_new <- beta
beta_0_new <- beta_0

beta_new <- stats::rnorm(1, beta[1], cprop_beta)
beta_0_new[1] <- stats::rnorm(1, beta_0[1], cprop_beta)
beta_0_new[2] <- stats::rnorm(1, beta_0[2], cprop_beta)

df_hist <- .dataframe_fun(Y_0, I_0, X_0, s, lambda_0, bp = bp_0, J)
df <- .dataframe_fun(Y, I, X, s, lambda, bp = bp, J) 

test_that("loglikelihood ratio for both beta/beta_0 works ", {

  expect_no_error(llikeli <- .llikelihood_ratio_beta(df, beta, beta_new))
  expect_no_error(llikeli_0 <- .llikelihood_ratio_beta(df_hist, beta_0, beta_0_new))
  expect_true(!is.na(llikeli))
  expect_true(!is.na(llikeli_0))
  
  X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
  xdpb <- X %*% beta
  xdpb_new <- X %*% beta_new
  
  llikelihood_ratio <- sum((xdpb_new  - xdpb) * df$I - ((df$Y - df$tstart) * df$lambda) * (exp(xdpb_new) - exp(xdpb)))
  
  expect_true(!is.na(llikelihood_ratio))
  expect_equal(llikelihood_ratio, llikeli)
  })

test_that("Beta updates are correct", {
  set.seed(2023)
  new_beta <- .beta_MH_RW(df, beta, bp, cprop_beta[1], count_beta = 0)
  expect_true(all(!is.na(new_beta)))
  
  count_beta <- beta*0
  set.seed(2023)
  for(k in 1:bp){
      
    beta_new <- beta
    beta_prop <- stats::rnorm(1, beta[k], cprop_beta[k])
    beta_new[k] <- beta_prop
      
    logacc <- .llikelihood_ratio_beta(df, beta, beta_new)
      
    if(logacc > log(stats::runif(1))){
      beta[k] <- beta_prop
      count_beta[k] <- count_beta[k] + 1
    }
      
  }
  
  expect_equal(beta, new_beta$beta)
  
  # check acceptance ratio
  count_beta <- 0
  for (ii in 1:100) {
    new_beta <- .beta_MH_RW(df, beta, bp, cprop_beta[1], count_beta)
    count_beta <- new_beta$count_beta
  }
  expect_gt(count_beta, 20)
  
})

test_that("Mean proposal (mu_prop)", {
  mu_prop <- .beta_mom(df, 1, beta, bp, cprop_beta[1])
  expect_true(!is.na(mu_prop))
})


test_that("lprop_density_beta are correct", {
  mu_prop <- 1
  set.seed(2023)
  beta_prop <- stats::rnorm(n = 1, mean = mu_prop, sd = cprop_beta[1])
  
  # Negative (log density)
  expect_true(!is.na(.lprop_density_beta(beta_prop, mu_prop, cprop_beta[1])))
  expect_lt(.lprop_density_beta(beta_prop, mu_prop, cprop_beta[1]), 0)
  l1 <- (-1 / (2 * cprop_beta[1]**2)) * (beta_prop - mu_prop)**2
  set.seed(2023)
  expect_equal(l1, .lprop_density_beta(beta_prop, mu_prop, cprop_beta[1]))
})

test_that("MALA proposal is correct", {
  set.seed(2023)
  new_beta <- .beta_MH_MALA(df, beta, bp, cprop_beta, count_beta = 0)
  expect_true(all(!is.na(new_beta)))
  
  count_beta <- beta*0
  set.seed(2023)
  for(k in 1:bp){
    beta_new <- beta
    
    mu_prop <- .beta_mom(df, k, beta, bp, cprop_beta[k])
    beta_prop <- rnorm(n = 1, mean = mu_prop, sd = cprop_beta[k])
    beta_new[k] <- beta_prop
    
    mu_old <- .beta_mom(df, k, beta_new, bp, cprop_beta)
    
    log_prop_ratio <- .lprop_density_beta(beta, mu_prop, cprop_beta[k]) - 
      .lprop_density_beta(beta_prop, mu_old, cprop_beta[k])
    target_ratio <- .llikelihood_ratio_beta(df, beta, beta_new)
    
    logacc <- target_ratio - log_prop_ratio  
    
    if(logacc > log(runif(1))) {
      beta[k] <- beta_prop
      count_beta[k] <- count_beta[k] + 1
    }
    
  }
  
  expect_equal(beta, new_beta$beta)
  
  # check acceptance ratio
  count_beta <- 0
  for (ii in 1:100) {
    new_beta <- .beta_MH_MALA(df, beta, bp, cprop_beta, count_beta[1])
    count_beta <- new_beta$count_beta
  }
  expect_gt(count_beta, 20)
})

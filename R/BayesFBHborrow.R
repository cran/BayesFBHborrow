## BayesFBHborrow
## the main user facing entrypoint 
BayesFBHborrow <- function(formula = NULL,
                           data,
                           data_hist = NULL,                           
                           subset,
                           na.action,
                           CntlOnly = FALSE,
                           model_choice = "mix",
                           tuning_parameters = NULL,
                           hyperparameters = NULL,
                           iter = 1000,
                           warmup_iter = 100,
                           refresh = 0,
                           verbose = FALSE,
                           max_grid = 2000,
                           standardise = TRUE,
                           G_compute = FALSE,
                           preprocess = TRUE, 
                           fnm_cnd_blh = "cnd_blh",
                           fnm_mgnl_haz0 = "mgnl_haz0",
                           fnm_mgnl_haz1 = "mgnl_haz1",
                           dbg=FALSE,
                           ...){

  borrow <- !missing(data_hist)
  if(borrow && CntlOnly) stop("Set the argument 'CntlOnly' to TRUE only when you want to fit the model on control data only")
  if (borrow){
    class(data) <- c("WBorrow", "data.frame")
  } else {
    class(data) <- c("NoBorrow", "data.frame")
  }

  UseMethod("BayesFBHborrow", data)
}

## BayesFBHborrow.WBorrow
## non user facing-- called via method dispatch
BayesFBHborrow.WBorrow <-
    function(formula = NULL,
             data,
             data_hist = NULL,             
             subset,
             na.action,
             CntlOnly = FALSE, 
             model_choice = "mix",
             tuning_parameters = NULL,
             hyperparameters = NULL,
             iter = 1000,
             warmup_iter = 100,
             refresh = 0,
             verbose = FALSE,
             max_grid = 2000,
             standardise = TRUE,
             G_compute = FALSE,
             preprocess = TRUE, 
             fnm_cnd_blh = "cnd_blh",
             fnm_mgnl_haz0 = "mgnl_haz0",
             fnm_mgnl_haz1 = "mgnl_haz1",
             dbg=FALSE,
             ...){
   ## when doing borrowing we have a current data set with a treatment indicator and controls only historical dataset
   if(dbg)browser()
  .call. <- mf <- match.call()
  .call.[[1]] <- as.name("BayesFBHborrow")

  r_sfx <- rand_hex_strng(6)
  fnm_cnd_blh = fnm_cnd_blh %,% "_" %,% r_sfx %,% "_wb.bin"
  fnm_mgnl_haz0 = fnm_mgnl_haz0 %,% "_" %,% r_sfx %,% "_wb.bin"
  fnm_mgnl_haz1 = fnm_mgnl_haz1 %,%  "_" %,%r_sfx %,% "_wb.bin"

  if(missing(CntlOnly)) CntlOnly <- FALSE
  if(!missing(formula))
  {
      ## formula interface with or without pre-processing as per user stipulation
      ## BFBHB.form.intrfc is the engine which does the preprocessing (or not)
      mf <- .call.
      mf[[1]] <- as.name("BFBHB.form.intrfc")

      mf$model_choice <- mf$tuning_parameters <- mf$hyperparameters <- 
          mf$iter <- mf$warmup_iter <- mf$refresh <- mf$verbose <- mf$max_grid <- 
              mf$standardise <- mf$G_compute <- mf$fnm_cnd_blh <- mf$fnm_mgnl_haz0 <- 
                  mf$fnm_mgnl_haz1 <- NULL
      mf <- eval(mf)
                 
      R <- mf$R
      X <- as.data.frame(mf$X)
      R_0 <- mf$R_0
      X_0 <- as.data.frame(mf$X_0)
      nms_fct <- mf$nms_fct
      nms_cnt <- mf$nms_cnt
      sd_X <- mf$sd_X
  } else {
      ## data only interface
      checkmate::assert_flag(verbose)
      checkmate::assert_flag(CntlOnly)
      checkmate::assert_data_frame(data, any.missing = FALSE, null.ok = FALSE)
      checkmate::assert_names(names(data), must.include = c("tte", "event"))
      checkmate::assert_data_frame(data[grepl("X",names(data))],any.missing=FALSE,min.cols=1*(!CntlOnly))
      checkmate::assert_data_frame(data_hist, any.missing = FALSE)
      checkmate::assert_names(names(data_hist), must.include = c("tte", "event"))

      ## if preprocess is set to TRUE, use a call to BFBHB.form.intrfc with a formula constructed from
      ## the expected naming format to do the preprocessing
      if(preprocess)
      {
          bfbhbfi.call <- as.call(expression(BFBHB.form.intrfc))         
          d_nms <- names(data)
          x_nms <- d_nms[grep("X_", d_nms)]
          if(length(x_nms)==0) x_nms <- "1"
          form <- as.formula("Surv(tte, event) ~ " %,% paste(x_nms, collapse=" + "))
          bfbhbfi.call$formula <- form
          bfbhbfi.call$data <- as.name("data")
          bfbhbfi.call$data_hist <- as.name("data_hist")
          if(!missing(subset)) bfbhbfi.call$subset <- as.name("subset")
          if(!missing(na.action)) bfbhbfi.call$na.action <- as.name("na.action")
          bfbhbfi.call$preprocess <- TRUE
          bfbhbfi <- eval(bfbhbfi.call)
                  
          R <- bfbhbfi$R
          X <- bfbhbfi$X
          R_0 <- bfbhbfi$R_0
          X_0 <- bfbhbfi$X_0
          nms_fct <- bfbhbfi$nms_fct
          nms_cnt <- bfbhbfi$nms_cnt
          sd_X <- bfbhbfi$sd_X                  
      } else {
          ## otherwise just pull the names according to the expected naming format directly from the data
          R <- as.matrix(data[,c("tte","event")])
          X <- data[grepl("X_", names(data))]
          R_0 <- data_hist[,c("tte","event")]
          if (length(data_hist[grepl("X_", names(data_hist))]) != 0) {
              checkmate::assert_data_frame(data_hist[grepl("X_", names(data_hist))], any.missing = FALSE, min.cols = 1)
              X_0 <- data_hist[grepl("X_", names(data_hist))]
          } else {
              X_0 <- NULL
          }
          nms_cnt <- NULL
      }
  }        
  if (verbose) {
    hyperparameters <- set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice)
    tuning_parameters <- set_tuning_parameters(tuning_parameters = tuning_parameters,
                                                borrow = TRUE, X = X, X_0 = X_0)
    message("Starting MCMC sampler")
  } else {
    suppressMessages(
      hyperparameters <-
        set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice)
    )
    suppressMessages(
      tuning_parameters <-
        set_tuning_parameters(tuning_parameters = tuning_parameters, borrow = TRUE, X = X, X_0 = X_0))
  }

  out <- GibbsMH_wborrow(Y = R[,1],
                         I = R[,2],
                         X = X,
                         Y_0 = R_0[,1],
                         I_0 = R_0[,2],
                         X_0 = X_0,
                         tuning_parameters = tuning_parameters,
                         hyperparameters = hyperparameters,
                         iter = iter,
                         warmup_iter = warmup_iter,
                         refresh = refresh,
                         max_grid = max_grid,
                         CntlOnly = CntlOnly,
                         standardise = standardise,
                         G_compute = G_compute,
                         fnm_cnd_blh = fnm_cnd_blh,
                         fnm_mgnl_haz0 = fnm_mgnl_haz0,
                         fnm_mgnl_haz1 = fnm_mgnl_haz1)

  out$beta_acc_ratio <- out$beta_move / iter
  if (!is.null(out$beta_0_move)) out$beta_0_acc_ratio <- out$beta_0_move / iter

  if(length(nms_cnt)>0)
  {
      out$samples[,nms_cnt] <- t(t(out$samples[,nms_cnt])/sd_X)
      out$samples[,nms_cnt %,% "_0"] <- t(t(out$samples[,nms_cnt %,% "_0"])/sd_X)
  }
        
  class(out) <- c("BayesFBHborrow", "list")
  out$call <- .call.
  return(out)
}

## Without borrowing method
## BayesFBHborrow.NoBorrow
BayesFBHborrow.NoBorrow <-
  function(formula = NULL,
             data,
             data_hist = NULL,             
             subset,
             na.action,
             CntlOnly = FALSE,
             model_choice = "mix",
             tuning_parameters = NULL,
             hyperparameters = NULL,
             iter = 1000,
             warmup_iter = 100,
             refresh = 0,
             verbose = FALSE,
             max_grid = 2000,
             standardise = TRUE,
             G_compute = FALSE,
             preprocess = TRUE, 
             fnm_cnd_blh = "cnd_blh",
             fnm_mgnl_haz0 = "mgnl_haz0",
             fnm_mgnl_haz1 = "mgnl_haz1",
             dbg=FALSE, 
             ...){

  if(dbg)browser()
  .call. <- match.call()
  .call.[[1]] <- as.name("BayesFBHborrow")

  r_sfx <- rand_hex_strng(6)
  fnm_cnd_blh = fnm_cnd_blh %,% "_" %,% r_sfx %,% "_nb.bin"
  fnm_mgnl_haz0 = fnm_mgnl_haz0 %,% "_" %,% r_sfx %,% "_nb.bin"
  fnm_mgnl_haz1 = fnm_mgnl_haz1 %,% "_" %,% r_sfx %,% "_nb.bin"

  if(!missing(formula))
  {
      ## formula interface with or without pre-processing as per user stipulation
      ## BFBHB.form.intrfc is the engine which does the preprocessing (or not)
      mf <- .call.
      mf[[1]] <- as.name("BFBHB.form.intrfc")
      
      mf$model_choice <- mf$tuning_parameters <- mf$hyperparameters <- 
          mf$iter <- mf$warmup_iter <- mf$refresh <- mf$verbose <- mf$max_grid <- 
              mf$standardise <- mf$G_compute <- mf$fnm_cnd_blh <- mf$fnm_mgnl_haz0 <- 
                  mf$fnm_mgnl_haz1 <- NULL
      mf <- eval(mf)
                 
      R <- mf$R
      X <- as.data.frame(mf$X)
      nms_fct <- mf$nms_fct
      nms_cnt <- mf$nms_cnt
      sd_X <- mf$sd_X
  } else {
      ## data only interface 
      checkmate::assert_data_frame(data, any.missing = FALSE, null.ok = FALSE)
      checkmate::assert_names(names(data), must.include = c("tte", "event"))
      checkmate::assert_flag(verbose)      
      checkmate::assert_flag(CntlOnly)
      checkmate::assert_data_frame(data[grepl("X", names(data))], any.missing = FALSE, min.cols = 1*(!CntlOnly))
      
      ## if preprocess is set to TRUE, use a call to BFBHB.form.intrfc with a formula constructed from
      ## the expected naming format to do the preprocessing
      if(preprocess)
      {
          bfbhbfi.call <- as.call(expression(BFBHB.form.intrfc))         
          d_nms <- names(data)
          x_nms <- d_nms[grep("X_", d_nms)]
          if(length(x_nms)==0) x_nms <- "1"
          form <- as.formula("Surv(tte, event) ~ " %,% paste(x_nms, collapse=" + "))
          bfbhbfi.call$formula <- form
          bfbhbfi.call$data <- as.name("data")
          if(!missing(subset)) bfbhbfi.call$subset <- as.name("subset")
          if(!missing(na.action)) bfbhbfi.call$na.action <- as.name("na.action")
          bfbhbfi.call$preprocess <- TRUE
          bfbhbfi <- eval(bfbhbfi.call)
                  
          R <- bfbhbfi$R
          X <- bfbhbfi$X
          nms_fct <- bfbhbfi$nms_fct
          nms_cnt <- bfbhbfi$nms_cnt
          sd_X <- bfbhbfi$sd_X                            
      } else {
          R <- as.matrix(data[,c("tte","event")])
          X <- data[grepl("X", names(data))]
      }
      nms_cnt <- NULL
  }

  if (verbose) {
    hyperparameters <- set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice)
    tuning_parameters <- set_tuning_parameters(tuning_parameters = tuning_parameters,
                                               borrow = FALSE, X = X, X_0 = NULL)
    message("Starting MCMC sampler")
  } else {
    suppressMessages(hyperparameters <- set_hyperparameters(hyperparameters = hyperparameters, model_choice = model_choice))
    suppressMessages(tuning_parameters <- set_tuning_parameters(tuning_parameters = tuning_parameters,
                                                                borrow = FALSE, X = X, X_0 = NULL))
  }

  out <- GibbsMH_noborrow(Y = R[,1],
                          I = R[,2],
                          X = X,
                          tuning_parameters = tuning_parameters,
                          hyperparameters = hyperparameters,
                          iter = iter,
                          warmup_iter = warmup_iter,
                          refresh = refresh,
                          max_grid = max_grid,
                          CntlOnly = CntlOnly, 
                          standardise = standardise,                          
                          G_compute = G_compute,
                          fnm_cnd_blh = fnm_cnd_blh,
                          fnm_mgnl_haz0 = fnm_mgnl_haz0,
                          fnm_mgnl_haz1 = fnm_mgnl_haz1)

  out$beta_acc_ratio <- out$beta_move / iter

  if(length(nms_cnt)>0) out$samples[,nms_cnt] <- t(t(out$samples[,nms_cnt])/sd_X)
        
  class(out) <- c("BayesFBHborrow", "list")
  out$call <- .call.
  return(out)
}

## The backend of the formula/data interface which will do covariate centering as explained in the documentation
## when preprocess=TRUE.  Called via a constructed formula for the data only interface, so this functionality is
## ubiquitous across the two interface methods
BFBHB.form.intrfc <-
function(formula, data, data_hist=NULL, subset, na.action, CntlOnly = FALSE, preprocess=TRUE, dbg=FALSE)
{
  if(dbg)browser()
  .call. <- match.call()
  xtrt.nm <- attr(terms(formula), "term.labels")[1]

  is.hist <- !missing(data_hist)
  n_cc <- nrow(data)
  n_hst <- 0
  
  if(is.hist) n_hst <- nrow(data_hist)
  
  mf <- as.call(expression(model.frame))        
  mf$formula <- formula
  mf$data <- as.name("data")
  mf$subset <- .call.$subset
  mf$na.action <- .call.$na.action
  
  mf <- eval(mf)
  tl <- attr(terms(mf), "term.labels")
  nc_X <- length(tl)
  is_X <- nc_X > 0
  is_2arm <- !CntlOnly
    
  R <- model.response(mf)

  ## apply the sum to 0 contrast to categorical variables
  Terms <- terms(mf)

  if(preprocess)
  {
      X_clss <- attr(Terms, "dataClasses")  
      nms_fct <- names(X_clss)[which(X_clss %in% c("factor", "character"))]
      if(length(nms_fct)>0)
      {
          con_lst <- centered_contrasts_lst(mf[,nms_fct,drop=FALSE])
          X <- model.matrix(Terms, data=mf, contrasts.arg=con_lst)[,-1]
      } else X <- model.matrix(Terms, data=mf)[,-1]
      
      ## which X other than the treatment variable (drop 1st one)
      ## are numeric? center and standardize those
      nms_cnt <- names(X_clss)[which(X_clss == "numeric")]
      if(is_2arm)nms_cnt <- nms_cnt[-1]
      if(length(nms_cnt)>0)
      {
          mu_X <- colMeans(X[,nms_cnt,drop=FALSE])
          sd_X <- apply(X[,nms_cnt,drop=FALSE], 2, FUN=var)^0.5
          X[,nms_cnt] <- t((t(X[,nms_cnt,drop=FALSE]) - mu_X)/sd_X)
      } else {
        nms_fct <- nms_cnt <- character(0)
        sd_X <- NULL
      }          
  } else {
        X <- model.matrix(Terms, data=mf)[,-1]
        nms_fct <- nms_cnt <- character(0)
        sd_X <- rep(1, ncol(X))
  }

  if(is.hist)
  {
      mf <- as.call(expression(model.frame))
      mf$formula <- drop_first_rhs(formula)
      mf$data <- as.name("data_hist")
      mf$subset <- .call.$subset
      mf$na.action <- .call.$na.action

      mf <- eval(mf)

      R_0 <- model.response(mf)

  ## apply the sum to 0 contrast to categorical variables
      Terms <- terms(mf)

      if(preprocess)
      {
          X_clss <- attr(Terms, "dataClasses")  
          nms_fct <- names(X_clss)[which(X_clss %in% c("factor", "character"))]
          if(length(nms_fct)>0)
          {
              con_lst <- centered_contrasts_lst(mf[,nms_fct,drop=FALSE])
              X_0 <- model.matrix(Terms, data=mf, contrasts.arg=con_lst)[,-1]
          } else X_0 <- model.matrix(Terms, data=mf)[,-1]
          
          ## which X other than the treatment variable (drop 1st one)
          ## are numeric? center and standardize those
          nms_cnt <- names(X_clss)[which(X_clss == "numeric")]
          if(length(nms_cnt)>0) X_0[,nms_cnt] <- t((t(X_0[,nms_cnt,drop=FALSE]) - mu_X)/sd_X)
      } else {
          X_0 <- model.matrix(Terms, data=mf)[,-1]
          nms_fct <- nms_cnt <- character(0)
          sd_X <- rep(1, ncol(X_0))
      }
  }
  
  out <- list(R=R, X=as.data.frame(X), nms_fct=nms_fct, nms_cnt=nms_cnt, sd_X=sd_X)
  if(is.hist)
  {
      out$R_0 <- R_0
      out$X_0 <- as.data.frame(X_0)
  }
  out
}

## Function which is passed to the model.matrix function which defines how contrasts are to be
## recoded into indicators and possibly centered.
centered_contrasts_lst <-
function(data)
{
  p <- ncol(data)
  CC_lst <- list()  
  for(k in 1:p)
  {
      x <- as.factor(data[,k])
      levs <- levels(x)
      d <- length(levs)
      n <- length(x)

      ## weighted: turn off
      ## Count occurrences of each level
      ## n_k <- table(x)
      ## n_ref <- n_k[levs[1]]
      ## a <- -n_k[setdiff(levs, levs[1])]/n_ref

      ## unweighted
      a <- rep(-1, length(levs)-1)
      
      CC <- rbind(a, diag(d-1))
      dimnames(CC) <- list(levs, levs[-1])
      CC_lst[[k]] <- CC
  }
  names(CC_lst) <- names(data)
  CC_lst
}

## function to standardise a continuous variable
stdrdz <- \(x)(x-mean(x))/var(x)^0.5

## uses the functionality of the rlang package to
## construct a new formula from an old one with the
## first regressor dropped
drop_first_rhs <- function(f) {
  lhs <- f_lhs(f)
  rhs <- f_rhs(f)
  env <- f_env(f)
  
  # Handle single variable case (e.g., Y ~ X_trt)
  if (is_symbol(rhs)) {
    return(new_formula(lhs, 1, env = env))
  }
  
  # Recursive helper to find and remove the leftmost leaf
  rm_leftmost <- function(x) {
    # If the left side of the '+' is not a call, it's the first variable
    if (!is_call(x[[2]], "+")) {
      return(x[[3]]) # Return only the right side of this '+'
    }
    
    # Otherwise, keep digging left
    new_left <- rm_leftmost(x[[2]])
    call2("+", new_left, x[[3]])
  }
  
  new_rhs <- rm_leftmost(rhs)
  new_formula(lhs, new_rhs, env = env)
}

## GibbsMH_wborrow
## The main MCMC/MH sampler function under the borrowing case 
## Arguements are the data for the current data Y, I, X,
## and historical trial data Y_0, I_0, X_0, 
## (time to event, event indicator covariates),
## and state variables. 
GibbsMH_wborrow <- function(Y,
                            I,
                            X,
                            Y_0,
                            I_0,
                            X_0,
                            tuning_parameters,
                            hyperparameters,
                            iter,
                            warmup_iter,
                            refresh,
                            max_grid,
                            CntlOnly,
                            standardise,
                            G_compute,
                            fnm_cnd_blh,
                            fnm_mgnl_haz0,
                            fnm_mgnl_haz1) {

  ### Initialize parameters ###
  # count accept for MH beta
  lambda_0_count <- 0
  lambda_count <- 0
  lambda_move <- 0
  lambda_0_move <- 0

  # proposal prior
  a_lambda <- tuning_parameters$a_lambda
  b_lambda <- tuning_parameters$b_lambda

  # hyperparameters
  a_sigma <- hyperparameters$a_sigma
  b_sigma <- hyperparameters$b_sigma
  clam <- hyperparameters$clam_smooth
  phi <- hyperparameters$phi

  stdz = 1.0
  if(standardise == TRUE){
    ## standardise time
    Y.raw <- Y
    Y_0.raw <- Y_0

    stdz <- sum(I) / sum(Y.raw)
    Y_0 <- Y_0.raw * stdz
    Y <- Y.raw * stdz
  }

  # Set initial values
  J <- phi
  sigma2 <- b_sigma / (a_sigma + 1)
  Y_mix <- Y[Y < max(Y_0)]
  quantiles <- quantile(Y_mix, probs = seq(0, 1, length.out = J + 2),
                        na.rm = TRUE, names = FALSE)
  s <- c(0, quantiles[1:J + 1], max(Y, Y_0))

  # s need to be unique
  if(length(unique(s)) < length(s)){
    J <- 1
    s <- c(0, max(Y, Y_0), max(Y))
  }

  ##Only works if X is a data.frame
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
  bp_0 <- ncol(X_0)
    
  df_curr <- dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp = ncol(X), J = J)

  if(bp_0 == 0){
    beta_0 <- NULL
  } else {
    df_hist <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0, bp = bp_0, J = J)
    glm.mom_0 <- glmFit(df_hist)
    beta_0 <- glm.mom_0$beta.mu
    #beta_0_count <- numeric(bp_0)
    ##### Testing: Multivariate proposal
    beta_0_count <- 0

    # prior as precision
    beta_0_prior <- hyperparameters$beta_0_prior
    beta_0_prior_prec <- diag(1/beta_0_prior, bp_0)
  }

  glm.mom <- glmFit(df_curr)
  beta <- glm.mom$beta.mu
  bp <- ncol(X)

  # Check for G-computaion
  if(bp == 1){
    G_check <- 0
  }else{
    G_check <- G_compute * 1
  }

  #beta_count <- numeric(bp)
  ##### Testing: Multivariate proposal
  beta_count <- 0
  beta_prior <- hyperparameters$beta_prior
  beta_prior_prec <- diag(1/beta_prior, bp)

  if(bp_0==0){
    df_hist <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0, bp = bp_0, J = J)
  }

  a_tau <- hyperparameters$a_tau
  b_tau <- hyperparameters$b_tau
  type <- hyperparameters$type
  Jmax <- hyperparameters$Jmax

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
  cprop_beta <- tuning_parameters$cprop_beta
  pi_b <- tuning_parameters$pi_b
  alpha <- tuning_parameters$alpha
  if (!is.null(X_0)) {
    cprop_beta_0 <- tuning_parameters$cprop_beta_0
  }

  maxSj <- min(max(Y), max(Y_0))

  ## treatment is always part of bp, with d pre-treatment covariates giving bp=d+1
    ## bp_0 is 0 without covariates and bp_0=d with d pre-treatment covariates
  nmsB <- dimnames(X)[[2]]  
  samples <- data.frame(matrix(NA, nrow = 0, ncol = bp + bp_0 + 3))   # 2d + 1 + 3
  if(bp > 0) colnames(samples)[1:bp] <- nmsB
  if(bp_0 > 0) {
    nmsB_0 <-   dimnames(X_0)[[2]] %,% "_0"
    colnames(samples)[(bp + 1):(bp + bp_0)] <- nmsB_0
  }
  colnames(samples)[bp + bp_0 + (1:3)] <- c("J", "mu", "sigma2")

  out_lambda <- data.frame(matrix(NA, nrow = 0, ncol = Jmax + 1))
  out_lambda_0 <- data.frame(matrix(NA, nrow = 0, ncol = Jmax + 1))
  out_s <- data.frame(matrix(NA, nrow = 0, ncol = Jmax + 2))

  if(type %in% c("uni", "mix")) {
    out_tau <- data.frame(matrix(NA, nrow = 0, ncol = Jmax +1))
  }else{
    out_tau <- data.frame(matrix(NA, nrow = 0, ncol = 1))
    colnames(out_tau) <- NULL
  }

  ## Define grid of time points
  time_grid <- seq(1e-8, max(Y), length.out =  max_grid)

  ## open binary file for writing out adjusted baseline hazard
  conn_cnd_blh <- file(fnm_cnd_blh, "wb")

  ## 1) ARRAY out_slam iter x max_grid(2000) ## conditional baseline hazards
  ## out_slam <-  data.frame(matrix(data = NA, nrow = 0, ncol =  max_grid))
  ## colnames(out_slam) <- time_grid

  # Array for marginal estimand from G-computation
  if(bp > 1 && !CntlOnly){
    nx <- nrow(X)
    trtAX <- cbind(rep(1, nx), as.matrix(X[,-1]))
    conAX <- as.matrix(X[,-1])

    ## If doing G-computation, open binary file for writing out
    ## conditional hazards in treatment (1) and control (0) arms
    conn_mgnl_haz1 <- file(fnm_mgnl_haz1, "wb")
    conn_mgnl_haz0 <- file(fnm_mgnl_haz0, "wb")

    ## 2) ARRAY lhaz_bar_trt iter x max_grid(2000) marginal treat hazard(t)
    ## lhaz_bar_trt <- data.frame(matrix(NA, nrow = 0, ncol = max_grid))
    ## colnames(lhaz_bar_trt) <- time_grid
    ## 3) ARRAY out_slam iter x max_grid(2000) marginal control hazard(t)
    ## lhaz_bar_con <- data.frame(matrix(NA, nrow = 0, ncol = max_grid))
    ## colnames(lhaz_bar_con) <- time_grid

    w_sxn_trt <- extraDistr::rdirichlet(iter, rep(1, nx)) # iter draws s x n
    w_sxn_con <- extraDistr::rdirichlet(iter, rep(1, nx))
  }

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

  j <- 1 # indx after warmup_iter

  for(i in 1:(iter + warmup_iter)) {
    if(i%%refresh == 0 && refresh != 0){message(mess[i])}

    ICAR <- ICAR_calc(s = s, J = J, clam = clam)
    Sigma_s <- ICAR$Sigma_s

    # 1. Conjugate posterior updates [Gibbs]
    mu <- mu_update(Sigma_s, lambda_0, sigma2, J)
    sigma2 <- sigma2_update(mu, lambda_0, Sigma_s, J, a_sigma, b_sigma)
    tau_all <- tau_update(lambda_0, lambda, J, s, a_tau, b_tau, c_tau, d_tau, p_0, type)
    tau <- tau_all$tau

    # 2. Update beta [MH NR]
    #beta_all <- beta_MH_NR(df = df_curr, beta, beta_prior_prec, bp, cprop_beta, beta_count)
    beta_all <- beta_MH_NR_vec(df_curr,
                                beta,
                                beta_prior_prec,
                                bp,
                                cprop_beta,
                                beta_count)
    beta <- beta_all$beta
    beta_count <- beta_all$beta_count

    if(bp_0 > 0) {
      #beta_0_all <- beta_MH_NR(df = df_hist, beta = beta_0, bp = bp_0,
      #                          cprop_beta = cprop_beta_0, beta_count = beta_0_count)
      beta_0_all <- beta_MH_NR_vec(df_hist,
                                    beta_0,
                                    beta_0_prior_prec,
                                    bp = bp_0,
                                    cprop_beta = cprop_beta_0,
                                    beta_count = beta_0_count)
      beta_0 <- beta_0_all$beta
      beta_0_count <- beta_0_all$beta_count
    }

    # 3. Update lambda_0, propose new lambda_0 from conditional [MH]
    #   conjugate posterior
    df_lambda_0 <- lambda_0_MH_cp(df_hist, Y_0, I_0, X_0, s,
                                   beta_0, mu, sigma2,  lambda, lambda_0, tau,
                                   bp_0, J, clam, a_lam = a_lambda, b_lam = b_lambda,
                                   lambda_0_count, lambda_0_move)
    lambda_0 <- df_lambda_0$lambda_0
    lambda_0_move <- df_lambda_0$lambda_0_move
    lambda_0_count <- df_lambda_0$lambda_0_count
    df_hist <- df_lambda_0$df_hist

    # 4. Update lambda, propose new lambda from conditional
    #   conjugate posterior [MH]
    df_lambda <- lambda_MH_cp(df_hist, df_curr, Y, I, X, s,
                               beta, beta_0, mu, sigma2, lambda, lambda_0, tau,
                               bp, bp_0, J, a_lam = a_lambda, b_lam = b_lambda, lambda_move,
                               lambda_count, alpha)
    lambda <- df_lambda$lambda
    lambda_move <- df_lambda$lambda_move
    lambda_count <- df_lambda$lambda_count
    df_curr <- df_lambda$df_curr

    # 5. Shuffle split point locations (accept/reject) [MH]
    if (J > 0) {
      swap_df <- shuffle_split_point_location(df_hist, df_curr, Y_0, I_0, X_0,
                                               lambda_0, beta_0, Y, I, X, lambda, beta, s, J, bp_0, bp,
                                               clam, maxSj)
      s <- swap_df$s
      Sigma_s <- swap_df$Sigma_s
      df_hist <- swap_df$df_hist
      df_curr <- swap_df$df_curr
    }

    # 6. Propose a birth/death of a split point via a reversible jump step [MH-Green]
    # This will update lambda_0, lambda, tau, s and J (via weighted mean) if accepted
    rjmcmc_out <- J_RJMCMC(df_hist, df_curr, Y, Y_0, I, I_0, X, X_0,
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

    if(i > warmup_iter){
      samples[j, 1:bp] <- beta

      if(bp_0 > 0) {
        samples[j, (bp + 1):(bp + bp_0)] <- beta_0
      }

      samples[j, bp + bp_0 + 1] <- J
      samples[j, bp + bp_0 + 2] <- mu
      samples[j, bp + bp_0 + 3] <- sigma2

      out_lambda[j, 1:(length(lambda))] <- lambda
      out_lambda_0[j, 1:(length(lambda_0))] <- lambda_0
      out_s[j, 1:(length(s))] <- s
      out_tau[j, 1:(length(tau))] <- tau

      #Grid of baseline hazards for shrunk estimate
      indx <- findInterval(time_grid, s, left.open = T)
      lambda_s <- lambda[indx]
      ## 1) ARRAY
      ## out_slam[j,] <- lambda_s
      writeBin(lambda_s, conn_cnd_blh)

      if(bp > 1 && !CntlOnly){
        #G-computation for marginal hazard ratio
        #draw from Dirichlet
        w_sxn_ts <- w_sxn_trt[j,]
        lhaz_up_trt <- log_haz_inroutine(trtAX,
                                          beta,
                                          lambda_s,
                                          time_grid,
                                          w_sxn_ts)
        ## 2) ARRAY
        ## lhaz_bar_trt[j,] <- lhaz_up_trt
        writeBin(lhaz_up_trt, conn_mgnl_haz1)

        w_sxn_cs <- w_sxn_con[j,]
        beta_c <- beta[-1]
        lhaz_up_con <- log_haz_inroutine(conAX,
                                          beta_c,
                                          lambda_s,
                                          time_grid,
                                          w_sxn_cs)
        ## 3) ARRAY
        ## lhaz_bar_con[j,] <- lhaz_up_con
        writeBin(lhaz_up_con, conn_mgnl_haz0)
      }

      # for
      j <- j + 1
    }

  }

  ## Close binary files
  if(bp > 1 && !CntlOnly){
      close(conn_mgnl_haz0)
      close(conn_mgnl_haz1)
      ##   out_list[["lhaz_bar_trt"]] <- lhaz_bar_trt
      ##   out_list[["lhaz_bar_con"]] <- lhaz_bar_con
  }
  close(conn_cnd_blh)
  ### END MCMC/MH loop


  ## Begin summary statistic block
  ## open binary files for reading
  if(file.exists(fnm_cnd_blh)) conn_cnd_blh <- file(fnm_cnd_blh, "rb")
  if(bp > 1 && !CntlOnly){
      if(file.exists(fnm_mgnl_haz1)) conn_mgnl_haz1 <- file(fnm_mgnl_haz1, "rb")
      if(file.exists(fnm_mgnl_haz0)) conn_mgnl_haz0 <- file(fnm_mgnl_haz0, "rb")
  }

  ## compute summary statistics
    ## matrix(...) is   max_grid x iter  matrix

  ## find TOS at 25%, 50%, 75% information time
  Y.I <- cbind(Y/stdz,I)
  Y.f <- Y.I[order(Y),]
  Y.f <- cbind(Y.f[,1],cumsum(Y.f[,2])/sum(Y.f[,2]))
  T.if <- c(sapply(0.25*(1:3), FUN=function(x, yf)yf[max(which(yf[,2]<=x)),1], yf=Y.f), max(Y)/stdz)

  ## 1. RESCALE HAZARD if standardise=TRUE
  ## 2. USE MEDIAN  in all cases

  ## case 1,2: X is just tretament OR  X is treatment + other covariates but G_check==0
  ## e.g.  G_compute==FALSE
    
  ## NOTE: in either case, when G_check==0 or 1, we can store the conditional hazard and
  ## survival functions in surv_dat. Thus this first block gets done in either case. 

  ## NOTE: in the following '0' means control group whereas above '0' means historical data
  ## so just keep that in mind
    
  ## 1st line: cumsum over time grid per row
  ## 2nd line: mean over iterations, per column
  ## 3rd line: ci over iterations, per column
  cnd_haz_0 <- matrix(readBin(conn_cnd_blh, what="double", n=max_grid*iter), max_grid, iter)*stdz
  med_cnd_haz_0 <- apply(cnd_haz_0, 1, FUN=median)
  CI_cnd_haz_0 <- t(apply(cnd_haz_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

  cnd_surv_0 <- exp(-apply(cnd_haz_0*DX(time_grid), 2, FUN=cumsum))
  med_cnd_surv_0 <- apply(cnd_surv_0, 1, FUN=median)
  CI_cnd_surv_0 <- t(apply(cnd_surv_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))
  
  cnd_haz_1 <- t(exp(samples[,nmsB[1]])*t(cnd_haz_0))
  med_cnd_haz_1 <- apply(cnd_haz_1, 1, FUN=median)
  CI_cnd_haz_1 <- t(apply(cnd_haz_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))
  
  cnd_surv_1 <- exp(-apply(cnd_haz_1*DX(time_grid), 2, FUN=cumsum))
  med_cnd_surv_1 <- apply(cnd_surv_1, 1, FUN=median)
  CI_cnd_surv_1 <- t(apply(cnd_surv_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

  t_dlt <- diff(time_grid/stdz)[1]*40

  ## smooth with ultra fine bandwidth (just larger than the grid spacing) to wipe out blips
  surv_dat <- data.frame(time = time_grid/stdz,
                           
                         med_cnd_haz_0 = med_cnd_haz_0,
                         L_cnd_haz_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_0[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_haz_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_0[,2], bandwidth=t_dlt, kernel="normal")$y,
                          
                         med_cnd_surv_0 = med_cnd_surv_0,
                         L_cnd_surv_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_0[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_surv_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_0[,2], bandwidth=t_dlt, kernel="normal")$y,
                         
                         med_cnd_haz_1 = med_cnd_haz_1,
                         L_cnd_haz_1 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_1[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_haz_1 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_1[,2], bandwidth=t_dlt, kernel="normal")$y,
                         
                         med_cnd_surv_1 = med_cnd_surv_1,
                         L_cnd_surv_1 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_1[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_surv_1 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_1[,2], bandwidth=t_dlt, kernel="normal")$y)
  MTE_dens_dfs <- NULL


  ## case 3: G_compute is TRUE so, necessarily, X is treatment + other covariantes
  ## We can add the 'mgnl' hazard and survival function to the (already existing) surv_dat data.frame
  ## We can compute the marginal treatment effect densities for information fractions 0.25*(1:4) and
  ## put them in MTE_dens_dfs
  if(bp > 1 && !CntlOnly)
  {
      log_mgnl_haz_0 <- matrix(readBin(conn_mgnl_haz0, what = "double", n = max_grid * iter), max_grid, iter) * stdz
      log_mgnl_haz_1 <- matrix(readBin(conn_mgnl_haz1, what = "double", n = max_grid * iter), max_grid, iter) * stdz

      ## marginal hazard from logged marginal hazard per sample
      mgnl_haz_0 <- exp(log_mgnl_haz_0)
      mgnl_haz_1 <- exp(log_mgnl_haz_1)

      ## lines 1-6: same as above for survival function
      med_mgnl_haz_0 <- apply(mgnl_haz_0, 1, FUN=median)
      CI_mgnl_haz_0 <- t(apply(mgnl_haz_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_surv_0 <- exp(-apply(mgnl_haz_0*DX(time_grid), 2, FUN=cumsum))
      med_mgnl_surv_0 <- apply(mgnl_surv_0, 1, FUN=median)
      CI_mgnl_surv_0 <- t(apply(mgnl_surv_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      ## lines 1-6: same as above for survival function
      med_mgnl_haz_1 <- apply(mgnl_haz_1, 1, FUN=median)
      CI_mgnl_haz_1 <- t(apply(mgnl_haz_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_surv_1 <- exp(-apply(mgnl_haz_1*DX(time_grid), 2, FUN=cumsum))
      med_mgnl_surv_1 <- apply(mgnl_surv_1, 1, FUN=median)
      CI_mgnl_surv_1 <- t(apply(mgnl_surv_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_trt_eff <- log_mgnl_haz_1 - log_mgnl_haz_0
      med_mgnl_trt_eff <- apply(mgnl_trt_eff, 1, FUN=median)
      CI_mgnl_trt_eff <- t(apply(mgnl_trt_eff, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))      

      surv_dat$med_mgnl_haz_0 <-  med_mgnl_haz_0
      surv_dat$L_mgnl_haz_0 <-  ksmooth(x = time_grid/stdz, y=CI_mgnl_haz_0[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_haz_0 <-  ksmooth(x = time_grid/stdz, y=CI_mgnl_haz_0[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_surv_0 <- med_mgnl_surv_0
      surv_dat$L_mgnl_surv_0 <- ksmooth(x = time_grid/stdz, y=CI_mgnl_surv_0[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_surv_0 <- ksmooth(x = time_grid/stdz, y=CI_mgnl_surv_0[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_haz_1 <-  med_mgnl_haz_1
      surv_dat$L_mgnl_haz_1 <-  ksmooth(x = time_grid/stdz, y=CI_mgnl_haz_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_haz_1 <-  ksmooth(x = time_grid/stdz, y=CI_mgnl_haz_1[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_surv_1 <- med_mgnl_surv_1
      surv_dat$L_mgnl_surv_1 <- ksmooth(x = time_grid/stdz, y=CI_mgnl_surv_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_surv_1 <- ksmooth(x = time_grid/stdz, y=CI_mgnl_surv_1[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_trt_eff <- med_mgnl_trt_eff
      surv_dat$L_mgnl_trt_eff <- ksmooth(x = time_grid/stdz, y=CI_mgnl_trt_eff[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_trt_eff <- ksmooth(x = time_grid/stdz, y=CI_mgnl_trt_eff[,2], bandwidth=t_dlt, kernel="normal")$y
      
      MTE_dens_dfs <- list()
      for (k in 1:4)
      {
          T.idx <- T.if[k]
          idx <- sum(time_grid/stdz <= T.idx)
          B <- mgnl_haz_1[idx,] - mgnl_haz_0[idx,]
          dens <- density(B)
          MTE_dens_dfs[[k]] <- data.frame(x=dens$x, y=dens$y)
      }
  }

  ## Close binary files
  if(bp > 1 && !CntlOnly){
      close(conn_mgnl_haz1)
      close(conn_mgnl_haz0)
  }
  close(conn_cnd_blh)
  ## End summary statistics block

  ## reorganize some of the output
  ## these are all samples of the state variables
  ## prior to snapping the hazards to the base grid
  ## the plot method can plot a trace on any given column of
  ## samples

  ## output QC
  ## drop any columns which contain all NA
  ## -- these will be at the end
  ## then out_s is the correct dimension
  ## name them s_0, ..., s_nSm1  
  s_nna <- apply(out_s, 2, \(x)sum(is.na(x)))
  out_s <- out_s[,(s_nna < iter)]
  nSm1 <- ncol(out_s) - 1
  
  ## drop any columns which contain all NA
  ## -- these will be at the end
  ## then out_lambda is the correct dimension
  ## name them lambda_1, ..., s_nL
  l_nna <- apply(out_lambda, 2, \(x)sum(is.na(x)))
  out_lambda <- out_lambda[,(l_nna < iter)]
  nL <- ncol(out_lambda)

  ## drop any columns which contain all NA
  ## -- these will be at the end
  ## then out_lambda_0 is the correct dimension
  ## name them lambda_0_1, ..., lambda_0_nL0
  l0_nna <- apply(out_lambda_0, 2, \(x)sum(is.na(x)))
  out_lambda_0 <- out_lambda_0[,(l0_nna < iter)]
  nL0 <- ncol(out_lambda_0)
    
  s.l.l0 <- as.data.frame(cbind(out_s, out_lambda, out_lambda_0))
  names(s.l.l0) <- c("s_" %,% (0:nSm1), "lambda_" %,% (1:nL), "lambda_0_" %,% (1:nL0))
  samples <- as.data.frame(samples)
  samples <- cbind(samples, s.l.l0)
  ## create output list

  if (type %in% c("uni", "mix")) {
      out_list <- list(
          Y = Y,
          I = I, 
          X = X,
          Y_0 = Y_0,
          I_0 = I_0,
          X_0 = X_0, 
          samples = samples,
          lambda = out_lambda,
          lambda_0 = out_lambda_0,
          s = out_s,
          tau = out_tau,
          lambda_0_move = lambda_0_move,
          lambda_move = lambda_move,
          beta_move = beta_count,
          surv_dat = surv_dat,
          MTE_dens_dfs = MTE_dens_dfs,
          stdz = stdz,
          T.if = T.if, 
          fnms=c(fnm_cnd_blh, fnm_mgnl_haz0, fnm_mgnl_haz1))
  } else {
      samples <- cbind(samples, out_tau)
      samples <- dplyr::rename_all(samples, dplyr::recode, out_tau = "tau")
      out_list <- list(
          Y = Y,
          I = I, 
          X = X,
          Y_0 = Y_0,
          I_0 = I_0,
          X_0 = X_0, 
          samples = samples,
          lambda = out_lambda,
          lambda_0 = out_lambda_0,
          s = out_s,
          lambda_0_move = lambda_0_move,
          lambda_move = lambda_move,
          beta_move = beta_count,
          surv_dat = surv_dat,
          MTE_dens_dfs = MTE_dens_dfs,
          stdz = stdz,
          T.if=T.if, 
          fnms = c(fnm_cnd_blh, fnm_mgnl_haz0, fnm_mgnl_haz1))
  }

  if(!is.null(beta_0)) {
    out_list$beta_0_move <- beta_0_count
  }
    
  class(out_list) <- c("BayesFBHborrow", "list")
  return(out_list)
}

## GibbsMH_wborrow
## The main MCMC/MH sampler function under the without
## borrowing case
## Arguements are the data for the current data Y, I, X, 
## (time to event, event indicator covariates),
## and state variables. 
GibbsMH_noborrow <- function(Y,
                             I,
                             X,
                             tuning_parameters,
                             hyperparameters,
                             iter,
                             warmup_iter,
                             refresh,
                             max_grid,
                             CntlOnly,
                             standardise,
                             G_compute,
                             fnm_cnd_blh,
                             fnm_mgnl_haz0,
                             fnm_mgnl_haz1){
    
  # Count accepts
  lambda_count <- 0
  lambda_move <- 0

  # proposal prior
  a_lambda <- tuning_parameters$a_lambda
  b_lambda <- tuning_parameters$b_lambda

  # hyperparameters
  a_sigma <- hyperparameters$a_sigma
  b_sigma <- hyperparameters$b_sigma
  clam_smooth <- hyperparameters$clam_smooth
  phi <- hyperparameters$phi
  Jmax <- hyperparameters$Jmax

  stdz <- 1.0
  if(standardise == TRUE){
    ## standardise time
    Y.raw <- Y
    stdz <- sum(I) / sum(Y.raw)
    Y <- Y.raw * stdz
  }

  bp <- ncol(X)
  # Set initial values
  J <- phi
  sigma2 <- b_sigma / (a_sigma + 1)
  quantiles <- quantile(Y, probs = seq(0, 1, length.out = J + 2),
                        na.rm = TRUE, names = FALSE)
  s <- c(0, quantiles[1:J + 1], max(Y))

  # s need to be unique
  if(length(unique(s)) < length(s)){
    J <- 1
    s <- c(0, max(Y)/ 2, max(Y))
  }


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

  beta_count <- 0
  if(bp == 0) {
    beta <- NULL
  } else {
    df_all <- dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda, bp =  ncol(X), J = J)
    glm.mom <- glmFit(df_all)
    beta <- glm.mom$beta.mu
    beta_count <- numeric(bp)

    # prior precision matrix
    beta_prior <- hyperparameters$beta_prior
    beta_prior_prec <- diag(1/beta_prior, bp)

  }

  if (bp == 0) {
    df_all <- dataframe_fun(Y = Y, I = I, X = NULL, s = s, lambda = lambda, bp =  0, J = J)
  }

  ## Check for G-computaion
  if(bp %in% c(0,1)){
    G_check <- 0
  }else{
    G_check <- G_compute * 1
  }

  ## Tuning parameters
  cprop_beta <- tuning_parameters$cprop_beta
  pi_b <- tuning_parameters$pi_b

  ## Output array
  ## J, mu, sigma2, beta
  nmsB <- dimnames(X)[[2]]  
  samples <- data.frame(matrix(NA, nrow = iter, ncol = bp + 3))
  if(bp > 0) colnames(samples)[1:bp] <- nmsB
  colnames(samples)[bp + (1:3)] <- c("J", "mu", "sigma2" )

  out_lambda <- data.frame(matrix(NA, nrow = iter, ncol = Jmax + 1))
  out_s <- data.frame(matrix(NA, nrow = iter, ncol = Jmax + 2))

  #Max number of grid points
  time_grid <- seq(1e-8, max(Y), length.out = max_grid)

  conn_cnd_blh <- file(fnm_cnd_blh, "wb")

  ## out_slam <-  data.frame(matrix(data = NA, nrow = 0, ncol = max_grid))
  ## colnames(out_slam) <- time_grid

  # Array for marginal estimand
  if(bp > 1 && !CntlOnly){
    nx <- nrow(X)
    trtAX <- cbind(rep(1, nx), as.matrix(X[,-1]))
    conAX <- as.matrix(X[,-1])

    conn_mgnl_haz1 <- file(fnm_mgnl_haz1, "wb")
    conn_mgnl_haz0 <- file(fnm_mgnl_haz0, "wb")

    ## 2) ARRAY lhaz_bar_trt iter x max_grid(2000) marginal treat hazard(t)
    ## lhaz_bar_trt <- data.frame(matrix(NA, nrow = 0, ncol = max_grid))
    ## colnames(lhaz_bar_trt) <- time_grid
    ## 3) ARRAY out_slam iter x max_grid(2000) marginal control hazard(t)
    ## lhaz_bar_con <- data.frame(matrix(NA, nrow = 0, ncol = max_grid))
    ## colnames(lhaz_bar_con) <- time_grid
    w_sxn_trt <- rdirichlet(iter, rep(1, nx)) # iter draws s x n
    w_sxn_con <- rdirichlet(iter, rep(1, nx))
  }


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

  j <- 1 # indx after warmup_iter
  for (i in 1:(iter + warmup_iter)) {
    if(i%%refresh == 0 && refresh != 0){message(mess[i])}

    ICAR <- ICAR_calc(s, J, clam_smooth)
    Sigma_s <- ICAR$Sigma_s

    mu <- mu_update(Sigma_s, lambda, sigma2, J)
    sigma2 <- sigma2_update(mu, lambda, Sigma_s, J, a_sigma, b_sigma)

    #Map lambda and introduce indicators.
    if (bp > 0) {
      beta_all <- beta_MH_NR(df = df_all, beta = beta, beta_prior_prec = beta_prior_prec,
                              bp =  bp, cprop_beta = cprop_beta, beta_count = beta_count)
      beta <- beta_all$beta
      beta_count <- beta_all$beta_count
    }

    #lambda, lambda and adjusted data frames
    dflambda <-lambda_0_MH_cp_NoBorrow(df_all, Y, I, X, s, beta, mu,
                                        sigma2, lambda,  bp, J, clam_smooth,
                                        a_lam = a_lambda, b_lam = b_lambda, lambda_count,
                                        lambda_move)

    lambda <- dflambda$lambda_0
    df_all <- dflambda$df
    lambda_count <- dflambda$lambda_0_count
    lambda_move <- dflambda$lambda_0_move

    #shuffle s
    if (J > 0) {
      swap_df <- shuffle_split_point_location_NoBorrow(df_all, Y, I,  X,
                                                        lambda, beta, s, J,  bp, clam_smooth)
      s <- swap_df$s
      Sigma_s <- swap_df$Sigma_s
      df_all <- swap_df$df_all
    }

    #J update
    rjmcmc_out <- J_RJMCMC_NoBorrow(df_all, Y, I, X, lambda, beta,
                                     mu, sigma2, s, J, Jmax,  bp, clam_smooth,
                                     phi, pi_b)
    J <- rjmcmc_out$J
    s <- rjmcmc_out$s
    lambda <- rjmcmc_out$lambda
    Sigma_s <- rjmcmc_out$Sigma_s
    df_all <- rjmcmc_out$df_all

    #J, mu, sigma2, beta, beta
    if(i > warmup_iter){
      if(bp > 0) samples[j, 1:bp] <- beta
      samples[j, bp + 1] <- J
      samples[j, bp + 2] <- mu
      samples[j, bp + 3] <- sigma2


      out_lambda[j, 1:(length(lambda))] <- lambda
      out_s[j, 1:(length(s))] <- s

      #Grid of baseline hazards for shrunk estimate
      indx <- findInterval(time_grid, s, left.open = T)
      lambda_s <- lambda[indx]

      ## 1) ARRAY
      ## out_slam[j,] <- lambda_s
      writeBin(lambda_s, conn_cnd_blh)

      if(bp > 1 && !CntlOnly){
        #G-computation output for marginal hazard ratio
        #draw from Dirichlet
        w_sxn_ts <- w_sxn_trt[j,]
        lhaz_up_trt <- log_haz_inroutine(trtAX,
                                          beta,
                                          lambda_s,
                                          time_grid,
                                          w_sxn_ts)
        ## 2) ARRAY
        ## lhaz_bar_trt[j,] <- lhaz_up_trt
        writeBin(lhaz_up_trt, conn_mgnl_haz1)

        w_sxn_cs <- w_sxn_trt[j,]
        beta_c <- beta[-1]
        lhaz_up_con <- log_haz_inroutine(conAX,
                                          beta_c,
                                          lambda_s,
                                          time_grid,
                                          w_sxn_cs)
        ## 3) ARRAY
        ## lhaz_bar_con[j,] <- lhaz_up_con
        writeBin(lhaz_up_con, conn_mgnl_haz0)
      }

      j <- j + 1
    }

  }

  if (bp > 1 && !CntlOnly){
      close(conn_mgnl_haz0)
      close(conn_mgnl_haz1)
      ## out_list[["lhaz_bar_trt"]] <- lhaz_bar_trt
      ## out_list[["lhaz_bar_con"]] <- lhaz_bar_con
  }
  close(conn_cnd_blh)

  ## Begin summary statistic block
  ## open binary files for reading
  if(file.exists(fnm_cnd_blh)) conn_cnd_bhl <- file(fnm_cnd_blh, "rb")
  if(bp > 1 && !CntlOnly){
      if(file.exists(fnm_mgnl_haz1)) conn_mgnl_haz1 <- file(fnm_mgnl_haz1, "rb")
      if(file.exists(fnm_mgnl_haz0)) conn_mgnl_haz0 <- file(fnm_mgnl_haz0, "rb")
  }

  ## find TOS at 25%, 50%, 75% information time
  Y.I <- cbind(Y/stdz,I)
  Y.f <- Y.I[order(Y),]
  Y.f <- cbind(Y.f[,1],cumsum(Y.f[,2])/sum(Y.f[,2]))
  T.if <- c(sapply(0.25*(1:3), FUN=function(x, yf)yf[max(which(yf[,2]<=x)),1], yf=Y.f), max(Y)/stdz)


  ## case 1,2: X is just tretament OR  X is treatment + other covariates but G_check==0
  ## e.g.  G_compute==FALSE
    
  ## NOTE: in either case, when G_check==0 or 1, we can store the conditional hazard and
  ## survival functions in surv_dat. Thus this first block gets done in either case. 

  ## NOTE: in the following '0' means control group whereas above '0' means historical data
  ## so just keep that in mind
    
  ## 1st line: cumsum over time grid per row
  ## 2nd line: mean over iterations, per column
  ## 3rd line: ci over iterations, per column
  cnd_haz_0 <- matrix(readBin(conn_cnd_blh, what="double", n=max_grid*iter), max_grid, iter)*stdz
  med_cnd_haz_0 <- apply(cnd_haz_0, 1, FUN=median)
  CI_cnd_haz_0 <- t(apply(cnd_haz_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

  cnd_surv_0 <- exp(-apply(cnd_haz_0*DX(time_grid), 2, FUN=cumsum))
  med_cnd_surv_0 <- apply(cnd_surv_0, 1, FUN=median)
  CI_cnd_surv_0 <- t(apply(cnd_surv_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

  t_dlt <- diff(time_grid/stdz)[1]*40

  surv_dat <- data.frame(time = time_grid/stdz,
                          
                         med_cnd_haz_0 = med_cnd_haz_0,
                         L_cnd_haz_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_0[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_haz_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_haz_0[,2], bandwidth=t_dlt, kernel="normal")$y,

                         med_cnd_surv_0 = med_cnd_surv_0,
                         L_cnd_surv_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_0[,1], bandwidth=t_dlt, kernel="normal")$y,
                         U_cnd_surv_0 = ksmooth(x=time_grid/stdz, y=CI_cnd_surv_0[,2], bandwidth=t_dlt, kernel="normal")$y)
                         
  if(!CntlOnly)
  {
      cnd_haz_1 <- t(exp(samples[,nmsB[1]])*t(cnd_haz_0))
      med_cnd_haz_1 <- apply(cnd_haz_1, 1, FUN=median)
      CI_cnd_haz_1 <- t(apply(cnd_haz_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))
  
      cnd_surv_1 <- exp(-apply(cnd_haz_1*DX(time_grid), 2, FUN=cumsum))
      med_cnd_surv_1 <- apply(cnd_surv_1, 1, FUN=median)
      CI_cnd_surv_1 <- t(apply(cnd_surv_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      surv_dat$med_cnd_haz_1 <- med_cnd_haz_1
      surv_dat$L_cnd_haz_1 <- ksmooth(x=time_grid/stdz, y=CI_cnd_haz_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_cnd_haz_1 <- ksmooth(x=time_grid/stdz, y=CI_cnd_haz_1[,2], bandwidth=t_dlt, kernel="normal")$y

      surv_dat$med_cnd_surv_1 <- med_cnd_surv_1
      surv_dat$L_cnd_surv_1 <- ksmooth(x=time_grid/stdz, y=CI_cnd_surv_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_cnd_surv_1 <- ksmooth(x=time_grid/stdz, y=CI_cnd_surv_1[,2], bandwidth=t_dlt, kernel="normal")$y
  }
    
  MTE_dens_dfs <- NULL

  ## case 3: G_compute is TRUE so, necessarily, X is treatment + other covariantes
  ## We can add the 'mgnl' hazard and survival function to the (already existing) surv_dat data.frame
  ## We can compute the marginal treatment effect densities for information fractions 0.25*(1:4) and
  ## put them in MTE_dens_dfs
  if(bp > 1 && !CntlOnly)
  {
      log_mgnl_haz_0 <- matrix(readBin(conn_mgnl_haz0, what = "double", n = max_grid * iter), max_grid, iter) * stdz
      log_mgnl_haz_1 <- matrix(readBin(conn_mgnl_haz1, what = "double", n = max_grid * iter), max_grid, iter) * stdz

      ## marginal hazard from logged marginal hazard per sample
      mgnl_haz_0 <- exp(log_mgnl_haz_0)
      mgnl_haz_1 <- exp(log_mgnl_haz_1)

      ## lines 1-6: same as above for survival function
      med_mgnl_haz_0 <- apply(mgnl_haz_0, 1, FUN=median)
      CI_mgnl_haz_0 <- t(apply(mgnl_haz_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_surv_0 <- exp(-apply(mgnl_haz_0*DX(time_grid), 2, FUN=cumsum))
      med_mgnl_surv_0 <- apply(mgnl_surv_0, 1, FUN=median)
      CI_mgnl_surv_0 <- t(apply(mgnl_surv_0, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      ## lines 1-6: same as above for survival function
      med_mgnl_haz_1 <- apply(mgnl_haz_1, 1, FUN=median)
      CI_mgnl_haz_1 <- t(apply(mgnl_haz_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_surv_1 <- exp(-apply(mgnl_haz_1*DX(time_grid), 2, FUN=cumsum))
      med_mgnl_surv_1 <- apply(mgnl_surv_1, 1, FUN=median)
      CI_mgnl_surv_1 <- t(apply(mgnl_surv_1, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))

      mgnl_trt_eff <- log_mgnl_haz_1 - log_mgnl_haz_0
      med_mgnl_trt_eff <- apply(mgnl_trt_eff, 1, FUN=median)
      CI_mgnl_trt_eff <- t(apply(mgnl_trt_eff, 1, FUN=\(x){i <- ci(x); c(i$CI_low, i$CI_high)}))      

      surv_dat$med_mgnl_haz_0 <-  med_mgnl_haz_0
      surv_dat$L_mgnl_haz_0 <-  ksmooth(x=time_grid/stdz, y=CI_mgnl_haz_0[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_haz_0 <-  ksmooth(x=time_grid/stdz, y=CI_mgnl_haz_0[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_surv_0 <- med_mgnl_surv_0
      surv_dat$L_mgnl_surv_0 <- ksmooth(x=time_grid/stdz, y=CI_mgnl_surv_0[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_surv_0 <- ksmooth(x=time_grid/stdz, y=CI_mgnl_surv_0[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_haz_1 <-  med_mgnl_haz_1
      surv_dat$L_mgnl_haz_1 <-  ksmooth(x=time_grid/stdz, y=CI_mgnl_haz_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_haz_1 <-  ksmooth(x=time_grid/stdz, y=CI_mgnl_haz_1[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_surv_1 <- med_mgnl_surv_1
      surv_dat$L_mgnl_surv_1 <- ksmooth(x=time_grid/stdz, y=CI_mgnl_surv_1[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_surv_1 <- ksmooth(x=time_grid/stdz, y=CI_mgnl_surv_1[,2], bandwidth=t_dlt, kernel="normal")$y
      
      surv_dat$med_mgnl_trt_eff <- med_mgnl_trt_eff
      surv_dat$L_mgnl_trt_eff <- ksmooth(x=time_grid/stdz, y=CI_mgnl_trt_eff[,1], bandwidth=t_dlt, kernel="normal")$y
      surv_dat$U_mgnl_trt_eff <- ksmooth(x=time_grid/stdz, y=CI_mgnl_trt_eff[,2], bandwidth=t_dlt, kernel="normal")$y
      
      MTE_dens_dfs <- list()
      for (k in 1:4)
      {
          T.idx <- T.if[k]
          idx <- sum(time_grid/stdz <= T.idx)
          B <- mgnl_haz_1[idx,] - mgnl_haz_0[idx,]
          dens <- density(B)
          MTE_dens_dfs[[k]] <- data.frame(x=dens$x, y=dens$y)
      }
  }

  ## Close binary files
  if(bp > 1 && !CntlOnly){
      close(conn_mgnl_haz1)
      close(conn_mgnl_haz0)
  }
  close(conn_cnd_bhl)

  ## output QC

  ## drop any columns which contain all NA
  ## -- these will be at the end
  ## then out_s is the correct dimension
  ## name them s_0, ..., s_nSm1  
  s_nna <- apply(out_s, 2, \(x)sum(is.na(x)))
  out_s <- out_s[,(s_nna < iter)]
  nSm1 <- ncol(out_s) - 1
  
  ## drop any columns which contain all NA
  ## -- these will be at the end
  ## then out_lambda is the correct dimension
  ## name them lambda_1, ..., s_nL
  l_nna <- apply(out_lambda, 2, \(x)sum(is.na(x)))
  out_lambda <- out_lambda[,(l_nna < iter)]
  nL <- ncol(out_lambda)

  s.l <- as.data.frame(cbind(out_s, out_lambda))
  names(s.l) <- c("s_" %,% (0:nSm1), "lambda_" %,% (1:nL))
  samples <- as.data.frame(samples)
  samples <- cbind(samples, s.l)
  ## create output list

  out_list <- list(
    Y = Y,
    I = I, 
    X = X,
    samples = samples,
    lambda = out_lambda,
    s = out_s,
    lambda_move = lambda_move,
    beta_move = beta_count,
    surv_dat = surv_dat,
    MTE_dens_dfs = MTE_dens_dfs,
    stdz = stdz,
    T.if=T.if, 
    fnms = c(fnm_cnd_blh, fnm_mgnl_haz0, fnm_mgnl_haz1))

  class(out_list) <- c("BayesFBHborrow", "list")
  return(out_list)
}

## birth_move
## Update for the number of split points for adding a splitpoint
birth_move <- function(U, sj, s_star, sjm1, x, j) {
  lxj <- log(x[j]) - (sj - s_star) / (sj - sjm1) * log((1 - U) / U)

  lxjp1 <- log(x[j]) + (s_star - sjm1) / (sj - sjm1) * log((1 - U) / U)

  x_prop <- append(x, exp(c(lxj, lxjp1)) , after = j)[-j]

  return(x_prop)
}

## death_move
## Update for the number of split points for removing a splitpoint
death_move <- function(sjp1, sj, sjm1, x, j) {
  lxj <- ((sj- sjm1) * log(x[j-1]) + (sjp1 - sj) * log(x[j])) / (sjp1 - sjm1)

  x_prop <-  append(x, exp(lxj), after = j)[-((j-1):j)]

  return(x_prop)
}

## ltau_dprior
## non user facing
## called by J_RJMCMC
ltau_dprior <- function(tau, a_tau, b_tau, c_tau = NULL, d_tau = NULL, p_0 = NULL, type) {
  if (type == "mix") {
    ldtau <- sum(log(p_0 * invgamma::dinvgamma(tau, shape = a_tau, rate = b_tau) +
                       (1 - p_0) * invgamma::dinvgamma(tau, shape = c_tau, rate = d_tau)))
  } else if (type == "uni") {
    ldtau <- sum(invgamma::dinvgamma(tau, shape = a_tau, rate = b_tau, log = T))
  }
}

## J_RJMCMC
## non user facing
## called by GibbsMH_wborrow
J_RJMCMC <- function(df_hist, df_curr, Y, Y_0, I, I_0,
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
  Sigma_s <- ICAR_calc(s, J, clam_smooth)$Sigma_s

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
    lambda0_prop <- birth_move(U = U3, sj = sup , s_star = s_star, sjm1 = slow, x = lambda_0, j = jlow)
    lambda_prop <-  birth_move(U = U2, sj = sup , s_star = s_star, sjm1 = slow, x = lambda, j = jlow)

    # Update data.frames and calculate llikelihood ratio
    df_hist_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    df_curr_prop <- dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda_prop, bp = bp, J = J_prop)

    llike_num <- log_likelihood(df_hist_prop, beta_0) + log_likelihood(df_curr_prop, beta)
    llike_den <- log_likelihood(df_hist, beta_0) + log_likelihood(df_curr, beta)

    # Calculate lpriors
    Sigma_s_prop <- ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s

    lprior_num <- mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J + 2), Sigma_s_prop * sigma2, log = T) +
      log(s_star - slow) + log(sup - s_star) + log(2 * J + 3) + log(2 * J + 2) +
      stats::dpois(J_prop, phi, log = T)

    lprior_den <- mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) +
      stats::dpois(J, phi, log = T) + log(sup - slow) + 2 * log(s_max)

    # Adjust for scalar tau if J == 0 and non-piecewise tau
    if (type %in% c("uni", "mix")) {
      tau_prop <- birth_move(U = U4, sj = sup , s_star = s_star, sjm1 = slow, x = tau, j = jlow)
      tau_star <- tau_prop[c(jlow, jup)]
      tau_curr <- tau[jlow]

      lprior_num <- lprior_num + mvtnorm::dmvnorm(log(lambda_prop), log(lambda0_prop), diag(tau_prop), log = T) +
        ltau_dprior(tau_star, a_tau, b_tau, c_tau, d_tau, p_0, type)
      lprior_den <- lprior_den + ltau_dprior(tau_curr, a_tau, b_tau, c_tau, d_tau, p_0, type)

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
    lambda0_prop <- death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda_0, j = j)
    lambda_prop <- death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda, j = j)

    # likelihood
    df_hist_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)
    df_curr_prop <- dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda_prop, bp = bp, J = J_prop)

    llike_num <- log_likelihood(df_hist_prop, beta_0) + log_likelihood(df_curr_prop, beta)
    llike_den <- log_likelihood(df_hist, beta_0) + log_likelihood(df_curr, beta)

    # lprior calculations
    Sigma_s_prop <- ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s

    lprior_num <- stats::dpois(J_prop, phi, log = T)  +
      mvtnorm::dmvnorm(log(lambda0_prop), rep(mu, J), Sigma_s_prop * sigma2, log = T) +
      2 * log(s_max) + log(sup - slow)

    lprior_den <- stats::dpois(J, phi, log = T)  +
      mvtnorm::dmvnorm(log(lambda_0), rep(mu, J + 1), Sigma_s * sigma2, log = T) +
      log(sj - slow) + log(sup - sj) + log(2 * J + 1) + log(2 * J)

    # Account for scalar / piecewise tau
    if (type %in% c("uni", "mix")) {
      tau_prop <- death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = tau, j = j)
      tau_star <- tau_prop[j-1]
      tau_curr <- tau[(j-1):j]

      lprior_num <- lprior_num + ltau_dprior(tau_star, a_tau, b_tau, c_tau, d_tau, p_0, type)
      lprior_den <- lprior_den + mvtnorm::dmvnorm(log(lambda), log(lambda_0), diag(tau), log = T) +
      ltau_dprior(tau_curr, a_tau, b_tau, c_tau, d_tau, p_0, type)
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

  return(list(J = J, s = s, lambda = lambda, lambda_0 = lambda_0,
              tau = tau, Sigma_s = Sigma_s, df_hist= df_hist, df_curr = df_curr))
}

## J_RJMCMC_NoBorrow
## non user facing
## called by GibbsMH_noborrow
J_RJMCMC_NoBorrow <- function(df, Y_0, I_0, X_0, lambda_0, beta_0, mu, sigma2,
                              s, J, Jmax, bp_0, clam_smooth, phi, pi_b) {
  sindx <- 1:(J + 2)
  Sigma_s <- ICAR_calc(s, J, clam_smooth)$Sigma_s

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
    lambda0_prop <- birth_move(U = U2, sj = sup , s_star = s_star, sjm1 = slow, x = lambda_0, j = jlow)

    ##Likelihood
    df_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)

    llike_num <- log_likelihood(df_prop, beta_0)
    llike_den <- log_likelihood(df, beta_0)

    ##Prior
    Sigma_s_prop <- ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s

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
    lambda0_prop <- death_move(sjp1 = sup, sj = sj, sjm1 = slow, x = lambda_0, j = j)

    ##Likelihood
    df_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda0_prop, bp = bp_0, J = J_prop)

    llike_num <- log_likelihood(df_prop, beta_0)
    llike_den <- log_likelihood(df, beta_0)

    ##Prior
    Sigma_s_prop <- ICAR_calc(s_prop, J_prop, clam_smooth)$Sigma_s

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

  return(list(J = J, s = s, lambda_0 = lambda_0, Sigma_s = Sigma_s, df_all = df))
}

## llikelihood_ratio_beta
## non user facing
## called by beta_XX functions
llikelihood_ratio_beta <- function(df, beta, beta_new) {

  X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
  xdpb <- X %*% beta
  xdpb_new <- X %*% beta_new

  llikelihood_ratio <- sum((xdpb_new  - xdpb) * df$I -
                             ((df$Y - df$tstart) * df$lambda) *
                             (exp(xdpb_new) - exp(xdpb)))

  return(llikelihood_ratio)
}

## beta_MH_RW
## non user facing
beta_MH_RW <- function(df, beta, bp, cprop_beta, beta_count) {
  for (k in 1:bp) {
    beta_new <- beta
    beta_prop <- stats::rnorm(1, beta[k], cprop_beta)
    beta_new[k] <- beta_prop

    logacc <- llikelihood_ratio_beta(df, beta, beta_new)

    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      beta_count[k] <- beta_count[k] + 1
    }

  }

  return(list(beta = beta, beta_count = beta_count))

}

## beta_mom
## non user facing
beta_mom <- function(df, k, beta, bp, cprop_beta) {
  X <- as.matrix(df[, paste0("X", 1:bp)])
  xdpb <- X %*% beta
  x <- X[, k]

  D1 <- sum(df$I * x - (df$Y - df$tstart) * df$lambda * x * exp(xdpb))
  mu_prop <- beta[k] + (cprop_beta[k]**2) / 2 * D1

  return(mu_prop)
}

## lprop_density_beta
## non user facing
lprop_density_beta <- function(beta_prop, mu, cprop_beta) {
  (-1 / (2 * cprop_beta**2)) * (beta_prop - mu)**2
}

## beta_MH_MALA
## non user facing
beta_MH_MALA <- function(df, beta, bp, cprop_beta, beta_count) {
  for(k in 1:bp){
    beta_new <- beta
    mu_prop <- beta_mom(df, k, beta, bp, cprop_beta)
    beta_prop <- stats::rnorm(n = 1, mean = mu_prop, sd = cprop_beta[k])
    beta_new[k] <- beta_prop

    mu_old <- beta_mom(df, k, beta_new, bp, cprop_beta)

    log_prop_ratio <- lprop_density_beta(beta[k], mu_old, cprop_beta[k]) - lprop_density_beta(beta_prop, mu_prop, cprop_beta[k])

    ## target_ratio <- llikelihood_ratio_beta(df, beta, beta_new, beta_prior_prec)
    target_ratio <- llikelihood_ratio_beta(df, beta, beta_new)

    logacc <- target_ratio - log_prop_ratio
    if(logacc > log(stats::runif(1))) {
      beta[k] <- beta_prop
      beta_count[k] <- beta_count[k] + 1
    }
  }

  return(list(beta = beta, beta_count = beta_count))

}

## glmFit
## non user facing
glmFit <- function(df){

  lenp <- length(df[1, grepl("X", names(df))])
  lab <- paste0("X", 1:lenp)
  splits <- length(unique(df$tstart))

  if(splits == 1){
    alllab <- paste(paste(lab, collapse= "+"), "+ offset(log(Y))")
    fmla <- stats::as.formula(paste("I ~", paste(paste(lab, collapse= "+"), "+ offset(log(Y))")))
  }else{
    df$tstart <- as.factor(df$tstart)
    alllab <- paste(paste(lab, collapse= "+"), "+ offset(log(Y))")
    fmla <- stats::as.formula(paste("I ~ tstart +", paste(paste(lab, collapse= "+"), "+ offset(log(Y))")))
  }

  #fit PWE model
  fit <- stats::glm(fmla, data = df)
  ss.x <- grepl("X", names(fit$coefficients))

  beta.mu <- fit$coefficients[ss.x]
  beta.vcov <- stats::vcov(fit)[ss.x, ss.x]

  return(list(beta.mu = beta.mu, beta.vcov = beta.vcov))

}

## beta.MH.RW.glm
## non user facing
beta.MH.RW.glm <- function(df, beta, beta_count, cprop_beta){

  glm.mom <- glmFit(df)
  beta.new <- beta
  cd2 <- cprop_beta**2 / length(beta)

  beta.prop <- as.vector(mvtnorm::rmvnorm(1, mean = beta, sigma = cd2 * glm.mom$beta.vcov))

  for (k in 1:length(beta.new)) {

    beta.new[k] <- beta.prop[k]
    logacc <- llikelihood_ratio_beta(df, beta, beta.new)

    if (logacc > log(stats::runif(1))) {
      beta <- beta.new
      beta_count[k] <- beta_count[k] + 1
    }

  }

  return(list(beta = beta, beta_count = beta_count))

}

## llike_ratio_beta_proper
## non user facing
llike_ratio_beta_proper <- function(df, beta, beta_new, beta_prior_prec){

  X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
  xdpb <- X %*% beta
  xdpb_new <- X %*% beta_new

  #prior dot products
  prior_dp <- t(beta) %*% beta_prior_prec %*% beta
  prior_dp_new <- t(beta_new) %*% beta_prior_prec %*% beta_new

  llikelihood_ratio <- sum((xdpb_new  - xdpb) * df$I -
                             ((df$Y - df$tstart) * df$lambda) *
                             (exp(xdpb_new) - exp(xdpb))) +
                       (1/2) * (prior_dp - prior_dp_new)

  return(llikelihood_ratio)
}

## lprop.dens.beta.NR
## non user facing
lprop.dens.beta.NR <- function(beta.prop, mu_old, var_old){
  ldens <- (-1 / (2 * var_old)) * (beta.prop - mu_old)**2
  return(ldens)
}

## beta_mom.NR.fun
## non user facing
beta_mom.NR.fun <- function(df, k, beta, beta_prior_prec, bp, cprop_beta) {

  bp <- length(beta)
  X <- as.matrix(df[, paste0("X", 1:bp)])
  xdpb <- X %*% beta
  x <- X[, k]
  betak <- beta[k]
  beta_prior_preck <- beta_prior_prec[k,k]

  D1 <- sum(df$I * x - (df$Y - df$tstart) * df$lambda * x * exp(xdpb)) - betak * beta_prior_preck
  D2 <- - sum(df$lambda * (df$Y - df$tstart) * x**2  * exp(xdpb)) - beta_prior_preck
  mu <- beta[k] - D1 / D2
  var <- -cprop_beta**2 / D2

  return(list(D1 = D1, D2 = D2, mu = mu, var = var))

}

## beta_MH_NR
## non user facing
beta_MH_NR <- function(df, beta, beta_prior_prec, bp, cprop_beta, beta_count){

  for(k in 1:bp){

    beta.new <- beta

    # proposal given beta
    mom.prop <- beta_mom.NR.fun(df, k, beta, beta_prior_prec, bp, cprop_beta)
    beta.prop <- stats::rnorm(n = 1, mean = mom.prop$mu, sd = sqrt(mom.prop$var))
    beta.new[k] <- beta.prop

    # mean and variance for proposal given new beta
    mom.old <- beta_mom.NR.fun(df, k, beta.new, beta_prior_prec, bp, cprop_beta)

    log_prop_ratio <- lprop.dens.beta.NR(beta[k], mom.old$mu,  mom.old$var) - lprop.dens.beta.NR(beta.prop,  mom.prop$mu,  mom.prop$var)

    ltarget_ratio <- llike_ratio_beta_proper(df, beta, beta.new, beta_prior_prec)

    logacc <- ltarget_ratio - log_prop_ratio

    if(logacc > log(stats::runif(1))){
      beta[k] <- beta.prop
      beta_count[k] <- beta_count[k] + 1
    }

  }

  return(list(beta = beta, beta_count = beta_count))

}


####################################################################################
##### Vector proposals
####################################################################################

## llikelihood_ratio_vbeta
## non user facing
llikelihood_ratio_vbeta <- function(df, beta, beta_new, beta_prior_prec) {

  X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
  xdpb <- X %*% beta
  xdpb_new <- X %*% beta_new

  #prior dot products
  if(length(beta) > 1){
    prior_dp <- t(beta) %*% beta_prior_prec %*% beta
    prior_dp_new <- t(beta_new) %*% beta_prior_prec %*% beta_new
  }else{
    prior_dp <- beta^2 * beta_prior_prec
    prior_dp_new <- beta_new^2 * beta_prior_prec
  }

  llikelihood_ratio <- sum((xdpb_new  - xdpb) * df$I -
                             ((df$Y - df$tstart) * df$lambda) *
                             (exp(xdpb_new) - exp(xdpb))) +
                        (1/2) * (prior_dp - prior_dp_new)

  return(llikelihood_ratio)
}

## beta_prop_vec
## non user facing
beta_prop_vec <- function(df,
                           beta,
                           beta_prior_prec,
                           cprop_beta,
                           bp){

  X <- as.matrix(df[, paste0("X", 1:bp)])
  eta <- as.numeric(X %*% beta)
  Aij <- (df$Y - df$tstart) * df$lambda
  w <- as.vector(Aij * exp(eta))

  # D1
  if (length(beta) > 1){
    P1 <- t(beta) %*%  beta_prior_prec
  }else{
    P1 <- beta * beta_prior_prec
  }
  D1 <- as.vector(colSums(df$I * X) - colSums((w) * X) - P1)

  # D2
  K <- crossprod(X, X * w) + beta_prior_prec # precision

  # Newton step: m = beta + Hessian %*% D1
  R <- chol(K)  # K = R' R
  Delta <- backsolve(R, forwardsolve(t(R), D1))    # K^{-1} D1 without inverting
  mu <- as.numeric(beta + Delta)

  # Sample y ~ N(0, K^{-1}) via solves
  z     <- rnorm(ncol(X))
  y     <- backsolve(R, z)                         # R y = z -> y = K^{-1/2} z

  # Proposal
  beta_prop <- mu + cprop_beta * y

  return(list(mu = mu,
              beta_prop = beta_prop,
              CholR = R))
}

## lprop_den_bvec
## non user facing
lprop_den_bvec <- function(beta_prop,
                            mu,
                            R,
                            cprop_beta){

  # efficient t(a) %*% B  %*% a
  d <- beta_prop - mu
  v <- R %*% d
  qf <- sum(v^2)
  ldet_K <-  2 * sum(log(diag(R)))
  logq <- -0.5 * (ldet_K + qf / (cprop_beta^2))

  return(logq)

}

## beta_MH_NR_vec
## non user facing
beta_MH_NR_vec <- function(df, beta, beta_prior_prec,
                            bp, cprop_beta, beta_count) {

  beta_new_prop <- beta_prop_vec(df, beta, beta_prior_prec, cprop_beta, bp)
  mu_prop <- beta_new_prop$mu
  beta_new <- beta_new_prop$beta_prop

  beta_old <- beta_prop_vec(df, beta_new, beta_prior_prec, cprop_beta, bp)

  log_prop_ratio <- lprop_den_bvec(beta,
                                    beta_old$mu,
                                    beta_old$CholR,
                                    cprop_beta) -
                    lprop_den_bvec(beta_new,
                                    mu_prop,
                                    beta_new_prop$CholR,
                                    cprop_beta)

  log_target_ratio <- llikelihood_ratio_vbeta(df, beta, beta_new, beta_prior_prec)

  logacc <- log_target_ratio  - log_prop_ratio

  if(logacc > log(stats::runif(1))) {
    beta <- beta_new
    beta_count <- beta_count + 1
  }

  return(list(beta = beta, beta_count = beta_count))

}

## lambda_conj_prop
## non user facing
lambda_conj_prop <- function(df, beta, j, bp,  alam = 0.01, blam = 0.01) {

  indx <- unique(df$tstart)
  df_ss <- df[df$tstart == indx[j],  ]

  if(!is.null(beta)) {
    X <- as.matrix(df_ss[, paste0("X", 1:bp)])
    xdpb <- X %*% beta
    rate_prop <- blam + sum((df_ss$Y - df_ss$tstart) * exp(xdpb))
  }else{
    rate_prop <- blam + sum((df_ss$Y - df_ss$tstart))
  }
  shape_prop <- alam + sum(df_ss$I)

  lambda_prop <- 0
  # if sum(df_ss) = 0, this causes lambda_prop --> 0, which will cause NA in logacc
  while (lambda_prop == 0) {
  lambda_prop <- stats:: rgamma(1, shape = shape_prop, rate = rate_prop)
  }

  return(list(lambda_prop = lambda_prop, shape_prop = shape_prop, rate_prop = rate_prop))

}

## llikelihood_ratio_lambda
## non user facing
llikelihood_ratio_lambda <- function(df, df_prop, beta) {
  if(!is.null(beta)) {
    X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
    xdpb <- X %*% beta

    llikelihood_ratio <- sum((log(df_prop$lambda) - log(df$lambda)) * df$I -
                         (df$Y - df$tstart) * (df_prop$lambda - df$lambda) * exp(xdpb))

  }else{
    llikelihood_ratio <- sum((log(df_prop$lambda) - log(df$lambda)) * df$I  -
                         (df$Y - df$tstart) * (df_prop$lambda - df$lambda))
  }
  return(llikelihood_ratio)
}

## nu_sigma_update
## non user facing
nu_sigma_update <- function(j, lambda_0, mu, sigma2, W, Q, J) {
    
  nu <- mu
  if(!is.null(dim(W))) nu <- mu + (lambda_0 - rep(mu, nrow(W))) %*% W[j,]

  if(J > 0) {
    sigma2j <- sigma2 * Q[j,j]
  }else{
    sigma2j <- sigma2 * 1
  }

  return(list(nu = nu, sigma2j = sigma2j))
}

## lgamma_ratio
## non user facing
lgamma_ratio <- function(x1, x2, shape, rate) {
  (shape - 1) * log(x1) - rate * x1 - (shape - 1) * log(x2) + rate * x2
}

## lambda_0_MH_cp
## non user facing
lambda_0_MH_cp <- function(df_hist, Y_0, I_0, X_0 = NULL, s, beta_0 = NULL,
                           mu, sigma2,  lambda, lambda_0, tau, bp_0 = 0, J,
                           clam, a_lam = 0.01, b_lam = 0.01, lambda_0_count = 0,
                           lambda_0_move = 0) {
  ICAR <- ICAR_calc(s, J, clam)
  Sigma_s <- ICAR$Sigma_s
  Q <- ICAR$Q
  W <- ICAR$W

  for (j in 1:(J + 1)) {

    lambda_0_new <- lambda_0
    lambda_0_prop_all <- lambda_conj_prop(df_hist, beta = beta_0, j, bp = bp_0, alam = a_lam,
                                           blam = b_lam)
    lambda_0_prop <- lambda_0_prop_all$lambda_prop

    lambda_0_new[j] <- lambda_0_prop
    shape_prop <- lambda_0_prop_all$shape_prop
    rate_prop <- lambda_0_prop_all$rate_prop

    df_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0_new, bp = bp_0, J = J)
    
    nu_sigma <-nu_sigma_update(j, lambda_0, mu, sigma2, W, Q, J)

    llikelihood_ratio <- llikelihood_ratio_lambda(df_hist, df_prop, beta_0)

    log_prop_ratio <- lgamma_ratio(x1 = lambda_0[j], x2 = lambda_0_prop, shape = shape_prop, rate = rate_prop)

    target_num <- stats::dnorm(log(lambda_0_prop), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                  log(lambda_0_prop)

    target_den <- stats::dnorm(log(lambda_0[j]), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                  log(lambda_0[j])

    if(length(tau) > 1) {
      target_num <- target_num +
        stats::dnorm(log(lambda[j]), log(lambda_0_prop), sqrt(tau[j]), log = T)

      target_den <- target_den +
        stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau[j]), log = T)
    }else{
      target_num <- target_num +
        stats::dnorm(log(lambda[j]), log(lambda_0_prop), sqrt(tau), log = T)

      target_den <- target_den +
        stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau), log = T)
    }

    logacc <- llikelihood_ratio + target_num - target_den + log_prop_ratio

    if(logacc > log(stats::runif(1))) {
      lambda_0 <- lambda_0_new
      df_hist <- df_prop
      lambda_0_move <- lambda_0_move + 1
    }

    lambda_0_count <- lambda_0_count + 1
  }

  return(list(lambda_0 = lambda_0, df_hist = df_hist, lambda_0_count = lambda_0_count, lambda_0_move = lambda_0_move))

}

## Lambda_0_MH_cp_NoBorrow
## non user facing
lambda_0_MH_cp_NoBorrow <- function(df_hist, Y_0, I_0, X_0 = NULL, s,
                                    beta_0 = NULL, mu, sigma2,  lambda_0,
                                    bp_0 = 0, J, clam, a_lam = 0.01,
                                    b_lam = 0.01, lambda_0_count = 0,
                                    lambda_0_move = 0) {
  ICAR <- ICAR_calc(s, J, clam)
  Sigma_s <- ICAR$Sigma_s
  Q <- ICAR$Q
  W <- ICAR$W

  for (j in 1:(J + 1)) {

    lambda_0_new <- lambda_0
    lambda_0_prop_all <- lambda_conj_prop(df_hist, beta = beta_0, j, bp = bp_0, alam = a_lam,
                                           blam = b_lam)
    lambda_0_prop <- lambda_0_prop_all$lambda_prop

    lambda_0_new[j] <- lambda_0_prop
    shape_prop <- lambda_0_prop_all$shape_prop
    rate_prop <- lambda_0_prop_all$rate_prop

    df_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s, lambda = lambda_0_new, bp = bp_0, J = J)

    nu_sigma <-nu_sigma_update(j, lambda_0, mu, sigma2, W, Q, J)

    llikelihood_ratio <- llikelihood_ratio_lambda(df_hist, df_prop, beta_0)
      log_prop_ratio <- stats::dgamma(lambda_0[j], shape = shape_prop, rate = rate_prop, log = T) -
                        stats::dgamma(lambda_0_prop, shape = shape_prop, rate = rate_prop, log = T)

      target_num <- log_likelihood(df_prop, beta_0) +
                      stats::dnorm(log(lambda_0_prop), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                      log(lambda_0_prop)

      target_den <- log_likelihood(df_hist, beta_0) +
                      stats::dnorm(log(lambda_0[j]), nu_sigma$nu, sqrt(nu_sigma$sigma2j),log = T) -
                      log(lambda_0[j])

      logacc <- target_num - target_den + log_prop_ratio

    if(logacc > log(stats::runif(1))) {
      lambda_0 <- lambda_0_new
      df_hist <- df_prop
      lambda_0_move <- lambda_0_move + 1
    }

      lambda_0_count <- lambda_0_count + 1
  }

  return(list(lambda_0 = lambda_0, df_hist = df_hist, lambda_0_count = lambda_0_count, lambda_0_move = lambda_0_move))

}

## lambda_MH_cp
## non user facing
lambda_MH_cp <- function(df_hist, df_curr, Y, I, X, s, beta, beta_0 = NULL, mu, sigma2, lambda, lambda_0, tau,
                         bp, bp_0 = 0, J, a_lam = 0.01, b_lam = 0.01, lambda_move = 0,
                         lambda_count = 0, alpha = 0.3) {

  for (j in 1:(J + 1)) {

    lambda_new <- lambda
    lambda_prop_cc <- lambda_conj_prop(df_curr, beta, j, bp = bp, alam = a_lam,
                                        blam = b_lam)
    cc_shape_prop <- lambda_prop_cc$shape_prop
    cc_rate_prop <- lambda_prop_cc$rate_prop

    lambda_prop_hist <- lambda_conj_prop(df_hist, beta_0, j, bp = bp_0, alam = a_lam,
                                          blam = b_lam)
    hist_shape_prop <- lambda_prop_hist$shape_prop
    hist_rate_prop <- lambda_prop_hist$rate_prop

    shape_prop <- a_lam + cc_shape_prop + alpha * hist_shape_prop
    rate_prop <- b_lam + cc_rate_prop + alpha * hist_rate_prop

    lambda_prop <- stats:: rgamma(1, shape = shape_prop, rate = rate_prop)
    lambda_new[j] <- lambda_prop

    df_prop <- dataframe_fun(Y = Y, I = I, X = X, s = s, lambda = lambda_new, bp = bp, J = J)

    llikelihood_ratio <- llikelihood_ratio_lambda(df_curr, df_prop, beta)
    log_prop <- lgamma_ratio(x1 = lambda[j], x2 = lambda_prop, shape = shape_prop, rate = rate_prop)

    target_num <- (- log(lambda_prop))
    target_den <- (- log(lambda[j]))

    # Adjust for non piecewise tau
    if(length(tau) > 1) {
      target_num <- target_num + stats::dnorm(log(lambda_prop), log(lambda_0[j]), sqrt(tau[j]), log = T)
      target_den <- target_den + stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau[j]), log = T)
    }else{
      target_num <- target_num + stats::dnorm(log(lambda_prop), log(lambda_0[j]), sqrt(tau), log = T)
      target_den <- target_den + stats::dnorm(log(lambda[j]), log(lambda_0[j]), sqrt(tau), log = T)
    }

    logacc <- llikelihood_ratio + target_num - target_den + log_prop

    if(logacc > log(stats::runif(1))) {
      lambda <- lambda_new
      df_curr <- df_prop
      lambda_move <- lambda_move + 1
    }

    lambda_count <- lambda_count + 1

  }

  return(list(lambda = lambda, df_curr = df_curr, lambda_count = lambda_count, lambda_move = lambda_move))

}

print.BayesFBHborrow <-
function(x, ...)
{
    ## print method is unaltered by borrowing/no borrowing
    ## so 4 possibilites
    ## no trtmt gp, no covs
    ## no trtmt gp, w covs
    ## trtd+cntls, no covs
    ## trtd+cntls, w/ covs

    x_call <- x$call
    
    G_chk <- FALSE
    if(!is.null(x_call$G_compute)) G_chk <- eval(x_call$G_compute)
    
    CntlOnly <- FALSE
    if(!is.null(x_call$CntlOnly)) CntlOnly <- eval(x_call$CntlOnly)

    ind_J <- which(names(x$samples)=="J")
    nc_X <- ind_J - 1
    is_X <- nc_X > 0
    b_nms <- NULL
    if(is_X) b_nms <- names(x$samples)[1:(ind_J-1)]

    is_2arm <- is_X && !CntlOnly
    is_Xadj <- nc_X > 1
    G_chk <- G_chk && is_Xadj
    
    ## do in all cases
    n_grid <- nrow(x$surv_dat)
    ind.t.if <- sapply(x$T.if, FUN=function(x, TT)max(which(TT<=x)), TT=x$surv_dat[,1])
    S <- x$surv_dat[ind.t.if,-1]
    p_nms <- c("0.25","0.50", "0.75", "1.00")
    t_nms <- round(x$T.if, 4)

    cat(sprintf("Call:\n"))
    print(x$call)
    cat('\n')
        

    ## controls only, w covs or w/o covs
    if(!is_2arm)
    {
        a_nms <- rep("C",4)
        if(!is_Xadj)
        {
            surv_hdng <- "Survival:\n"
            surv_df <- data.frame(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")])                                     
            surv_df <- cbind(a_nms, p_nms, t_nms, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))
        }
        if(is_Xadj)
        {
            beta_reps <- x$samples[,b_nms,drop=FALSE]
            beta_ests_cis <- t(apply(beta_reps, 2, FUN=\(x){ci_x <- ci(x); c(logHR=median(x), 'exp(logHR)'=exp(median(x)), `CrIn.L`=ci_x$CI_low, `CrIn.H`=ci_x$CI_high)}))
            dimnames(beta_ests_cis)[[1]] <- b_nms
            
            surv_hdng <- "Conditional Survival:\n"
            surv_df <- data.frame(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")])                                     
            surv_df <- cbind(a_nms, p_nms, t_nms, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))

            cat(sprintf("Conditional Treatment Effect and Regression Coefficients:\n"))
            print(beta_ests_cis)
            cat('\n')
        }
    }
    ## trtd+cntls, w covs or w/o covs
    if(is_2arm)
    {
        a_nms <- c(rep("C", 4), rep("I", 4))
        p_nms2 <- rep(p_nms, 2)
        t_nms2 <- rep(t_nms, 2)
        if(!G_chk)
        {
            ## prepare table of conditional coefficient estimates
            beta_reps <- x$samples[,b_nms,drop=FALSE]
            beta_ests_cis <- t(apply(beta_reps, 2, FUN=\(x){ci_x <- ci(x); c(logHR=median(x), 'exp(logHR)'=exp(median(x)), `CrIn.L`=ci_x$CI_low, `CrIn.H`=ci_x$CI_high)}))
            dimnames(beta_ests_cis)[[1]] <- b_nms
            
            surv_hdng <- "Conditional Survival:\n"
            surv_df <- data.frame(rbind(as.matrix(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")]), 
                                        as.matrix(S[,c("med_cnd_surv_1","L_cnd_surv_1","U_cnd_surv_1")])))
            surv_df <- cbind(a_nms, p_nms2, t_nms2, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))
        }
        if(G_chk && is_Xadj)
        {
            surv_hdng <- "Marginal Survival:\n" 
            
            surv_df <- data.frame(rbind(as.matrix(S[,c("med_mgnl_surv_0","L_mgnl_surv_0","U_mgnl_surv_0")]), 
                                        as.matrix(S[,c("med_mgnl_surv_1","L_mgnl_surv_1","U_mgnl_surv_1")])))
            surv_df <- cbind(a_nms, p_nms2, t_nms2, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))
            
            mte <- cbind(S[,"med_mgnl_trt_eff"], exp(S[,"med_mgnl_trt_eff"]), S[, c("L_mgnl_trt_eff", "U_mgnl_trt_eff")])
            MTE <- cbind(p_nms, t_nms, mte)
            dimnames(MTE) <- list(1:nrow(MTE), c("InfFrac", "Time", "MTE", "exp(MTE)", "MTE_L", "MTE_U"))
        }
        
        if(!G_chk)
        {
            cat(sprintf("Conditional Treatment Effect and Regression Coefficients:\n"))
            print(beta_ests_cis)
            cat('\n')
        }
        if(G_chk && is_Xadj)
        {
            cat(sprintf("Marginal Treatment Effect at Landmark Times:\n"))
            print(MTE)
            cat('\n')
        }
    }    

    cat(sprintf(surv_hdng))
    print(surv_df, row.names=FALSE)
    
    invisible(x)
}

## summary.BayesFBHborrow
summary.BayesFBHborrow <-
function(object, ...)
{
    ## summary method is the same as print method except
    ##      - no print out
    ##      - anything that is made inside is attached to 'x'.
    ##      - 'x' is returned

    x_call <- x$call
    
    G_chk <- FALSE
    if(!is.null(x_call$G_compute)) G_chk <- eval(x_call$G_compute)
    
    CntlOnly <- FALSE
    if(!is.null(x_call$CntlOnly)) CntlOnly <- eval(x_call$CntlOnly)

    ind_J <- which(names(x$samples)=="J")
    nc_X <- ind_J - 1
    is_X <- nc_X > 0
    b_nms <- NULL
    if(is_X) b_nms <- names(x$samples)[1:(ind_J-1)]

    is_2arm <- is_X && !CntlOnly
    is_Xadj <- nc_X > 1
    G_chk <- G_chk && is_Xadj
    
    ## do in all cases
    n_grid <- nrow(x$surv_dat)
    ind.t.if <- sapply(x$T.if, FUN=function(x, TT)max(which(TT<=x)), TT=x$surv_dat[,1])
    S <- x$surv_dat[ind.t.if,-1]
    p_nms <- c("0.25","0.50", "0.75", "1.00")
    t_nms <- round(x$T.if, 4)

    ## controls only, w covs or w/o covs
    if(!is_2arm)
    {
        a_nms <- rep("C",4)
        if(!is_Xadj)
        {
            surv_hdng <- "Survival:\n"
            surv_df <- data.frame(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")])                                     
            surv_df <- cbind(a_nms, p_nms, t_nms, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))
        }
        if(is_Xadj)
        {
            beta_reps <- x$samples[,b_nms,drop=FALSE]
            beta_ests_cis <- t(apply(beta_reps, 2, FUN=\(x){ci_x <- ci(x); c(logHR=median(x), 'exp(logHR)'=exp(median(x)), `CrIn.L`=ci_x$CI_low, `CrIn.H`=ci_x$CI_high)}))
            dimnames(beta_ests_cis)[[1]] <- b_nms
            
            surv_hdng <- "Conditional Survival:\n"
            surv_df <- data.frame(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")])                                     
            surv_df <- cbind(a_nms, p_nms, t_nms, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))

            x$coefficients <- beta_ests_cis
        }
        x$surv_summary <- surv_df
    }
    ## trtd+cntls, w covs or w/o covs
    if(is_2arm)
    {
        a_nms <- c(rep("C", 4), rep("I", 4))
        p_nms2 <- rep(p_nms, 2)
        t_nms2 <- rep(t_nms, 2)
        if(!G_chk)
        {
            ## prepare table of conditional coefficient estimates
            beta_reps <- x$samples[,b_nms,drop=FALSE]
            beta_ests_cis <- t(apply(beta_reps, 2, FUN=\(x){ci_x <- ci(x); c(logHR=median(x), 'exp(logHR)'=exp(median(x)), `CrIn.L`=ci_x$CI_low, `CrIn.H`=ci_x$CI_high)}))
            dimnames(beta_ests_cis)[[1]] <- b_nms
            
            surv_hdng <- "Conditional Survival:\n"
            surv_df <- data.frame(rbind(as.matrix(S[,c("med_cnd_surv_0","L_cnd_surv_0","U_cnd_surv_0")]), 
                                        as.matrix(S[,c("med_cnd_surv_1","L_cnd_surv_1","U_cnd_surv_1")])))
            surv_df <- cbind(a_nms, p_nms2, t_nms2, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))

            x$coefficients <- beta_ests_cis
        }
        if(G_chk && is_Xadj)
        {
            surv_hdng <- "Marginal Survival:\n" 
            
            surv_df <- data.frame(rbind(as.matrix(S[,c("med_mgnl_surv_0","L_mgnl_surv_0","U_mgnl_surv_0")]), 
                                        as.matrix(S[,c("med_mgnl_surv_1","L_mgnl_surv_1","U_mgnl_surv_1")])))
            surv_df <- cbind(a_nms, p_nms2, t_nms2, surv_df)
            dimnames(surv_df) <- list(1:nrow(surv_df), c("Arm", "InfFrac", "Time", "Survival", "Surv_L", "Surv_U"))
            
            mte <- cbind(S[,"med_mgnl_trt_eff"], exp(S[,"med_mgnl_trt_eff"]), S[, c("L_mgnl_trt_eff", "U_mgnl_trt_eff")])
            MTE <- cbind(p_nms, t_nms, mte)
            dimnames(MTE) <- list(1:nrow(MTE), c("InfFrac", "Time", "MTE", "exp(MTE)", "MTE_L", "MTE_U"))

            x$coefficients <- MTE
        }
        x$surv_summary <- surv_df
    }    

    x
}

## coef.BayesFBHborrow
## needs a help file -- called via method dispatch from coef
coef.BayesFBHborrow <-
function(object, ...) {
    return(summary(object)$coefficients)
}
coefficients <- coef

## plot.BayesFBHborrow
## this needs a help file
## called via method dispatch from plot()
plot.BayesFBHborrow <-
function(x, type=c("survival", "hazard", "trace", "TrtEff"), col=NULL, InfFrac=NULL, ylim=NULL, ...)
{
    x_call <- x$call
    minmax <- \(x, lim)pmax(pmin(x,lim[2]), lim[1])
    
    if(missing(ylim)) ylim=c(-Inf, Inf)
    G_chk <- FALSE
    if(!is.null(x_call$G_compute)) G_chk <- eval(x_call$G_compute)
    
    CntlOnly <- FALSE
    if(!is.null(x_call$CntlOnly)) CntlOnly <- eval(x_call$CntlOnly)

    ind_J <- which(names(x$samples)=="J")
    nc_X <- ind_J - 1
    is_X <- nc_X > 0
    b_nms <- NULL
    if(is_X) b_nms <- names(x$samples)[1:(ind_J-1)]

    is_2arm <- is_X && !CntlOnly
    is_Xadj <- nc_X > 1
    G_chk <- G_chk && is_Xadj

    if(missing(type)) type <- "survival"

    if(type=="TrtEff" && !is_2arm)
    {
        cat(sprintf("Call:"))
        print(x$call)
        stop("There is no treatment effect in a controls only data set!")
    }
    
    if(!(type %in% c("TrtEff","trace")))
    {
        am <- c("cnd","mgnl")[1 + G_chk]
        sh <- c(`survival`="surv", `hazard`="haz")[type]
        ti <- x$surv_dat$time
        n_grid <- nrow(x$surv_dat)
         
        y0 <- x$surv_dat[["med_" %,% am %,% "_" %,% sh %,% "_0"]]
        y0.min <- x$surv_dat[["L_" %,% am %,% "_" %,% sh %,% "_0"]]
        y0.max <- x$surv_dat[["U_" %,% am %,% "_" %,% sh %,% "_0"]]
        
        if(!is_2arm)
        {
            hdng.pfx <- c("", "Conditional ")[1 + is_Xadj]
            
            DAT <- data.frame(time=ti, y=minmax(y0,ylim), y.min=minmax(y0.min,ylim), y.max=minmax(y0.max,ylim))
            
            p <- ggplot(data=DAT) + geom_ribbon(aes(x=time, ymin=y.min, ymax=y.max), fill="#A6CEE3", alpha=0.40) +
                   geom_line(aes(x=time, y=y), color="#1F78B4") + xlab("time") + ylab(hdng.pfx %,% type)
        }    
        if(is_2arm)
        {
            hdng.pfx <- c("Conditional ", "Marginal ")[1 + G_chk]
            
            y1 <- x$surv_dat[["med_" %,% am %,% "_" %,% sh %,% "_1"]]
            y1.min <- x$surv_dat[["L_" %,% am %,% "_" %,% sh %,% "_1"]]
            y1.max <- x$surv_dat[["U_" %,% am %,% "_" %,% sh %,% "_1"]]
        
        
            DAT <- data.frame(time=c(ti, ti), y=minmax(c(y0,y1),ylim), y.min=minmax(c(y0.min, y1.min),ylim),
                              y.max=minmax(c(y0.max, y1.max),ylim), group=as.factor(rep(c("C","I"), each=n_grid)))
            
            p <- ggplot(data=DAT) + geom_ribbon(aes(x=time, ymin=y.min, ymax=y.max, fill=group), alpha=0.40) +
                   scale_fill_manual(values=c("C"="#A6CEE3", "I"="#B2DF8A")) +
                   scale_color_manual(values=c("C"="#1F78B4", "I"="#33A02C")) +
                   geom_line(aes(x=time, y=y, color=group)) + xlab("time") + ylab(hdng.pfx %,% type )
        }    
    }
    if(type=="TrtEff")
    {
        if(!G_chk)
        {
            DAT <- x$samples[,b_nms[1],drop=FALSE]
            p <- ggplot(data=DAT, aes(x=.data[[b_nms[1]]])) + geom_density(fill="#A6CEE3", alpha=0.4) +
                xlab("Cdtl Trt Eff")
        } else {
            if(missing(InfFrac)){
                InfFrac <- 1
            }else{
                if(!(InfFrac %in% (0.25*(1:4))))
                    stop("'InfFrac' must be 0.25, 0.5, 0.75 or 1")
            }
            dens_df <- x$MTE_dens_dfs[[4*InfFrac]]
            
            xl <- "Mgnl Trt Eff"
            if(InfFrac < 1) xl <- xl %,% ", InfFrac=" %,% InfFrac
            p <- ggplot(dens_df, aes(x = x, y = y)) +
                geom_ribbon(aes(ymin = 0, ymax = y), fill = "#A6CEE3", alpha = 0.4) +  #
                geom_line(color = "black", linewidth = 1/2) + xlab(xl) + ylab("density") 
        }
    }
    if(type=="trace")
    {
        if(missing(col)) col <- "J"
        if(is.null(x$call$iter)){
            n.iter <- formals(BayesFBHborrow)$iter
        } else n.iter <- x$call$iter
        DAT <- data.frame(iter=1:n.iter, y=x$samples[,col])
        p <- ggplot(data=DAT) + geom_line(aes(x=iter, y=y)) + ylab(col)
    }
    p
}

Combine <-
function(p1, p2, lgnd=NULL)
{
    is.cnd.1 <- (class(p1$layers[[1]]$stat)[1] == "StatDensity")
    tmp <- (attributes(p1)$labels$y=="density")
    is.mgnl.1 <-  ifelse(length(tmp)==0, FALSE, tmp)
    
    is.cnd.2 <- (class(p2$layers[[1]]$stat)[1] == "StatDensity")
    tmp <- (attributes(p2)$labels$y=="density")
    is.mgnl.2 <- ifelse(length(tmp)==0, FALSE, tmp)
    
    is.dens.1 <- is.cnd.1 + is.mgnl.1
    is.dens.2 <- is.cnd.2 + is.mgnl.2

    is.bad <- (is.dens.1 + is.dens.2) == 1
    is.dens <- (is.dens.1 + is.dens.2) == 2
    
    if(is.bad) stop("Combine is for plots of similar type")
                    
    if(!is.dens)
    {        
        ## Extract y values from the data used to create the plots
        yvals1 <- c(ggplot_build(p1)$data[[1]]$y,
                    ggplot_build(p1)$data[[1]]$ymin,
                    ggplot_build(p1)$data[[1]]$ymax)
        
        yvals2 <- c(ggplot_build(p2)$data[[1]]$y,
                    ggplot_build(p2)$data[[1]]$ymin,
                    ggplot_build(p2)$data[[1]]$ymax)
        
        ymin <- min(c(yvals1, yvals2), na.rm = TRUE)
        ymax <- max(c(yvals1, yvals2), na.rm = TRUE) * 1.05 # Add 5% buffer
        
        p1 <- p1 + coord_cartesian(ylim = c(ymin, ymax))
        p2 <- p2 + coord_cartesian(ylim = c(ymin, ymax))

        p2 <- p2 + scale_y_continuous(position = "right")
    
        plt <-
        (p1 | p2) +  
            plot_layout(guides = "collect", axes = "collect") &  
            theme(
                plot.margin = margin(0, 1, 0, 1),
                legend.position = "bottom")
    } else
    {
        p.lst <- list(p1, p2)
        dens.lst <- list()
        rng.12 <- NULL
        is.cnd <- c(is.cnd.1, is.cnd.2)
        is.same <- sum(is.cnd) %in% c(0,2)
        for(k in 1:2)
        {
            if(is.cnd[k])
            {
                tmp <- density(attributes(p.lst[[k]])$data[[1]])
                Group <- 1*(1-is.same) + k*is.same
                dens.lst[[k]] <- data.frame(x=tmp$x, y=tmp$y, Group=Group)
            } else
            {
                Group <- 2*(1-is.same) + k*is.same               
                dens.lst[[k]] <- cbind(attributes(p.lst[[k]])$data, Group=Group)
            }
            rng.12 <- cbind(rng.12, range(dens.lst[[k]]$x))
        }
        ord <- order(rng.12[1,])
        rng.12 <- rng.12[,ord]        
        rng <- range(rng.12)

        ## quad-extend the leftmost density to the rhs
        pnts <- dens.lst[[ord[1]]]
        p0 <- as.matrix(pnts[nrow(pnts),])[-3]
        p1 <- c(rng[2], 0)
        m <- (p1[2]-p0[2])/(p1[1]-p0[1])
        qdr <- try(quadextend(p0, p1), silent=TRUE)
        if(class(qdr)[1]=="try-error") qdr <- c(0, m, p1[2] - p0[1]/m)
        .x. <- seq(from=p0[1], to=p1[1], len=100)
        dens.lst[[ord[1]]] <- rbind(dens.lst[[ord[1]]], data.frame(x=.x., y=pmax(qdr[1]*.x.^2 + qdr[2]*.x. + qdr[3],0),
                                                                   Group=dens.lst[[ord[1]]]$Group[1]))
                                    
        ## quad-extend the rightmost density to the lhs
        pnts <- dens.lst[[ord[2]]]
        p0 <- as.matrix(pnts[1,])
        p1 <- c(rng[1], 0)
        qdr <- try(quadextend(p0, p1), silent=TRUE)
        if(class(qdr)[1]=="try-error") qdr <- c(0, m, p1[2] - p0[1]/m)
        .x. <- seq(from=p0[1], to=p1[1], len=100)
        dens.lst[[ord[2]]] <- rbind(data.frame(x=.x., y=pmax(qdr[1]*.x.^2 + qdr[2]*.x. + qdr[3],0),
                                               Group=dens.lst[[ord[2]]]$Group[1]), dens.lst[[ord[2]]])

        dens.df <- rbind(dens.lst[[ord[1]]], dens.lst[[ord[2]]])
        if(missing(lgnd)) lgnd <- c("Cndl", "Mgnl")
        dens.df$Group <- lgnd[dens.df$Group]
        cols <- c("#A6CEE3", "#B2DF8A")
        names(cols) <- lgnd
        plt <- 
        ggplot(dens.df, aes(x = x, y = y, group = Group, fill = Group)) +
            geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.4) +
            geom_line(color = "black", linewidth = 0.5) +
            xlab("Treatment Effect") + ylab("Posterior Density") +
            scale_fill_manual(values = cols) +
            theme(legend.position = "bottom")
    }
    plt
}

quadextend <-
function(p0, p1)
{
    ## connect p0 to p1 with a parabola subject to y'(p1)=0
    ## a*x^2 + b*x + c
    ## a*x0^2 + b*x0 + c = y0
    ## a*x1^2 + b*x1 + c = y1
    ## 2*a*x1 + b         = 0
    x0 <- p0[1]; y0 <- p0[2]
    x1 <- p1[1]; y1 <- p1[2]
    M <- matrix(c(x0^2, x1^2, 2*x1, x0, x1, 1, 1, 1, 0), 3, 3)
    b <- c(y0, y1, 0)
    solve(M)%*%b
}

update.BayesFBHborrow <-
function(object,
         formula = NULL,
         data,
         data_hist = NULL,                           
         subset,
         na.action,
         CntlOnly = FALSE,
         model_choice = "mix",
         tuning_parameters = NULL,
         hyperparameters = NULL,
         iter = 1000,
         warmup_iter = 100,
         refresh = 0,
         verbose = FALSE,
         max_grid = 2000,
         standardise = TRUE,
         G_compute = FALSE,
         preprocess = TRUE, 
         fnm_cnd_blh = "cnd_blh",
         fnm_mgnl_haz0 = "mgnl_haz0",
         fnm_mgnl_haz1 = "mgnl_haz1",
         dbg=FALSE,
         ...){
    
    .call. <- um <- match.call()
    m <- object$call
    um <- as.list(um[-2])
    nms.um <- names(um)[-1]
    for(e in nms.um) m[[e]] <- um[[e]]
    object$call <- m

    if((length(nms.um)==1) && (nms.um[1] == "G_compute")) ans <- object
    else ans <- eval(object$call, sys.parent())
    ans
}

sigma2_update <- function(mu, lambda_0, Sigma_s, J, a_sigma, b_sigma) {
  a_post <- a_sigma + (J + 1) / 2

  one <- rep(1, J + 1)
  cp <- t(mu * one - log(lambda_0)) %*% solve(Sigma_s) %*% (mu * one - log(lambda_0))
  b_post <- b_sigma + cp / 2

  sigma2 <- invgamma::rinvgamma(1, shape = a_post, rate = b_post)
}

mu_update <- function(Sigma_s, lambda_0, sigma2, J) {
  one <- rep(1, J + 1)
  mu_num <- (t(one) %*% solve(Sigma_s) %*% log(lambda_0))
  mu_den <- t(one) %*% solve(Sigma_s) %*% one
  mu_mu <- mu_num / mu_den

  mu_var <- sigma2 / mu_den

  mu <- stats::rnorm(1, mean = mu_mu, sd = sqrt(mu_var))
}

## old.ICAR_calc
## the old one (see below)
old.ICAR_calc <- function(s, J, clam) {

  W <- matrix(rep(0,(J + 1) * (J + 1)), nrow = J + 1)
  Q <- matrix(rep(0,(J + 1) * (J + 1)), nrow = J + 1)

  interval_length <- diff(s[!(is.na(s))])

  if (J < 2) {
    if (J == 1) {

      W[1, 2] <- clam * (interval_length[1] + interval_length[2]) / (2 * interval_length[1] + interval_length[2])
      W[2, 1] <- clam * (interval_length[2]+ interval_length[1]) / (interval_length[1] + 2 * interval_length[2])
      Q[1, 1] <- 2 / (2 * interval_length[1] + interval_length[2])
      Q[2, 2] <- 2 / (interval_length[1] + 2 * interval_length[2])
      Sigma_s <- solve(diag(J + 1) - W) %*% Q

    } else {

      Sigma_s <- as.matrix(1)

    }
  } else {

    for (j in 2:J) {

      W[j, j + 1] <- clam * (interval_length[j] + interval_length[j + 1]) / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])
      W[j, j-1] <- clam * (interval_length[j] + interval_length[j-1]) / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])
      Q[j, j] <- 2 / (interval_length[j-1] + 2 * interval_length[j] + interval_length[j + 1])

    }

    Q[j + 1, j + 1] <- 2 / (interval_length[J] + 2 * interval_length[j + 1])
    Q[1, 1] <- 2 / (2 * interval_length[1] + interval_length[2])
    W[1, 2] <- clam * (interval_length[1] + interval_length[2]) / (2 * interval_length[1] + interval_length[2])
    W[j + 1, J] <- clam * (interval_length[j + 1] + interval_length[J]) / (interval_length[J]+ 2 * interval_length[j + 1])

    Sigma_s <- solve(diag(j + 1) - W) %*% Q

  }

  return(list(Sigma_s = Sigma_s, W = W, Q = Q))

}

## ICAR_calc
## non user facing
ICAR_calc <- function(s, J, clam){

    interval_length <- diff(s[!(is.na(s))])
    if(J == 0)
    {
       Sigma_s <- as.matrix(1)
       W <- 0; Q <- 0
    }
    else
    {
        il <- c(0, interval_length, 0)

        wL <- clam*(il[-((J+2):(J+3))] + il[-c(1,J+3)])/(il[-((J+2):(J+3))] + 2*il[-c(1,J+3)] + il[-(1:2)])
        wU <- clam*(il[-c(1,J+3) ] + il[-( 1:2)]  )/(il[-((J+2):(J+3))] + 2*il[-c(1,J+3)] + il[-(1:2)])

        W <- matrix(rep(0,(J + 1) * (J + 1)), nrow = J + 1)
        W[outer((1:(J+2)), 2:(J+2), "==")[-(J+2),]] <- wL[-1]
        W[outer((1:(J+2)), 1:(J+1), "==")[-1,]] <- wU[-(J+1)]

        Q <- diag(2/(il[-((J+2):(J+3))] + 2*il[-c(1,J+3)] + il[-(1:2)]))

        ## cat(sprintf("I think it bombs %g %g %g", J, dim(W)[1], dim(Q)[1]))
        Sigma_s <- qr.solve(diag(J+1)-W, Q)
        ## cat(sprintf(" here!\n"))
    }
    return(list(Sigma_s = Sigma_s, W = W, Q = Q))

}

## tau_update
## non user facing
tau_update <- function(lambda_0, lambda, J, s,
                    a_tau, b_tau, c_tau = NULL, d_tau = NULL,
                    p_0 = NULL, type) {

  sq_diff <- (log(lambda_0) - log(lambda))**2

  if (type == "mix") {

    #compute prob on the log scale
    lw_0_num <- a_tau * log(b_tau) + lgamma(a_tau + 0.5)
    lw_0_den <- (0.5 + a_tau) * log((sq_diff / 2) + b_tau) + lgamma(a_tau)
    lw_1_num <- c_tau * log(d_tau) + lgamma(c_tau + 0.5)
    lw_1_den <- (0.5 + c_tau) * log((sq_diff / 2) + d_tau) + lgamma(c_tau)

    p_0_new <- log(p_0) +  lw_0_num - lw_0_den
    p_1_new <- log(1 - p_0) +  lw_1_num - lw_1_den

    probability_mat <- cbind(p_0_new, p_1_new)

    # normalize with log sum exp trick - avoid overflow
    p_new <- apply(probability_mat, 1, normalize_prob)

    # sample mixture
    comp <- apply(p_new, 2, sample, x = 1:2, size = 1, replace = F)

    # hyperparameters
    ac  <- matrix(rep(c(0.5 + a_tau, 0.5 + c_tau), J + 1), nrow = J + 1, byrow = T)
    bd <- cbind(sq_diff / 2 + b_tau, sq_diff / 2 + d_tau)

    call <- cbind(1:(J + 1), comp)
    tau <- invgamma::rinvgamma(n = J + 1, shape = ac[call], rate = bd[call])

  }else if (type == "uni") {

    shape_tau <- 0.5 + a_tau
    rate_tau <- sq_diff / 2 + b_tau

    # placeholder
    p_new <- 1

    tau <- invgamma::rinvgamma(n = J + 1, shape = shape_tau, rate = rate_tau)


  } else {

    sq_diff_all <- sum(sq_diff)

    # compute prob on the log scale
    lw_0_num <- a_tau * log(b_tau) + lgamma(a_tau + (J + 1) / 2)
    lw_0_den <- ((J + 1) / 2 + a_tau) * log((sq_diff_all / 2) + b_tau) + lgamma(a_tau)
    lw_1_num <- c_tau * log(d_tau) + lgamma(c_tau + (J + 1) / 2)
    lw_1_den <- ((J + 1) / 2 + c_tau) * log((sq_diff_all / 2) + d_tau) + lgamma(c_tau)

    p_0_new <- log(p_0) +  lw_0_num - lw_0_den
    p_1_new <- log(1 - p_0) +  lw_1_num - lw_1_den

    probability_mat <- cbind(p_0_new, p_1_new)

    # normalize with log sum exp trick - avoid overflow
    p_new <- apply(probability_mat, 1, normalize_prob)

    # sample mixture
    mix <- sample(x = 1:2, size = 1, replace = F, prob = p_new)

    # hyperparameters
    ac  <- c((J + 1) / 2 + a_tau, (J + 1) / 2 + c_tau)
    bd <- c(sq_diff_all / 2 + b_tau, sq_diff_all / 2 + d_tau)

    tau <- invgamma::rinvgamma(n = 1, shape = ac[mix], rate =  bd[mix])

  }

  return(list(tau = tau, p_new = p_new))

}

## shuffle_split_point_location
## non user facing
shuffle_split_point_location <- function(df_hist, df_curr, Y_0, I_0, X_0, lambda_0,
                                    beta_0, Y, I, X, lambda, beta, s, J, bp_0,
                                    bp, clam_smooth, maxSj) {
  Sigma_s <- ICAR_calc(s, J, clam_smooth)$Sigma_s

  #star individual proposals
  #prop vector of proposals and current
  for (j in 1:J) {
    if (j == J) {
      s_star <- stats::runif(1, min = s[j], max = maxSj)
    } else {
      s_star <- stats::runif(1, min = s[j], max = s[j + 2])
    }
    s_prop <- s
    s_prop[j + 1] <- s_star

    # Create new data.frames for s_prop
    df_hist_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda_0, bp = bp_0, J = J)
    df_curr_prop <- dataframe_fun(Y = Y, I = I, X = X, s = s_prop, lambda = lambda, bp = bp, J = J)

    # Update ICAR
    Sigma_s_prop <- ICAR_calc(s_prop, J, clam_smooth)$Sigma_s

    # Probability of accepting
    if (j == J) {
      lprior_num <- log(maxSj - s_star) + log(s_star - s[j])
      lprior_den <- log(maxSj - s[j + 1]) + log(s[j + 1] - s[j])
    } else {
      lprior_num <- log(s[j + 2] - s_star) + log(s_star - s[j])
      lprior_den <- log(s[j + 2] - s[j + 1]) + log(s[j + 1] - s[j])
    }

    llike_num <- log_likelihood(df_hist_prop, beta_0) + log_likelihood(df_curr_prop, beta)
    llike_den <- log_likelihood(df_hist, beta_0) + log_likelihood(df_curr, beta)

    # Acceptance ratio
    logacc <- llike_num - llike_den + lprior_num - lprior_den

    if (logacc > log(stats::runif(1))) {
      Sigma_s <- Sigma_s_prop
      df_hist <- df_hist_prop
      df_curr <- df_curr_prop
      s <- s_prop
    }

  }

  return(list(s = s, Sigma_s = Sigma_s, df_hist = df_hist, df_curr = df_curr))

}

## shuffle_split_point_location_NoBorrow
## non user facing
shuffle_split_point_location_NoBorrow <- function(df, Y_0, I_0, X_0,
                     lambda_0, beta_0, s, J,
                     bp_0, clam_smooth) {

  Sigma_s <- ICAR_calc(s, J, clam_smooth)$Sigma_s

  #star individual proposals
  #prop vector of proposals and current

  for (j in 1:J) {

    s_star <- stats::runif(1, min = s[j], max = s[j + 2])
    s_prop <- s
    s_prop[j + 1] <- s_star

    ##like
    df_prop <- dataframe_fun(Y = Y_0, I = I_0, X = X_0, s = s_prop, lambda = lambda_0, bp = bp_0, J = J)

    #ICAR
    Sigma_s_prop <- ICAR_calc(s_prop, J, clam_smooth)$Sigma_s

    #Prob of accepting
    lprior_num <- log(s[j + 2] - s_star) + log(s_star - s[j])
    lprior_denom <- log(s[j + 2] - s[j + 1]) + log(s[j + 1] - s[j])

    llike_num <- log_likelihood(df_prop, beta_0)
    llike_den <- log_likelihood(df, beta_0)

    #Prob
    logacc <- llike_num - llike_den + lprior_num - lprior_denom

    if (logacc > log(stats::runif(1))) {
      Sigma_s <- Sigma_s_prop
      df <- df_prop
      s <- s_prop
    }

  }

  return(list(s = s, Sigma_s = Sigma_s, df_all = df))
}

## input_check
## non user facing
input_check <- function(Y, Y_0, X, X_0, tuning_parameters, initial_values = NULL, hyperparameters) {

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

dataframe_fun <- function(Y, I, X, s, lambda, bp, J) {
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

logsumexp <- function(x) {
  c <- max(x)
  p <- c + log(sum(exp(x - c)))
}

normalize_prob <- function(x) {
  exp(x - logsumexp(x))
}

log_likelihood <- function(df, beta) {
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

    return(list(events = events, time_exposed = time_exposed,
                num_at_risk = num_at_risk, num_cnsr = cnsr))

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

    return(list(events_c = events_c, events_trt = events_trt,
                time_c = time_exposed_c, time_trt = time_exposed_trt,
                num_at_risk_c = num_at_risk_c, num_at_risk_trt = num_at_risk_trt))
  }

}

init_lambda_hyperparameters <- function(group_data,
                                        s,
                                        w = 0.5) {
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

  return(list(shape = shape, rate = rate, t2 = t2))
}

set_tuning_parameters <- function(tuning_parameters = NULL,
                                   borrow,
                                   X,
                                   X_0 = NULL) {
  tuning_parameters_out <- tuning_parameters

  if (borrow) {
    n_beta = ncol(X)
    defaults <- list(cprop_beta = 1.1,
                     a_lambda = 0.01,
                     b_lambda = 0.01,
                     pi_b = 0.5,
                     alpha = 0.4)

    if (!is.null(X_0)) {
      defaults$cprop_beta_0 <- 1.1
    }

    for (key in names(defaults)) {
      if (!key %in% names(tuning_parameters)) {
        tuning_parameters_out[[key]] <- defaults[[key]]
      }
    }

  } else {
    n_beta = ncol(X)
    defaults <- list(cprop_beta = 1.1,
                     a_lambda = 0.01,
                     b_lambda = 0.01,
                     pi_b = 0.5)

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

set_hyperparameters <- function(hyperparameters = NULL, model_choice) {
  hyperparameters_out <- hyperparameters
  if (model_choice == "mix") {
    defaults <- list(beta_prior = 10^2,
                     beta_0_prior = 10^2,
                     a_tau = 1,
                     b_tau = 0.001,
                     c_tau = 1,
                     d_tau = 5,
                     type = "mix",
                     p_0 = 0.8,
                     a_sigma = 1,
                     b_sigma = 1,
                     phi = 3,
                     clam_smooth = 0.8,
                     Jmax = 5)

    for (key in names(defaults)) {
      if (!key %in% names(hyperparameters)) {
        hyperparameters_out[[key]] <- defaults[[key]]
      }
    }

  } else if (model_choice == "all") {
    defaults <- list(beta_prior = 10^2,
                     beta_0_prior = 10^2,
                     a_tau = 1,
                     b_tau = 0.001,
                     c_tau = 1,
                     d_tau = 5,
                     type = "all",
                     p_0 = 0.8,
                     a_sigma = 1,
                     b_sigma = 1,
                     phi = 3,
                     clam_smooth = 0.8,
                     Jmax = 5)

    for (key in names(defaults)) {
      if (!key %in% names(hyperparameters)) {
        hyperparameters_out[[key]] <- defaults[[key]]
      }
    }

  } else if (model_choice == "uni") {
    defaults <- list(beta_prior = 10^2,
                     beta_0_prior = 10^2,
                     a_tau = 1,
                     b_tau = 0.001,
                     type = "uni",
                     a_sigma = 1,
                     b_sigma = 1,
                     phi = 3,
                     clam_smooth = 0.8,
                     Jmax = 5)

  for (key in names(defaults)) {
    if (!key %in% names(hyperparameters)) {
      hyperparameters_out[[key]] <- defaults[[key]]
      }
  }

  } else if (model_choice == "no_borrow") {
    defaults <- list(beta_prior = 10^2,
                     a_sigma = 1,
                     b_sigma = 1,
                     phi = 3,
                     clam_smooth = 0.8,
                     Jmax = 5)

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

log_haz_inroutine <- function(AX,
                              beta_s,
                              lambda_s,
                              time_grid,
                              w_sxn){

  xbeta <- AX %*% beta_s
  hz_mult <- as.numeric(exp(drop(xbeta))) #numeric vector
  cw <- as.numeric(-lambda_s * time_grid)
  cumhazard_s <-  tcrossprod(hz_mult, cw)
  post_surv <- exp(cumhazard_s)

  # w_sxn is a draw from the Dirichlet
  post_surv_adj <- post_surv * w_sxn

  # Average over the covariates for the survival function
  S_bar_s <- apply(post_surv_adj, 2, sum)
  eps <- 1e-15
  S_bar_clamped <- pmin(pmax(S_bar_s, eps), 1 - eps)

  # Log hazard per sample 1 x time grid
  log_hazard_marg <-  log(-log(S_bar_clamped))    # n x G

  return(log_hazard_marg)
}

deviance_BIC <- function(df, lambda, beta) {

  if(!is.null(beta)) {
    X <- as.matrix(df[, substr(colnames(df), 1, 1) == "X"])
    xdpb <- X %*% beta

    llikelihood <- sum(log(df$lambda) * xdpb * df$I -
                  (df$Y - df$tstart) * df$lambda * exp(xdpb))
  }else{
    llikelihood <- sum(log(df$lambda) * df$I -
                   (df$Y - df$tstart) * df$lambda )
  }

  dev <- -2 * llikelihood
  return(dev)
}

mean_for_bins <- function(s_breaks, # s from RJMCMC
                          haz_vals, # baseline hazards
                          bin_spec  # approx matrix for BDIC
){

  S0 <- s_breaks[-length(s_breaks)]
  S1 <- s_breaks[-1]
  B0 <- bin_spec[-length(bin_spec)]
  B1 <- bin_spec[-1]

  # Build min-right and max-left as proper I x J matrices
  # computing [x, y) [x, max-left and  y), min right
  min_right <- outer(S1, B1, function(a, b) pmin(a, b))
  max_left  <- outer(S0, B0, function(a, b) pmax(a, b))

  # Ensure matrix type
  min_right <- matrix(min_right, nrow = length(S1), ncol = length(B1))
  max_left  <- matrix(max_left,  nrow = length(S0), ncol = length(B0))

  # Overlap matrix (I x J)
  overlap <- apply(min_right - max_left, 2, function(x) pmax(x,0))

  # Numerator and denominator
  num <- drop(matrix(haz_vals, nrow = 1) %*% overlap)
  den <- colSums(overlap)

  means <- ifelse(den > 0, num / den, NA_real_)

  return(means)
  ## plug means into matrix
}

dataframe_fun_BIC <- function(Y,
                               I,
                               X,
                               bins_DIC,
                               bp) {

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
                                  data = df_like, cut = bins_DIC)

  return(df_split)
}

BIC_dataframe_lambda <- function(df_approx,
                                  lambda,
                                  bins_DIC){

  lam_mat <- as.data.frame(cbind(bins_DIC[-length(bins_DIC)], lambda))
  colnames(lam_mat)[1] <- "tstart"

  df_all <- merge(df_approx, lam_mat, by = "tstart", all.x = T)

  return(df_all)

}

p0tp <- function(xi, b_tau, d_tau, p_0) {
  (1 + ((1-p_0) * d_tau) / (p_0 *  b_tau) * ((xi**2 + b_tau * 2) / (xi**2 + d_tau *2))**(1.5))**(-1) - 0.5
}

p0pa <- function(xi, b_tau, d_tau){

  (1 +(b_tau/d_tau) * ((xi^2 + 2 * b_tau) / (xi^2 + 2 * d_tau))^(-3/2))^(-1)

}

## xifinder
xifinder <- function(b_tau, d_tau, p_0){
##  uniroot(p0j_tp, lower = 0, upper = 2, b_tau = b_tau, d_tau = d_tau, p_0 = p_0)$root
  uniroot(p0tp, lower = 0, upper = 2, b_tau = b_tau, d_tau = d_tau, p_0 = p_0)$root
}

## read hazard function samples
read_haz_mcmc_smpls <-
function(obj)
{
  cnd_haz_0 <- mgnl_haz_0 <- mgnl_haz_1 <- NULL
  fmls <- formals(BayesFBHborrow)
  cll <- obj$call
  iter <- ifelse(!is.null(cll$iter), cll$iter, fmls$iter)
  max_grid <- ifelse(!is.null(cll$max_grid), cll$max_grid, fmls$max_grid)

  stdz <- obj$stdz

  fnm_cnd_blh <- obj$fnms[1]
  fnm_mgnl_haz0 <- obj$fnms[2]
  fnm_mgnl_haz1 <- obj$fnms[3]

  if (file.exists(fnm_cnd_blh))
  {
    conn_cnd_blh <- file(fnm_cnd_blh, "rb")
    cnd_haz_0 <- matrix(readBin(conn_cnd_blh, what = "double", n = max_grid * iter), max_grid, iter) * stdz
    close(conn_cnd_blh)
  }
  if (file.exists(fnm_mgnl_haz0))
  {
    conn_mgnl_haz0 <- file(fnm_mgnl_haz0, "rb")
    mgnl_haz_0 <- matrix(readBin(conn_mgnl_haz0, what = "double", n = max_grid * iter), max_grid, iter) * stdz
    close(conn_mgnl_haz0)
  }
  if (file.exists(fnm_mgnl_haz1))
  {
    conn_mgnl_haz1 <- file(fnm_mgnl_haz1, "rb")
    mgnl_haz_1 <- matrix(readBin(conn_mgnl_haz1, what = "double", n = max_grid * iter), max_grid, iter) * stdz
    close(conn_mgnl_haz1)
  }
  list(cnd_haz_0=cnd_haz_0, mgnl_haz_0=mgnl_haz_0, mgnl_haz_1=mgnl_haz_1)
}

## data simulation function
genBFBHBdat <-
function(n_cc_1, n_cc_0, n_hst, B_trt, B_x_cc, int_cc, B_x_hst, int_hst, X_fact_levs=NULL, shape=1, t_er, t_fin, dbg=FALSE)
{
    if(dbg)browser()
    ## number of factor variables
    p_f <- length(X_fact_levs)

    ## number of coefficients for factor variables
    bp_f <- sum(X_fact_levs) - p_f

    ## total number of coefficients
    bp <- length(B_x_cc)

    ## number of numeric variables = number of coefficients for num vars
    p_n <- bp - bp_f

    ## total number of variables
    p <- p_n + p_f
    
    DAT.X.n <- data.frame(matrix(rnorm((n_cc_0 + n_cc_1)*p_n), n_cc_0 + n_cc_1, p_n))
    DAT.X.f <- NULL
    if(p_f > 0) DAT.X.f <- data.frame(sapply(X_fact_levs, \(x, n)factor(sample(letters[1:x], n, replace=TRUE), levels=letters[1:x]), n=n_cc_0+n_cc_1))
    DAT.X <- cbind(DAT.X.n, DAT.X.f)
    names(DAT.X) <- "X_0" %,%  (1:p)
    form <- "~" %,% paste("X_0" %,% (1:p), collapse="+")
    mf <- as.call(expression(model.frame))
    mf$formula <- form
    mf$data <- as.name("DAT.X")
    mf <- eval(mf)
    Terms <- terms(mf)
    
    X <- model.matrix(Terms, data=mf)[,-1]
    if(p_n > 0)
    {
        mu_X <- colMeans(X[, 1:p_n,drop=FALSE])
        sd_X <- apply(X[, 1:p_n,drop=FALSE], 2, FUN=var)^0.5
        X[, 1:p_n] <- t((t(X[, 1:p_n,drop=FALSE]) - mu_X)/sd_X)  
    }
    arm <- c(rep(0,n_cc_0), rep(1,n_cc_1))

    TT <- rweibull(n=n_cc_0 + n_cc_1, shape=shape, scale=exp(-(arm*B_trt + X%*%B_x_cc + int_cc)/shape))
    CNS <- t_fin - t_er*runif(n_cc_0 + n_cc_1)
    TOS <- pmin(TT, CNS)
    D <- 1*(TOS==TT)
    DAT_cc <- data.frame(tte=TOS, event=D, X_trt=arm)
    DAT_cc <- cbind(DAT_cc, DAT.X)
    
    DAT.X0.n <- data.frame(matrix(rnorm(n_hst*p_n), n_hst, p_n))
    DAT.X0.f <- NULL
    if(p_f > 0) DAT.X0.f <- data.frame(sapply(X_fact_levs, \(x, n)factor(sample(letters[1:x], n, replace=TRUE), levels=letters[1:x]), n=n_hst))
    DAT.X0 <- cbind(DAT.X0.n, DAT.X0.f)
    names(DAT.X0) <- "X_0" %,%  (1:p)
    form <- "~" %,% paste("X_0" %,% (1:p), collapse="+")
    mf0 <- as.call(expression(model.frame))
    mf0$formula <- form
    mf0$data <- as.name("DAT.X0")
    mf0 <- eval(mf0)
    Terms0 <- terms(mf0)

    X0 <- model.matrix(Terms0, data=mf0)[,-1]
    if(p_n > 0) X0[, 1:p_n] <- t((t(X0[, 1:p_n,drop=FALSE]) - mu_X)/sd_X)  

    TT0 <- rweibull(n=n_hst, shape=shape, scale=exp(-(X0%*%B_x_hst + int_hst)/shape))
    CNS0 <- t_fin - t_er*runif(n_hst)
    TOS0 <- pmin(TT0, CNS0)
    D0 <- 1*(TOS0==TT0)
    DAT_hst <- data.frame(tte=TOS0, event=D0)
    DAT_hst <- cbind(DAT_hst, DAT.X0)

    names(DAT_cc)[-(1:3)] <- names(DAT_hst)[-(1:2)] <- "X_" %,% add.zeros(1:p, max(p,10))
    list(DAT_cc=DAT_cc, DAT_hst=DAT_hst)
}

add.zeros <-
function (x, B)
{
    do.one <- function(x, B) {
        delta <- floor(logb(B, 10)) - floor(logb(x, 10))
        z <- ""
        if (delta > 0)
            z <- paste(rep("0", delta), collapse = "")
        z %,% x
    }
    if (length(x) == 1)
        ans <- do.one(x, B)
    if (length(x) > 1)
        ans <- sapply(x, FUN = do.one, B = B)
    ans
}

## paste operator
"%,%" <- function(a,b)paste0(a,b)

## first element and diff of x
DX <- function(x)c(x[1],diff(x))

rand_hex_strng <-
function(m)paste(sample(c(as.character(0:9), letters[1:6]), m, replace=TRUE), collapse="")

## bullshit repellent
if (getRversion() >= "2.15.1")  utils::globalVariables(c("x", "y", "y.min", "y.max", "group", "iter"))
##  Prints package name and version on load
.onAttach <- function(libname, pkgname)
{
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    msg <- paste(pkgname, ver)
    msg <- msg %,% "\n" %,% "Please cite this package in your work; see citation(\"BayesFBHborrow\") \n" %,%
                              "or toBibtex(citation(\"BayesFBHborrow\")) \n"
    packageStartupMessage(msg)
}

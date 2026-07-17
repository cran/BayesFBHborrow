### R code from vignette source 'Using-BayesFBHborrow.Rnw'

###################################################
### code chunk number 1: sim-exmpls
###################################################
########################################################################################################
##   SET do_examples_option:  uncomment a line
##
     do_examples_option <- character(0) # don't run the examples, there are saved versions in RData sets
##   do_examples_option <- 1            # do just the simulated data examples
##   do_examples_option <- 2            # do just the GBCS with borrowing examples
##   do_examples_option <- 3            # do options 1 & 2
##   do_examples_option <- 4            # do just the GBCS no borrowing examples
##   do_examples_option <- 5            # do everything

## chunk 1

options(prompt = "R> ", continue = "   + ")

## chunk 2
library(BayesFBHborrow)
library(dplyr)
library(kableExtra)
library(survival)
library(survminer)
library(patchwork)
"%,%" <- paste0

set.seed(round(exp(1)*1e8))

## chunk 3
hyperparameters_sim <-
    list(beta_prior = 10^2,
         beta_0_prior = 10^2,
         a_tau = 1,
         b_tau = 0.001,
         c_tau = 1,
         d_tau = 5,
         p_0 = 0.8,
         a_sigma = 1,
         b_sigma = 1,
         clam_smooth = 0.8,
         phi = 3,
         Jmax = 5)

## chunk 4
tuning_parameters_sim <-
    list(cprop_beta = 1.35,
         cprop_beta_0 = 1.35,
         a_lambda = 0.01,
         b_lambda = 0.01,
         pi_b = 0.5,
         alpha = 0.4)

## chunk 5
o.pre <- objects()

n_cc_1 <- 200
n_cc_0 <- 100
n_hst <- 100
shape <- 2
B_trt <- log(0.55)
X_fact_levs <- 3
B_x_cc <- B_x_hst <- c(-0.3,0.5,0.25,-0.5)
int_cc <- int_hst <- -log(3)

o.post <- objects()
o.sv <- setdiff(o.post, o.pre)
params <- list()
for(o in o.sv) params[[o]] <- get(o)

dat_lst <- 
    genBFBHBdat(n_cc_1=n_cc_1, n_cc_0=n_cc_0, n_hst=n_hst,
                B_trt=B_trt,
                B_x_cc=B_x_cc, 
                B_x_hst=B_x_hst,
                int_cc=int_cc, int_hst=int_hst, shape=shape,
                t_er=0.5, t_fin=1.5, X_fact_levs=X_fact_levs)

DAT_cc <- dat_lst$DAT_cc
DAT_hst <- dat_lst$DAT_hst



## chunk 6
do_sim_exmpls <- FALSE
## raw <- Sys.getenv("clparam1")
## raw <- commandArgs(trailingOnly=TRUE)
raw <- do_examples_option
is_missing <- is.character(raw) && length(raw) == 0
if(!is_missing) do_sim_exmpls <- (as.numeric(raw) %in% c(1,3,5))
sim_exmpls_fnms <- suppressWarnings(system("ls -t *bfbhb-sim-exmpls.rda", intern=TRUE, ignore.stderr=TRUE))[1]
do_sim_exmpls <- do_sim_exmpls || (length(sim_exmpls_fnms)==0)

if(do_sim_exmpls)
{
    ## chunk 6
    ## msg_conn <- textConnection("captured_msgs", "w", local = TRUE)
    ## sink(msg_conn, type = "message")

    fit_sim_hst_cndl <-
        BayesFBHborrow(Surv(tte, event)~X_01+X_02+X_03,
                       data = DAT_hst, CntlOnly = TRUE, 
                       tuning_parameters = tuning_parameters_sim,
                       hyperparameters = hyperparameters_sim,
                       iter = 6000, warmup_iter = 2000, refresh = 2000,
                       verbose = TRUE)


    fit_sim_nb_cndl <-
        BayesFBHborrow(Surv(tte, event)~X_trt+X_01+X_02+X_03,
                       data = DAT_cc, 
                       tuning_parameters = tuning_parameters_sim,
                       hyperparameters = hyperparameters_sim,
                       iter = 6000, warmup_iter = 2000, refresh = 2000,
                       verbose = TRUE)

    fit_sim_wb_cndl <- BayesFBHborrow(Surv(tte, event)~X_trt+X_01+X_02+X_03,
                                      data = DAT_cc, data_hist = DAT_hst,
                                      model_choice = 'mix',
                                      tuning_parameters = tuning_parameters_sim,
                                      hyperparameters = hyperparameters_sim,
                                      iter = 6000, warmup_iter = 2000, refresh = 2000,
                                      verbose = TRUE)

    date.stamp <- format(Sys.time(), "%Y-%m-%d-%H:%M:%S")

    fnm <- "bfbhb-sim-exmpls.rda"
    save(list=c("params", "DAT_cc", "DAT_hst", "fit_sim_hst_cndl", "fit_sim_nb_cndl", "fit_sim_wb_cndl"), file=fnm)
}
if(!do_sim_exmpls)
{
    do_sim_exmpls <- do_gbcs_nb <- do_gbcs_nexmpls <- FALSE
    load(sim_exmpls_fnms[1])
}


###################################################
### code chunk number 2: gbcs-exmpls
###################################################
## chunk 1

xifinder(b_tau = 0.001, d_tau = 5, p_0 = 0.5)

hyperparameters_gbcs <- list(beta_prior = 10^2,
                             beta_0_prior = 10^2,
                             a_tau = 1,
                             b_tau = 0.001,
                             c_tau = 1,
                             d_tau = 5,
                             p_0 = 0.5,
                             a_sigma = 1,
                             b_sigma = 1,
                             clam_smooth = 0.8,
                             phi = 3,
                             Jmax = 5)

tuning_parameters_gbcs <- list(cprop_beta = 1.17, 
                               cprop_beta_0 = 1.21, 
                               a_lambda = 0.5,
                               b_lambda = 0.5,
                               pi_b = 0.5,
                               alpha = 0.4)

##     hst,trt  trt + 1 - hst_flg + 1
##     1   0    1
##     0   0    2
##     0   1    3 
data(gbcsCS, package="condSURV")
gbcs_full <- gbcsCS
gbcs_full$diagdateb <- as.Date(format(as.character(gbcs_full$diagdateb),
                                      format="%d-%m-%Y"), format="%d-%m-%Y")
gbcs_full$tamoxifen <- gbcs_full$hormone - 1    
gbcs_full$menopause <- gbcs_full$menopause - 1
gbcs_full$grade <- factor(gbcs_full$grade, levels=as.character(1:3))   

gbcs_full$hst_flg <- 1*with(gbcs_full, diagdateb < median(diagdateb))

grp_lvls <- c("Hist Cntl","Curr Cntl","Curr Trt")
gbcs_full$group <-
    factor(with(gbcs_full,
                grp_lvls[tamoxifen + (1 - hst_flg) + 1]),
           levels=grp_lvls)

gbcs_curr <- gbcs_full[gbcs_full$hst_flg==0,]
gbcs_hist <- gbcs_full[with(gbcs_full, (hst_flg==1) & (tamoxifen==0)),]

do_gbcs_nb <- FALSE
## raw <- Sys.getenv("clparam1")
raw <- do_examples_option 
is_missing <- is.character(raw) && length(raw) == 0
if(!is_missing) do_gbcs_nb <- (as.numeric(raw) %in% c(4,5))
gbcs_nb_fnms <- suppressWarnings(system("ls -t *bfbhb-gbcs-nb-exmpl.rda", intern=TRUE, ignore.stderr=TRUE))[1]
do_gbcs_nb <- do_gbcs_nb || (length(gbcs_nb_fnms)==0)

if(do_gbcs_nb)
{
    fit_gbcs_nb_cndl <- 
    BayesFBHborrow(Surv(rectime, censrec)~tamoxifen + menopause + size + grade,
                   data=gbcs_curr,
                   tuning_parameters=tuning_parameters_gbcs,
                   hyperparameters=hyperparameters_gbcs,
                   iter=6000, warmup_iter=2000, refresh=2000,
                   verbose=TRUE,
                   max_grid=2000,
                   standardise=TRUE)
    
    date.stamp <- format(Sys.time(), "%Y-%m-%d-%H:%M:%S")

    fnm <- "bfbhb-gbcs-nb-exmpl.rda"
    save(list=c("gbcs_full","gbcs_curr", "gbcs_hist", "fit_gbcs_nb_cndl"), file=fnm)
}
if(!do_gbcs_nb) 
{
    do_sim_exmpls <- do_gbcs_nb <- do_gbcs_nexmpls <- FALSE
    load(gbcs_nb_fnms[1])
}

do_gbcs_exmpls <- FALSE
## raw <- Sys.getenv("clparam1")
raw <- do_examples_option
is_missing <- is.character(raw) && length(raw) == 0
if(!is_missing) do_gbcs_exmpls <- (as.numeric(raw) %in% c(2,3,5))
gbcs_exmpls_fnms <- suppressWarnings(system("ls -t *bfbhb-gbcs-exmpls.rda", intern=TRUE, ignore.stderr=TRUE))[1]
do_gbcs_exmpls <- do_gbcs_exmpls || (length(gbcs_exmpls_fnms)==0)

if(do_gbcs_exmpls)
{
    test <- 
    BayesFBHborrow(Surv(rectime, censrec) ~ tamoxifen + menopause + size +  grade,
                   data=gbcs_curr, data_hist=gbcs_hist,
                   model_choice="mix",
                   tuning_parameters=tuning_parameters_gbcs,
                   hyperparameters=hyperparameters_gbcs,
                   iter=500, warmup_iter=100, refresh=0, 
                   max_grid=2000,
                   standardise=TRUE)
    
    fit_gbcs_wb_cndl <- 
        BayesFBHborrow(Surv(rectime, censrec)~tamoxifen + menopause + size + grade,
                       data=gbcs_curr, data_hist=gbcs_hist,
                       model_choice="mix",
                       tuning_parameters=tuning_parameters_gbcs,
                       hyperparameters=hyperparameters_gbcs,
                       iter=6000, warmup_iter=2000, refresh=2000,
                       verbose=TRUE,
                       max_grid=2000,
                       standardise=TRUE)

##   borrowing confounds it 
##    fit_gbcs_wb_mgnl <- update(fit_gbcs_wb_cndl, G_compute=TRUE)

    date.stamp <- format(Sys.time(), "%Y-%m-%d-%H:%M:%S")

    fnm <- "bfbhb-gbcs-exmpls.rda"
    save(list=c("gbcs_full","gbcs_curr", "gbcs_hist", "test", "fit_gbcs_wb_cndl"), file=fnm)
}
if(!do_gbcs_exmpls)
{
    do_sim_exmpls <- do_gbcs_nb <- do_gbcs_nexmpls <- FALSE
    load(gbcs_exmpls_fnms[1])
}


###################################################
### code chunk number 3: define-params-1 (eval = FALSE)
###################################################
## n_cc_1 <- 200
## n_cc_0 <- 100
## n_hst <- 100


###################################################
### code chunk number 4: define-params-2 (eval = FALSE)
###################################################
## shape <- 2
## B_trt <- log(0.55)


###################################################
### code chunk number 5: define-params-3 (eval = FALSE)
###################################################
## X_fact_levs <- 3
## B_x_cc <- B_x_hst <- c(-0.3,0.5,0.25,-0.5)
## int_cc <- int_hst <- -log(3)


###################################################
### code chunk number 6: genBFBHBdat (eval = FALSE)
###################################################
## dat_lst <- 
##     genBFBHBdat(n_cc_1=n_cc_1, n_cc_0=n_cc_0, n_hst=n_hst,
##                 B_trt=B_trt,
##                 B_x_cc=B_x_cc, 
##                 B_x_hst=B_x_hst,
##                 int_cc=int_cc, int_hst=int_hst, shape=shape,
##                 t_er=0.5, t_fin=1.5, X_fact_levs=X_fact_levs)
## 
## DAT_cc <- dat_lst$DAT_cc
## DAT_hst <- dat_lst$DAT_hst


###################################################
### code chunk number 7: xifinder
###################################################
xifinder(b_tau = 0.001, d_tau = 5, p_0 = 0.8)


###################################################
### code chunk number 8: show-hyperparams (eval = FALSE)
###################################################
##     hyperparameters_sim <-
##        list(beta_prior = 10^2,
##             beta_0_prior = 10^2,
##             a_tau = 1,
##             b_tau = 0.001,
##             c_tau = 1,
##             d_tau = 5,
##             p_0 = 0.8,
##             a_sigma = 1,
##             b_sigma = 1,
##             clam_smooth = 0.8,
##             phi = 3,
##             Jmax = 5)


###################################################
### code chunk number 9: tuningparam (eval = FALSE)
###################################################
## tuning_parameters_sim <-
##     list(cprop_beta = 1.35,
##          cprop_beta_0 = 1.35,
##          a_lambda = 0.01,
##          b_lambda = 0.01,
##          pi_b = 0.5,
##          alpha = 0.4)


###################################################
### code chunk number 10: fit_hst_cndl (eval = FALSE)
###################################################
## fit_sim_hst_cndl <-
##     BayesFBHborrow(formula=Surv(tte, event)~X_01+X_02+X_03,
##                    data = DAT_hst, CntlOnly = TRUE,
##                    hyperparameters = hyperparameters_sim,
##                    tuning_parameters = tuning_parameters_sim,                
##                    warmup_iter = 2000, iter = 6000,
##                    refresh = 2000, verbose = TRUE)


###################################################
### code chunk number 11: fit_sim_nb-no-borrow (eval = FALSE)
###################################################
## fit_sim_nb_cndl <-
##     BayesFBHborrow(formula=Surv(tte, event)~X_trt+X_01+X_02+X_03,
##                    data = DAT_cc,
##                    hyperparameters = hyperparameters_sim,
##                    tuning_parameters = tuning_parameters_sim,                
##                    warmup_iter = 2000, iter = 6000,
##                    refresh = 2000, verbose = TRUE)


###################################################
### code chunk number 12: fit_sim_nb-no-borrow (eval = FALSE)
###################################################
## fit_sim_wb_cndl <-
##     BayesFBHborrow(formula=Surv(tte, event)~X_trt+X_01+X_02+X_03,
##                    data = DAT_cc, data_hist = DAT_hst,
##                    model_choice="mix",
##                    hyperparameters = hyperparameters_sim,
##                    tuning_parameters = tuning_parameters_sim,                
##                    warmup_iter = 2000, iter = 6000,
##                    refresh = 2000, verbose = TRUE)


###################################################
### code chunk number 13: prt_fit_sim_wb_cndl_coef
###################################################
trteff_sim_wb_cndl <- coef(fit_sim_wb_cndl)[1,1]
fit_sim_wb_cndl_coef_kbl <-
    kable(coef(fit_sim_wb_cndl), digits=4, format="latex", booktabs=TRUE, linesep="", table.envir=NULL)
fit_sim_wb_cndl_coef_kbl


###################################################
### code chunk number 14: mk_sim_cndl_surv_tbl
###################################################
fit_sim_wb_cndl_surv_kbl <-
    kable(summary(fit_sim_wb_cndl)$surv_summary, digits=4, format="latex", booktabs=TRUE, linesep="", table.envir=NULL)


###################################################
### code chunk number 15: prt_sim_cdnl_surv_tbl
###################################################
  fit_sim_wb_cndl_surv_kbl <-
    kable(summary(fit_sim_wb_cndl)$surv_summary, digits=4, format="latex", booktabs=TRUE, linesep="", table.envir=NULL)
  fit_sim_wb_cndl_surv_kbl


###################################################
### code chunk number 16: mgnl_estmd
###################################################
mgnl_cntrst_sim_wb <- update(fit_sim_wb_cndl, G_compute=TRUE)


###################################################
### code chunk number 17: mgnl_estmd
###################################################
mgnl_cntrst_sim_wb_txt <- coef(mgnl_cntrst_sim_wb)[4,3]
mgnl_cntrst_sim_wb_kbl <- kable(coef(mgnl_cntrst_sim_wb), digits=4, format="latex", booktabs=TRUE, 
                                linesep="", table.envir=NULL)
mgnl_cntrst_sim_wb_kbl


###################################################
### code chunk number 18: show-read-call (eval = FALSE)
###################################################
## surv_dat_full <- read_haz_mcmc_smpls(fit_sim_wb_cndl)


###################################################
### code chunk number 19: mk-haz-plots-1
###################################################
p_cndl_haz_hst <- plot(fit_sim_hst_cndl, type="hazard")


###################################################
### code chunk number 20: mk-haz-plots-2
###################################################
p_cndl_haz_nb <- plot(fit_sim_nb_cndl, type="hazard")
p_cndl_haz_wb <- plot(fit_sim_wb_cndl, type="hazard")


###################################################
### code chunk number 21: mk-haz-plots-3
###################################################
p_cndl_haz_cmb <- Combine(p_cndl_haz_nb, p_cndl_haz_wb)


###################################################
### code chunk number 22: plot_cndl_haz_nb_wb
###################################################
p_cndl_haz_cmb


###################################################
### code chunk number 23: gbcs_recode (eval = FALSE)
###################################################
## data(gbcsCS, package="condSURV")
## gbcs_full <- gbcsCS
## gbcs_full$diagdateb <- as.Date(format(as.character(gbcs_full$diagdateb),
##                                       format="%d-%m-%Y"), format="%d-%m-%Y")
## gbcs_full$tamoxifen <- gbcs_full$hormone - 1
## gbcs_full$menopause <- gbcs_full$menopause - 1
## gbcs_full$grade <- factor(gbcs_full$grade, levels=as.character(1:3))
## 
## gbcs_full$hst_flg <- 1*with(gbcs_full, diagdateb < median(diagdateb))
## 
## grp_lvls <- c("Hist Cntl","Curr Cntl","Curr Trt")
## gbcs_full$group <-
##   factor(with(gbcs_full,
##               grp_lvls[tamoxifen + (1 - hst_flg) + 1]),
##          levels=grp_lvls)
## 
## gbcs_curr <- gbcs_full[gbcs_full$hst_flg==0,]
## gbcs_hist <- gbcs_full[with(gbcs_full, (hst_flg==1) & (tamoxifen==0)),]


###################################################
### code chunk number 24: gbcs-show
###################################################
gbcs_kbl <- kable(head(gbcs_curr[,c("id","diagdateb","rectime","censrec","tamoxifen","menopause","estrg_recp","size","grade")]),
                  digits=4, format="latex", booktabs=TRUE, linesep="", row.names=FALSE, table.envir=NULL)
gbcs_kbl


###################################################
### code chunk number 25: gbcs_KM_mk
###################################################
grp_lvls <- c("Hist Cntl","Curr Cntl","Curr Trt")
fit_sf <- survfit(Surv(rectime, censrec)~group, data=gbcs_full)
KM_dat <- data.frame(time=fit_sf$time, surv=fit_sf$surv)
n_strat <- fit_sf$strata
KM_dat$group <-
    c(rep(grp_lvls[1], n_strat[1]), rep(grp_lvls[2], n_strat[2]),
      rep(grp_lvls[3], n_strat[3]))
p_KM <- ggplot(data=KM_dat) + geom_step(aes(time, surv, group=group, color=group))


###################################################
### code chunk number 26: ggcs_KM_do
###################################################
p_KM


###################################################
### code chunk number 27: xi-gbcs
###################################################
xifinder(b_tau = 0.001, d_tau = 5, p_0 = 0.5)


###################################################
### code chunk number 28: hyperpars (eval = FALSE)
###################################################
## hyperparameters_gbcs <- list(beta_prior = 10^2,
##                              beta_0_prior = 10^2,
##                              a_tau = 1,
##                              b_tau = 0.001,
##                              c_tau = 1,
##                              d_tau = 5,
##                              p_0 = 0.5,
##                              a_sigma = 1,
##                              b_sigma = 1,
##                              clam_smooth = 0.8,
##                              phi = 3,
##                              Jmax = 5)


###################################################
### code chunk number 29: tune (eval = FALSE)
###################################################
## tuning_parameters_gbcs <- list(cprop_beta = 1.17,
##                                cprop_beta_0 = 1.21,
##                                a_lambda = 0.5,
##                                b_lambda = 0.5,
##                                pi_b = 0.5,
##                                alpha = 0.4)


###################################################
### code chunk number 30: test (eval = FALSE)
###################################################
## test <- 
## BayesFBHborrow(Surv(rectime, censrec) ~ tamoxifen + menopause + size + grade,
##                data=gbcs_curr, data_hist=gbcs_hist,
##                model_choice="mix",
##                tuning_parameters=tuning_parameters_gbcs,
##                hyperparameters=hyperparameters_gbcs,
##                iter=500, warmup_iter=100, refresh=0)


###################################################
### code chunk number 31: gbcs_fit_nb (eval = FALSE)
###################################################
## fit_gbcs_nb_cndl <- 
##     BayesFBHborrow(Surv(rectime, censrec)~tamoxifen + menopause + size + grade,
##                    data=gbcs_curr,
##                    tuning_parameters=tuning_parameters_gbcs,
##                    hyperparameters=hyperparameters_gbcs,
##                    iter=6000, warmup_iter=2000, refresh=2000,
##                    verbose=TRUE)


###################################################
### code chunk number 32: gbcs_fit_nb (eval = FALSE)
###################################################
## fit_gbcs_wb_cndl <- 
##   BayesFBHborrow(Surv(rectime, censrec)~ tamoxifen + menopause + size + grade,
##                  data=gbcs_curr, data_hist=gbcs_hist,
##                  model_choice="mix",
##                  tuning_parameters=tuning_parameters_gbcs,
##                  hyperparameters=hyperparameters_gbcs,
##                  iter=6000, warmup_iter=2000, refresh=2000,
##                  verbose=TRUE,
##                  max_grid=2000,
##                  standardise=TRUE)
## 
##     


###################################################
### code chunk number 33: trace
###################################################
p_gbcs_wb_trace <- plot(fit_gbcs_wb_cndl, type="trace", col="tamoxifen")


###################################################
### code chunk number 34: fit_gbcs_wb_cndl_coef
###################################################
trteff_gbcs_wb_cndl_all <- coef(fit_gbcs_wb_cndl)[1,]
trteff_gbcs_wb_cndl <- trteff_gbcs_wb_cndl_all[1]
trteff_gbcs_wb_cndl_L <- trteff_gbcs_wb_cndl_all[3]
trteff_gbcs_wb_cndl_U <- trteff_gbcs_wb_cndl_all[4]

fit_gbcs_wb_cndl_coef_kbl <-
    kable(coef(fit_gbcs_wb_cndl), digits=4, format="latex", booktabs=TRUE, linesep="", table.envir=NULL)
fit_gbcs_wb_cndl_coef_kbl


###################################################
### code chunk number 35: mk-plots
###################################################
p_gbcs_cnd_nb_haz <- plot(fit_gbcs_nb_cndl, type="hazard", ylim=c(0,0.0010))
p_gbcs_cnd_wb_haz <- plot(fit_gbcs_wb_cndl, type="hazard", ylim=c(0,0.0010))
p_gbcs_cnd_haz_cmb <- Combine(p_gbcs_cnd_nb_haz, p_gbcs_cnd_wb_haz)

p_gbcs_cnd_nb_TrtEff <- plot(fit_gbcs_nb_cndl, type="TrtEff")
p_gbcs_cnd_wb_TrtEff <- plot(fit_gbcs_wb_cndl, type="TrtEff")
p_gbcs_cnd_TrtEff_cmb <- Combine(p_gbcs_cnd_nb_TrtEff, p_gbcs_cnd_wb_TrtEff,
                                 lgnd=c("NB","WB"))


###################################################
### code chunk number 36: p_gbcs_cnd_haz_plot
###################################################
p_gbcs_cnd_haz_cmb


###################################################
### code chunk number 37: p_gbcs_cnd_TrtEff_plot
###################################################
p_gbcs_cnd_TrtEff_cmb


###################################################
### code chunk number 38: plot_cndl_haz_hst
###################################################
p_cndl_haz_hst


###################################################
### code chunk number 39: trace
###################################################
p_gbcs_wb_trace

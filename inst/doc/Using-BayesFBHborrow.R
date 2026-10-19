### R code from vignette source 'Using-BayesFBHborrow.Rnw'

###################################################
### code chunk number 1: production
###################################################
## Here we specify whether this is a production run or validation run. Production runs use 6000 sampler iterations and
## 2000 warm-up iterations, whereas validation runs use 100 and 50, respectively. A production run takes about an hour
## and saves the results whereas a validation run checks for previously saved results to load and is very quick.
##
## run_type is one of "test", "load saved", or "production"
##
##   where:
##
##    test: runs the sampler for only a bit and doesn't save or load anything
##          just using the minimal runs everywhere
##
##    load saved: does a test as above on a minimal set of sampler runs
##                but then loads all the required objects from a previously
##                saved run
##
##    production: runs the sampler for a good while, using these results
##                everywhere in the paper, and saves the results
##
run_type <- "load saved"
################################################################################################


test <- run_type == "test"
load_saved <- run_type == "load saved"
production <- (run_type == "production")

if(production) {
  iter <- 6000
  warmup_iter <- 2000
  iter_tst <- 100
  warmup_iter_tst <- 50
} else {
  iter <- iter_tst <- 10
  warmup_iter <- warmup_iter_tst <- 5
}


###################################################
### code chunk number 2: options
###################################################
options(prompt = "R> ", continue = "   + ")


###################################################
### code chunk number 3: load-libs
###################################################
set.seed(123456789)
library(BayesFBHborrow)
library(dplyr)
library(kableExtra)
library(survival)
library(survminer)
library(patchwork)
"%,%" <- paste0


###################################################
### code chunk number 4: define-params-1
###################################################
n_cc_1 <- 300
n_cc_0 <- 100
n_hst <- 200


###################################################
### code chunk number 5: define-params-2
###################################################
shape <- 2
B_trt <- log(0.55)


###################################################
### code chunk number 6: define-params-3
###################################################
X_fact_levs <- 3
B_x_cc <- B_x_hst <- c(-0.3,0.5,0.25,-0.5)
int_cc <- int_hst <- -log(3)


###################################################
### code chunk number 7: collect_params
###################################################
par.objs <- setdiff(objects(), "%,%")
params <- list()
for (o in par.objs) params[[o]] <- get(o)


###################################################
### code chunk number 8: genBFBHBdat
###################################################
dat_lst <- 
    genBFBHBdat(n_cc_1=n_cc_1, n_cc_0=n_cc_0, n_hst=n_hst,
                B_trt=B_trt,
                B_x_cc=B_x_cc, 
                B_x_hst=B_x_hst,
                int_cc=int_cc, int_hst=int_hst, shape=shape,
                t_er=0.5, t_fin=1.5, X_fact_levs=X_fact_levs)

DAT_cc <- dat_lst$DAT_cc
DAT_hst <- dat_lst$DAT_hst


###################################################
### code chunk number 9: xifinder
###################################################
xifinder(b_tau = 0.001, d_tau = 5, p_0 = 0.8)


###################################################
### code chunk number 10: show-hyperparams
###################################################
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


###################################################
### code chunk number 11: tuningparam
###################################################
tuning_parameters_sim <-
    list(cprop_beta = 1.35,
         cprop_beta_0 = 1.35,
         a_lambda = 0.01,
         b_lambda = 0.01,
         pi_b = 0.5,
         alpha = 0.4)


###################################################
### code chunk number 12: fit_hst_cndl
###################################################
fit_sim_hst_cndl <-
    BayesFBHborrow(formula=Surv(tte, event)~X_01+X_02+X_03,
                   data = DAT_hst, CntlOnly = TRUE,
                   hyperparameters = hyperparameters_sim,
           tuning_parameters = tuning_parameters_sim,
           warmup_iter = warmup_iter, iter = iter,
           refresh = warmup_iter, verbose = TRUE)


###################################################
### code chunk number 13: fit_sim_nb-no-borrow
###################################################
fit_sim_nb_cndl <-
    BayesFBHborrow(formula=Surv(tte, event)~X_trt+X_01+X_02+X_03,
                   data = DAT_cc,
                   hyperparameters = hyperparameters_sim,
           tuning_parameters = tuning_parameters_sim,
           warmup_iter = warmup_iter, iter = iter,
           refresh = warmup_iter, verbose = TRUE)


###################################################
### code chunk number 14: fit_sim_wb-mix
###################################################
fit_sim_wb_cndl <-
    BayesFBHborrow(formula=Surv(tte, event)~X_trt+X_01+X_02+X_03,
                   data = DAT_cc, data_hist = DAT_hst,
                   model_choice="mix",
                   hyperparameters = hyperparameters_sim,
                   tuning_parameters = tuning_parameters_sim,
                   warmup_iter = warmup_iter, iter = iter,
                   refresh = warmup_iter, verbose = TRUE)


###################################################
### code chunk number 15: save-sim-data-rslts
###################################################
if(production)
{
  save(list=c("params", "DAT_cc", "DAT_hst", "hyperparameters_sim", "tuning_parameters_sim", "fit_sim_hst_cndl",
              "fit_sim_nb_cndl", "fit_sim_wb_cndl"), file="rdata/bfbhb-sim-exmpls.rda")
}
if(load_saved)
{
  sim_exmpls_fnms <- suppressWarnings(system("ls -t rdata/*bfbhb-sim-exmpls.rda", intern=TRUE, ignore.stderr=TRUE))[1]
  if(!is.na(sim_exmpls_fnms)) load(sim_exmpls_fnms[1])
}


###################################################
### code chunk number 16: prt_fit_sim_wb_cndl_coef
###################################################
trteff_sim_wb_cndl <- coef(fit_sim_wb_cndl)["X_trt", "logHR"]
fit_sim_wb_cndl_coef_kbl <-
    kable(coef(fit_sim_wb_cndl), digits=4, format="latex", booktabs=TRUE,
          linesep="", table.envir=NULL)
fit_sim_wb_cndl_coef_kbl


###################################################
### code chunk number 17: mk_sim_cndl_surv_tbl
###################################################
fit_sim_wb_cndl_surv_kbl <-
    kable(summary(fit_sim_wb_cndl)$surv_summary, digits=4, format="latex",
          booktabs=TRUE, linesep="", table.envir=NULL)


###################################################
### code chunk number 18: prt_sim_cdnl_surv_tbl
###################################################
  fit_sim_wb_cndl_surv_kbl <-
      kable(summary(fit_sim_wb_cndl)$surv_summary, digits=4, format="latex",
            booktabs=TRUE, linesep="", table.envir=NULL)
  fit_sim_wb_cndl_surv_kbl


###################################################
### code chunk number 19: mgnl_estmd
###################################################
mgnl_cntrst_sim_wb <- update(fit_sim_wb_cndl, G_compute=TRUE)


###################################################
### code chunk number 20: mgnl_estmd
###################################################
mgnl_cntrst_sim_wb_MTE <- MTE(mgnl_cntrst_sim_wb)
mgnl_cntrst_sim_wb_value <- mgnl_cntrst_sim_wb_MTE[4, "MTE"]
mgnl_cntrst_sim_wb_kbl <- kable(MTE(mgnl_cntrst_sim_wb), digits=4,
                                format="latex", booktabs=TRUE, linesep="",
                                table.envir=NULL)
mgnl_cntrst_sim_wb_kbl


###################################################
### code chunk number 21: show-read-call (eval = FALSE)
###################################################
## surv_dat_full <- read_haz_mcmc_smpls(fit_sim_wb_cndl)


###################################################
### code chunk number 22: setup
###################################################
dir.create("figures", showWarnings = FALSE, recursive = TRUE)


###################################################
### code chunk number 23: mk-haz-plots-1
###################################################
p_cndl_haz_hst <- plot(fit_sim_hst_cndl, type="hazard")


###################################################
### code chunk number 24: mk-haz-plots-2
###################################################
p_cndl_haz_nb <- plot(fit_sim_nb_cndl, type="hazard")
p_cndl_haz_wb <- plot(fit_sim_wb_cndl, type="hazard")


###################################################
### code chunk number 25: mk-haz-plots-3
###################################################
p_cndl_haz_cmb <- Combine(p_cndl_haz_nb, p_cndl_haz_wb)


###################################################
### code chunk number 26: plot_cndl_haz_nb_wb
###################################################
p_cndl_haz_cmb


###################################################
### code chunk number 27: gbcs_recode
###################################################
data(gbcsCS, package="condSURV")
gbcs_full <- gbcsCS
## to years
gbcs_full$rectime <- gbcs_full$rectime/365.25
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


###################################################
### code chunk number 28: gbcs-show
###################################################
gbcs_kbl <- kable(head(gbcs_curr[,c("id","diagdateb","rectime","censrec","tamoxifen","menopause","estrg_recp","size","grade")]),
                  digits=4, format="latex", booktabs=TRUE, linesep="",
                  row.names=FALSE, table.envir=NULL)
gbcs_kbl


###################################################
### code chunk number 29: gbcs_KM_mk
###################################################
grp_lvls <- c("Hist Cntl","Curr Cntl","Curr Trt")
fit_sf <- survfit(Surv(rectime, censrec)~group, data=gbcs_full)
KM_dat <- data.frame(time=fit_sf$time, surv=fit_sf$surv)
n_strat <- fit_sf$strata
KM_dat$group <-
    c(rep(grp_lvls[1], n_strat[1]), rep(grp_lvls[2], n_strat[2]),
      rep(grp_lvls[3], n_strat[3]))
p_KM <- ggplot(data=KM_dat) + geom_step(aes(time, surv, group=group,
                                            color=group))


###################################################
### code chunk number 30: ggcs_KM_do
###################################################
p_KM


###################################################
### code chunk number 31: xi-gbcs
###################################################
xifinder(b_tau = 0.001, d_tau = 5, p_0 = 0.5)


###################################################
### code chunk number 32: hyperpars
###################################################
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


###################################################
### code chunk number 33: tune
###################################################
tuning_parameters_gbcs <- list(cprop_beta = 1.17,
                               cprop_beta_0 = 1.21,
                               a_lambda = 0.5,
                               b_lambda = 0.5,
                               pi_b = 0.5,
                               alpha = 0.4)


###################################################
### code chunk number 34: test
###################################################
test <- 
BayesFBHborrow(Surv(rectime, censrec) ~ tamoxifen + menopause + size + grade,
               data=gbcs_curr, data_hist=gbcs_hist,
               model_choice="mix",
               tuning_parameters=tuning_parameters_gbcs,
               hyperparameters=hyperparameters_gbcs,
               iter=iter_tst, warmup_iter=warmup_iter_tst, refresh=0)


###################################################
### code chunk number 35: gbcs_fit_nb
###################################################
fit_gbcs_nb_cndl <- 
    BayesFBHborrow(Surv(rectime, censrec)~tamoxifen + menopause + size + grade,
                   data=gbcs_curr,
                   tuning_parameters=tuning_parameters_gbcs,
                   hyperparameters=hyperparameters_gbcs,
                   iter=iter, warmup_iter=warmup_iter,
                   refresh=warmup_iter,
                   verbose=TRUE)


###################################################
### code chunk number 36: save-gbcs-nb-rslts
###################################################

if(production)
{
  save(list=c("gbcs_full","gbcs_curr", "gbcs_hist", "hyperparameters_gbcs", "tuning_parameters_gbcs",
             "fit_gbcs_nb_cndl"), file="rdata/bfbhb-gbcs-nb-exmpl.rda")
} 
if(load_saved) {
  gbcs_nb_fnms <- suppressWarnings(system("ls -t rdata/*bfbhb-gbcs-nb-exmpl.rda", intern=TRUE, ignore.stderr=TRUE))[1]
  if(!is.na(gbcs_nb_fnms)) load(gbcs_nb_fnms[1])
}


###################################################
### code chunk number 37: gbcs_fit_wb
###################################################
fit_gbcs_wb_cndl <- 
  BayesFBHborrow(Surv(rectime, censrec)~ tamoxifen + menopause + size + grade,
                 data=gbcs_curr, data_hist=gbcs_hist,
                 model_choice="mix",
                 tuning_parameters=tuning_parameters_gbcs,
                 hyperparameters=hyperparameters_gbcs,
                 iter=iter, warmup_iter=warmup_iter,
                 refresh=warmup_iter,
                 verbose=TRUE,
                 max_grid=2000,
                 standardise=TRUE)


###################################################
### code chunk number 38: save-gbcs-data-rslts
###################################################
if(production){
  save(list=c("gbcs_full","gbcs_curr", "gbcs_hist", "hyperparameters_gbcs", "tuning_parameters_gbcs",
              "test", "fit_gbcs_wb_cndl"), file="rdata/bfbhb-gbcs-exmpls.rda")
} else {
  gbcs_exmpls_fnms <- suppressWarnings(system("ls -t rdata/*bfbhb-gbcs-exmpls.rda", intern=TRUE, ignore.stderr=TRUE))[1]
  if(!is.na(gbcs_exmpls_fnms)) load(gbcs_exmpls_fnms[1])
}


###################################################
### code chunk number 39: trace
###################################################
p_gbcs_wb_trace <- plot(fit_gbcs_wb_cndl, type="trace", col="tamoxifen")


###################################################
### code chunk number 40: fit_gbcs_wb_cndl_coef
###################################################
trteff_gbcs_wb_cndl_all <- coef(fit_gbcs_wb_cndl)[1,]
trteff_gbcs_wb_cndl <- unname(trteff_gbcs_wb_cndl_all["logHR"])
trteff_gbcs_wb_cndl_L <- unname(trteff_gbcs_wb_cndl_all["CrIn.L"])
trteff_gbcs_wb_cndl_U <- unname(trteff_gbcs_wb_cndl_all["CrIn.H"])
trteff_gbcs_wb_cndl_all <- coef(fit_gbcs_wb_cndl)[1,]
fit_gbcs_wb_cndl_coef_kbl <-
    kable(coef(fit_gbcs_wb_cndl), digits=4, format="latex", booktabs=TRUE, linesep="", table.envir=NULL)
fit_gbcs_wb_cndl_coef_kbl


###################################################
### code chunk number 41: mk-plots
###################################################
p_gbcs_cndl_nb_haz <- plot(fit_gbcs_nb_cndl, type="hazard", ylim=c(0,0.0010))
p_gbcs_cndl_wb_haz <- plot(fit_gbcs_wb_cndl, type="hazard", ylim=c(0,0.0010))
p_gbcs_cndl_haz_cmb <- Combine(p_gbcs_cndl_nb_haz, p_gbcs_cndl_wb_haz)

p_gbcs_cndl_nb_TrtEff <- plot(fit_gbcs_nb_cndl, type="TrtEff")
p_gbcs_cndl_wb_TrtEff <- plot(fit_gbcs_wb_cndl, type="TrtEff")
p_gbcs_cndl_TrtEff_cmb <- Combine(p_gbcs_cndl_nb_TrtEff, p_gbcs_cndl_wb_TrtEff,
                                 lgnd=c("NB","WB"))


###################################################
### code chunk number 42: p_gbcs_cndl_haz_plot
###################################################
p_gbcs_cndl_haz_cmb


###################################################
### code chunk number 43: p_gbcs_cndl_TrtEff_plot
###################################################
p_gbcs_cndl_TrtEff_cmb


###################################################
### code chunk number 44: plot_cndl_haz_hst
###################################################
p_gbcs_cndl_haz_cmb


###################################################
### code chunk number 45: trace
###################################################
p_gbcs_wb_trace

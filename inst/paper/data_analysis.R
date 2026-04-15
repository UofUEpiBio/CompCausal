###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(matrixStats)
library(tableone)

source(paste0(getwd(), "/HelperFunction.R"))
source(paste0(getwd(), "/singleindexmodelfunctions.R"))
source(paste0(getwd(), "/est_s_t_y.R"))
source(paste0(getwd(), "/est_exchange.R"))
source(paste0(getwd(), "/SIDR_Ravinew.R"))
source(paste0(getwd(), "/SIDRnew.R"))
source(paste0(getwd(), "/SensIAT_sim_outcome_modeler_mave.R"))

###############
## load data ##
###############

df <- readRDS(paste0(getwd(), "/TOIB3/TOIB_aggregated.rds"))
R <- df$isRCT
t <- df$isTropical
Y <- df$womac_12m
M <- ifelse(is.na(df$womac_12m), 0, 1)
n= length(t)

## X
X <- df %>% dplyr::select(c(age, womac_bq, expectationb, ChronicPainb))
## complete cases
complete_indx <- complete.cases(X)
X <- X[complete_indx, ]
R <- R[complete_indx]
t <- t[complete_indx]
Y <- Y[complete_indx]
M <- M[complete_indx]
n <- length(t)

data <- data.frame(cbind(Y, M, R, t))
colnames(data) <- c("Y", "M", "R", "t")

## relevel covariates
X$expectationb <- ifelse(X$expectationb %in% c(0, 1), "Much/A little worse",
                         ifelse(X$expectationb==2, "About the same", "A little/much better"))
X$expectationb <- relevel(factor(X$expectationb), ref="Much/A little worse")
X$ChronicPainb <- ifelse(X$ChronicPainb %in% c(0, 1), "1-2", "3-4")
X$ChronicPainb <- relevel(factor(X$ChronicPainb), ref="1-2")

#############
## table 1 ##
#############

dat <- cbind(R, t, X)
tab <- CreateTableOne(data=dat, strata=c("R", "t"))
print(tab, showAllLevels = TRUE)
tab_all <- CreateTableOne(data=dat, strata="R")
print(tab_all, showAllLevels = TRUE)

################
## estimation ##
################

gamma1 <- seq(-2, 2, length=51)
gamma0 <- seq(-2, 2, length=51)

fit_one_analysis <- function(trt_val, fold, seed, gamma) {
  # try mave first
  out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                     t=t, trt=trt_val, gamma=gamma, fold=fold,
                     IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                     single_index_method="norm1coef", use_mave=TRUE, seed=seed),
             silent = TRUE)
  
  if (inherits(out, "try-error")) {
    # fallback cumSIR
    out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                       t=t, trt=trt_val, gamma=gamma, fold=fold,
                       IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                       single_index_method="norm1coef", use_mave=FALSE, seed=seed),
               silent = TRUE)
  }
  
  if (inherits(out, "try-error") || is.null(out)) {
    return(NULL)
  }else{
    return(out)
  }
}

topical_WOMAC_12m <- fit_one_analysis(trt_val=1, fold=5, seed=2000056, gamma=gamma1)
oral_WOMAC_12m    <- fit_one_analysis(trt_val=0, fold=5, seed=2000056, gamma=gamma0)


##########################
## parametric bootstrap ##
##########################

nboot <- 500
set.seed(2000019)
seed_list <- sample(1:1000000, nboot)

## empty containers
est_topical_boot <- est_R0_topical_boot <- est_R1_topical_boot <- matrix(0, nrow=nboot, ncol=51)
var_topical_boot <- var_R0_topical_boot <- var_R1_topical_boot <- matrix(0, nrow=nboot, ncol=51)
est_oral_boot <- est_R0_oral_boot <- est_R1_oral_boot <- matrix(0, nrow=nboot, ncol=51)
var_oral_boot <- var_R0_oral_boot <- var_R1_oral_boot <- matrix(0, nrow=nboot, ncol=51)

diff_boot <- expand.grid(gamma_topical=gamma1, gamma_oral=gamma0, boot=1:nboot)
diff_boot$est <- diff_boot$est_R0 <- diff_boot$est_R1 <- 0
diff_boot$var <- diff_boot$var_R0 <- diff_boot$var_R1 <- 0

## fit models
X_with_t <- cbind(as.factor(t), X)
colnames(X_with_t)[1] <- "treatment"
gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
gam.var.M <- paste(gam.variables(X_with_t), collapse = "+") ## gam variables for missing data model
index.var.Y <- paste(single.index.variables(X), collapse = "+")
X_adjust <- model.matrix(as.formula(paste("~", index.var.Y)), data = X)[,-1]
X_adjust_scale <- scale(X_adjust)

X_adjust_scale_t_R0 <- X_adjust_scale[which(t==1 & R==0), ]
X_adjust_scale_t_R1 <- X_adjust_scale[which(t==1 & R==1), ]
X_adjust_scale_t0_R0 <- X_adjust_scale[which(t==0 & R==0), ]
X_adjust_scale_t0_R1 <- X_adjust_scale[which(t==0 & R==1), ]
M_t_R0 <- M[which(t==1 & R==0)]
M_t_R1 <- M[which(t==1 & R==1)]
M_t0_R0 <- M[which(t==0 & R==0)]
M_t0_R1 <- M[which(t==0 & R==1)]
Y_t_R0 <- Y[which(t==1 & R==0)]
Y_t_R1 <- Y[which(t==1 & R==1)]
Y_t0_R0 <- Y[which(t==0 & R==0)]
Y_t0_R1 <- Y[which(t==0 & R==1)]

X_R0 <- X[which(R==0), ]
t_R0 <- t[which(R==0)]
X_R1 <- X[which(R==1), ]
t_R1 <- t[which(R==1)]

X_with_t_R0 <- X_with_t[which(R==0), ]
X_with_t_R1 <- X_with_t[which(R==1), ]
M_R0 <- M[which(R==0)]
M_R1 <- M[which(R==1)]

g.fit <- mgcv::gam(as.formula(paste("R ~", gam.var)), data=X, family=binomial) 
t_R0.fit <- mgcv::gam(as.formula(paste("t_R0 ~", gam.var)), data=X_R0, family=binomial) ## treatment model
t_R1.fit <- mgcv::gam(as.formula(paste("t_R1 ~", gam.var)), data=X_R1, family=binomial) ## treatment model
M_R0.fit <- mgcv::gam(as.formula(paste("M_R0 ~", gam.var.M)), data=X_with_t_R0, family=binomial) ## missing data model
M_R1.fit <- mgcv::gam(as.formula(paste("M_R1 ~", gam.var.M)), data=X_with_t_R1, family=binomial) ## missing data model

fit_t_R1_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R1[which(M_t_R1==1), ],
                                                            Y = Y_t_R1[which(M_t_R1==1)],
                                                            ids=1:length(Y_t_R1[which(M_t_R1==1)]), 
                                                            kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                            use_mave=TRUE), silent = TRUE)
if (inherits(fit_t_R1_h, "try-error")) {
  fit_t_R1_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R1[which(M_t_R1==1), ],
                                                              Y = Y_t_R1[which(M_t_R1==1)],
                                                              ids=1:length(Y_t_R1[which(M_t_R1==1)]), 
                                                              kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                              use_mave=FALSE), silent = TRUE)
}

fit_t_R0_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R0[which(M_t_R0==1), ],
                                                            Y = Y_t_R0[which(M_t_R0==1)],
                                                            ids=1:length(Y_t_R0[which(M_t_R0==1)]), 
                                                            kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                            use_mave=TRUE), silent = TRUE)
if (inherits(fit_t_R0_h, "try-error")) {
  fit_t_R0_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R0[which(M_t_R0==1), ],
                                                              Y = Y_t_R0[which(M_t_R0==1)],
                                                              ids=1:length(Y_t_R0[which(M_t_R0==1)]), 
                                                              kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                              use_mave=FALSE), silent = TRUE)
}

fit_t0_R1_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t0_R1[which(M_t0_R1==1), ],
                                                             Y = Y_t0_R1[which(M_t0_R1==1)],
                                                             ids=1:length(Y_t0_R1[which(M_t0_R1==1)]), 
                                                             kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                             use_mave=TRUE), silent = TRUE)
if (inherits(fit_t0_R1_h, "try-error")) {
  fit_t0_R1_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t0_R1[which(M_t0_R1==1), ],
                                                               Y = Y_t0_R1[which(M_t0_R1==1)],
                                                               ids=1:length(Y_t0_R1[which(M_t0_R1==1)]), 
                                                               kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                               use_mave=FALSE), silent = TRUE)
}

fit_t0_R0_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t0_R0[which(M_t0_R0==1), ],
                                                             Y = Y_t0_R0[which(M_t0_R0==1)],
                                                             ids=1:length(Y_t0_R0[which(M_t0_R0==1)]), 
                                                             kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                             use_mave=TRUE), silent = TRUE)
if (inherits(fit_t0_R0_h, "try-error")) {
  fit_t0_R0_h  <- try(fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t0_R0[which(M_t0_R0==1), ],
                                                               Y = Y_t0_R0[which(M_t0_R0==1)],
                                                               ids=1:length(Y_t0_R0[which(M_t0_R0==1)]), 
                                                               kernel="dnorm", bw.selection="ise", bw.method="optim", 
                                                               use_mave=FALSE), silent = TRUE)
}

X_t_R0_beta_t_R0 <- as.vector(X_adjust_scale_t_R0[which(M_t_R0==1), ] %*% fit_t_R0_h$coef)
X_t_R1_beta_t_R1 <- as.vector(X_adjust_scale_t_R1[which(M_t_R1==1), ] %*% fit_t_R1_h$coef)
X_t0_R0_beta_t0_R0 <- as.vector(X_adjust_scale_t0_R0[which(M_t0_R0==1), ] %*% fit_t0_R0_h$coef)
X_t0_R1_beta_t0_R1 <- as.vector(X_adjust_scale_t0_R1[which(M_t0_R1==1), ] %*% fit_t0_R1_h$coef)
y_t_R0 = sort(unique(Y_t_R0[which(M_t_R0==1)]))    
ny_t_R0 = length(y_t_R0) 
y_t_R1 = sort(unique(Y_t_R1[which(M_t_R1==1)]))    
ny_t_R1 = length(y_t_R1) 
y_t0_R0 = sort(unique(Y_t0_R0[which(M_t0_R0==1)]))    
ny_t0_R0 = length(y_t0_R0) 
y_t0_R1 = sort(unique(Y_t0_R1[which(M_t0_R1==1)]))    
ny_t0_R1 = length(y_t0_R1) 

## parametric bootstrap
for(i in 1:nboot){
  
  ## simulate data
  boot.indices <- sample(1:n, size=n, replace = T)
  X.new <- X[boot.indices,]
  X_adjust_scale.new <- X_adjust_scale[boot.indices, ]
  ## new R
  prob.R <- stats::predict(g.fit, newdata=X.new, type="response") 
  R.new <- rbinom(n=n, size=1, prob=prob.R)
  ## new T
  T.new_R1 <- rbinom(n=n,size=1,prob=0.5)
  prob.T_R0 <- stats::predict(t_R0.fit, newdata=X.new, type="response")
  T.new_R0 <- rbinom(n=n,size=1,prob=prob.T_R0)
  T.new <- T.new_R1*R.new+T.new_R0*(1-R.new)
  
  ## new Y
  X.new_beta_t_R0 <- as.vector(X_adjust_scale.new %*% fit_t_R0_h$coef)
  X.new_beta_t_R1 <- as.vector(X_adjust_scale.new %*% fit_t_R1_h$coef)
  X.new_beta_t0_R0 <- as.vector(X_adjust_scale.new %*% fit_t0_R0_h$coef)
  X.new_beta_t0_R1 <- as.vector(X_adjust_scale.new %*% fit_t0_R1_h$coef)
  
  F_X_t_R0 <- NW_new(Xb=X_t_R0_beta_t_R0, Y=Y_t_R0[which(M_t_R0==1)], 
                     xb=X.new_beta_t_R0, y=y_t_R0, h=fit_t_R0_h$bandwidth, 
                     kernel = "dnorm")
  F_X_t_R1 <- NW_new(Xb=X_t_R1_beta_t_R1, Y=Y_t_R1[which(M_t_R1==1)], 
                     xb=X.new_beta_t_R1, y=y_t_R1, h=fit_t_R1_h$bandwidth, 
                     kernel = "dnorm")
  F_X_t0_R0 <- NW_new(Xb=X_t0_R0_beta_t0_R0, Y=Y_t0_R0[which(M_t0_R0==1)], 
                      xb=X.new_beta_t0_R0, y=y_t0_R0, h=fit_t0_R0_h$bandwidth, 
                      kernel = "dnorm")
  F_X_t0_R1 <- NW_new(Xb=X_t0_R1_beta_t0_R1, Y=Y_t0_R1[which(M_t0_R1==1)], 
                      xb=X.new_beta_t0_R1, y=y_t0_R1, h=fit_t0_R1_h$bandwidth, 
                      kernel = "dnorm")
  
  i1 = which(apply(F_X_t_R0==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X.new_beta_t_R0[i1], X.new_beta_t_R0[-i1], FUN = "-")), 1, which.min)
  F_X_t_R0[i1, ] <- F_X_t_R0[-i1, ][i1.closest, ]
  
  i1 = which(apply(F_X_t_R1==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X.new_beta_t_R1[i1], X.new_beta_t_R1[-i1], FUN = "-")), 1, which.min)
  F_X_t_R1[i1, ] <- F_X_t_R1[-i1, ][i1.closest, ]
  
  i1 = which(apply(F_X_t0_R0==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X.new_beta_t0_R0[i1], X.new_beta_t0_R0[-i1], FUN = "-")), 1, which.min)
  F_X_t0_R0[i1, ] <- F_X_t0_R0[-i1, ][i1.closest, ]
  
  i1 = which(apply(F_X_t0_R1==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X.new_beta_t0_R1[i1], X.new_beta_t0_R1[-i1], FUN = "-")), 1, which.min)
  F_X_t0_R1[i1, ] <- F_X_t0_R1[-i1, ][i1.closest, ]
  
  Y.new_t_R0 <- apply(F_X_t_R0, 1, function(x){
    ecdf <- runif(1,0,1)
    y_t_R0[min(which(ecdf<=x))]
  })
  Y.new_t_R1 <- apply(F_X_t_R1, 1, function(x){
    ecdf <- runif(1,0,1)
    y_t_R1[min(which(ecdf<=x))]
  })
  Y.new_t0_R0 <- apply(F_X_t0_R0, 1, function(x){
    ecdf <- runif(1,0,1)
    y_t0_R0[min(which(ecdf<=x))]
  })
  Y.new_t0_R1 <- apply(F_X_t0_R1, 1, function(x){
    ecdf <- runif(1,0,1)
    y_t0_R1[min(which(ecdf<=x))]
  })
  Y.new <- Y.new_t_R0*T.new*(1-R.new)+Y.new_t_R1*T.new*R.new+
    Y.new_t0_R0*(1-T.new)*(1-R.new)+Y.new_t0_R1*(1-T.new)*R.new
  
  ## new M
  X_with_t.new <- cbind(T.new, X.new)
  colnames(X_with_t.new)[1] <- "treatment"
  prob.M_R0 <- stats::predict(M_R0.fit, newdata=X_with_t.new, type="response")
  prob.M_R1 <- stats::predict(M_R1.fit, newdata=X_with_t.new, type="response")
  M.new_R0 <- rbinom(n=n,size=1,prob=prob.M_R0)
  M.new_R1 <- rbinom(n=n,size=1,prob=prob.M_R1)
  M.new <- M.new_R0*(1-R.new)+M.new_R1*R.new
  
  Y.new[M.new==0] <- NA
  
  data.new <- data.frame(cbind(Y.new, M.new, R.new, T.new, X.new))
  colnames(data.new) <- c("Y", "M", "R", "t", colnames(X.new))  
  
  ## run algorithm
  fit_one <- function(trt_val) {
    # try mave first
    out <- try(est_psi(Y=data.new$Y, M=data.new$M, R=data.new$R, X=X.new,
                       t=data.new$t, trt=trt_val, gamma=gamma1, fold=5,
                       IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm",
                       method="optim", single_index_method="norm1coef",
                       use_mave=TRUE, seed=seed_list[i]),
               silent = TRUE)
    
    if (inherits(out, "try-error")) {
      # fallback cumSIR
      out <- try(est_psi(Y=data.new$Y, M=data.new$M, R=data.new$R, X=X.new,
                         t=data.new$t, trt=trt_val, gamma=seq(-2, 2, length=51), fold=5,
                         IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm",
                         method="optim", single_index_method="norm1coef",
                         use_mave=FALSE, seed=seed_list[i]),
                 silent = TRUE)
    }
    
    if (inherits(out, "try-error") || is.null(out)) {
      return(NULL)
    }else{
      return(out)
    }
  }
  
  topical_vals <- fit_one(1)
  oral_vals    <- fit_one(0)
  
  if(is.null(topical_vals)){
    est_topical_boot[i, ] <- est_R0_topical_boot[i, ] <- est_R1_topical_boot[i, ] <- NA
    var_topical_boot[i, ] <- var_R0_topical_boot[i, ] <- var_R1_topical_boot[i, ] <- NA
  }else{
    est_topical_boot[i, ] <- topical_vals$est_trunc
    est_R0_topical_boot[i, ] <- topical_vals$est_trunc_R0
    est_R1_topical_boot[i, ] <- topical_vals$est_trunc_R1
    var_topical_boot[i, ] <- topical_vals$var_trunc
    var_R0_topical_boot[i, ] <- topical_vals$var_trunc_R0
    var_R1_topical_boot[i, ] <- topical_vals$var_trunc_R1

  }
  
  if(is.null(oral_vals)){
    est_oral_boot[i, ] <- est_R0_oral_boot[i, ] <- est_R1_oral_boot[i, ] <- NA
    var_oral_boot[i, ] <- var_R0_oral_boot[i, ] <- var_R1_oral_boot[i, ] <- NA
  }else{
    est_oral_boot[i, ] <- oral_vals$est_trunc
    est_R0_oral_boot[i, ] <- oral_vals$est_trunc_R0
    est_R1_oral_boot[i, ] <- oral_vals$est_trunc_R1
    var_oral_boot[i, ] <- oral_vals$var_trunc
    var_R0_oral_boot[i, ] <- oral_vals$var_trunc_R0
    var_R1_oral_boot[i, ] <- oral_vals$var_trunc_R1
  }

  
  ## Precompute once
  fold_list <- c(rep(1, 112), rep(c(2, 3, 4), each=113), rep(5, 112))
  fold_idx <- split(seq_along(fold_list),fold_list)
  
  ## If n is constant
  scale_n <- 1 / n
  
  topical_IF_mat <- t(do.call(rbind, topical_vals$IF_trunc))
  topical_IF_R0_mat <- t(do.call(rbind, topical_vals$IF_trunc_R0))
  topical_IF_R1_mat <- t(do.call(rbind, topical_vals$IF_trunc_R1))
  
  for (g_0 in seq_along(gamma0)) {
    
    indx_boot <- which(diff_boot$gamma_oral==gamma0[g_0]&diff_boot$boot==i)
    
    diff_boot$est[indx_boot]    <- topical_vals$est_trunc - oral_vals$est_trunc[g_0]
    diff_boot$est_R0[indx_boot] <- topical_vals$est_trunc_R0 - oral_vals$est_trunc_R0[g_0]
    diff_boot$est_R1[indx_boot] <- topical_vals$est_trunc_R1 - oral_vals$est_trunc_R1[g_0]
    
    IF_diff    <- topical_IF_mat - oral_vals$IF_trunc[[g_0]]
    IF_R0_diff <- topical_IF_R0_mat - oral_vals$IF_trunc_R0[[g_0]]
    IF_R1_diff <- topical_IF_R1_mat - oral_vals$IF_trunc_R1[[g_0]]
    
    var_temp    <- vapply(fold_idx, function(id) colVars(IF_diff[id, ]),    numeric(51))
    var_R0_temp <- vapply(fold_idx, function(id) colVars(IF_R0_diff[id, ]), numeric(51))
    var_R1_temp <- vapply(fold_idx, function(id) colVars(IF_R1_diff[id, ]), numeric(51))
    
    diff_boot$var[indx_boot]    <- rowMeans(var_temp)    * scale_n
    diff_boot$var_R0[indx_boot] <- rowMeans(var_R0_temp) * scale_n
    diff_boot$var_R1[indx_boot] <- rowMeans(var_R1_temp) * scale_n
    
  }
}

abs.topical <- abs(est_topical_boot-matrix(topical_WOMAC_12m$est_trunc, 
                                           ncol=length(gamma1), 
                                           nrow=nboot, byrow=T))/sqrt(var_topical_boot)
abs_R1.topical <- abs(est_R1_topical_boot-matrix(topical_WOMAC_12m$est_trunc_R1, 
                                                 ncol=length(gamma1), 
                                                 nrow=nboot, byrow=T))/sqrt(var_R1_topical_boot)
abs_R0.topical <- abs(est_R0_topical_boot-matrix(topical_WOMAC_12m$est_trunc_R0, 
                                                 ncol=length(gamma1), 
                                                 nrow=nboot, byrow=T))/sqrt(var_R0_topical_boot)

abs.oral <- abs(est_oral_boot-matrix(oral_WOMAC_12m$est_trunc, 
                                     ncol=length(gamma0), 
                                     nrow=nboot, byrow=T))/sqrt(var_oral_boot)
abs_R1.oral <- abs(est_R1_oral_boot-matrix(oral_WOMAC_12m$est_trunc_R1, 
                                           ncol=length(gamma1), 
                                           nrow=nboot, byrow=T))/sqrt(var_R1_oral_boot)
abs_R0.oral <- abs(est_R0_oral_boot-matrix(oral_WOMAC_12m$est_trunc_R0, 
                                           ncol=length(gamma1), 
                                           nrow=nboot, byrow=T))/sqrt(var_R0_oral_boot)

t.topical <- apply(abs.topical, 2, function(x){quantile(x, probs=0.95, na.rm=T)})
t_R1.topical <- apply(abs_R1.topical, 2, function(x){quantile(x, probs=0.95, na.rm=T)})
t_R0.topical <- apply(abs_R0.topical, 2, function(x){quantile(x, probs=0.95, na.rm=T)})

t.oral <- apply(abs.oral, 2, function(x){quantile(x, probs=0.95, na.rm=T)})
t_R1.oral <- apply(abs_R1.oral, 2, function(x){quantile(x, probs=0.95, na.rm=T)})
t_R0.oral <- apply(abs_R0.oral, 2, function(x){quantile(x, probs=0.95, na.rm=T)})

diff_t <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_t_R1 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_t_R0 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))

for(g_1 in 1:length(gamma1)){
  for(g_0 in 1:length(gamma0)){
    
    indx_t <- which(diff_boot$gamma_topical==gamma1[g_1]&diff_boot$gamma_oral==gamma0[g_0])
    
    abs_t <- abs(diff_boot$est[indx_t]-(topical_WOMAC_12m$est_trunc[g_1]-oral_WOMAC_12m$est_trunc[g_0]))/sqrt(diff_boot$var[indx_t])
    abs_t_R1 <- abs(diff_boot$est_R1[indx_t]-(topical_WOMAC_12m$est_trunc_R1[g_1]-oral_WOMAC_12m$est_trunc_R1[g_0]))/sqrt(diff_boot$var_R1[indx_t])
    abs_t_R0 <- abs(diff_boot$est_R0[indx_t]-(topical_WOMAC_12m$est_trunc_R0[g_1]-oral_WOMAC_12m$est_trunc_R0[g_0]))/sqrt(diff_boot$var_R0[indx_t])
    
    diff_t[g_1, g_0] <- quantile(abs_t, probs=0.95, na.rm=T)
    diff_t_R1[g_1, g_0] <- quantile(abs_t_R1, probs=0.95, na.rm=T)
    diff_t_R0[g_1, g_0] <- quantile(abs_t_R0, probs=0.95, na.rm=T)
    
  }
}


#####################
## estimation plot ##
#####################

r_topical <- data.frame(est=topical_WOMAC_12m$est_trunc, var=topical_WOMAC_12m$var_trunc, 
                        est_R0=topical_WOMAC_12m$est_trunc_R0, var_R0=topical_WOMAC_12m$var_trunc_R0, 
                        est_R1=topical_WOMAC_12m$est_trunc_R1, var_R1=topical_WOMAC_12m$var_trunc_R1, 
                        t=t.topical, t_R0=t_R0.topical, t_R1=t_R1.topical)

r_oral <- data.frame(est=oral_WOMAC_12m$est_trunc, var=oral_WOMAC_12m$var_trunc, 
                     est_R0=oral_WOMAC_12m$est_trunc_R0, var_R0=oral_WOMAC_12m$var_trunc_R0, 
                     est_R1=oral_WOMAC_12m$est_trunc_R1, var_R1=oral_WOMAC_12m$var_trunc_R1, 
                     t=t.oral, t_R0=t_R0.oral, t_R1=t_R1.oral)

est_plot <- data.frame(est=c(r_topical$est, r_oral$est, 
                             r_topical$est_R0, r_oral$est_R0, 
                             r_topical$est_R1, r_oral$est_R1), 
                       treatment=c(rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0)), 
                                   rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0)), 
                                   rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0))), 
                       lowerCI=c(r_topical$est-r_topical$t*sqrt(r_topical$var), 
                                 r_oral$est-r_oral$t*sqrt(r_oral$var), 
                                 r_topical$est_R0-r_topical$t_R0*sqrt(r_topical$var_R0), 
                                 r_oral$est_R0-r_oral$t_R0*sqrt(r_oral$var_R0), 
                                 r_topical$est_R1-r_topical$t_R1*sqrt(r_topical$var_R1), 
                                 r_oral$est_R1-r_oral$t_R1*sqrt(r_oral$var_R1)), 
                       upperCI=c(r_topical$est+r_topical$t*sqrt(r_topical$var), 
                                 r_oral$est+r_oral$t*sqrt(r_oral$var), 
                                 r_topical$est_R0+r_topical$t_R0*sqrt(r_topical$var_R0), 
                                 r_oral$est_R0+r_oral$t_R0*sqrt(r_oral$var_R0), 
                                 r_topical$est_R1+r_topical$t_R1*sqrt(r_topical$var_R1), 
                                 r_oral$est_R1+r_oral$t_R1*sqrt(r_oral$var_R1)), 
                       gamma=c(gamma1, gamma0, gamma1, gamma0, gamma1, gamma0), 
                       group=c(rep("CC", each=length(gamma1)+length(gamma0)), 
                               rep("PPS", each=length(gamma1)+length(gamma0)), 
                               rep("RCT", each=length(gamma1)+length(gamma0))))
est_plot$treatment <- factor(est_plot$treatment, levels=c("Topical", "Oral"))
est_plot$group <- factor(est_plot$group, levels=c("CC", "PPS", "RCT"))

## plot for estimates
plot1 <- ggplot(est_plot[which(est_plot$treatment=="Topical"), ], aes(x=gamma, y=est, color=group))+
  geom_line(linewidth = 1.2)+
  geom_line(aes(x=gamma, y=lowerCI), alpha=0.3, linewidth = 1.2, linetype = "dashed")+
  geom_line(aes(x=gamma, y=upperCI), alpha=0.3, linewidth = 1.2, linetype = "dashed")+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_color_manual(values=c("black", "blue", "green4"))+
  scale_x_continuous(breaks=seq(-2, 2, 1), limits = c(-2.05, 2.1)) +
  scale_y_continuous(breaks=seq(30, 55, 5), limits = c(30, 55)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))
plot1_f <-   plot1+
  annotate("text", x = 1.9, y = 40, label = as.character(expression(E*"["*Y(1)*"]")), parse = TRUE, size = 4)+
  annotate("text", x = 1.9, y = 42.8, label = as.character(expression(E*"["*Y(1)~'|'~R==0*"]")), parse = TRUE, size = 4, color="blue")+
  annotate("text", x = 1.9, y = 36.6, label = as.character(expression(E*"["*Y(1)~'|'~R==1*"]")), parse = TRUE, size = 4, color="green4")

plot2 <- plot1 + est_plot[which(est_plot$treatment=="Oral"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("text", x = 1.9, y = 41, label = as.character(expression(E*"["*Y(0)*"]")), parse = TRUE, size = 4)+
  annotate("text", x = 1.9, y = 43.8, label = as.character(expression(E*"["*Y(0)~'|'~R==0*"]")), parse = TRUE, size = 4, color="blue")+
  annotate("text", x = 1.9, y = 36.8, label = as.character(expression(E*"["*Y(0)~'|'~R==1*"]")), parse = TRUE, size = 4, color="green4")


png("est_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()

####################
## Counterfactual ##
####################

Counter_topical <- counterfactual(Y=Y, Y0=X$womac_bq, M=M, R=R, X=X, t=t, trt=1, 
                                  gamma=gamma1, est=topical_WOMAC_12m$est_trunc, 
                                  est_R0=topical_WOMAC_12m$est_trunc_R0, 
                                  est_R1=topical_WOMAC_12m$est_trunc_R1)

Counter_oral <- counterfactual(Y=Y, Y0=X$womac_bq, M=M, R=R, X=X, t=t, trt=0, 
                               gamma=gamma0, est=oral_WOMAC_12m$est_trunc, 
                               est_R0=oral_WOMAC_12m$est_trunc_R0, 
                               est_R1=oral_WOMAC_12m$est_trunc_R1)

Counter.plot.dat <- tibble(Data=c(Counter_topical$EY_t0, Counter_oral$EY_t0, Counter_topical$EY_t, Counter_oral$EY_t, 
                                  Counter_topical$EY_R0_t0, Counter_oral$EY_R0_t0, Counter_topical$EY_R0_t, Counter_oral$EY_R0_t, 
                                  Counter_topical$EY_R1_t0, Counter_oral$EY_R1_t0, Counter_topical$EY_R1_t, Counter_oral$EY_R1_t), 
                           Gamma=c(gamma1, gamma0, gamma1, gamma0, gamma1, gamma0, gamma1, gamma0, 
                                   gamma1, gamma0, gamma1, gamma0), 
                           Treatment=factor(c(rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0))), 
                                            levels=c("Topical", "Oral")), 
                           Estimate=factor(c(rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0)), 
                                             rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0)), 
                                             rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0))), 
                                           levels=c("Observed", "Induced")), 
                           group=factor(c(rep("CC", each=2*length(gamma1)+2*length(gamma0)), 
                                          rep("PPS", each=2*length(gamma1)+2*length(gamma0)), 
                                          rep("RCT", each=2*length(gamma1)+2*length(gamma0))), 
                                        levels=c("CC", "PPS", "RCT")))

Counter.plot.data_split <- split(Counter.plot.dat, f=Counter.plot.dat$Treatment)


Counter.plot1 <- ggplot(data=Counter.plot.data_split$Topical, 
                        aes(x=Gamma, y=Data, color=group, linetype = Estimate)) +
  geom_line() +
  scale_color_manual(values=c("black", "blue", "green4"))+
  scale_x_continuous(breaks=seq(-2, 2, 1), limits = c(-2.1, 2.1)) +
  scale_y_continuous(breaks=seq(25, 55, 5), limits = c(25, 55)) +
  facet_wrap(~Treatment, scales = "fixed")+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months") +
  theme_bw()+
  scale_linetype_manual(values= c("dashed", "solid"), guide = "none")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 18, face="bold"),
        axis.text.y = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 18, face="bold"),
        strip.text.x = element_text(size = 16, face="bold"))

Counter.plot1_f <- Counter.plot1+
  annotate("segment", x = 0, xend = 0.5, y = 34,  yend = 34, linetype = "dashed", color="black")+
  annotate("text",  x = 0.7, y = 34, label = as.character(expression(E*"["*Y(1)~'|'~T==1*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = 0, xend = 0.5, y = 32.5,  yend = 32.5, linetype = "solid", color="black")+
  annotate("text",  x = 0.7, y = 32.5, label = as.character(expression(E*"["*Y(1)~'|'~T==0*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = 0, xend = 0.5, y = 31,  yend = 31, linetype = "dashed", color="blue")+
  annotate("text",  x = 0.7, y = 31, label = as.character(expression(E*"["*Y(1)~'|'~T==1* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = 0, xend = 0.5, y = 29.5,  yend = 29.5, linetype = "solid", color="blue")+
  annotate("text",  x = 0.7, y = 29.5, label = as.character(expression(E*"["*Y(1)~'|'~T==0* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = 0, xend = 0.5, y = 28,  yend = 28, linetype = "dashed", color="green4")+
  annotate("text",  x = 0.7, y = 28, label = as.character(expression(E*"["*Y(1)~'|'~T==1* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")+
  annotate("segment", x = 0, xend = 0.5, y = 26.5,  yend = 26.5, linetype = "solid", color="green4")+
  annotate("text",  x = 0.7, y = 26.5, label = as.character(expression(E*"["*Y(1)~'|'~T==0* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")

Counter.plot2 <- Counter.plot1 + Counter.plot.data_split$Oral+
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x = 0, xend = 0.5, y = 34,  yend = 34, linetype = "dashed", color="black")+
  annotate("text",  x = 0.7, y = 34, label = as.character(expression(E*"["*Y(0)~'|'~T==0*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = 0, xend = 0.5, y = 32.5,  yend = 32.5, linetype = "solid", color="black")+
  annotate("text",  x = 0.7, y = 32.5, label = as.character(expression(E*"["*Y(0)~'|'~T==1*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = 0, xend = 0.5, y = 31,  yend = 31, linetype = "dashed", color="blue")+
  annotate("text",  x = 0.7, y = 31, label = as.character(expression(E*"["*Y(0)~'|'~T==0* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = 0, xend = 0.5, y = 29.5,  yend = 29.5, linetype = "solid", color="blue")+
  annotate("text",  x = 0.7, y = 29.5, label = as.character(expression(E*"["*Y(0)~'|'~T==1* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = 0, xend = 0.5, y = 28,  yend = 28, linetype = "dashed", color="green4")+
  annotate("text",  x = 0.7, y = 28, label = as.character(expression(E*"["*Y(0)~'|'~T==0* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")+
  annotate("segment", x = 0, xend = 0.5, y = 26.5,  yend = 26.5, linetype = "solid", color="green4")+
  annotate("text",  x = 0.7, y = 26.5, label = as.character(expression(E*"["*Y(0)~'|'~T==1* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")

png("Counter_Plot.png", width = 10, height = 7, units = 'in',res = 600)
grid.arrange(Counter.plot1_f,Counter.plot2, ncol=2)
dev.off()

#######################
## Treatment effects ##
#######################

diff <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_R0 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_R1 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var_R0 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var_R1 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))

fold_list <- c(rep(1, 112), rep(c(2, 3, 4), each=113), rep(5, 112))
fold_idx <- split(seq_along(fold_list),fold_list)
topical_IF_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc))
topical_IF_R0_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc_R0))
topical_IF_R1_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc_R1))
scale_n <- 1/n

for(g_0 in seq_along(gamma0)){
  
  diff[, g_0] <- topical_WOMAC_12m$est_trunc-oral_WOMAC_12m$est_trunc[g_0]
  diff_R0[, g_0] <- topical_WOMAC_12m$est_trunc_R0-oral_WOMAC_12m$est_trunc_R0[g_0]
  diff_R1[, g_0] <- topical_WOMAC_12m$est_trunc_R1-oral_WOMAC_12m$est_trunc_R1[g_0]
  
  IF_diff    <- topical_IF_mat - oral_WOMAC_12m$IF_trunc[[g_0]]
  IF_R0_diff <- topical_IF_R0_mat - oral_WOMAC_12m$IF_trunc_R0[[g_0]]
  IF_R1_diff <- topical_IF_R1_mat - oral_WOMAC_12m$IF_trunc_R1[[g_0]]
  
  var_temp    <- vapply(fold_idx, function(id) colVars(IF_diff[id, ]),    numeric(51))
  var_R0_temp <- vapply(fold_idx, function(id) colVars(IF_R0_diff[id, ]), numeric(51))
  var_R1_temp <- vapply(fold_idx, function(id) colVars(IF_R1_diff[id, ]), numeric(51))
  
  var[, g_0]    <- rowMeans(var_temp)    * scale_n
  var_R0[, g_0] <- rowMeans(var_R0_temp) * scale_n
  var_R1[, g_0] <- rowMeans(var_R1_temp) * scale_n
  
}

lowerCI <- diff-diff_t*sqrt(var)
lowerCI_R0 <- diff_R0-diff_t_R0*sqrt(var_R0)
lowerCI_R1 <- diff_R1-diff_t_R1*sqrt(var_R1)

upperCI <- diff+diff_t*sqrt(var)
upperCI_R0 <- diff_R0+diff_t_R0*sqrt(var_R0)
upperCI_R1 <- diff_R1+diff_t_R1*sqrt(var_R1)

## contour plot
png("contour_est_s_t_y.png", width = 10, height = 10, units = 'in',res = 300)
par(cex.lab = 2, cex.axis = 1.5, mgp = c(2.5, 1, 0))
contour(gamma1, gamma0, diff, nlevels = 15 , xlab = substitute(paste(bold(gamma[1]))), ylab = substitute(paste(bold(gamma[0]))), 
        lty="solid", lwd=1, col="black", labcex = 1.5, axes = F)
axis(1, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
axis(2, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
contour(gamma1, gamma0, diff, levels = 0, lwd = 3, add = T, col = "red", drawlabels = FALSE)
contour(gamma1, gamma0, lowerCI, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
contour(gamma1, gamma0, upperCI, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
text(x=0, y=0, labels="Inconclusive", cex=4, col="grey")
#rect(xleft=-2, ybottom=-2, xright=2, ytop=2, lty = 2)
#points(0, 0, pch = 3, col = "black", cex = 2, lwd = 2)
dev.off()

png("contour_est_R0_s_t_y.png", width = 10, height = 10, units = 'in',res = 300)
par(cex.lab = 2, cex.axis = 1.5, mgp = c(2.5, 1, 0))
contour(gamma1, gamma0, diff_R0, nlevels = 15 , xlab = substitute(paste(bold(gamma[1]))), ylab = substitute(paste(bold(gamma[0]))), 
        lty="solid", lwd=1, col="black", labcex = 1.5, axes = F)
axis(1, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
axis(2, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
contour(gamma1, gamma0, diff_R0, levels = 0, lwd = 3, add = T, col = "red", drawlabels = FALSE)
contour(gamma1, gamma0, lowerCI_R0, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
contour(gamma1, gamma0, upperCI_R0, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
text(x=0, y=0, labels="Inconclusive", cex=4, col="grey")
text(x=1.7, y=-1.95, labels="Favors Topical NSAIDs", cex=1, col="grey")
text(x=-1.75, y=1.95, labels="Favors Oral NSAIDs", cex=1, col="grey")
arrows(x0=-1.75, y0=1.75, x1=-1.9, y1=1.85, col="steelblue", lwd=2, length=0.15)
arrows(x0=1.7, y0=-1.7, x1=1.9, y1=-1.85, col="steelblue", lwd=2, length=0.15)
#rect(xleft=-2, ybottom=-2, xright=2, ytop=2, lty = 2)
#points(0, 0, pch = 3, col = "black", cex = 2, lwd = 2)
dev.off()

####################
## Exchangability ##
####################

gamma1_extend <- seq(-2.5, 2.5, length=51)
gamma0_extend <- seq(-2.5, 2.5, length=51)

fit_one_analysis <- function(trt_val, fold, seed, gamma) {
  # try mave first
  out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                     t=t, trt=trt_val, gamma=gamma, fold=fold,
                     IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim",
                     single_index_method="norm1coef", use_mave=TRUE, seed=seed),
             silent = TRUE)
  
  if (inherits(out, "try-error")) {
    # fallback cumSIR
    out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                       t=t, trt=trt_val, gamma=gamma, fold=fold,
                       IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim",
                       single_index_method="norm1coef", use_mave=FALSE, seed=seed),
               silent = TRUE)
  }
  
  if (inherits(out, "try-error") || is.null(out)) {
    return(NULL)
  }else{
    return(out)
  }
}

topical_WOMAC_12m_extend <- fit_one_analysis(trt_val=1, fold=5, seed=2000056, gamma=gamma1_extend)
oral_WOMAC_12m_extend <- fit_one_analysis(trt_val=0, fold=5, seed=2000056, gamma=gamma0_extend)


topical_WOMAC_12m_exchange <- est_psi_exchange(Y=Y, M=M, R=R, X=X, t=t, trt=1, 
                                               gamma=gamma1_extend, fold=5, seed=2000056, IF_output=TRUE, 
                                               simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                                               single_index_method="norm1coef", use_mave=TRUE)
oral_WOMAC_12m_exchange <- est_psi_exchange(Y=Y, M=M, R=R,  X=X, t=t, trt=0, 
                                            gamma=gamma0_extend, fold=5, seed=2000056, IF_output=TRUE, 
                                            simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                                            single_index_method="norm1coef", use_mave=TRUE)

est_plot_combine <- data.frame(est=c(topical_WOMAC_12m_extend$est_trunc, oral_WOMAC_12m_extend$est_trunc, 
                                     topical_WOMAC_12m_extend$est_trunc_R0, oral_WOMAC_12m_extend$est_trunc_R0), 
                               treatment=c(rep("Topical", each=length(gamma1_extend)), rep("Oral", each=length(gamma0_extend)), 
                                           rep("Topical", each=length(gamma1_extend)), rep("Oral", each=length(gamma0_extend))), 
                               gamma=c(gamma1_extend, gamma0_extend, gamma1_extend, gamma0_extend), 
                               group=c(rep("CC", each=length(gamma1_extend)+length(gamma0_extend)), rep("PPS", each=length(gamma1_extend)+length(gamma0_extend))))
est_plot_combine$treatment <- factor(est_plot_combine$treatment, levels=c("Topical", "Oral"))
est_plot_combine$group <- factor(est_plot_combine$group, levels=c("CC", "PPS"))

plot1 <- ggplot(est_plot_combine[which(est_plot_combine$treatment=="Topical"&est_plot_combine$group=="CC"), ], aes(x=gamma, y=est))+
  geom_line(linewidth = 1.2)+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_x_continuous(breaks=seq(-2.5, 2.5, 1), limits = c(-2.5, 2.5)) +
  scale_y_continuous(breaks=seq(33, 48, 5), limits = c(33, 48)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))

plot1_f <- plot1+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which(gamma1_extend==0)], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which(gamma1_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-0.5))], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-1))], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[1]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[1]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[1]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

plot2 <- plot1 + est_plot_combine[which(est_plot_combine$treatment=="Oral"&est_plot_combine$group=="CC"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which(gamma0_extend==0)], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which(gamma0_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-0.5))], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-1))], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[0]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[0]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[0]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

png("combine_exchange_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()

plot1 <- ggplot(est_plot_combine[which(est_plot_combine$treatment=="Topical"&est_plot_combine$group=="PPS"), ], aes(x=gamma, y=est))+
  geom_line(linewidth = 1.2)+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_x_continuous(breaks=seq(-2.5, 2.5, 1), limits = c(-2.5, 2.5)) +
  scale_y_continuous(breaks=seq(33, 48, 5), limits = c(33, 48)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))

plot1_f <- plot1+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which(gamma1_extend==0)], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which(gamma1_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-0.5))], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-1))], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[1]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[1]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[1]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

plot2 <- plot1 + est_plot_combine[which(est_plot_combine$treatment=="Oral"&est_plot_combine$group=="PPS"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which(gamma0_extend==0)], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which(gamma0_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-0.5))], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-1))], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[0]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[0]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[0]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

png("combine_exchange_R0_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()





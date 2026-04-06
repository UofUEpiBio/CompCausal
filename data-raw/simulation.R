###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(betareg)

source(paste0(getwd(), "/HelperFunction.R"))
source(paste0(getwd(), "/singleindexmodelfunctions.R"))
source(paste0(getwd(), "/est_s_t_y.R"))
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
## relevel covariates
X$expectationb <- ifelse(X$expectationb %in% c(0, 1), "Much/A little worse",
                         ifelse(X$expectationb==2, "About the same", "A little/much better"))
X$expectationb <- relevel(factor(X$expectationb), ref="Much/A little worse")
X$ChronicPainb <- ifelse(X$ChronicPainb %in% c(0, 1), "1-2", "3-4")
X$ChronicPainb <- relevel(factor(X$ChronicPainb), ref="1-2")  


###########
## truth ##
###########

truth_topical <- truth_beta(Y=Y, Y0=X$womac_bq, M=M, R=R, t=t, X=X, trt=1, gamma=seq(-2, 2, by=0.5))
truth_oral <- truth_beta(Y=Y, Y0=X$womac_bq, M=M, R=R, t=t, X=X, trt=0, gamma=seq(-2, 2, by=0.5))

################
## fit models ##
################

# rescale outcome to (0,1)
eps <- 1e-6
Y01 <- Y/100
# nudge any 0/1 slightly inside (only needed if any exact 0 or 20)
Y01[Y01 <= 0] <- eps
Y01[Y01 >= 1] <- 1 - eps

## set up covariates list and design matrix for modeling
X_with_t <- cbind(as.factor(t), X)
colnames(X_with_t)[1] <- "treatment"
gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
gam.var.M <- paste(gam.variables(X_with_t), collapse = "+") ## gam variables for missing data model
index.var.Y <- paste(single.index.variables(X), collapse = "+")

X_with_t_R0 <- X_with_t[which(R==0), ]
X_with_t_R1 <- X_with_t[which(R==1), ]
X_R0 <- X[which(R==0), ]
t_R0 <- t[which(R==0)]
M_R0 <- M[which(R==0)]
M_R1 <- M[which(R==1)]

## fit models 
g.fit <- mgcv::gam(as.formula(paste("R ~", gam.var)), data=cbind(R, X), family=binomial) 

t_R0.fit <- mgcv::gam(as.formula(paste("t_R0 ~", gam.var)), data=X_R0, family=binomial) ## treatment model

Y_t_R0.fit <- betareg(as.formula(paste("Y01 ~", index.var.Y)), data=cbind(Y01, X)[which(t==1 & R==0 & M==1), ])
Y_t0_R0.fit <- betareg(as.formula(paste("Y01 ~", index.var.Y)), data=cbind(Y01, X)[which(t==0 & R==0 & M==1), ])
Y_t_R1.fit <- betareg(as.formula(paste("Y01 ~", index.var.Y)), data=cbind(Y01, X)[which(t==1 & R==1 & M==1), ])
Y_t0_R1.fit <- betareg(as.formula(paste("Y01 ~", index.var.Y)), data=cbind(Y01, X)[which(t==0 & R==1 & M==1), ])

M_R0.fit <- mgcv::gam(as.formula(paste("M_R0 ~", gam.var.M)), data=X_with_t_R0, family=binomial) ## missing data model
M_R1.fit <- mgcv::gam(as.formula(paste("M_R1 ~", gam.var.M)), data=X_with_t_R1, family=binomial) ## missing data model

###################
## simulate data ##
###################

set.seed(10000898)

## get random sample from data
indices <- sample.int(n, size=sim.size, replace=T)
X.new <- X[indices,]

## new R
prob.R <- predict(g.fit, newdata=X.new, type="response") 
R.new <- rbinom(n=sim.size,size=1,prob=prob.R)

## new T
T.new_R1 <- rbinom(n=sim.size,size=1,prob=0.5)
prob.T_R0 <- predict(t_R0.fit, newdata=X.new, type="response")
T.new_R0 <- rbinom(n=sim.size,size=1,prob=prob.T_R0)
T.new <- T.new_R1*R.new+T.new_R0*(1-R.new)

## new Y
mu_t_R0 <- predict(Y_t_R0.fit, newdata=X.new, type="response")
mu_t0_R0 <- predict(Y_t0_R0.fit, newdata=X.new, type="response")
mu_t_R1 <- predict(Y_t_R1.fit, newdata=X.new, type="response")
mu_t0_R1 <- predict(Y_t0_R1.fit, newdata=X.new, type="response")

phi_t_R0 <- Y_t_R0.fit$coefficients$precision
phi_t0_R0 <- Y_t0_R0.fit$coefficients$precision
phi_t_R1 <- Y_t_R1.fit$coefficients$precision
phi_t0_R1 <- Y_t0_R1.fit$coefficients$precision

Y.new_t_R0 <- 100*mapply(rbeta, n=1, shape1=mu_t_R0*phi_t_R0, shape2=(1-mu_t_R0)*phi_t_R0)
Y.new_t0_R0 <- 100*mapply(rbeta, n=1, shape1=mu_t0_R0*phi_t0_R0, shape2=(1-mu_t0_R0)*phi_t0_R0)
Y.new_t_R1 <- 100*mapply(rbeta, n=1, shape1=mu_t_R1*phi_t_R1, shape2=(1-mu_t_R1)*phi_t_R1)
Y.new_t0_R1 <- 100*mapply(rbeta, n=1, shape1=mu_t0_R1*phi_t0_R1, shape2=(1-mu_t0_R1)*phi_t0_R1)
Y.new <- Y.new_t_R0*T.new*(1-R.new)+Y.new_t_R1*T.new*R.new+
  Y.new_t0_R0*(1-T.new)*(1-R.new)+Y.new_t0_R1*(1-T.new)*R.new

## new M
X_with_t.new <- cbind(T.new, X.new)
colnames(X_with_t.new)[1] <- "treatment"
prob.M_R0 <- predict(M_R0.fit, newdata=X_with_t.new, type="response")
prob.M_R1 <- predict(M_R1.fit, newdata=X_with_t.new, type="response")
M.new_R0 <- rbinom(n=sim.size,size=1,prob=prob.M_R0)
M.new_R1 <- rbinom(n=sim.size,size=1,prob=prob.M_R1)
M.new <- M.new_R0*(1-R.new)+M.new_R1*R.new

Y.new[M.new==0] <- NA

data <- data.frame(cbind(Y.new, M.new, R.new, T.new, X.new))
colnames(data) <- c("Y", "M", "R", "t", colnames(X.new))

################
## estimation ##
################

rand_seed <- 18890
fold <- 5

X = data.frame(dplyr::select(data, c(age, womac_bq, expectationb, ChronicPainb)))

topical_WOMAC_12m_sim <- fit_one_analysis(trt_val=1, fold=fold, seed=rand_seed, gamma=seq(-2, 2, by=0.5), IF_output=FALSE,
                                          single_index_method="norm1coef", method="optim", kernel="dnorm")
oral_WOMAC_12m_sim    <- fit_one_analysis(trt_val=0, fold=fold, seed=rand_seed, gamma=seq(-2, 2, by=0.5), IF_output=FALSE,
                                          single_index_method="norm1coef", method="optim", kernel="dnorm")






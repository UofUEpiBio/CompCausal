###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)

source(paste0(getwd(), "HelperFunction.R"))
source(paste0(getwd(), "singleindexmodelfunctions.R"))
source(paste0(getwd(), "est_s_t_y.R"))
source(paste0(getwd(), "SIDR_Ravinew.R"))
source(paste0(getwd(), "SIDRnew.R"))
source(paste0(getwd(), "SensIAT_sim_outcome_modeler_mave.R"))

nboot_B <- 10
nboot_C <- 10

gamma_length <- length(seq(-5, 5, by=1))

################
## simulation ##
################

## parametric bootstrap
load(paste0("/uufs/chpc.utah.edu/common/home/u6070035/CCS/simData/imputed", imputed, "/n", sim.size, "/sim", j, ".rds"))
  X_sim = data.frame(dplyr::select(data, c(age, pain_bq, expectationb, ChronicPainb)))
  
Boot.Est <- foreach(i=1:nrow(grid), .combine=rbind, .packages=package_list, .errorhandling = "pass") %dopar% {
  
  ## simulate outer bootstrap sample
  data.new <- pboot(data=data, X_sim=X_sim, sim.size=sim.size, seed=seed_B[grid$b[i]])
  X.new <- data.frame(dplyr::select(data.new, c(age, pain_bq, expectationb, ChronicPainb)))
  if(grid$c[i]==1){
    ## outer estimator
    topical_vals <- fit_one(data=data.new, X=X.new, trt_val=1)
    oral_vals    <- fit_one(data=data.new, X=X.new, trt_val=0)
  }else{
    topical_vals <- rep(NA_real_, 6*gamma_length)
    oral_vals <- rep(NA_real_, 6*gamma_length)
  }

  ## simulate inner bootstrap sample
  data.new2 <- pboot(data=data.new, X_sim=X.new, sim.size=sim.size, seed=NULL)
  X.new2 <- data.frame(dplyr::select(data.new2, c(age, pain_bq, expectationb, ChronicPainb)))
  ## inner estimator
  temp_topical <- fit_one(data=data.new2, X=X.new2, trt_val=1)
  temp_oral <- fit_one(data=data.new2, X=X.new2, trt_val=0)
  
  ## return result
  c(topical_vals, oral_vals, temp_topical, temp_oral)
  

load(paste0("/uufs/chpc.utah.edu/common/home/u6070035/CCS/simResult/imputed", imputed, "/n", sim.size, "/", kernel, "_", single_index_method, "_fold5/sim", j, ".RData"))
Q_b_topical <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 133:143]<=matrix(topical_WOMAC_12m$est_trunc, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})
Q_b_R1_topical <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 144:154]<=matrix(topical_WOMAC_12m$est_trunc_R1, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})
Q_b_R0_topical <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 155:165]<=matrix(topical_WOMAC_12m$est_trunc_R0, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})
Q_b_oral <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 199:209]<=matrix(oral_WOMAC_12m$est_trunc, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})
Q_b_R1_oral <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 210:220]<=matrix(oral_WOMAC_12m$est_trunc_R1, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})
Q_b_R0_oral <- sapply(1:nboot_B, function(x){
  loc <- which(grid$b==x)
  colSums(Boot.Est[loc, 221:231]<=matrix(oral_WOMAC_12m$est_trunc_R0, 
                                         nrow=nboot_C, ncol=gamma_length, 
                                         byrow=T), na.rm=T)/nboot_C
})

t_bc_topical <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 133:143]-Boot.Est[loc_b, 1:11])/sqrt(Boot.Est[x, 166:176])
})
t_b_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 1:11]-topical_WOMAC_12m$est_trunc)/sqrt(Boot.Est[loc_b, 34:44])
})
Q_t_b_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_topical[loc_b, ]<=matrix(t_b_topical[x, ], 
                                        nrow=nboot_C, ncol=gamma_length, 
                                        byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_topical[loc_b, ])<=matrix(abs(t_b_topical[x, ]), 
                                        nrow=nboot_C, ncol=gamma_length, 
                                        byrow=T), na.rm=T)/nboot_C
})
t_bc_R1_topical <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 144:154]-Boot.Est[loc_b, 12:22])/sqrt(Boot.Est[x, 177:187])
})
t_b_R1_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 12:22]-topical_WOMAC_12m$est_trunc_R1)/sqrt(Boot.Est[loc_b, 45:55])
})
Q_t_b_R1_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_R1_topical[loc_b, ]<=matrix(t_b_R1_topical[x, ], 
                                        nrow=nboot_C, ncol=gamma_length, 
                                        byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_R1_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_R1_topical[loc_b, ])<=matrix(abs(t_b_R1_topical[x, ]), 
                                             nrow=nboot_C, ncol=gamma_length, 
                                             byrow=T), na.rm=T)/nboot_C
})
t_bc_R0_topical <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 155:165]-Boot.Est[loc_b, 23:33])/sqrt(Boot.Est[x, 188:198])
})
t_b_R0_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 23:33]-topical_WOMAC_12m$est_trunc_R0)/sqrt(Boot.Est[loc_b, 56:66])
})
Q_t_b_R0_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_R0_topical[loc_b, ]<=matrix(t_b_R0_topical[x, ], 
                                           nrow=nboot_C, ncol=gamma_length, 
                                           byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_R0_topical <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_R0_topical[loc_b, ])<=matrix(abs(t_b_R0_topical[x, ]), 
                                                nrow=nboot_C, ncol=gamma_length, 
                                                byrow=T), na.rm=T)/nboot_C
})
t_bc_oral <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 199:209]-Boot.Est[loc_b, 67:77])/sqrt(Boot.Est[x, 232:242])
})
t_b_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 67:77]-oral_WOMAC_12m$est_trunc)/sqrt(Boot.Est[loc_b, 100:110])
})
Q_t_b_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_oral[loc_b, ]<=matrix(t_b_oral[x, ], 
                                     nrow=nboot_C, ncol=gamma_length, 
                                     byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_oral[loc_b, ])<=matrix(abs(t_b_oral[x, ]), 
                                          nrow=nboot_C, ncol=gamma_length, 
                                          byrow=T), na.rm=T)/nboot_C
})
t_bc_R1_oral <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 210:220]-Boot.Est[loc_b, 78:88])/sqrt(Boot.Est[x, 243:253])
})
t_b_R1_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 78:88]-oral_WOMAC_12m$est_trunc_R1)/sqrt(Boot.Est[loc_b, 111:121])
})
Q_t_b_R1_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_R1_oral[loc_b, ]<=matrix(t_b_R1_oral[x, ], 
                                        nrow=nboot_C, ncol=gamma_length, 
                                        byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_R1_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_R1_oral[loc_b, ])<=matrix(abs(t_b_R1_oral[x, ]), 
                                             nrow=nboot_C, ncol=gamma_length, 
                                             byrow=T), na.rm=T)/nboot_C
})
t_bc_R0_oral <- sapply(1:nrow(grid), function(x){
  loc_b <- which(grid$b==grid$b[x]&grid$c==1)
  (Boot.Est[x, 221:231]-Boot.Est[loc_b, 89:99])/sqrt(Boot.Est[x, 254:264])
})
t_b_R0_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 89:99]-oral_WOMAC_12m$est_trunc_R0)/sqrt(Boot.Est[loc_b, 122:132])
})
Q_t_b_R0_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(t_bc_R0_oral[loc_b, ]<=matrix(t_b_R0_oral[x, ], 
                                        nrow=nboot_C, ncol=gamma_length, 
                                        byrow=T), na.rm=T)/nboot_C
})
abs_Q_t_b_R0_oral <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  colSums(abs(t_bc_R0_oral[loc_b, ])<=matrix(abs(t_b_R0_oral[x, ]), 
                                             nrow=nboot_C, ncol=gamma_length, 
                                             byrow=T), na.rm=T)/nboot_C
})
t_b_sd_topical <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 1:11]-topical_WOMAC_12m$est_trunc)/apply(Boot.Est[loc, 133:143], 2, function(x){sd(x, na.rm=TRUE)})
})
t_b_R1_sd_topical <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 12:22]-topical_WOMAC_12m$est_trunc_R1)/apply(Boot.Est[loc, 144:154], 2, function(x){sd(x, na.rm=TRUE)})
})
t_b_R0_sd_topical <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 23:33]-topical_WOMAC_12m$est_trunc_R0)/apply(Boot.Est[loc, 155:165], 2, function(x){sd(x, na.rm=TRUE)})
})
t_b_sd_oral <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 67:77]-oral_WOMAC_12m$est_trunc)/apply(Boot.Est[loc, 199:209], 2, function(x){sd(x, na.rm=TRUE)})
})
t_b_R1_sd_oral <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 78:88]-oral_WOMAC_12m$est_trunc_R1)/apply(Boot.Est[loc, 210:220], 2, function(x){sd(x, na.rm=TRUE)})
})
t_b_R0_sd_oral <- sapply(1:nboot_B, function(x){ 
  loc <- which(grid$b==x)
  loc_b <- which(grid$b==x&grid$c==1)
  (Boot.Est[loc_b, 89:99]-oral_WOMAC_12m$est_trunc_R0)/apply(Boot.Est[loc, 221:231], 2, function(x){sd(x, na.rm=TRUE)})
})

est_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 1:11]
})
est_R1_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 12:22]
})
est_R0_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 23:33]
})
var_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 34:44]
})
var_R1_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 45:55]
})
var_R0_topical_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 56:66]
})
est_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 67:77]
})
est_R1_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 78:88]
})
est_R0_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 89:99]
})
var_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 100:110]
})
var_R1_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 111:121]
})
var_R0_oral_boot <- sapply(1:nboot_B, function(x){
  loc_b <- which(grid$b==x&grid$c==1)
  Boot.Est[loc_b, 122:132]
})

## bootstrap result
save(est_topical_boot, est_R1_topical_boot, est_R0_topical_boot, 
     est_oral_boot, est_R1_oral_boot, est_R0_oral_boot, 
     var_topical_boot, var_R1_topical_boot, var_R0_topical_boot, 
     var_oral_boot, var_R1_oral_boot, var_R0_oral_boot, 
     Q_b_topical, Q_b_R1_topical, Q_b_R0_topical, 
     Q_b_oral, Q_b_R1_oral, Q_b_R0_oral, 
     t_b_topical, t_b_R1_topical, t_b_R0_topical, 
     t_b_oral, t_b_R1_oral, t_b_R0_oral, 
     Q_t_b_topical, Q_t_b_R1_topical, Q_t_b_R0_topical, 
     Q_t_b_oral, Q_t_b_R1_oral, Q_t_b_R0_oral, 
     abs_Q_t_b_topical, abs_Q_t_b_R1_topical, abs_Q_t_b_R0_topical, 
     abs_Q_t_b_oral, abs_Q_t_b_R1_oral, abs_Q_t_b_R0_oral, 
     t_b_sd_topical, t_b_R1_sd_topical, t_b_R0_sd_topical, 
     t_b_sd_oral, t_b_R1_sd_oral, t_b_R0_sd_oral,
     file=paste0("/uufs/chpc.utah.edu/common/home/u6070035/CCS/simResult/imputed", imputed, "/double_pboot_n", sim.size, "/", kernel, "_", single_index_method, "/sim", j, ".RData"))
  print(j)
  
}

stopCluster(cl)





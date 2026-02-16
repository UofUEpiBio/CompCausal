###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)

source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/HelperFunction.R")
source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/singleindexmodelfunctions.R")
source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/est_s_t_y.R")
source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/SIDR_Ravinew.R")
source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/SIDRnew.R")
source("/uufs/chpc.utah.edu/common/home/u6070035/CCS/code/SensIAT_sim_outcome_modeler_mave.R")

ntasks <- Sys.getenv("SLURM_NTASKS")
if (ntasks == '') {
  ntasks <- 4
} else {
  ntasks <- strtoi(ntasks) }
cat("This script use ", ntasks, " cores\n")

package_list <- c("splines", "tidyverse", "gridExtra")

nboot_B <- 500

################
## simulation ##
################

cl <- makeCluster(ntasks)
registerDoParallel(cl)

for (j in (40*(index-1)+1):(40*index)){
  
Boot.Est <- foreach(i=1:nboot_B, .combine=rbind, .packages=package_list, .errorhandling = "pass") %dopar% {
    
  gamma_length <- length(seq(-5, 5, by=1))
  nboot_C <- 100
  
  load(paste0("/uufs/chpc.utah.edu/common/home/u6070035/CCS/simData/imputed", imputed, "/n", sim.size, "/sim", j, ".rds"))
  X_sim = data.frame(dplyr::select(data, c(age, pain_bq, expectationb, ChronicPainb)))
  
  ## simulate outer bootstrap sample
  boot_b_temp <- pboot(data=data, X_sim=X_sim, sim.size=sim.size, seed=NULL)
  data.new <- boot_b_temp$data.new
  
  X.new <- data.frame(dplyr::select(data.new, c(age, pain_bq, expectationb, ChronicPainb)))
  topical_vals <- fit_one(data=data.new, X=X.new, trt_val=1, coef_g.fit=boot_b_temp$coef_g.fit, 
                          coef_t_R0.fit=boot_b_temp$coef_t_R0.fit, coef_t_R1.fit=boot_b_temp$coef_t_R1.fit, 
                          coef_M_R0.fit=boot_b_temp$coef_M_R0.fit, coef_M_R1.fit=boot_b_temp$coef_M_R1.fit)
  oral_vals    <- fit_one(data=data.new, X=X.new, trt_val=0, coef_g.fit=boot_b_temp$coef_g.fit, 
                          coef_t_R0.fit=boot_b_temp$coef_t_R0.fit, coef_t_R1.fit=boot_b_temp$coef_t_R1.fit, 
                          coef_M_R0.fit=boot_b_temp$coef_M_R0.fit, coef_M_R1.fit=boot_b_temp$coef_M_R1.fit)
  boot_b_fit <- pboot_model(data=data.new, X_sim=X.new)
  
  load(paste0("/uufs/chpc.utah.edu/common/home/u6070035/CCS/simResult/imputed", imputed, "/n", sim.size, "/", kernel, "_", single_index_method, "_fold5/sim", j, ".RData"))
  
  t_b_topical <- (topical_vals[1:11]-topical_WOMAC_12m$est_trunc)/sqrt(topical_vals[34:44])
  t_b_R1_topical <- (topical_vals[12:22]-topical_WOMAC_12m$est_trunc_R1)/sqrt(topical_vals[45:55])
  t_b_R0_topical <- (topical_vals[23:33]-topical_WOMAC_12m$est_trunc_R0)/sqrt(topical_vals[56:66])
  t_b_oral <- (oral_vals[1:11]-oral_WOMAC_12m$est_trunc)/sqrt(oral_vals[34:44])
  t_b_R1_oral <- (oral_vals[12:22]-oral_WOMAC_12m$est_trunc_R1)/sqrt(oral_vals[45:55])
  t_b_R0_oral <- (oral_vals[23:33]-oral_WOMAC_12m$est_trunc_R0)/sqrt(oral_vals[56:66])
  
  ## empty containers
  Q_b_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  Q_b_R1_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  Q_b_R0_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  Q_b_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  Q_b_R1_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_R1_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_R0_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_R1_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  t_bc_R0_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_R1_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_R0_topical <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_R1_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  est_R0_oral <- matrix(NA, nrow=nboot_C, ncol=gamma_length)
  
  for(j in 1:nboot_C){
    
    data.new2 <- pboot_sim(data=data.new, X_sim=X.new, sim.size=sim.size, seed=NULL, 
                           fit_t_R0_h=boot_b_fit$fit_t_R0_h, fit_t_R1_h=boot_b_fit$fit_t_R1_h, 
                           fit_t0_R0_h=boot_b_fit$fit_t0_R0_h, fit_t0_R1_h=boot_b_fit$fit_t0_R1_h, 
                           g.fit=boot_b_fit$g.fit, t_R0.fit=boot_b_fit$t_R0.fit, M_R0.fit=boot_b_fit$M_R0.fit, 
                           M_R1.fit=boot_b_fit$M_R1.fit)
    X.new2 <- data.frame(dplyr::select(data.new2, c(age, pain_bq, expectationb, ChronicPainb)))
    ## inner estimator
    temp_topical <- fit_one(data=data.new2, X=X.new2, trt_val=1, coef_g.fit=coef(boot_b_fit$g.fit), 
                            coef_t_R0.fit=coef(boot_b_fit$t_R0.fit), coef_t_R1.fit=coef(boot_b_fit$t_R1.fit),
                            coef_M_R0.fit=coef(boot_b_fit$M_R0.fit), coef_M_R1.fit=coef(boot_b_fit$M_R1.fit))
    temp_oral <- fit_one(data=data.new2, X=X.new2, trt_val=0, coef_g.fit=coef(boot_b_fit$g.fit), 
                         coef_t_R0.fit=coef(boot_b_fit$t_R0.fit), coef_t_R1.fit=coef(boot_b_fit$t_R1.fit),
                         coef_M_R0.fit=coef(boot_b_fit$M_R0.fit), coef_M_R1.fit=coef(boot_b_fit$M_R1.fit))
    
    Q_b_topical[j, ] <- as.numeric(temp_topical[1:11]<=topical_WOMAC_12m$est_trunc)
    Q_b_R1_topical[j, ] <- as.numeric(temp_topical[12:22]<=topical_WOMAC_12m$est_trunc_R1)
    Q_b_R0_topical[j, ] <- as.numeric(temp_topical[23:33]<=topical_WOMAC_12m$est_trunc_R0)
    Q_b_oral[j, ] <- as.numeric(temp_oral[1:11]<=oral_WOMAC_12m$est_trunc)
    Q_b_R1_oral[j, ] <- as.numeric(temp_oral[12:22]<=oral_WOMAC_12m$est_trunc_R1)
    Q_b_R0_oral[j, ] <- as.numeric(temp_oral[23:33]<=oral_WOMAC_12m$est_trunc_R0)
    
    t_bc_topical[j, ] <- (temp_topical[1:11]-topical_vals[1:11])/sqrt(temp_topical[34:44])
    t_bc_R1_topical[j, ] <- (temp_topical[12:22]-topical_vals[12:22])/sqrt(temp_topical[45:55])
    t_bc_R0_topical[j, ] <- (temp_topical[23:33]-topical_vals[23:33])/sqrt(temp_topical[56:66])
    t_bc_oral[j, ] <- (temp_oral[1:11]-oral_vals[1:11])/sqrt(temp_topical[34:44])
    t_bc_R1_oral[j, ] <- (temp_oral[12:22]-oral_vals[12:22])/sqrt(temp_oral[45:55])
    t_bc_R0_oral[j, ] <- (temp_oral[23:33]-oral_vals[23:33])/sqrt(temp_oral[56:66])
    
    est_topical[j, ] <- temp_topical[1:11]
    est_R1_topical[j, ] <- temp_topical[12:22]
    est_R0_topical[j, ] <- temp_topical[23:33]
    est_oral[j, ] <- temp_oral[1:11]
    est_R1_oral[j, ] <- temp_oral[12:22]
    est_R0_oral[j, ] <- temp_oral[23:33]
    
  }
  
  r_Q_b_topical <- colMeans(Q_b_topical, na.rm=T)
  r_Q_b_R1_topical <- colMeans(Q_b_R1_topical, na.rm=T)
  r_Q_b_R0_topical <- colMeans(Q_b_R0_topical, na.rm=T)
  r_Q_b_oral <- colMeans(Q_b_oral, na.rm=T)
  r_Q_b_R1_oral <- colMeans(Q_b_R1_oral, na.rm=T)
  r_Q_b_R0_oral <- colMeans(Q_b_R0_oral, na.rm=T)
  Q_t_b_topical <- colMeans(t_bc_topical<=matrix(t_b_topical, nrow=nboot_C, ncol=gamma_length, 
                                                 byrow=T), na.rm=T)
  abs_Q_t_b_topical <- colMeans(abs(t_bc_topical)<=matrix(abs(t_b_topical), nrow=nboot_C, ncol=gamma_length, 
                                                          byrow=T), na.rm=T)
  Q_t_b_R1_topical <- colMeans(t_bc_R1_topical<=matrix(t_b_R1_topical, nrow=nboot_C, ncol=gamma_length, 
                                                       byrow=T), na.rm=T)
  abs_Q_t_b_R1_topical <- colMeans(abs(t_bc_R1_topical)<=matrix(abs(t_b_R1_topical), nrow=nboot_C, ncol=gamma_length, 
                                                                byrow=T), na.rm=T)
  Q_t_b_R0_topical <- colMeans(t_bc_R0_topical<=matrix(t_b_R0_topical, nrow=nboot_C, ncol=gamma_length, 
                                                       byrow=T), na.rm=T)
  abs_Q_t_b_R0_topical <- colMeans(abs(t_bc_R0_topical)<=matrix(abs(t_b_R0_topical), nrow=nboot_C, ncol=gamma_length, 
                                                                byrow=T), na.rm=T)
  Q_t_b_oral <- colMeans(t_bc_oral<=matrix(t_b_oral, nrow=nboot_C, ncol=gamma_length, 
                                           byrow=T), na.rm=T)
  abs_Q_t_b_oral <- colMeans(abs(t_bc_oral)<=matrix(abs(t_b_oral), nrow=nboot_C, ncol=gamma_length, 
                                                    byrow=T), na.rm=T)
  Q_t_b_R1_oral <- colMeans(t_bc_R1_oral<=matrix(t_b_R1_oral, nrow=nboot_C, ncol=gamma_length, 
                                                 byrow=T), na.rm=T)
  abs_Q_t_b_R1_oral <- colMeans(abs(t_bc_R1_oral)<=matrix(abs(t_b_R1_oral), nrow=nboot_C, ncol=gamma_length, 
                                                          byrow=T), na.rm=T)
  Q_t_b_R0_oral <- colMeans(t_bc_R0_oral<=matrix(t_b_R0_oral, nrow=nboot_C, ncol=gamma_length, 
                                                 byrow=T), na.rm=T)
  abs_Q_t_b_R0_oral <- colMeans(abs(t_bc_R0_oral)<=matrix(abs(t_b_R0_oral), nrow=nboot_C, ncol=gamma_length, 
                                                          byrow=T), na.rm=T)
  t_b_sd_topical <- (topical_vals[1:11]-topical_WOMAC_12m$est_trunc)/apply(est_topical, 2, function(x){sd(x, na.rm=TRUE)})
  t_b_R1_sd_topical <- (topical_vals[12:22]-topical_WOMAC_12m$est_trunc_R1)/apply(est_R1_topical, 2, function(x){sd(x, na.rm=TRUE)})
  t_b_R0_sd_topical <- (topical_vals[23:33]-topical_WOMAC_12m$est_trunc_R0)/apply(est_R0_topical, 2, function(x){sd(x, na.rm=TRUE)})
  t_b_sd_oral <- (oral_vals[1:11]-oral_WOMAC_12m$est_trunc)/apply(est_oral, 2, function(x){sd(x, na.rm=TRUE)})
  t_b_R1_sd_oral <- (oral_vals[12:22]-oral_WOMAC_12m$est_trunc_R1)/apply(est_R1_oral, 2, function(x){sd(x, na.rm=TRUE)})
  t_b_R0_sd_oral <- (oral_vals[23:33]-oral_WOMAC_12m$est_trunc_R0)/apply(est_R0_oral, 2, function(x){sd(x, na.rm=TRUE)})

return(c(topical_vals, oral_vals, r_Q_b_topical, r_Q_b_R1_topical, r_Q_b_R0_topical, 
         r_Q_b_oral, r_Q_b_R1_oral, r_Q_b_R0_oral, Q_t_b_topical, Q_t_b_R1_topical, Q_t_b_R0_topical, 
         Q_t_b_oral, Q_t_b_R1_oral, Q_t_b_R0_oral, abs_Q_t_b_topical, abs_Q_t_b_R1_topical, abs_Q_t_b_R0_topical, 
         abs_Q_t_b_oral, abs_Q_t_b_R1_oral, abs_Q_t_b_R0_oral, t_b_sd_topical, t_b_R1_sd_topical, t_b_R0_sd_topical, 
         t_b_sd_oral, t_b_R1_sd_oral, t_b_R0_sd_oral))

}

est_topical_boot <- Boot.Est[, 1:11]
est_R1_topical_boot <- Boot.Est[, 12:22]
est_R0_topical_boot <- Boot.Est[, 23:33]
var_topical_boot <- Boot.Est[, 34:44]
var_R1_topical_boot <- Boot.Est[, 45:55]
var_R0_topical_boot <- Boot.Est[, 56:66]

est_oral_boot <- Boot.Est[, 67:77]
est_R1_oral_boot <- Boot.Est[, 78:88]
est_R0_oral_boot <- Boot.Est[, 89:99]
var_oral_boot <- Boot.Est[, 100:110]
var_R1_oral_boot <- Boot.Est[, 111:121]
var_R0_oral_boot <- Boot.Est[, 122:132]

Q_b_topical <- Boot.Est[, 133:143]
Q_b_R1_topical <- Boot.Est[, 144:154]
Q_b_R0_topical <- Boot.Est[, 155:165]
Q_b_oral <- Boot.Est[, 166:176]
Q_b_R1_oral <- Boot.Est[, 177:187]
Q_b_R0_oral <- Boot.Est[, 188:198]

Q_t_b_topical <- Boot.Est[, 199:209]
Q_t_b_R1_topical <- Boot.Est[, 210:220]
Q_t_b_R0_topical <- Boot.Est[, 221:231]
Q_t_b_oral <- Boot.Est[, 232:242]
Q_t_b_R1_oral <- Boot.Est[, 243:253]
Q_t_b_R0_oral <- Boot.Est[, 254:264]

abs_Q_t_b_topical <- Boot.Est[, 265:275]
abs_Q_t_b_R1_topical <- Boot.Est[, 276:286]
abs_Q_t_b_R0_topical <- Boot.Est[, 287:297]
abs_Q_t_b_oral <- Boot.Est[, 298:308]
abs_Q_t_b_R1_oral <- Boot.Est[, 309:319]
abs_Q_t_b_R0_oral <- Boot.Est[, 320:330]

t_b_sd_topical <- Boot.Est[, 331:341]
t_b_R1_sd_topical <- Boot.Est[, 342:352]
t_b_R0_sd_topical <- Boot.Est[, 353:363]
t_b_sd_oral <- Boot.Est[, 364:374]
t_b_R1_sd_oral <- Boot.Est[, 375:385]
t_b_R0_sd_oral <- Boot.Est[, 386:396]

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







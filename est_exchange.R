#' Create containers for cross-fitted estimator
#' 
#' This helper function initializes an environment to hold influence-function
#' vectors, point estimates, variance estimates, and other intermediate 
#' quantities for the cross-fitted estimator. 
#' The structure is designed to facilitate storage and retrieval of these 
#' objects across multiple folds and sensitivity parameter values.
#' @noRd
est_exchange_create_containers <- function(gamma, fold) {
  
  e <- new.env(parent = emptyenv())
  e$IF               <- vector(mode = "list", length = length(gamma))
  e$IF_R0            <- vector(mode = "list", length = length(gamma))
  e$est_temp         <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R0_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_temp         <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R0_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$IF_diff          <- vector(mode = "list", length = length(gamma))
  e$IF_R0_diff       <- vector(mode = "list", length = length(gamma))
  e$est_temp_diff    <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R0_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_temp_diff    <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R0_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$vector_R0        <- c()
  e$id_list          <- c()
  
  structure(
    e,
    class = c("est_container", class(e))
  )
}




#' One-step, split sample estimator for E[Y(t)], E[Y(t)|R=0], 
#'   and E[Y(t)-Y0], E[Y(t)-Y0|R=0]

est_psi_exchange <- function(Y, M, R, X, t, trt, gamma, fold, seed, IF_output, 
                             simple_trunc, quant, kernel, method="optim", single_index_method, 
                             use_mave=TRUE, s_t_y=NULL){
  
  if(is.null(s_t_y)){
    s_t_y <- function(y){pnorm((y-60)/25)}
  }
  
  n <- length(t)
  Y[is.na(Y)] <- 0
  
  trt.ind <- as.numeric(t==trt) # create treatment indicator variable
  
  ## set up covariates list and design matrix for modeling
  X_with_T <- cbind(as.factor(t), X)
  colnames(X_with_T)[1] <- "treatment"
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model
  gam.var.M <- paste(gam.variables(X_with_T), collapse = "+") ## gam variables for missing data model
  index.var.Y <- single.index.variables(X)
  X_adjust <- model.matrix(as.formula(paste("~", paste(index.var.Y, collapse = "+"))), data = X)[,-1]
  X_adjust_scale <- scale(X_adjust)
  
  ## empty containers
  containers <- est_exchange_create_containers(gamma, fold)
  if (!simple_trunc){
    containers_trunc <- est_exchange_create_containers(gamma, fold)
  }
  
  ## cross fit
  if(!is.null(seed)){set.seed(seed)}
  indx <- sample(1:n)
  fold_list <- split(1:n, indx %% fold)
  fold_nk_list <- vector(length=fold)
  
  ## containers
  pain_bq_reordered <- c()
  pain_bq_reordered_R0 <- c()
  pain_bq_reordered_R1 <- c()
  pi_R1_l <- c()
  g_l <- c()
  eta_T_R0_l <- c()
  eta_T_R1_l <- c()
  fold_index_pi_R1_l <- c()
  fold_index_g_l <- c()
  fold_index_eta_T_R0_l <- c()
  fold_index_eta_T_R1_l <- c()
  fold_index_pain <- c()
  
  ## compute weights across folds
  for (k in 1:fold){
    ## out-of-fold data
    if(fold>1){
      out_fold_id_list <- setdiff(1:n, fold_list[[k]])
    }else{
      out_fold_id_list <- 1:n
    }
    M_out_fold <- M[out_fold_id_list]
    R_out_fold <- R[out_fold_id_list]
    t_out_fold <- t[out_fold_id_list]
    trt.ind_out_fold <- trt.ind[out_fold_id_list]
    X_out_fold <- X[out_fold_id_list, ]
    X_with_T_out_fold <- X_with_T[out_fold_id_list, ]
    
    trt.ind_out_fold_R0 <- trt.ind_out_fold[which(R_out_fold==0)]
    X_out_fold_R0 <- X_out_fold[which(R_out_fold==0), ]
    t_out_fold_R1 <- t_out_fold[which(R_out_fold==1)]
    X_out_fold_R1 <- X_out_fold[which(R_out_fold==1), ]
    
    M_out_fold_R0 <- M_out_fold[which(R_out_fold==0)]
    X_with_T_out_fold_R0 <- X_with_T_out_fold[which(R_out_fold==0), ]
    M_out_fold_R1 <- M_out_fold[which(R_out_fold==1)]
    X_with_T_out_fold_R1 <- X_with_T_out_fold[which(R_out_fold==1), ]
    
    ## in-fold data
    nk_in_fold <- length(fold_list[[k]])
    fold_nk_list[k] <- nk_in_fold
    R_in_fold <- R[fold_list[[k]]]
    X_in_fold <- X[fold_list[[k]], ]
    X_with_T_in_fold <- X_with_T[fold_list[[k]], ]
    
    X_in_fold_R1 <- X_in_fold[which(R_in_fold==1), ]
    X_in_fold_R0 <- X_in_fold[which(R_in_fold==0), ]
    
    X_with_T_in_fold_R0 <- X_with_T_in_fold[which(R_in_fold==0), ]
    X_with_T_in_fold_R1 <- X_with_T_in_fold[which(R_in_fold==1), ]
    
    ## fit models
    t_R1.fit <- mgcv::gam(as.formula(paste("t_out_fold_R1 ~", gam.var)), data=X_out_fold_R1, family=binomial) ## treatment model
    g.fit <- mgcv::gam(as.formula(paste("R_out_fold ~", gam.var)), data=X_out_fold, family=binomial)
    
    M_R0.fit <- mgcv::gam(as.formula(paste("M_out_fold_R0 ~", gam.var.M)), data=X_with_T_out_fold_R0, family=binomial) ## missing data model
    M_R1.fit <- mgcv::gam(as.formula(paste("M_out_fold_R1 ~", gam.var.M)), data=X_with_T_out_fold_R1, family=binomial) ## missing data model
    
    prop.R1 <- mean(R_in_fold)
    
    ## get predictions for pi
    if(trt==1){
      pi_R1 <- predict(t_R1.fit, newdata=X_in_fold_R1, type="response")  
    }else{
      pi_R1 <- 1-predict(t_R1.fit, newdata=X_in_fold_R1, type="response")  
    }
    pi_R1_l <- c(pi_R1_l, pi_R1)
    fold_index_pi_R1_l <- c(fold_index_pi_R1_l, rep(k, length(pi_R1)))
    
    ## get predictions for g
    g1 <- predict(g.fit, newdata=X_in_fold_R1, type="response")
    g_l <- c(g_l, g1)
    fold_index_g_l <- c(fold_index_g_l, rep(k, length(g1)))
    
    ## get predictions for eta
    eta_T_R0 <- predict(M_R0.fit, newdata=X_with_T_in_fold_R0, type="response")
    eta_T_R1 <- predict(M_R1.fit, newdata=X_with_T_in_fold_R1, type="response")
    eta_T_R0_l <- c(eta_T_R0_l, eta_T_R0)
    eta_T_R1_l <- c(eta_T_R1_l, eta_T_R1)
    fold_index_eta_T_R0_l <- c(fold_index_eta_T_R0_l, rep(k, length(eta_T_R0)))
    fold_index_eta_T_R1_l <- c(fold_index_eta_T_R1_l, rep(k, length(eta_T_R1)))
    
    pain_bq_temp <- c(X_in_fold_R0$womac_bq, X_in_fold_R1$womac_bq)
    pain_bq_R0_temp <- c(X_in_fold_R0$womac_bq, rep(0, length(X_in_fold_R1$womac_bq)))/(1-prop.R1)
    pain_bq_R1_temp <- c(rep(0, length(X_in_fold_R0$womac_bq)), X_in_fold_R1$womac_bq)/prop.R1
    pain_bq_reordered <- c(pain_bq_reordered, pain_bq_temp)
    pain_bq_reordered_R0 <- c(pain_bq_reordered_R0, pain_bq_R0_temp)
    pain_bq_reordered_R1 <- c(pain_bq_reordered_R1, pain_bq_R1_temp)
    fold_index_pain <- c(fold_index_pain, rep(k, length(pain_bq_temp)))
  }
  
  ## weight truncation
  if(!simple_trunc){
    pi_R1_l <- pmin(pmax(pi_R1_l,0.01),0.99)
    g_l <- pmin(pmax(g_l,0.01),0.99)
    eta_T_R0_l <- pmin(pmax(eta_T_R0_l,0.01),0.99)
    eta_T_R1_l <- pmin(pmax(eta_T_R1_l,0.01),0.99)
  }
  pi_R1_weight_l <- 1/pi_R1_l
  g_weight_l <- (1-g_l)/g_l
  eta_T_R0_weight_l <- 1/eta_T_R0_l
  eta_T_R1_weight_l <- 1/eta_T_R1_l
  
  if(simple_trunc){
    pi_R1_weight_l[which(pi_R1_weight_l >= quantile(pi_R1_weight_l, probs = quant))] <- quantile(pi_R1_weight_l, probs = quant)
    g_weight_l[which(g_weight_l>=quantile(g_weight_l, probs=quant))] <- quantile(g_weight_l, probs=quant)
    eta_T_R0_weight_l[which(eta_T_R0_weight_l>=quantile(eta_T_R0_weight_l, probs = quant))] <- quantile(eta_T_R0_weight_l, probs = quant)
    eta_T_R1_weight_l[which(eta_T_R1_weight_l>=quantile(eta_T_R1_weight_l, probs = quant))] <- quantile(eta_T_R1_weight_l, probs = quant)
  }
  
  ## compute outcome model and influence function for each fold
  for (k in 1:fold){
    
    ## out-of-fold data
    if(fold>1){
      out_fold_id_list <- setdiff(1:n, fold_list[[k]])
    }else{
      out_fold_id_list <- 1:n
    }
    Y_out_fold <- Y[out_fold_id_list]
    M_out_fold <- M[out_fold_id_list]
    R_out_fold <- R[out_fold_id_list]
    t_out_fold <- t[out_fold_id_list]
    X_out_fold_adjust_scale <- X_adjust_scale[out_fold_id_list, ]
    
    Y_out_fold_t_R1 <- Y_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    X_out_fold_t_R1_adjust_scale <- X_out_fold_adjust_scale[which(t_out_fold==trt & R_out_fold==1), ]
    M_out_fold_t_R1 <- M_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    
    ## in-fold data
    nk_in_fold <- length(fold_list[[k]])
    Y_in_fold <- Y[fold_list[[k]]]
    M_in_fold <- M[fold_list[[k]]]
    R_in_fold <- R[fold_list[[k]]]
    t_in_fold <- t[fold_list[[k]]]
    trt.ind_in_fold <- trt.ind[fold_list[[k]]]
    X_in_fold_adjust_scale <- X_adjust_scale[fold_list[[k]], ]
    
    X_in_fold_R0_adjust_scale <- X_in_fold_adjust_scale[which(R_in_fold==0), ]
    X_in_fold_R1_adjust_scale <- X_in_fold_adjust_scale[which(R_in_fold==1), ]
    
    M_in_fold_R0 <- M_in_fold[which(R_in_fold==0)]
    M_in_fold_R1 <- M_in_fold[which(R_in_fold==1)]
    
    Y_in_fold_R1 <- Y_in_fold[which(R_in_fold==1)]
    
    trt.ind_in_fold_R1 <- trt.ind_in_fold[which(R_in_fold==1)]
    
    ## get the weights in fold k
    pi_R1_weight <- pi_R1_weight_l[which(fold_index_pi_R1_l==k)]
    g_weight <- g_weight_l[which(fold_index_g_l==k)]
    eta_T_R0_weight <- eta_T_R0_weight_l[which(fold_index_eta_T_R0_l==k)]
    eta_T_R1_weight <- eta_T_R1_weight_l[which(fold_index_eta_T_R1_l==k)]
    
    ## get baseline pain score in fold k
    pain_bq_temp <- pain_bq_reordered[which(fold_index_pain==k)]
    pain_bq_R0_temp <- pain_bq_reordered_R0[which(fold_index_pain==k)]
    pain_bq_R1_temp <- pain_bq_reordered_R1[which(fold_index_pain==k)]
    
    ## fit models
    if(single_index_method=="fixed_bandwidth"){
      if(use_mave){
        requireNamespace('MAVE', quietly = TRUE)
        SDR_t_R1 <- coef(MAVE::mave.compute(X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                            Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], max.dim = 1), 1)
      }else{
        SDR_t_R1 <- cumuSIR(X = X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                            Y = Y_out_fold_t_R1[which(M_out_fold_t_R1==1)])
      }
      fit_t_R1_h <- SIDRnew_fixed_bandwidth(X = X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                            Y = Y_out_fold_t_R1[which(M_out_fold_t_R1==1)],
                                            initial = SDR_t_R1,
                                            kernel = kernel,
                                            method = method,
                                            ids = 1:length(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)]))
    }else if(single_index_method=="fixed_coef"){
      if(use_mave){
        requireNamespace('MAVE', quietly = TRUE)
        SDR_t_R1 <- coef(MAVE::mave.compute(X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                            Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], max.dim = 1), 1)
      }else{
        SDR_t_R1 <- cumuSIR(X=X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ], 
                            Y=Y_out_fold_t_R1[which(M_out_fold_t_R1==1)])
      }
      fit_t_R1_h <- SIDR_Ravinew(X = X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                 Y = Y_out_fold_t_R1[which(M_out_fold_t_R1==1)],
                                 initial=c(1,as.vector(SDR_t_R1[-1]/SDR_t_R1[1])), 
                                 kernel = kernel,
                                 method = method, 
                                 index_ID=1:length(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)]))
    }else if(single_index_method=="norm1coef"){
      fit_t_R1_h  <- fit_SensIAT_single_index_norm1coef_model(X = X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                                              Y = Y_out_fold_t_R1[which(M_out_fold_t_R1==1)],
                                                              ids=1:length(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)]), 
                                                              kernel=kernel, bw.selection="ise", bw.method="optim", use_mave=use_mave)
    }
    
    ## get prediction for X'beta
    X_in_fold_R0_beta_t_R1 <- as.vector(X_in_fold_R0_adjust_scale %*% fit_t_R1_h$coef)
    X_in_fold_R1_beta_t_R1 <- as.vector(X_in_fold_R1_adjust_scale %*% fit_t_R1_h$coef)
    X_out_fold_t_R1_beta_t_R1 <- as.vector(X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ] %*% fit_t_R1_h$coef)
    
    ## compute CDF
    y_t_R1 = sort(unique(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)]))    
    ny_t_R1 = length(y_t_R1) 
    F_t_R1_X_R1 <- NW_new(Xb=X_out_fold_t_R1_beta_t_R1, Y=Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], 
                          xb=X_in_fold_R1_beta_t_R1, y=y_t_R1, h=fit_t_R1_h$bandwidth, 
                          kernel = kernel)
    F_t_R1_X_R0 <- NW_new(Xb=X_out_fold_t_R1_beta_t_R1, Y=Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], 
                          xb=X_in_fold_R0_beta_t_R1, y=y_t_R1, h=fit_t_R1_h$bandwidth, 
                          kernel = kernel)
    
    #fix up cases where CDF is improper by finding closest people with closest X'beta
    i1 = which(apply(F_t_R1_X_R1==0,1,prod)==1)
    i1.closest <- apply(abs(outer(X_in_fold_R1_beta_t_R1[i1], X_in_fold_R1_beta_t_R1[-i1], FUN = "-")), 1, which.min)
    F_t_R1_X_R1[i1, ] <- F_t_R1_X_R1[-i1, ][i1.closest, ]
    
    i1 = which(apply(F_t_R1_X_R0==0,1,prod)==1)
    i1.closest <- apply(abs(outer(X_in_fold_R0_beta_t_R1[i1], X_in_fold_R0_beta_t_R1[-i1], FUN = "-")), 1, which.min)
    F_t_R1_X_R0[i1, ] <- F_t_R1_X_R0[-i1, ][i1.closest, ]
    
    #compute PDF
    dF_t_R1_X_R1 <- F_t_R1_X_R1[, -1, drop = FALSE]-F_t_R1_X_R1[,-ny_t_R1, drop = FALSE]
    dF_t_R1_X_R0 <- F_t_R1_X_R0[, -1, drop = FALSE]-F_t_R1_X_R0[,-ny_t_R1, drop = FALSE]
    
    ## conditional expectation of Y given R=1, T=t and X, for I(R=1)
    mu_Y_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (y_t_R1[-1]+y_t_R1[-ny_t_R1])/2)
    
    ## P(R=1)
    prop.R1 <- mean(R_in_fold)
    
    ## Compute influence function for each gamma_t
    for (g in 1:length(gamma)){
      
      mu_Yexp_t_R1_X_R0 <- c(dF_t_R1_X_R0 %*% (y_t_R1[-1]*exp(gamma[g]*pnorm((y_t_R1[-1]-60)/25))+
                                                 y_t_R1[-ny_t_R1]*exp(gamma[g]*pnorm((y_t_R1[-ny_t_R1]-60)/25)))/2)
      
      mu_exp_t_R1_X_R0 <- c(dF_t_R1_X_R0 %*% (exp(gamma[g]*pnorm((y_t_R1[-1]-60)/25))+
                                                exp(gamma[g]*pnorm((y_t_R1[-ny_t_R1]-60)/25)))/2)
      
      mu_Yexp_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (y_t_R1[-1]*exp(gamma[g]*pnorm((y_t_R1[-1]-60)/25))+
                                                 y_t_R1[-ny_t_R1]*exp(gamma[g]*pnorm((y_t_R1[-ny_t_R1]-60)/25)))/2)
      
      mu_exp_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (exp(gamma[g]*pnorm((y_t_R1[-1]-60)/25))+
                                                exp(gamma[g]*pnorm((y_t_R1[-ny_t_R1]-60)/25)))/2)
      
      ## IF+psi within each fold
      if_temp <- c(M_in_fold_R1*eta_T_R1_weight*trt.ind_in_fold_R1*g_weight*exp(gamma[g]*pnorm((Y_in_fold_R1-60)/25))/mu_exp_t_R1_X_R1*
                     (Y_in_fold_R1-mu_Yexp_t_R1_X_R1/mu_exp_t_R1_X_R1)+ 
                     M_in_fold_R1*eta_T_R1_weight*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1), 
                   M_in_fold_R0*eta_T_R0_weight*mu_Yexp_t_R1_X_R0/mu_exp_t_R1_X_R0)+
        c((1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1, 
          (1-M_in_fold_R0*eta_T_R0_weight)*mu_Yexp_t_R1_X_R0/mu_exp_t_R1_X_R0)
      if_temp_diff <- if_temp-pain_bq_temp
      
      if_R0_temp <- c(M_in_fold_R1*eta_T_R1_weight/(1-prop.R1)*trt.ind_in_fold_R1*g_weight*exp(gamma[g]*pnorm((Y_in_fold_R1-60)/25))/mu_exp_t_R1_X_R1*
                        (Y_in_fold_R1-mu_Yexp_t_R1_X_R1/mu_exp_t_R1_X_R1), 
                      M_in_fold_R0*eta_T_R0_weight/(1-prop.R1)*mu_Yexp_t_R1_X_R0/mu_exp_t_R1_X_R0)+
        c(rep(0, length(M_in_fold_R1)), 
          (1-M_in_fold_R0*eta_T_R0_weight)/(1-prop.R1)*mu_Yexp_t_R1_X_R0/mu_exp_t_R1_X_R0)
      if_R0_temp_diff <- if_R0_temp-pain_bq_R0_temp
      
      vector_mean_R0_temp <- c(rep(0, length(M_in_fold_R1)), rep(1/(1-prop.R1), length(M_in_fold_R0)))*mean(if_R0_temp)
      vector_mean_R0_temp_diff <- c(rep(0, length(M_in_fold_R1)), rep(1/(1-prop.R1), length(M_in_fold_R0)))*mean(if_R0_temp_diff)
      
      if_R1_temp <- c(rep(0, length(M_in_fold_R0)), 
                      M_in_fold_R1*eta_T_R1_weight/prop.R1*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
        c(rep(0, length(M_in_fold_R0)), (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1/prop.R1)
      if_R1_temp_diff <- if_R1_temp-pain_bq_R1_temp
      
      IF[[g]] <- c(IF[[g]], if_temp)
      IF_R0[[g]] <- c(IF_R0[[g]], if_R0_temp)
      IF_diff[[g]] <- c(IF_diff[[g]], if_temp_diff)
      IF_R0_diff[[g]] <- c(IF_R0_diff[[g]], if_R0_temp_diff)
      
      est_temp[k, g] <- mean(if_temp)
      est_R0_temp[k, g] <- mean(if_R0_temp)
      var_temp[k, g] <- sum((if_temp - mean(if_temp))^2)/(nk_in_fold-1)
      var_R0_temp[k, g] <- sum((if_R0_temp - vector_mean_R0_temp)^2)/(nk_in_fold-1)
      est_temp_diff[k, g] <- mean(if_temp_diff)
      est_R0_temp_diff[k, g] <- mean(if_R0_temp_diff)
      var_temp_diff[k, g] <- sum((if_temp_diff - mean(if_temp_diff))^2)/(nk_in_fold-1)
      var_R0_temp_diff[k, g] <- sum((if_R0_temp_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
      
      if(!simple_trunc){
        
        if (eq(max(abs(if_R0_temp)),z=if_R0_temp)>0) {
          trunc <- max(abs(if_R0_temp))
        } else {
          trunc <- uniroot(eq, z=if_R0_temp, interval=c(0.01, max(abs(if_R0_temp))))$root
        }
        if_R0_temp_trunc <- pmin(abs(if_R0_temp),trunc) * sign(if_R0_temp)
        IF_R0_trunc[[g]] <- c(IF_R0_trunc[[g]], if_R0_temp_trunc)
        
        if (eq(max(abs(if_R0_temp_diff)),z=if_R0_temp_diff)>0) {
          trunc <- max(abs(if_R0_temp_diff))
        } else {
          trunc <- uniroot(eq, z=if_R0_temp_diff, interval=c(0.01, max(abs(if_R0_temp_diff))))$root
        }
        if_R0_temp_trunc_diff <- pmin(abs(if_R0_temp_diff),trunc) * sign(if_R0_temp_diff)
        IF_R0_trunc_diff[[g]] <- c(IF_R0_trunc_diff[[g]], if_R0_temp_trunc_diff)
        
        vector_mean_R0_temp <- c(rep(0, length(M_in_fold_R1)), rep(1/(1-prop.R1), length(M_in_fold_R0)))*mean(if_R0_temp_trunc)
        vector_mean_R0_temp_diff <- c(rep(0, length(M_in_fold_R1)), rep(1/(1-prop.R1), length(M_in_fold_R0)))*mean(if_R0_temp_trunc_diff)
        
        if (eq(max(abs(if_R1_temp)),z=if_R1_temp)>0) {
          trunc <- max(abs(if_R1_temp))
        } else {
          trunc <- uniroot(eq, z=if_R1_temp, interval=c(0.01, max(abs(if_R1_temp))))$root
        }
        if_R1_temp_trunc <- pmin(abs(if_R1_temp),trunc) *  sign(if_R1_temp)
        
        if (eq(max(abs(if_R1_temp_diff)),z=if_R1_temp_diff)>0) {
          trunc <- max(abs(if_R1_temp_diff))
        } else {
          trunc <- uniroot(eq, z=if_R1_temp_diff, interval=c(0.01, max(abs(if_R1_temp_diff))))$root
        }
        if_R1_temp_trunc_diff <- pmin(abs(if_R1_temp_diff),trunc) * sign(if_R1_temp_diff)
        
        if_temp_trunc <- if_R1_temp_trunc*prop.R1+if_R0_temp_trunc*(1-prop.R1)
        IF_trunc[[g]] <- c(IF_trunc[[g]], if_temp_trunc)
        if_temp_trunc_diff <- if_R1_temp_trunc_diff*prop.R1+if_R0_temp_trunc_diff*(1-prop.R1)
        IF_trunc_diff[[g]] <- c(IF_trunc_diff[[g]], if_temp_trunc_diff)
        
        est_temp_trunc[k, g] <- mean(if_temp_trunc)
        est_R0_temp_trunc[k, g] <- mean(if_R0_temp_trunc)
        var_temp_trunc[k, g] <- sum((if_temp_trunc - mean(if_temp_trunc))^2)/(nk_in_fold-1)
        var_R0_temp_trunc[k, g] <- sum((if_R0_temp_trunc - vector_mean_R0_temp)^2)/(nk_in_fold-1)
        
        est_temp_trunc_diff[k, g] <- mean(if_temp_trunc_diff)
        est_R0_temp_trunc_diff[k, g] <- mean(if_R0_temp_trunc_diff)
        var_temp_trunc_diff[k, g] <- sum((if_temp_trunc_diff - mean(if_temp_trunc_diff))^2)/(nk_in_fold-1)
        var_R0_temp_trunc_diff[k, g] <- sum((if_R0_temp_trunc_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
        
      }
    }
  }
  
  ## output containers
  r_est <- colMeans(est_temp)
  r_est_R0 <- colMeans(est_R0_temp)
  r_var <- colSums(var_temp)/(fold*n)
  r_var_R0 <- colSums(var_R0_temp)/(fold*n)
  r_lowerCI <- r_est-qnorm(0.975)*sqrt(r_var)
  r_lowerCI_R0 <- r_est_R0-qnorm(0.975)*sqrt(r_var_R0)
  r_upperCI <- r_est+qnorm(0.975)*sqrt(r_var) 
  r_upperCI_R0 <- r_est_R0+qnorm(0.975)*sqrt(r_var_R0) 
  r_est_diff <- colMeans(est_temp_diff)
  r_est_R0_diff <- colMeans(est_R0_temp_diff)
  r_var_diff <- colSums(var_temp_diff)/(fold*n)
  r_var_R0_diff <- colSums(var_R0_temp_diff)/(fold*n)
  r_lowerCI_diff <- r_est_diff-qnorm(0.975)*sqrt(r_var_diff)
  r_lowerCI_R0_diff <- r_est_R0_diff-qnorm(0.975)*sqrt(r_var_R0_diff)
  r_upperCI_diff <- r_est_diff+qnorm(0.975)*sqrt(r_var_diff) 
  r_upperCI_R0_diff <- r_est_R0_diff+qnorm(0.975)*sqrt(r_var_R0_diff)
  if(!simple_trunc){
    r_est_trunc <- colMeans(est_temp_trunc)
    r_est_trunc_R0 <- colMeans(est_R0_temp_trunc)
    r_var_trunc <- colSums(var_temp_trunc)/(fold*n)
    r_var_trunc_R0 <- colSums(var_R0_temp_trunc)/(fold*n)
    r_lowerCI_trunc <- r_est_trunc-qnorm(0.975)*sqrt(r_var_trunc)
    r_lowerCI_trunc_R0 <- r_est_trunc_R0-qnorm(0.975)*sqrt(r_var_trunc_R0)
    r_upperCI_trunc <- r_est_trunc+qnorm(0.975)*sqrt(r_var_trunc) 
    r_upperCI_trunc_R0 <- r_est_trunc_R0+qnorm(0.975)*sqrt(r_var_trunc_R0) 
    r_est_trunc_diff <- colMeans(est_temp_trunc_diff)
    r_est_trunc_R0_diff <- colMeans(est_R0_temp_trunc_diff)
    r_var_trunc_diff <- colSums(var_temp_trunc_diff)/(fold*n)
    r_var_trunc_R0_diff <- colSums(var_R0_temp_trunc_diff)/(fold*n)
    r_lowerCI_trunc_diff <- r_est_trunc_diff-qnorm(0.975)*sqrt(r_var_trunc_diff)
    r_lowerCI_trunc_R0_diff <- r_est_trunc_R0_diff-qnorm(0.975)*sqrt(r_var_trunc_R0_diff)
    r_upperCI_trunc_diff <- r_est_trunc_diff+qnorm(0.975)*sqrt(r_var_trunc_diff) 
    r_upperCI_trunc_R0_diff <- r_est_trunc_R0_diff+qnorm(0.975)*sqrt(r_var_trunc_R0_diff) 
    
  }
  
  ## output final results
  if(IF_output){
    if(simple_trunc){
      result <- list(est=r_est, est_R0=r_est_R0, 
                     var=r_var, var_R0=r_var_R0, 
                     lowerCI=r_lowerCI, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R0_diff=r_var_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     IF=IF, IF_R0=IF_R0, 
                     IF_diff=IF_diff, IF_R0_diff=IF_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0)
    }else{
      result <- list(est=r_est, est_R0=r_est_R0, 
                     est_trunc=r_est_trunc, est_trunc_R0=r_est_trunc_R0,
                     var=r_var, var_R0=r_var_R0,
                     var_trunc=r_var_trunc, var_trunc_R0=r_var_trunc_R0,
                     lowerCI=r_lowerCI, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R0_diff=r_est_R0_diff,
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     var_diff=r_var_diff, var_R0_diff=r_var_R0_diff,
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                     upperCI_diff=r_upperCI_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff, 
                     upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff, 
                     IF=IF, IF_R0=IF_R0,
                     IF_trunc=IF_trunc, IF_trunc_R0=IF_R0_trunc,
                     IF_diff=IF_diff, IF_R0_diff=IF_R0_diff,
                     IF_trunc_diff=IF_trunc_diff, IF_trunc_R0_diff=IF_R0_trunc_diff,
                     pain_bq_reordered=pain_bq_reordered, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0)
    }
  }else{
    if(simple_trunc){
      result <- list(est=r_est, est_R0=r_est_R0, 
                     var=r_var, var_R0=r_var_R0, 
                     lowerCI=r_lowerCI, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R0_diff=r_var_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0)
    }else{
      result <- list(est=r_est, est_R0=r_est_R0, 
                     est_trunc=r_est_trunc, est_trunc_R0=r_est_trunc_R0,
                     var=r_var, var_R0=r_var_R0,
                     var_trunc=r_var_trunc, var_trunc_R0=r_var_trunc_R0,
                     lowerCI=r_lowerCI, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R0_diff=r_est_R0_diff,
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     var_diff=r_var_diff, var_R0_diff=r_var_R0_diff,
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                     upperCI_diff=r_upperCI_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff, 
                     upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0)
    }
  }
  
}


est_psi <- function(Y, M, R, X, t, trt, gamma, fold, seed, IF_output, 
                    simple_trunc, quant, kernel, method, single_index_method, 
                    use_mave=TRUE, coef_g.fit=NULL, coef_t_R0.fit=NULL, coef_t_R1.fit=NULL, 
                    coef_M_R0.fit=NULL, coef_M_R1.fit=NULL){
  
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
  IF <- vector(mode = "list", length = length(gamma))
  IF_R1 <- vector(mode = "list", length = length(gamma))
  IF_R0 <- vector(mode = "list", length = length(gamma))
  est_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  est_R1_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  est_R0_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  var_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  var_R1_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  var_R0_temp <- matrix(0, nrow=fold, ncol=length(gamma))
  IF_diff <- vector(mode = "list", length = length(gamma))
  IF_R1_diff <- vector(mode = "list", length = length(gamma))
  IF_R0_diff <- vector(mode = "list", length = length(gamma))
  est_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  est_R1_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  est_R0_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  var_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  var_R1_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  var_R0_temp_diff <- matrix(0, nrow=fold, ncol=length(gamma))  
  vector_R1 <- vector_R0 <- c()
  id_list <- c()
  if(!simple_trunc){
    IF_trunc <- vector(mode = "list", length = length(gamma))
    IF_R1_trunc <- vector(mode = "list", length = length(gamma))
    IF_R0_trunc <- vector(mode = "list", length = length(gamma))
    est_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    est_R1_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    est_R0_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    var_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    var_R1_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    var_R0_temp_trunc <- matrix(0, nrow=fold, ncol=length(gamma))
    IF_trunc_diff <- vector(mode = "list", length = length(gamma))
    IF_R1_trunc_diff <- vector(mode = "list", length = length(gamma))
    IF_R0_trunc_diff <- vector(mode = "list", length = length(gamma))
    est_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
    est_R1_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
    est_R0_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
    var_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
    var_R1_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
    var_R0_temp_trunc_diff <- matrix(0, nrow=fold, ncol=length(gamma))
  
  }
  
  ## cross fit
  if(!is.null(seed)){set.seed(seed)}
  indx <- sample(1:n)
  fold_list <- split(1:n, indx %% fold)
  fold_nk_list <- vector(length=fold)
  
  ## containers
  pain_bq_reordered <- c()
  pain_bq_reordered_R1 <- c()
  pain_bq_reordered_R0 <- c()
  pi_R0_l <- c()
  pi_R1_l <- c()
  eta_t_R0_l <- c()
  eta_t0_R0_l <- c()
  eta_T_R1_l <- c()
  fold_index_pi_R0_l <- c()
  fold_index_pi_R1_l <- c()
  fold_index_eta_t_R0_l <- c()
  fold_index_eta_t0_R0_l <- c()
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
    
    t_out_fold_R0 <- t_out_fold[which(R_out_fold==0)]
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
    t_in_fold <- t[fold_list[[k]]]
    X_in_fold <- X[fold_list[[k]], ]
    X_with_T_in_fold <- X_with_T[fold_list[[k]], ]
    
    X_in_fold_t_R0 <- X_in_fold[which(t_in_fold==trt & R_in_fold==0), ]
    X_in_fold_t0_R0 <- X_in_fold[which(t_in_fold!=trt & R_in_fold==0), ]
    X_in_fold_R1 <- X_in_fold[which(R_in_fold==1), ]
    
    X_with_T_in_fold_t_R0 <- X_with_T_in_fold[which(t_in_fold==trt & R_in_fold==0), ]
    X_with_T_in_fold_t0_R0 <- X_with_T_in_fold[which(t_in_fold!=trt & R_in_fold==0), ]
    X_with_T_in_fold_R1 <- X_with_T_in_fold[which(R_in_fold==1), ]
    
    ## fit models
    t_R0.fit <- mgcv::gam(as.formula(paste("t_out_fold_R0 ~", gam.var)), data=X_out_fold_R0, 
                          family=binomial, start=coef_t_R0.fit) ## treatment model
    t_R1.fit <- mgcv::gam(as.formula(paste("t_out_fold_R1 ~", gam.var)), data=X_out_fold_R1, 
                          family=binomial, start=coef_t_R1.fit) ## treatment model
    M_R0.fit <- mgcv::gam(as.formula(paste("M_out_fold_R0 ~", gam.var.M)), data=X_with_T_out_fold_R0, 
                          family=binomial, start=coef_M_R0.fit) ## missing data model
    M_R1.fit <- mgcv::gam(as.formula(paste("M_out_fold_R1 ~", gam.var.M)), data=X_with_T_out_fold_R1, 
                          family=binomial, start=coef_M_R1.fit) ## missing data model
    
    prop.R1 <- mean(R_in_fold)
    
    ## get predictions for pi
    if(trt==1){
      pi_R0 <- predict(t_R0.fit, newdata=X_in_fold_t_R0, type="response") 
      pi_R1 <- predict(t_R1.fit, newdata=X_in_fold_R1, type="response")  
    }else{
      pi_R0 <- 1-predict(t_R0.fit, newdata=X_in_fold_t_R0, type="response") 
      pi_R1 <- 1-predict(t_R1.fit, newdata=X_in_fold_R1, type="response")  
    }
    pi_R0_l <- c(pi_R0_l, pi_R0)
    pi_R1_l <- c(pi_R1_l, pi_R1)
    fold_index_pi_R0_l <- c(fold_index_pi_R0_l, rep(k, length(pi_R0)))
    fold_index_pi_R1_l <- c(fold_index_pi_R1_l, rep(k, length(pi_R1)))
    
    ## get predictions for eta
    eta_t_R0 <- predict(M_R0.fit, newdata=mutate(X_with_T_in_fold_t_R0, treatment=trt), type="response")
    eta_t0_R0 <- predict(M_R0.fit, newdata=mutate(X_with_T_in_fold_t0_R0, treatment=1-trt), type="response")
    eta_T_R1 <- predict(M_R1.fit, newdata=X_with_T_in_fold_R1, type="response")
    eta_t_R0_l <- c(eta_t_R0_l, eta_t_R0)
    eta_t0_R0_l <- c(eta_t0_R0_l, eta_t0_R0)
    eta_T_R1_l <- c(eta_T_R1_l, eta_T_R1)
    fold_index_eta_t_R0_l <- c(fold_index_eta_t_R0_l, rep(k, length(eta_t_R0)))
    fold_index_eta_t0_R0_l <- c(fold_index_eta_t0_R0_l, rep(k, length(eta_t0_R0)))
    fold_index_eta_T_R1_l <- c(fold_index_eta_T_R1_l, rep(k, length(eta_T_R1)))
    
    if(trt==1){
      pain_bq_temp <- c(X_in_fold_t_R0$pain_bq, X_in_fold_t0_R0$pain_bq, X_in_fold_R1$pain_bq)
      pain_bq_R1_temp <- c(rep(0, length(X_in_fold_t_R0$pain_bq)+length(X_in_fold_t0_R0$pain_bq)), X_in_fold_R1$pain_bq)/prop.R1
      pain_bq_R0_temp <- c(X_in_fold_t_R0$pain_bq, X_in_fold_t0_R0$pain_bq, rep(0, length(X_in_fold_R1$pain_bq)))/(1-prop.R1)
    }else{
      pain_bq_temp <- c(X_in_fold_t0_R0$pain_bq, X_in_fold_t_R0$pain_bq, X_in_fold_R1$pain_bq)
      pain_bq_R1_temp <- c(rep(0, length(X_in_fold_t0_R0$pain_bq)+length(X_in_fold_t_R0$pain_bq)), X_in_fold_R1$pain_bq)/prop.R1
      pain_bq_R0_temp <- c(X_in_fold_t0_R0$pain_bq, X_in_fold_t_R0$pain_bq, rep(0, length(X_in_fold_R1$pain_bq)))/(1-prop.R1)
    }
    pain_bq_reordered <- c(pain_bq_reordered, pain_bq_temp)
    pain_bq_reordered_R1 <- c(pain_bq_reordered_R1, pain_bq_R1_temp)
    pain_bq_reordered_R0 <- c(pain_bq_reordered_R0, pain_bq_R0_temp)
    fold_index_pain <- c(fold_index_pain, rep(k, length(pain_bq_temp)))
  }
  
  ## weight truncation
  if(!simple_trunc){
  pi_R0_l <- pmin(pmax(pi_R0_l,0.01),0.99)
  pi_R1_l <- pmin(pmax(pi_R1_l,0.01),0.99)
  eta_t_R0_l <- pmin(pmax(eta_t_R0_l,0.01),0.99)
  eta_t0_R0_l <- pmin(pmax(eta_t0_R0_l,0.01),0.99)
  eta_T_R1_l <- pmin(pmax(eta_T_R1_l,0.01),0.99)
  }
  pi_R0_weight_l <- 1/pi_R0_l
  pi_R1_weight_l <- 1/pi_R1_l
  eta_t_R0_weight_l <- 1/eta_t_R0_l
  eta_t0_R0_weight_l <- 1/eta_t0_R0_l
  eta_T_R1_weight_l <- 1/eta_T_R1_l
  if(simple_trunc){
  pi_R0_weight_l[which(pi_R0_weight_l >= quantile(pi_R0_weight_l, probs = quant))] <- quantile(pi_R0_weight_l, probs = quant)
  pi_R1_weight_l[which(pi_R1_weight_l >= quantile(pi_R1_weight_l, probs = quant))] <- quantile(pi_R1_weight_l, probs = quant)
  quant_eta_T_R0_weight <- quantile(c(eta_t_R0_weight_l, eta_t0_R0_weight_l), probs = quant)
  eta_t_R0_weight_l[which(eta_t_R0_weight_l>=quant_eta_T_R0_weight)] <- quant_eta_T_R0_weight
  eta_t0_R0_weight_l[which(eta_t0_R0_weight_l>=quant_eta_T_R0_weight)] <- quant_eta_T_R0_weight
  eta_T_R1_weight_l[which(eta_T_R1_weight_l >= quantile(eta_T_R1_weight_l, probs = quant))] <- quantile(eta_T_R1_weight_l, probs = quant)
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
    X_out_fold <- X[out_fold_id_list, ]
    X_out_fold_adjust_scale <- X_adjust_scale[out_fold_id_list, ]

    Y_out_fold_t_R0 <- Y_out_fold[which(t_out_fold==trt & R_out_fold==0)]
    Y_out_fold_t_R1 <- Y_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    X_out_fold_t_R0_adjust_scale <- X_out_fold_adjust_scale[which(t_out_fold==trt & R_out_fold==0), ]
    X_out_fold_t_R1_adjust_scale <- X_out_fold_adjust_scale[which(t_out_fold==trt & R_out_fold==1), ]
    M_out_fold_t_R0 <- M_out_fold[which(t_out_fold==trt & R_out_fold==0)]
    M_out_fold_t_R1 <- M_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    #X_out_fold_t_R0 <- X_out_fold[which(t_out_fold==trt & R_out_fold==0), ]
    #X_out_fold_t_R1 <- X_out_fold[which(t_out_fold==trt & R_out_fold==1), ]
    #X_out_fold_t_R0_adjust_scale <- model.matrix(as.formula(paste("~", paste(index.var.Y, collapse = "+"))), data = X)[,-1]

    ## in-fold data
    nk_in_fold <- length(fold_list[[k]])
    Y_in_fold <- Y[fold_list[[k]]]
    M_in_fold <- M[fold_list[[k]]]
    R_in_fold <- R[fold_list[[k]]]
    t_in_fold <- t[fold_list[[k]]]
    trt.ind_in_fold <- trt.ind[fold_list[[k]]]
    X_in_fold_adjust_scale <- X_adjust_scale[fold_list[[k]], ]

    X_in_fold_t_R0_adjust_scale <- X_in_fold_adjust_scale[which(t_in_fold==trt & R_in_fold==0), ]
    X_in_fold_t0_R0_adjust_scale <- X_in_fold_adjust_scale[which(t_in_fold!=trt & R_in_fold==0), ]
    X_in_fold_R1_adjust_scale <- X_in_fold_adjust_scale[which(R_in_fold==1), ]
    
    M_in_fold_t_R0 <- M_in_fold[which(t_in_fold==trt & R_in_fold==0)]
    M_in_fold_t0_R0 <- M_in_fold[which(t_in_fold!=trt & R_in_fold==0)]
    M_in_fold_R1 <- M_in_fold[which(R_in_fold==1)]

    Y_in_fold_t_R0 <- Y_in_fold[which(t_in_fold==trt & R_in_fold==0)]
    Y_in_fold_R1 <- Y_in_fold[which(R_in_fold==1)]
    
    trt.ind_in_fold_R1 <- trt.ind_in_fold[which(R_in_fold==1)]
    
    id_list <- c(id_list, fold_list[[k]])
    
    ## get the weights in fold k
    pi_R0_weight <- pi_R0_weight_l[which(fold_index_pi_R0_l==k)]
    pi_R1_weight <- pi_R1_weight_l[which(fold_index_pi_R1_l==k)]
    eta_t_R0_weight <- eta_t_R0_weight_l[which(fold_index_eta_t_R0_l==k)]
    eta_t0_R0_weight <- eta_t0_R0_weight_l[which(fold_index_eta_t0_R0_l==k)]
    eta_T_R1_weight <- eta_T_R1_weight_l[which(fold_index_eta_T_R1_l==k)]
    
    ## get baseline pain score in fold k
    pain_bq_temp <- pain_bq_reordered[which(fold_index_pain==k)]
    pain_bq_R1_temp <- pain_bq_reordered_R1[which(fold_index_pain==k)]
    pain_bq_R0_temp <- pain_bq_reordered_R0[which(fold_index_pain==k)]
    
    ## fit models
    if(single_index_method=="fixed_bandwidth"){
      if(use_mave){
        requireNamespace('MAVE', quietly = TRUE)
        SDR_t_R0 <- coef(MAVE::mave.compute(X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                                            Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], max.dim = 1), 1)
      }else{
        SDR_t_R0 <- cumuSIR(X = X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                            Y = Y_out_fold_t_R0[which(M_out_fold_t_R0==1)])
      }
      fit_t_R0_h <- SIDRnew_fixed_bandwidth(X = X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                                            Y = Y_out_fold_t_R0[which(M_out_fold_t_R0==1)],
                                            initial = SDR_t_R0,
                                            kernel = kernel,
                                            method = method,
                                            ids = 1:length(Y_out_fold_t_R0[which(M_out_fold_t_R0==1)]))

    }else if(single_index_method=="fixed_coef"){
    if(use_mave){
      requireNamespace('MAVE', quietly = TRUE)
      SDR_t_R0 <- coef(MAVE::mave.compute(X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                                          Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], max.dim = 1), 1)
    }else{
      SDR_t_R0 <- cumuSIR(X = X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                          Y = Y_out_fold_t_R0[which(M_out_fold_t_R0==1)])
    }
    fit_t_R0_h <- SIDR_Ravinew(X = X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                               Y = Y_out_fold_t_R0[which(M_out_fold_t_R0==1)],
                               initial=c(1,as.vector(SDR_t_R0[-1]/SDR_t_R0[1])), 
                               kernel = kernel,
                               method = method, 
                               index_ID=1:length(Y_out_fold_t_R0[which(M_out_fold_t_R0==1)]))
    }else if(single_index_method=="norm1coef"){
    fit_t_R0_h  <- fit_SensIAT_single_index_norm1coef_model(X = X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ],
                                                            Y = Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], 
                                                            ids=1:length(Y_out_fold_t_R0[which(M_out_fold_t_R0==1)]), 
                                                            kernel=kernel, bw.selection="ise", bw.method="optim", use_mave=use_mave)
    }
    
    if(single_index_method=="fixed_bandwidth"){
    if(use_mave){
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
        SDR_t_R1 <- coef(MAVE::mave.compute(X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                                            Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], max.dim = 1), 1)
      }else{
        SDR_t_R1 <- cumuSIR(X = X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ],
                            Y = Y_out_fold_t_R1[which(M_out_fold_t_R1==1)])
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
    X_in_fold_t_R0_beta_t_R0 <- as.vector(X_in_fold_t_R0_adjust_scale %*% fit_t_R0_h$coef)
    X_in_fold_t0_R0_beta_t_R0 <- as.vector(X_in_fold_t0_R0_adjust_scale %*% fit_t_R0_h$coef)
    X_in_fold_R1_beta_t_R1 <- as.vector(X_in_fold_R1_adjust_scale %*% fit_t_R1_h$coef)
    X_out_fold_t_R0_beta_t_R0 <- as.vector(X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ] %*% fit_t_R0_h$coef)
    X_out_fold_t_R1_beta_t_R1 <- as.vector(X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ] %*% fit_t_R1_h$coef)
    
    ## compute CDF
    y_t_R0 = sort(unique(Y_out_fold_t_R0[which(M_out_fold_t_R0==1)]))    
    ny_t_R0 = length(y_t_R0) 
    F_t_R0_X_t_R0 <- NW_new(Xb=X_out_fold_t_R0_beta_t_R0, Y=Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], 
                            xb=X_in_fold_t_R0_beta_t_R0, y=y_t_R0, h=fit_t_R0_h$bandwidth, 
                            kernel = kernel)
    F_t_R0_X_t0_R0 <- NW_new(Xb=X_out_fold_t_R0_beta_t_R0, Y=Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], 
                             xb=X_in_fold_t0_R0_beta_t_R0, y=y_t_R0, h=fit_t_R0_h$bandwidth, 
                             kernel = kernel)
    
    y_t_R1 = sort(unique(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)]))    
    ny_t_R1 = length(y_t_R1) 
    F_t_R1_X_R1 <- NW_new(Xb=X_out_fold_t_R1_beta_t_R1, Y=Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], 
                          xb=X_in_fold_R1_beta_t_R1, y=y_t_R1, h=fit_t_R1_h$bandwidth, 
                          kernel = kernel)
    
    #fix up cases where CDF is improper by finding closest people with closest X'beta
    i1 = which(apply(F_t_R0_X_t_R0==0,1,prod)==1)
    i1.closest <- apply(abs(outer(X_in_fold_t_R0_beta_t_R0[i1], X_in_fold_t_R0_beta_t_R0[-i1], FUN = "-")), 1, which.min)
    F_t_R0_X_t_R0[i1, ] <- F_t_R0_X_t_R0[-i1, ,drop = FALSE][i1.closest, ]

    i1 = which(apply(F_t_R0_X_t0_R0==0,1,prod)==1)
    i1.closest <- apply(abs(outer(X_in_fold_t0_R0_beta_t_R0[i1], X_in_fold_t0_R0_beta_t_R0[-i1], FUN = "-")), 1, which.min)
    F_t_R0_X_t0_R0[i1, ] <- F_t_R0_X_t0_R0[-i1, ,drop = FALSE][i1.closest, ]

    i1 = which(apply(F_t_R1_X_R1==0,1,prod)==1)
    i1.closest <- apply(abs(outer(X_in_fold_R1_beta_t_R1[i1], X_in_fold_R1_beta_t_R1[-i1], FUN = "-")), 1, which.min)
    F_t_R1_X_R1[i1, ] <- F_t_R1_X_R1[-i1, ,drop = FALSE][i1.closest, ]

    #compute PDF
    dF_t_R0_X_t_R0 <- F_t_R0_X_t_R0-cbind(rep(0, dim(F_t_R0_X_t_R0)[1]), F_t_R0_X_t_R0[,-ny_t_R0, drop = FALSE])
    dF_t_R0_X_t0_R0 <- F_t_R0_X_t0_R0-cbind(rep(0, dim(F_t_R0_X_t0_R0)[1]), F_t_R0_X_t0_R0[,-ny_t_R0, drop = FALSE])
    dF_t_R1_X_R1 <- F_t_R1_X_R1-cbind(rep(0, dim(F_t_R1_X_R1)[1]), F_t_R1_X_R1[,-ny_t_R1, drop = FALSE])
    
    ## conditional expectation of Y given R=0, T=t and X, for I(R=0, T=t)
    mu_Y_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0+c(0, y_t_R0[-ny_t_R0]))/2)
    
    ## conditional expectation of Y given R=1, T=t and X, for I(R=1)
    mu_Y_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (y_t_R1+c(0, y_t_R1[-ny_t_R1]))/2)
    
    ## P(R=1)
    prop.R1 <- mean(R_in_fold)
    
    ## vector_R1 and vector_R0
    vector_R1 <- c(vector_R1, c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                                rep(1/prop.R1, length(M_in_fold_R1))))
    vector_R0 <- c(vector_R0, c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                                rep(0, length(M_in_fold_R1))))
    
    ## Compute influence function for each gamma_t
    for (g in 1:length(gamma)){
      
    ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=t)
    mu_Yexp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0*exp(gamma[g]*pnorm((y_t_R0-60)/25))+
                                                     c(0, y_t_R0[-ny_t_R0])*exp(gamma[g]*c(0, pnorm((y_t_R0-60)/25)[-ny_t_R0])))/2)
      
    ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=t)
    mu_exp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (exp(gamma[g]*pnorm((y_t_R0-60)/25))+
                                                  exp(gamma[g]*c(0, pnorm((y_t_R0-60)/25)[-ny_t_R0])))/2)
    
    ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=1-t)
    mu_Yexp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (y_t_R0*exp(gamma[g]*pnorm((y_t_R0-60)/25))+
                                                     c(0, y_t_R0[-ny_t_R0])*exp(gamma[g]*c(0, pnorm((y_t_R0-60)/25)[-ny_t_R0])))/2)
    
    ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=1-t)
    mu_exp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (exp(gamma[g]*pnorm((y_t_R0-60)/25))+
                                                    exp(gamma[g]*c(0, pnorm((y_t_R0-60)/25)[-ny_t_R0])))/2)
    
    ## IF+psi within each fold
    if(trt==1){
      if_temp <- c(M_in_fold_t_R0*eta_t_R0_weight*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*pnorm((Y_in_fold_t_R0-60)/25))/mu_exp_t_R0_X_t_R0*
                                              (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                   M_in_fold_t0_R0*eta_t0_R0_weight*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                   M_in_fold_R1*eta_T_R1_weight*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
        c((1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0, 
          (1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0,
          (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1)
      if_temp_diff <- if_temp-pain_bq_temp
      
      if_R0_temp <- c(M_in_fold_t_R0*eta_t_R0_weight/(1-prop.R1)*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*pnorm((Y_in_fold_t_R0-60)/25))/mu_exp_t_R0_X_t_R0*
                                                 (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                      M_in_fold_t0_R0*eta_t0_R0_weight/(1-prop.R1)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                      rep(0, length(M_in_fold_R1)))+
        c((1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0/(1-prop.R1), 
          (1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/(mu_exp_t_R0_X_t0_R0*(1-prop.R1)),
          rep(0, length(M_in_fold_R1)))
      if_R0_temp_diff <- if_R0_temp-pain_bq_R0_temp
    }else{
      if_temp <- c(M_in_fold_t0_R0*eta_t0_R0_weight*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                   M_in_fold_t_R0*eta_t_R0_weight*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*pnorm((Y_in_fold_t_R0-60)/25))/mu_exp_t_R0_X_t_R0*
                                              (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                   M_in_fold_R1*eta_T_R1_weight*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
        c((1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0,
          (1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0, 
          (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1)
      if_temp_diff <- if_temp-pain_bq_temp
      
      if_R0_temp  <- c(M_in_fold_t0_R0*eta_t0_R0_weight/(1-prop.R1)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                         M_in_fold_t_R0*eta_t_R0_weight/(1-prop.R1)*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*pnorm((Y_in_fold_t_R0-60)/25))/mu_exp_t_R0_X_t_R0*
                                                    (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                         rep(0, length(M_in_fold_R1)))+
        c((1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/(mu_exp_t_R0_X_t0_R0*(1-prop.R1)),
          (1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0/(1-prop.R1), 
          rep(0, length(M_in_fold_R1)))
      if_R0_temp_diff <- if_R0_temp-pain_bq_R0_temp
    }
    
    vector_mean_R0_temp <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                             rep(0, length(M_in_fold_R1)))*mean(if_R0_temp)
    vector_mean_R0_temp_diff <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                                  rep(0, length(M_in_fold_R1)))*mean(if_R0_temp_diff)
    
    if_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                    M_in_fold_R1*eta_T_R1_weight/prop.R1*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
      c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1/prop.R1)
    if_R1_temp_diff <- if_R1_temp-pain_bq_R1_temp
    
    vector_mean_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                             rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp)
    vector_mean_R1_temp_diff <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                                  rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp_diff)
    
    IF[[g]] <- c(IF[[g]], if_temp)
    IF_R1[[g]] <- c(IF_R1[[g]], if_R1_temp)
    IF_R0[[g]] <- c(IF_R0[[g]], if_R0_temp)
    IF_diff[[g]] <- c(IF_diff[[g]], if_temp_diff)
    IF_R1_diff[[g]] <- c(IF_R1_diff[[g]], if_R1_temp_diff)
    IF_R0_diff[[g]] <- c(IF_R0_diff[[g]], if_R0_temp_diff)
    
    est_temp[k, g] <- mean(if_temp)
    est_R1_temp[k, g] <- mean(if_R1_temp)
    est_R0_temp[k, g] <- mean(if_R0_temp)
    var_temp[k, g] <- sum((if_temp - mean(if_temp))^2)/(nk_in_fold-1)
    var_R1_temp[k, g] <- sum((if_R1_temp - vector_mean_R1_temp)^2)/(nk_in_fold-1)
    var_R0_temp[k, g] <- sum((if_R0_temp - vector_mean_R0_temp)^2)/(nk_in_fold-1)
    est_temp_diff[k, g] <- mean(if_temp_diff)
    est_R1_temp_diff[k, g] <- mean(if_R1_temp_diff)
    est_R0_temp_diff[k, g] <- mean(if_R0_temp_diff)
    var_temp_diff[k, g] <- sum((if_temp_diff - mean(if_temp_diff))^2)/(nk_in_fold-1)
    var_R1_temp_diff[k, g] <- sum((if_R1_temp_diff - vector_mean_R1_temp_diff)^2)/(nk_in_fold-1)
    var_R0_temp_diff[k, g] <- sum((if_R0_temp_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
    
    if(!simple_trunc){
      
      if (eq(max(abs(if_R1_temp)),z=if_R1_temp)>0) {
        trunc <- max(abs(if_R1_temp))
      } else {
        trunc <- uniroot(eq, z=if_R1_temp, interval=c(0.01, max(abs(if_R1_temp))))$root
      }
      if_R1_temp_trunc <- pmin(abs(if_R1_temp),trunc) *  sign(if_R1_temp)
      IF_R1_trunc[[g]] <- c(IF_R1_trunc[[g]], if_R1_temp_trunc)
      
      if (eq(max(abs(if_R1_temp_diff)),z=if_R1_temp_diff)>0) {
        trunc <- max(abs(if_R1_temp_diff))
      } else {
        trunc <- uniroot(eq, z=if_R1_temp_diff, interval=c(0.01, max(abs(if_R1_temp_diff))))$root
      }
      if_R1_temp_trunc_diff <- pmin(abs(if_R1_temp_diff),trunc) * sign(if_R1_temp_diff)
      IF_R1_trunc_diff[[g]] <- c(IF_R1_trunc_diff[[g]], if_R1_temp_trunc_diff)
      
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
      
      if_temp_trunc <- if_R1_temp_trunc*prop.R1+if_R0_temp_trunc*(1-prop.R1)
      IF_trunc[[g]] <- c(IF_trunc[[g]], if_temp_trunc)
      if_temp_trunc_diff <- if_R1_temp_trunc_diff*prop.R1+if_R0_temp_trunc_diff*(1-prop.R1)
      IF_trunc_diff[[g]] <- c(IF_trunc_diff[[g]], if_temp_trunc_diff)
      
      vector_mean_R0_temp <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                               rep(0, length(M_in_fold_R1)))*mean(if_R0_temp_trunc)
      vector_mean_R0_temp_diff <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                                    rep(0, length(M_in_fold_R1)))*mean(if_R0_temp_trunc_diff)
      
      vector_mean_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                               rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp_trunc)
      vector_mean_R1_temp_diff <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                                    rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp_trunc_diff)
      # 
      est_temp_trunc[k, g] <- mean(if_temp_trunc)
      est_R1_temp_trunc[k, g] <- mean(if_R1_temp_trunc)
      est_R0_temp_trunc[k, g] <- mean(if_R0_temp_trunc)
      var_temp_trunc[k, g] <- sum((if_temp_trunc - mean(if_temp_trunc))^2)/(nk_in_fold-1)
      var_R1_temp_trunc[k, g] <- sum((if_R1_temp_trunc - vector_mean_R1_temp)^2)/(nk_in_fold-1)
      var_R0_temp_trunc[k, g] <- sum((if_R0_temp_trunc - vector_mean_R0_temp)^2)/(nk_in_fold-1)
      
      est_temp_trunc_diff[k, g] <- mean(if_temp_trunc_diff)
      est_R1_temp_trunc_diff[k, g] <- mean(if_R1_temp_trunc_diff)
      est_R0_temp_trunc_diff[k, g] <- mean(if_R0_temp_trunc_diff)
      var_temp_trunc_diff[k, g] <- sum((if_temp_trunc_diff - mean(if_temp_trunc_diff))^2)/(nk_in_fold-1)
      var_R1_temp_trunc_diff[k, g] <- sum((if_R1_temp_trunc_diff - vector_mean_R1_temp_diff)^2)/(nk_in_fold-1)
      var_R0_temp_trunc_diff[k, g] <- sum((if_R0_temp_trunc_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
      
    }
    }
  }
  
  ## output containers
  r_est <- colMeans(est_temp)
  r_est_R1 <- colMeans(est_R1_temp)
  r_est_R0 <- colMeans(est_R0_temp)
  r_var <- colSums(var_temp)/(fold*n)
  r_var_R1 <- colSums(var_R1_temp)/(fold*n)
  r_var_R0 <- colSums(var_R0_temp)/(fold*n)
  r_lowerCI <- r_est-qnorm(0.975)*sqrt(r_var)
  r_lowerCI_R1 <- r_est_R1-qnorm(0.975)*sqrt(r_var_R1)
  r_lowerCI_R0 <- r_est_R0-qnorm(0.975)*sqrt(r_var_R0)
  r_upperCI <- r_est+qnorm(0.975)*sqrt(r_var) 
  r_upperCI_R1 <- r_est_R1+qnorm(0.975)*sqrt(r_var_R1) 
  r_upperCI_R0 <- r_est_R0+qnorm(0.975)*sqrt(r_var_R0) 
  r_est_diff <- colMeans(est_temp_diff)
  r_est_R1_diff <- colMeans(est_R1_temp_diff)
  r_est_R0_diff <- colMeans(est_R0_temp_diff)
  r_var_diff <- colSums(var_temp_diff)/(fold*n)
  r_var_R1_diff <- colSums(var_R1_temp_diff)/(fold*n)
  r_var_R0_diff <- colSums(var_R0_temp_diff)/(fold*n)
  r_lowerCI_diff <- r_est_diff-qnorm(0.975)*sqrt(r_var_diff)
  r_lowerCI_R1_diff <- r_est_R1_diff-qnorm(0.975)*sqrt(r_var_R1_diff)
  r_lowerCI_R0_diff <- r_est_R0_diff-qnorm(0.975)*sqrt(r_var_R0_diff)
  r_upperCI_diff <- r_est_diff+qnorm(0.975)*sqrt(r_var_diff) 
  r_upperCI_R1_diff <- r_est_R1_diff+qnorm(0.975)*sqrt(r_var_R1_diff) 
  r_upperCI_R0_diff <- r_est_R0_diff+qnorm(0.975)*sqrt(r_var_R0_diff)
  if(!simple_trunc){
    r_est_trunc <- colMeans(est_temp_trunc)
    r_est_trunc_R1 <- colMeans(est_R1_temp_trunc)
    r_est_trunc_R0 <- colMeans(est_R0_temp_trunc)
    r_var_trunc <- colSums(var_temp_trunc)/(fold*n)
    r_var_trunc_R1 <- colSums(var_R1_temp_trunc)/(fold*n)
    r_var_trunc_R0 <- colSums(var_R0_temp_trunc)/(fold*n)
    r_lowerCI_trunc <- r_est_trunc-qnorm(0.975)*sqrt(r_var_trunc)
    r_lowerCI_trunc_R1 <- r_est_trunc_R1-qnorm(0.975)*sqrt(r_var_trunc_R1)
    r_lowerCI_trunc_R0 <- r_est_trunc_R0-qnorm(0.975)*sqrt(r_var_trunc_R0)
    r_upperCI_trunc <- r_est_trunc+qnorm(0.975)*sqrt(r_var_trunc) 
    r_upperCI_trunc_R1 <- r_est_trunc_R1+qnorm(0.975)*sqrt(r_var_trunc_R1) 
    r_upperCI_trunc_R0 <- r_est_trunc_R0+qnorm(0.975)*sqrt(r_var_trunc_R0) 
    r_est_trunc_diff <- colMeans(est_temp_trunc_diff)
    r_est_trunc_R1_diff <- colMeans(est_R1_temp_trunc_diff)
    r_est_trunc_R0_diff <- colMeans(est_R0_temp_trunc_diff)
    r_var_trunc_diff <- colSums(var_temp_trunc_diff)/(fold*n)
    r_var_trunc_R1_diff <- colSums(var_R1_temp_trunc_diff)/(fold*n)
    r_var_trunc_R0_diff <- colSums(var_R0_temp_trunc_diff)/(fold*n)
    r_lowerCI_trunc_diff <- r_est_trunc_diff-qnorm(0.975)*sqrt(r_var_trunc_diff)
    r_lowerCI_trunc_R1_diff <- r_est_trunc_R1_diff-qnorm(0.975)*sqrt(r_var_trunc_R1_diff)
    r_lowerCI_trunc_R0_diff <- r_est_trunc_R0_diff-qnorm(0.975)*sqrt(r_var_trunc_R0_diff)
    r_upperCI_trunc_diff <- r_est_trunc_diff+qnorm(0.975)*sqrt(r_var_trunc_diff) 
    r_upperCI_trunc_R1_diff <- r_est_trunc_R1_diff+qnorm(0.975)*sqrt(r_var_trunc_R1_diff) 
    r_upperCI_trunc_R0_diff <- r_est_trunc_R0_diff+qnorm(0.975)*sqrt(r_var_trunc_R0_diff) 

  }
  
  ## another variance estimation & skewness & multiplier bootstrap: empty container
  # r_var_2 <- r_var_R1_2 <- r_var_R0_2 <- c()
  # r_var_diff_2 <- r_var_R1_diff_2 <- r_var_R0_diff_2 <- c()
  # skew <- skew_R1 <- skew_R0 <- c()
  # skew_diff <- skew_R1_diff <- skew_R0_diff <- c()
  # if(!simple_trunc){
  #   r_var_trunc_2 <- r_var_trunc_R1_2 <- r_var_trunc_R0_2 <- c()
  #   r_var_trunc_diff_2 <- r_var_trunc_R1_diff_2 <- r_var_trunc_R0_diff_2 <- c()
  #   skew_trunc <- skew_trunc_R1 <- skew_trunc_R0 <- c()
  #   skew_trunc_diff <- skew_trunc_R1_diff <- skew_trunc_R0_diff <- c()
  #   
  #   B <- 2000
  #   q <- q_R1 <- q_R0 <- c()
  #   q_diff <- q_R1_diff <- q_R0_diff <- c()
  # }

  ## another variance estimation & skewness & multiplier bootstrap
  # for (g in 1:length(gamma)){
  #   
  #   r_var_2 <- c(r_var_2, sum((IF[[g]]-mean(IF[[g]]))^2)/(n*(n-1)))
  #   r_var_R1_2 <- c(r_var_R1_2, sum((IF_R1[[g]]-vector_R1*mean(IF_R1[[g]]))^2)/(n*(n-1)))
  #   r_var_R0_2 <- c(r_var_R0_2, sum((IF_R0[[g]]-vector_R0*mean(IF_R0[[g]]))^2)/(n*(n-1)))
  #   r_var_diff_2 <- c(r_var_diff_2, sum((IF_diff[[g]]-mean(IF_diff[[g]]))^2)/(n*(n-1)))
  #   r_var_R1_diff_2 <- c(r_var_R1_diff_2, sum((IF_R1_diff[[g]]-vector_R1*mean(IF_R1_diff[[g]]))^2)/(n*(n-1)))
  #   r_var_R0_diff_2 <- c(r_var_R0_diff_2, sum((IF_R0_diff[[g]]-vector_R0*mean(IF_R0_diff[[g]]))^2)/(n*(n-1)))
  #   
  #   skew <- c(skew, mean(IF[[g]]^3)/(mean(IF[[g]]^2)^(3/2)))
  #   skew_R1 <- c(skew_R1, mean(IF_R1[[g]]^3)/(mean(IF_R1[[g]]^2)^(3/2)))
  #   skew_R0 <- c(skew_R0, mean(IF_R0[[g]]^3)/(mean(IF_R0[[g]]^2)^(3/2)))
  #   skew_diff <- c(skew_diff, mean(IF_diff[[g]]^3)/(mean(IF_diff[[g]]^2)^(3/2)))
  #   skew_R1_diff <- c(skew_R1_diff, mean(IF_R1_diff[[g]]^3)/(mean(IF_R1_diff[[g]]^2)^(3/2)))
  #   skew_R0_diff <- c(skew_R0_diff, mean(IF_R0_diff[[g]]^3)/(mean(IF_R0_diff[[g]]^2)^(3/2)))
  #   
  #   if(!simple_trunc){
  #     
  #   r_var_trunc_2 <- c(r_var_trunc_2, sum((IF_trunc[[g]]-mean(IF_trunc[[g]]))^2)/(n*(n-1)))
  #   r_var_trunc_R1_2 <- c(r_var_trunc_R1_2, sum((IF_R1_trunc[[g]]-vector_R1*mean(IF_R1_trunc[[g]]))^2)/(n*(n-1)))
  #   r_var_trunc_R0_2 <- c(r_var_trunc_R0_2, sum((IF_R0_trunc[[g]]-vector_R0*mean(IF_R0_trunc[[g]]))^2)/(n*(n-1)))
  #   r_var_trunc_diff_2 <- c(r_var_trunc_diff_2, sum((IF_trunc_diff[[g]]-mean(IF_trunc_diff[[g]]))^2)/(n*(n-1)))
  #   r_var_trunc_R1_diff_2 <- c(r_var_trunc_R1_diff_2, sum((IF_R1_trunc_diff[[g]]-vector_R1*mean(IF_R1_trunc_diff[[g]]))^2)/(n*(n-1)))
  #   r_var_trunc_R0_diff_2 <- c(r_var_trunc_R0_diff_2, sum((IF_R0_trunc_diff[[g]]-vector_R0*mean(IF_R0_trunc_diff[[g]]))^2)/(n*(n-1)))
  #   
  #   skew_trunc <- c(skew_trunc, mean(IF_trunc[[g]]^3)/(mean(IF_trunc[[g]]^2)^(3/2)))
  #   skew_trunc_R1 <- c(skew_trunc_R1, mean(IF_R1_trunc[[g]]^3)/(mean(IF_R1_trunc[[g]]^2)^(3/2)))
  #   skew_trunc_R0 <- c(skew_trunc_R0, mean(IF_R0_trunc[[g]]^3)/(mean(IF_R0_trunc[[g]]^2)^(3/2)))
  #   skew_trunc_diff <- c(skew_trunc_diff, mean(IF_trunc_diff[[g]]^3)/(mean(IF_trunc_diff[[g]]^2)^(3/2)))
  #   skew_trunc_R1_diff <- c(skew_trunc_R1_diff, mean(IF_R1_trunc_diff[[g]]^3)/(mean(IF_R1_trunc_diff[[g]]^2)^(3/2)))
  #   skew_trunc_R0_diff <- c(skew_trunc_R0_diff, mean(IF_R0_trunc_diff[[g]]^3)/(mean(IF_R0_trunc_diff[[g]]^2)^(3/2)))
  #   
  #   wi <- matrix(sample(c(-1, 1), n*B, replace=TRUE), nrow=B, ncol=n)
  #   est_star <- r_est_trunc[g]+rowMeans(wi*(IF_trunc[[g]]-mean(IF[[g]])))
  #   est_R1_star <- r_est_trunc_R1[g]+rowMeans(wi*(IF_R1_trunc[[g]]-vector_R1*mean(IF_R1_trunc[[g]])))
  #   est_R0_star <- r_est_trunc_R0[g]+rowMeans(wi*(IF_R0_trunc[[g]]-vector_R0*mean(IF_R0_trunc[[g]])))
  #   est_diff_star <- r_est_trunc_diff[g]+rowMeans(wi*(IF_trunc_diff[[g]]-mean(IF_trunc_diff[[g]])))
  #   est_R1_diff_star <- r_est_trunc_R1_diff[g]+rowMeans(wi*(IF_R1_trunc_diff[[g]]-vector_R1*mean(IF_R1_trunc_diff[[g]])))
  #   est_R0_diff_star <- r_est_trunc_R0_diff[g]+rowMeans(wi*(IF_R0_trunc_diff[[g]]-vector_R0*mean(IF_R0_trunc_diff[[g]])))
  #   
  #   t_star <- (est_star-r_est_trunc[g])/sqrt(r_var_trunc[g])
  #   t_R1_star <- (est_R1_star-r_est_trunc_R1[g])/sqrt(r_var_trunc_R1[g])
  #   t_R0_star <- (est_R0_star-r_est_trunc_R0[g])/sqrt(r_var_trunc_R0[g])
  #   t_diff_star <- (est_diff_star-r_est_trunc_diff[g])/sqrt(r_var_trunc_diff[g])
  #   t_R1_diff_star <- (est_R1_diff_star-r_est_trunc_R1_diff[g])/sqrt(r_var_trunc_R1_diff[g])
  #   t_R0_diff_star <- (est_R0_diff_star-r_est_trunc_R0_diff[g])/sqrt(r_var_trunc_R0_diff[g])
  #   
  #   q <- c(q, quantile(abs(t_star), 1-0.05))
  #   q_R1 <- c(q_R1, quantile(abs(t_R1_star), 1-0.05))
  #   q_R0 <- c(q_R0, quantile(abs(t_R0_star), 1-0.05))
  #   q_diff <- c(q_diff, quantile(abs(t_diff_star), 1-0.05))
  #   q_R1_diff <- c(q_R1_diff, quantile(abs(t_R1_diff_star), 1-0.05))
  #   q_R0_diff <- c(q_R0_diff, quantile(abs(t_R0_diff_star), 1-0.05))
  #   
  #   }
  # }
  
  ## one-step jackknife
  # jack_est_trunc <- c()
  # jack_est_trunc_R0 <- c()
  # jack_est_trunc_R1 <- c()
  # for (g in 1:length(gamma)){
  #   jack_est_trunc <- cbind(jack_est_trunc, r_est_trunc[g]-IF_trunc[[g]]/(n-1))
  #   jack_est_trunc_R0 <- cbind(jack_est_trunc_R0, r_est_trunc_R0[g]-IF_R0_trunc[[g]]/(n-1))
  #   jack_est_trunc_R1 <- cbind(jack_est_trunc_R1, r_est_trunc_R1[g]-IF_R1_trunc[[g]]/(n-1))
  # }
  # jack_var_trunc <- (n-1)/n*colSums((jack_est_trunc-matrix(colMeans(jack_est_trunc), nrow=n, ncol=length(gamma), byrow=T))^2)
  # jack_var_trunc_R0 <- (n-1)/n*colSums((jack_est_trunc_R0-matrix(colMeans(jack_est_trunc_R0), nrow=n, ncol=length(gamma), byrow=T))^2)
  # jack_var_trunc_R1 <- (n-1)/n*colSums((jack_est_trunc_R1-matrix(colMeans(jack_est_trunc_R1), nrow=n, ncol=length(gamma), byrow=T))^2)
  # 
  ## output final results
  if(IF_output){
    if(simple_trunc){
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0, 
                     # var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2, 
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     # var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2, 
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     IF=IF, IF_R1=IF_R1, IF_R0=IF_R0, 
                     IF_diff=IF_diff, IF_R1_diff=IF_R1_diff, IF_R0_diff=IF_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     id_list=id_list)
    }else{
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     est_mat=est_temp, est_R1_mat=est_R1_temp, est_R0_mat=est_R0_temp, 
                     est_trunc=r_est_trunc, est_trunc_R1=r_est_trunc_R1, est_trunc_R0=r_est_trunc_R0,
                     est_trunc_mat=est_temp_trunc, est_trunc_R1_mat=est_R1_temp_trunc, est_trunc_R0_mat=est_R0_temp_trunc, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0,
                     # var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2, 
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                     # var_trunc_2=r_var_trunc_2, var_trunc_R1_2=r_var_trunc_R1_2, var_trunc_R0_2=r_var_trunc_R0_2, 
                     # skew_trunc=skew_trunc, skew_trunc_R1=skew_trunc_R1, skew_trunc_R0=skew_trunc_R0,
                     # jack_var_trunc=jack_var_trunc, jack_var_trunc_R0=jack_var_trunc_R0, jack_var_trunc_R1=jack_var_trunc_R1, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                     est_diff_mat=est_temp_diff, est_R1_diff_mat=est_R1_temp_diff, est_R0_diff_mat=est_R0_temp_diff, 
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     est_trunc_diff_mat=est_temp_trunc_diff, est_R1_trunc_diff_mat=est_R1_temp_trunc_diff, est_R0_trunc_diff_mat=est_R0_temp_trunc_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                     # var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2, 
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     # var_trunc_diff_2=r_var_trunc_diff_2, var_trunc_R1_diff_2=r_var_trunc_R1_diff_2, var_trunc_R0_diff_2=r_var_trunc_R0_diff_2, 
                     # skew_trunc_diff=skew_trunc_diff, skew_trunc_R1_diff=skew_trunc_R1_diff, skew_trunc_R0_diff=skew_trunc_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R1_diff=r_lowerCI_trunc_R1_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff, 
                     upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R1_diff=r_upperCI_trunc_R1_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff, 
                     IF=IF, IF_R1=IF_R1, IF_R0=IF_R0,
                     IF_trunc=IF_trunc, IF_trunc_R1=IF_R1_trunc, IF_trunc_R0=IF_R0_trunc,
                     IF_diff=IF_diff, IF_R1_diff=IF_R1_diff, IF_R0_diff=IF_R0_diff,
                     IF_trunc_diff=IF_trunc_diff, IF_trunc_R1_diff=IF_R1_trunc_diff, IF_trunc_R0_diff=IF_R0_trunc_diff,
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1,
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     id_list=id_list)
                     # q=q, q_R1=q_R1, q_R0=q_R0, q_diff=q_diff, q_R1_diff=q_R1_diff, q_R0_diff=q_R0_diff)
    }
  }else{
    if(simple_trunc){
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0, 
                     # var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2, 
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     # var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2, 
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     id_list=id_list)
    }else{
        result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0,
                       est_mat=est_temp, est_R1_mat=est_R1_temp, est_R0_mat=est_R0_temp, 
                       est_trunc=r_est_trunc, est_trunc_R1=r_est_trunc_R1, est_trunc_R0=r_est_trunc_R0,
                       est_trunc_mat=est_temp_trunc, est_trunc_R1_mat=est_R1_temp_trunc, est_trunc_R0_mat=est_R0_temp_trunc, 
                       var=r_var, var_R1=r_var_R1, var_R0=r_var_R0,
                       # var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2, 
                       # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                       var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                       # var_trunc_2=r_var_trunc_2, var_trunc_R1_2=r_var_trunc_R1_2, var_trunc_R0_2=r_var_trunc_R0_2, 
                       # skew_trunc=skew_trunc, skew_trunc_R1=skew_trunc_R1, skew_trunc_R0=skew_trunc_R0,
                       # jack_var_trunc=jack_var_trunc, jack_var_trunc_R0=jack_var_trunc_R0, jack_var_trunc_R1=jack_var_trunc_R1, 
                       lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                       upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                       lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                       upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                       est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                       est_diff_mat=est_temp_diff, est_R1_diff_mat=est_R1_temp_diff, est_R0_diff_mat=est_R0_temp_diff, 
                       est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                       est_trunc_diff_mat=est_temp_trunc_diff, est_R1_trunc_diff_mat=est_R1_temp_trunc_diff, est_R0_trunc_diff_mat=est_R0_temp_trunc_diff, 
                       var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                       # var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2, 
                       # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                       var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                       # var_trunc_diff_2=r_var_trunc_diff_2, var_trunc_R1_diff_2=r_var_trunc_R1_diff_2, var_trunc_R0_diff_2=r_var_trunc_R0_diff_2, 
                       # skew_trunc_diff=skew_trunc_diff, skew_trunc_R1_diff=skew_trunc_R1_diff, skew_trunc_R0_diff=skew_trunc_R0_diff, 
                       lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                       upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff,
                       lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R1_diff=r_lowerCI_trunc_R1_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff,
                       upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R1_diff=r_upperCI_trunc_R1_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff,
                       pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1,
                       pain_bq_reordered_R0=pain_bq_reordered_R0,
                       id_list=id_list)
                       # q=q, q_R1=q_R1, q_R0=q_R0, q_diff=q_diff, q_R1_diff=q_R1_diff, q_R0_diff=q_R0_diff)
    }
  }
  
}



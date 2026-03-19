#' Create containers for cross-fitted estimator
#' 
#' This helper function initializes an environment to hold influence-function
#' vectors, point estimates, variance estimates, and other intermediate 
#' quantities for the cross-fitted estimator. 
#' The structure is designed to facilitate storage and retrieval of these 
#' objects across multiple folds and sensitivity parameter values.
#' @noRd
est_s_t_y_create_containers <- function(gamma, fold) {
  
  e <- new.env(parent = emptyenv())
  e$IF               <- vector(mode = "list", length = length(gamma))
  e$IF_R1            <- vector(mode = "list", length = length(gamma))
  e$IF_R0            <- vector(mode = "list", length = length(gamma))
  e$est_temp         <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R1_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R0_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_temp         <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R1_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R0_temp      <- matrix(0, nrow = fold, ncol = length(gamma))
  e$IF_diff          <- vector(mode = "list", length = length(gamma))
  e$IF_R1_diff       <- vector(mode = "list", length = length(gamma))
  e$IF_R0_diff       <- vector(mode = "list", length = length(gamma))
  e$est_temp_diff    <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R1_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$est_R0_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_temp_diff    <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R1_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$var_R0_temp_diff <- matrix(0, nrow = fold, ncol = length(gamma))
  e$vector_R1        <- c()
  e$vector_R0        <- c()
  e$id_list          <- c()
  
  structure(
    e,
    class = c("est_container", class(e))
  )
}

#' Data adaptive influence function truncation
#' 
#' This helper function perform influence function truncation by (Wang et al., 2021)
#' @noRd
IF_trunc_func <- function(IF){
  
  if (eq(max(abs(IF)),z=IF)>0) {
    trunc <- max(abs(IF))
  } else {
    trunc <- uniroot(eq, z=IF, interval=c(0.01, max(abs(IF))))$root
  }
  IF_trunc <- pmin(abs(IF),trunc) *  sign(IF)
  return(IF_trunc)
  
}

#' Single index models
#' 
#' @noRd
SIM <- function(X, Y, kernel, method, single_index_method, use_mave){
  if(single_index_method=="fixed_bandwidth"){
    if(use_mave){
      requireNamespace('MAVE', quietly = TRUE)
      SDR <- coef(MAVE::mave.compute(X, Y, max.dim = 1), 1)
    }else{
      SDR <- cumuSIR(X=X, Y=Y)
    }
    fit_h <- SIDRnew_fixed_bandwidth(X=X, Y=Y, initial = SDR, kernel = kernel, method=method, ids = 1:length(Y))
    
  }else if(single_index_method=="fixed_coef"){
    if(use_mave){
      requireNamespace('MAVE', quietly = TRUE)
      SDR <- coef(MAVE::mave.compute(X, Y, max.dim = 1), 1)
    }else{
      SDR <- cumuSIR(X=X, Y=Y)
    }
    fit_h <- SIDR_Ravinew(X=X, Y=Y, initial=c(1,as.vector(SDR[-1]/SDR[1])), kernel=kernel, method=method, index_ID=1:length(Y))
  }else if(single_index_method=="norm1coef"){
    fit_h  <- fit_SensIAT_single_index_norm1coef_model(X=X, Y=Y,  ids=1:length(Y), kernel=kernel, bw.selection="ise", 
                                                       bw.method="optim", use_mave=use_mave)
  }
  return(fit_h)
}


#' One-step, split sample estimator for E[Y(t)], E[Y(t)|R=1], E[Y(t)|R=0], 
#'   and E[Y(t)-Y0], E[Y(t)-Y0|R=1], E[Y(t)-Y0|R=0]
#'   
#' Estimates study-specific and overall outcome means (and difference from baseline)
#' using cross-fitting with single index models (SIMs) and
#' nuisance models for treatment and outcome missingness (`mgcv::gam`). The function
#' also computes influence-function-based variances, confidence intervals, and
#' optional truncated influence-function diagnostics.
#'
#' @param Y Numeric outcome vector. Missing values are internally replaced with
#'   `0` prior to model fitting.
#' @param M Binary indicator for observed outcome (`1` = observed, `0` =
#'   missing).
#' @param R Binary group indicator used to stratify nuisance and outcome models.
#' @param X Data frame or matrix of baseline covariates.
#' @param t Treatment assignment vector.
#' @param trt Treatment level for which the target estimand is computed.
#' @param gamma Numeric vector of sensitivity parameters.
#' @param fold Number of cross-fitting folds.
#' @param seed Optional integer random seed for fold assignment. Use `NULL` to
#'   leave RNG state unchanged.
#' @param IF_output Logical; if `TRUE`, include influence-function vectors in
#'   the returned list.
#' @param simple_trunc Logical; if `TRUE`, apply quantile truncation to inverse
#'   probability weights. If `FALSE`, apply IF truncation diagnostics.
#' @param quant Numeric in `(0, 1)` used as the upper quantile for simple weight
#'   truncation when `simple_trunc = TRUE`.
#' @param kernel Characters; Kernel used for SIMs. `K2_Biweight` for Epanechnikov kernel, 
#'   `dnorm` for Gaussian kernel. 
#' @param single_index_method Characters; Three implementations for SIMs: `fixed_bandwidth` 
#'    for setting bandwidth to 1, `fixed_coef` for setting the first coefficient to 1, and `norm1coef`
#'    for setting the norm of coefficients to 1. 
#' @param method Characters; Optimization method used for SIMs. Choices are: `optim`, `nlminb`, `nmk`. 
#'    Note that method is set to `optim` if single_index_method=`norm1coef`. 
#' @param use_mave Logical; if `TRUE`, use Minimum Average Variance Estimation (MAVE) method for initial
#'    coefficients value for SIMs. If `FALSE`, use sliced inverse regression. Default is `TRUE`. 
#' @param s_t_y A function of Y in the exponential tilting model. If NULL, s_t_y is set to pnorm((y-60)/25). 
#' @param coef_g.fit Optional starting values for a treatment model; currently
#'   retained for interface compatibility.
#' @param coef_t_R0.fit Optional starting coefficients for treatment model fit
#'   in `t=trt` and `R = 0`stratum.
#' @param coef_t_R1.fit Optional starting coefficients for treatment model fit
#'   in `t=trt` and `R = 1` stratum.
#' @param coef_M_R0.fit Optional starting coefficients for missingness model fit
#'   in `R = 0` stratum.
#' @param coef_M_R1.fit Optional starting coefficients for missingness model fit
#'   in `R = 1` stratum.
#'
#' @return A named list of estimates and uncertainty summaries for each value in
#'   `gamma`. Core elements include point estimates (`est`, `est_R1`, `est_R0`), variance
#'   estimates (`var`, `var_R1`, `var_R0`), and confidence interval bounds (`lowerCI*`, `upperCI*`). 
#'   Additional components depend on `simple_trunc` and `IF_output`:
#'   \itemize{
#'   \item `simple_trunc = TRUE`: returns quantile-weight-truncated summaries only.
#'   \item `simple_trunc = FALSE`: additionally returns truncated summaries and
#'   truncated IF objects when requested.
#'   \item `IF_output = TRUE`: includes influence-function lists (`IF*`) and,
#'   when relevant, truncated IF lists (`IF_trunc*`).
#'   }
#'
#' @details
#' The procedure uses sample-splitting/cross-fitting to reduce overfitting bias
#' in nuisance estimation. Outcome regressions are fit with single index models
#' and then integrated over estimated conditional distributions to obtain
#' conditional means and sensitivity-adjusted moments.
#'
#' @examples
#' # out <- est_psi(Y, M, R, X, t, trt = 1, gamma = c(0, 0.5),
#' #                fold = 5, seed = 1, IF_output = FALSE,
#' #                simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
#' #                single_index_method="norm1coef", method="optim")
est_psi <- function(Y, M, R, X, t, trt, gamma, fold, seed, IF_output, 
                    simple_trunc, quant, kernel, single_index_method, method="optim", 
                    use_mave=TRUE, s_t_y=NULL, coef_g.fit=NULL, coef_t_R0.fit=NULL, 
                    coef_t_R1.fit=NULL, coef_M_R0.fit=NULL, coef_M_R1.fit=NULL){
  
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
  containers <- est_s_t_y_create_containers(gamma, fold)
  if (!simple_trunc){
    containers_trunc <- est_s_t_y_create_containers(gamma, fold)
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
    
    trt.ind_out_fold_R0 <- trt.ind_out_fold[which(R_out_fold==0)]
    X_out_fold_R0 <- X_out_fold[which(R_out_fold==0), ]
    trt.ind_out_fold_R1 <- trt.ind_out_fold[which(R_out_fold==1)]
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
    
    containers$id_list <- c(containers$id_list, fold_list[[k]])
    
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
    fit_t_R0_h <- SIM(X=X_out_fold_t_R0_adjust_scale[which(M_out_fold_t_R0==1), ], Y=Y_out_fold_t_R0[which(M_out_fold_t_R0==1)], 
                      kernel=kernel, method=method, single_index_method=single_index_method, use_mave=use_mave)
    fit_t_R1_h <- SIM(X=X_out_fold_t_R1_adjust_scale[which(M_out_fold_t_R1==1), ], Y=Y_out_fold_t_R1[which(M_out_fold_t_R1==1)], 
                      kernel=kernel, method=method, single_index_method=single_index_method, use_mave=use_mave)
    
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
    dF_t_R0_X_t_R0 <- F_t_R0_X_t_R0[, -1, drop = FALSE]-F_t_R0_X_t_R0[,-ny_t_R0, drop = FALSE]
    dF_t_R0_X_t0_R0 <- F_t_R0_X_t0_R0[, -1, drop = FALSE]-F_t_R0_X_t0_R0[,-ny_t_R0, drop = FALSE]
    dF_t_R1_X_R1 <- F_t_R1_X_R1[, -1, drop = FALSE]-F_t_R1_X_R1[,-ny_t_R1, drop = FALSE]
    
    ## conditional expectation of Y given R=0, T=t and X, for I(R=0, T=t)
    mu_Y_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0+c(0, y_t_R0[-ny_t_R0]))/2)
    
    ## conditional expectation of Y given R=1, T=t and X, for I(R=1)
    mu_Y_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (y_t_R1+c(0, y_t_R1[-ny_t_R1]))/2)
    
    ## P(R=1)
    prop.R1 <- mean(R_in_fold)
    
    ## containers$vector_R1 and containers$vector_R0
    ind_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), rep(1/prop.R1, length(M_in_fold_R1)))
    ind_R0_temp <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),  rep(0, length(M_in_fold_R1)))
    containers$vector_R1 <- c(containers$vector_R1, ind_R1_temp)
    containers$vector_R0 <- c(containers$vector_R0, ind_R0_temp)
    
    ## Compute influence function for each gamma_t
    for (g in 1:length(gamma)){
    
      ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=t)
      mu_Yexp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0[-1]*exp(gamma[g]*s_t_y(y_t_R0[-1]))+
                                                     y_t_R0[-ny_t_R0]*exp(gamma[g]*s_t_y(y_t_R0[-ny_t_R0])))/2)
      
      ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=t)
      mu_exp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (exp(gamma[g]*s_t_y(y_t_R0[-1]))+
                                                    exp(gamma[g]*s_t_y(y_t_R0[-ny_t_R0])))/2)
      
      ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=1-t)
      mu_Yexp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (y_t_R0[-1]*exp(gamma[g]*s_t_y(y_t_R0[-1]))+
                                                       y_t_R0[-ny_t_R0]*exp(gamma[g]*s_t_y(y_t_R0[-ny_t_R0])))/2)
      
      ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=1-t)
      mu_exp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (exp(gamma[g]*s_t_y(y_t_R0[-1]))+
                                                      exp(gamma[g]*s_t_y(y_t_R0[-ny_t_R0])))/2)

    ## IF+psi within each fold
    if(trt==1){
      if_temp <- c(M_in_fold_t_R0*eta_t_R0_weight*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*s_t_y(Y_in_fold_t_R0))/mu_exp_t_R0_X_t_R0*
                                              (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                   M_in_fold_t0_R0*eta_t0_R0_weight*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                   M_in_fold_R1*eta_T_R1_weight*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
        c((1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0, 
          (1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0,
          (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1)
      if_temp_diff <- if_temp-pain_bq_temp
      
      if_R0_temp <- c(M_in_fold_t_R0*eta_t_R0_weight/(1-prop.R1)*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*s_t_y(Y_in_fold_t_R0))/mu_exp_t_R0_X_t_R0*
                                                 (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                      M_in_fold_t0_R0*eta_t0_R0_weight/(1-prop.R1)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                      rep(0, length(M_in_fold_R1)))+
        c((1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0/(1-prop.R1), 
          (1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/(mu_exp_t_R0_X_t0_R0*(1-prop.R1)),
          rep(0, length(M_in_fold_R1)))
      if_R0_temp_diff <- if_R0_temp-pain_bq_R0_temp
    }else{
      if_temp <- c(M_in_fold_t0_R0*eta_t0_R0_weight*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                   M_in_fold_t_R0*eta_t_R0_weight*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*s_t_y(Y_in_fold_t_R0))/mu_exp_t_R0_X_t_R0*
                                              (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                   M_in_fold_R1*eta_T_R1_weight*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
        c((1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0,
          (1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0, 
          (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1)
      if_temp_diff <- if_temp-pain_bq_temp
      
      if_R0_temp  <- c(M_in_fold_t0_R0*eta_t0_R0_weight/(1-prop.R1)*mu_Yexp_t_R0_X_t0_R0/mu_exp_t_R0_X_t0_R0, 
                         M_in_fold_t_R0*eta_t_R0_weight/(1-prop.R1)*(Y_in_fold_t_R0+(pi_R0_weight-1)*exp(gamma[g]*s_t_y(Y_in_fold_t_R0))/mu_exp_t_R0_X_t_R0*
                                                    (Y_in_fold_t_R0-mu_Yexp_t_R0_X_t_R0/mu_exp_t_R0_X_t_R0)), 
                         rep(0, length(M_in_fold_R1)))+
        c((1-M_in_fold_t0_R0*eta_t0_R0_weight)*mu_Yexp_t_R0_X_t0_R0/(mu_exp_t_R0_X_t0_R0*(1-prop.R1)),
          (1-M_in_fold_t_R0*eta_t_R0_weight)*mu_Y_t_R0_X_t_R0/(1-prop.R1), 
          rep(0, length(M_in_fold_R1)))
      if_R0_temp_diff <- if_R0_temp-pain_bq_R0_temp
    }
    
    if_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                    M_in_fold_R1*eta_T_R1_weight/prop.R1*(trt.ind_in_fold_R1*pi_R1_weight*(Y_in_fold_R1-mu_Y_t_R1_X_R1)+mu_Y_t_R1_X_R1))+
      c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), (1-M_in_fold_R1*eta_T_R1_weight)*mu_Y_t_R1_X_R1/prop.R1)
    if_R1_temp_diff <- if_R1_temp-pain_bq_R1_temp
    
    containers$IF[[g]] <- c(containers$IF[[g]], if_temp)
    containers$IF_R1[[g]] <- c(containers$IF_R1[[g]], if_R1_temp)
    containers$IF_R0[[g]] <- c(containers$IF_R0[[g]], if_R0_temp)
    containers$IF_diff[[g]] <- c(containers$IF_diff[[g]], if_temp_diff)
    containers$IF_R1_diff[[g]] <- c(containers$IF_R1_diff[[g]], if_R1_temp_diff)
    containers$IF_R0_diff[[g]] <- c(containers$IF_R0_diff[[g]], if_R0_temp_diff)
    
    containers$est_temp[k, g] <- mean(if_temp)
    containers$est_R1_temp[k, g] <- mean(if_R1_temp)
    containers$est_R0_temp[k, g] <- mean(if_R0_temp)
    containers$var_temp[k, g] <- sum((if_temp - mean(if_temp))^2)/(nk_in_fold-1)
    containers$var_R1_temp[k, g] <- sum((if_R1_temp - ind_R1_temp*mean(if_R1_temp))^2)/(nk_in_fold-1)
    containers$var_R0_temp[k, g] <- sum((if_R0_temp - ind_R0_temp*mean(if_R0_temp))^2)/(nk_in_fold-1)
    containers$est_temp_diff[k, g] <- mean(if_temp_diff)
    containers$est_R1_temp_diff[k, g] <- mean(if_R1_temp_diff)
    containers$est_R0_temp_diff[k, g] <- mean(if_R0_temp_diff)
    containers$var_temp_diff[k, g] <- sum((if_temp_diff - mean(if_temp_diff))^2)/(nk_in_fold-1)
    containers$var_R1_temp_diff[k, g] <- sum((if_R1_temp_diff - ind_R1_temp*mean(if_R1_temp_diff))^2)/(nk_in_fold-1)
    containers$var_R0_temp_diff[k, g] <- sum((if_R0_temp_diff - ind_R0_temp*mean(if_R0_temp_diff))^2)/(nk_in_fold-1)
    
    if(!simple_trunc){
      
      if_R1_temp_trunc <- IF_trunc_func(if_R1_temp)
      containers_trunc$IF_R1[[g]] <- c(containers_trunc$IF_R1[[g]], if_R1_temp_trunc)
      
      if_R1_temp_trunc_diff <- IF_trunc_func(if_R1_temp_diff)
      containers_trunc$IF_R1_diff[[g]] <- c(containers_trunc$IF_R1_diff[[g]], if_R1_temp_trunc_diff)
      
      if_R0_temp_trunc <- IF_trunc_func(if_R0_temp)
      containers_trunc$IF_R0[[g]] <- c(containers_trunc$IF_R0[[g]], if_R0_temp_trunc)
      
      if_R0_temp_trunc_diff <- IF_trunc_func(if_R0_temp_diff)
      containers_trunc$IF_R0_diff[[g]] <- c(containers_trunc$IF_R0_diff[[g]], if_R0_temp_trunc_diff)
      
      if_temp_trunc <- if_R1_temp_trunc*prop.R1+if_R0_temp_trunc*(1-prop.R1)
      containers_trunc$IF[[g]] <- c(containers_trunc$IF[[g]], if_temp_trunc)
      if_temp_trunc_diff <- if_R1_temp_trunc_diff*prop.R1+if_R0_temp_trunc_diff*(1-prop.R1)
      containers_trunc$IF_diff[[g]] <- c(containers_trunc$IF_diff[[g]], if_temp_trunc_diff)

      containers_trunc$est_temp[k, g] <- mean(if_temp_trunc)
      containers_trunc$est_R1_temp[k, g] <- mean(if_R1_temp_trunc)
      containers_trunc$est_R0_temp[k, g] <- mean(if_R0_temp_trunc)
      containers_trunc$var_temp[k, g] <- sum((if_temp_trunc - mean(if_temp_trunc))^2)/(nk_in_fold-1)
      containers_trunc$var_R1_temp[k, g] <- sum((if_R1_temp_trunc - ind_R1_temp*mean(if_R1_temp_trunc))^2)/(nk_in_fold-1)
      containers_trunc$var_R0_temp[k, g] <- sum((if_R0_temp_trunc - ind_R0_temp*mean(if_R0_temp_trunc))^2)/(nk_in_fold-1)
      
      containers_trunc$est_temp_diff[k, g] <- mean(if_temp_trunc_diff)
      containers_trunc$est_R1_temp_diff[k, g] <- mean(if_R1_temp_trunc_diff)
      containers_trunc$est_R0_temp_diff[k, g] <- mean(if_R0_temp_trunc_diff)
      containers_trunc$var_temp_diff[k, g] <- sum((if_temp_trunc_diff - mean(if_temp_trunc_diff))^2)/(nk_in_fold-1)
      containers_trunc$var_R1_temp_diff[k, g] <- sum((if_R1_temp_trunc_diff - ind_R1_temp*mean(if_R1_temp_trunc_diff))^2)/(nk_in_fold-1)
      containers_trunc$var_R0_temp_diff[k, g] <- sum((if_R0_temp_trunc_diff - ind_R0_temp*mean(if_R0_temp_trunc_diff))^2)/(nk_in_fold-1)
      
    }
    }
  }
  
  ## output containers
  r_est <- colMeans(containers$est_temp)
  r_est_R1 <- colMeans(containers$est_R1_temp)
  r_est_R0 <- colMeans(containers$est_R0_temp)
  r_var <- colSums(containers$var_temp)/(fold*n)
  r_var_R1 <- colSums(containers$var_R1_temp)/(fold*n)
  r_var_R0 <- colSums(containers$var_R0_temp)/(fold*n)
  r_lowerCI <- r_est-qnorm(0.975)*sqrt(r_var)
  r_lowerCI_R1 <- r_est_R1-qnorm(0.975)*sqrt(r_var_R1)
  r_lowerCI_R0 <- r_est_R0-qnorm(0.975)*sqrt(r_var_R0)
  r_upperCI <- r_est+qnorm(0.975)*sqrt(r_var) 
  r_upperCI_R1 <- r_est_R1+qnorm(0.975)*sqrt(r_var_R1) 
  r_upperCI_R0 <- r_est_R0+qnorm(0.975)*sqrt(r_var_R0) 
  r_est_diff <- colMeans(containers$est_temp_diff)
  r_est_R1_diff <- colMeans(containers$est_R1_temp_diff)
  r_est_R0_diff <- colMeans(containers$est_R0_temp_diff)
  r_var_diff <- colSums(containers$var_temp_diff)/(fold*n)
  r_var_R1_diff <- colSums(containers$var_R1_temp_diff)/(fold*n)
  r_var_R0_diff <- colSums(containers$var_R0_temp_diff)/(fold*n)
  r_lowerCI_diff <- r_est_diff-qnorm(0.975)*sqrt(r_var_diff)
  r_lowerCI_R1_diff <- r_est_R1_diff-qnorm(0.975)*sqrt(r_var_R1_diff)
  r_lowerCI_R0_diff <- r_est_R0_diff-qnorm(0.975)*sqrt(r_var_R0_diff)
  r_upperCI_diff <- r_est_diff+qnorm(0.975)*sqrt(r_var_diff) 
  r_upperCI_R1_diff <- r_est_R1_diff+qnorm(0.975)*sqrt(r_var_R1_diff) 
  r_upperCI_R0_diff <- r_est_R0_diff+qnorm(0.975)*sqrt(r_var_R0_diff)
  if(!simple_trunc){
    r_est_trunc <- colMeans(containers_trunc$est_temp)
    r_est_trunc_R1 <- colMeans(containers_trunc$est_R1_temp)
    r_est_trunc_R0 <- colMeans(containers_trunc$est_R0_temp)
    r_var_trunc <- colSums(containers_trunc$var_temp)/(fold*n)
    r_var_trunc_R1 <- colSums(containers_trunc$var_R1_temp)/(fold*n)
    r_var_trunc_R0 <- colSums(containers_trunc$var_R0_temp)/(fold*n)
    r_lowerCI_trunc <- r_est_trunc-qnorm(0.975)*sqrt(r_var_trunc)
    r_lowerCI_trunc_R1 <- r_est_trunc_R1-qnorm(0.975)*sqrt(r_var_trunc_R1)
    r_lowerCI_trunc_R0 <- r_est_trunc_R0-qnorm(0.975)*sqrt(r_var_trunc_R0)
    r_upperCI_trunc <- r_est_trunc+qnorm(0.975)*sqrt(r_var_trunc) 
    r_upperCI_trunc_R1 <- r_est_trunc_R1+qnorm(0.975)*sqrt(r_var_trunc_R1) 
    r_upperCI_trunc_R0 <- r_est_trunc_R0+qnorm(0.975)*sqrt(r_var_trunc_R0) 
    r_est_trunc_diff <- colMeans(containers_trunc$est_temp_diff)
    r_est_trunc_R1_diff <- colMeans(containers_trunc$est_R1_temp_diff)
    r_est_trunc_R0_diff <- colMeans(containers_trunc$est_R0_temp_diff)
    r_var_trunc_diff <- colSums(containers_trunc$var_temp_diff)/(fold*n)
    r_var_trunc_R1_diff <- colSums(containers_trunc$var_R1_temp_diff)/(fold*n)
    r_var_trunc_R0_diff <- colSums(containers_trunc$var_R0_temp_diff)/(fold*n)
    r_lowerCI_trunc_diff <- r_est_trunc_diff-qnorm(0.975)*sqrt(r_var_trunc_diff)
    r_lowerCI_trunc_R1_diff <- r_est_trunc_R1_diff-qnorm(0.975)*sqrt(r_var_trunc_R1_diff)
    r_lowerCI_trunc_R0_diff <- r_est_trunc_R0_diff-qnorm(0.975)*sqrt(r_var_trunc_R0_diff)
    r_upperCI_trunc_diff <- r_est_trunc_diff+qnorm(0.975)*sqrt(r_var_trunc_diff) 
    r_upperCI_trunc_R1_diff <- r_est_trunc_R1_diff+qnorm(0.975)*sqrt(r_var_trunc_R1_diff) 
    r_upperCI_trunc_R0_diff <- r_est_trunc_R0_diff+qnorm(0.975)*sqrt(r_var_trunc_R0_diff) 
    
  }
  

  ## output final results
  if(IF_output){
    if(simple_trunc){
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     containers$IF=containers$IF, containers$IF_R1=containers$IF_R1, containers$IF_R0=containers$IF_R0, 
                     containers$IF_diff=containers$IF_diff, containers$IF_R1_diff=containers$IF_R1_diff, containers$IF_R0_diff=containers$IF_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     containers$id_list=containers$id_list)
    }else{
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     est_mat=containers$est_temp, est_R1_mat=containers$est_R1_temp, est_R0_mat=containers$est_R0_temp, 
                     est_trunc=r_est_trunc, est_trunc_R1=r_est_trunc_R1, est_trunc_R0=r_est_trunc_R0,
                     est_trunc_mat=containers_trunc$est_temp, est_trunc_R1_mat=containers_trunc$est_R1_temp, est_trunc_R0_mat=containers_trunc$est_R0_temp, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0,
                     var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                     est_diff_mat=containers$est_temp_diff, est_R1_diff_mat=containers$est_R1_temp_diff, est_R0_diff_mat=containers$est_R0_temp_diff, 
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     est_trunc_diff_mat=containers_trunc$est_temp_diff, est_R1_trunc_diff_mat=containers_trunc$est_R1_temp_diff, est_R0_trunc_diff_mat=containers_trunc$est_R0_temp_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R1_diff=r_lowerCI_trunc_R1_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff, 
                     upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R1_diff=r_upperCI_trunc_R1_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff, 
                     containers$IF=containers$IF, containers$IF_R1=containers$IF_R1, containers$IF_R0=containers$IF_R0,
                     containers_trunc$IF=containers_trunc$IF, containers_trunc$IF_R1=containers_trunc$IF_R1, containers_trunc$IF_R0=containers_trunc$IF_R0,
                     containers$IF_diff=containers$IF_diff, containers$IF_R1_diff=containers$IF_R1_diff, containers$IF_R0_diff=containers$IF_R0_diff,
                     containers_trunc$IF_diff=containers_trunc$IF_diff, containers_trunc$IF_R1_diff=containers_trunc$IF_R1_diff, containers_trunc$IF_R0_diff=containers_trunc$IF_R0_diff,
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1,
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     containers$id_list=containers$id_list)
    }
  }else{
    if(simple_trunc){
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff, 
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1, 
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     containers$id_list=containers$id_list)
    }else{
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     est_mat=containers$est_temp, est_R1_mat=containers$est_R1_temp, est_R0_mat=containers$est_R0_temp, 
                     est_trunc=r_est_trunc, est_trunc_R1=r_est_trunc_R1, est_trunc_R0=r_est_trunc_R0,
                     est_trunc_mat=containers_trunc$est_temp, est_trunc_R1_mat=containers_trunc$est_R1_temp, est_trunc_R0_mat=containers_trunc$est_R0_temp, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0,
                     var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                     est_diff_mat=containers$est_temp_diff, est_R1_diff_mat=containers$est_R1_temp_diff, est_R0_diff_mat=containers$est_R0_temp_diff, 
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     est_trunc_diff_mat=containers_trunc$est_temp_diff, est_R1_trunc_diff_mat=containers_trunc$est_R1_temp_diff, est_R0_trunc_diff_mat=containers_trunc$est_R0_temp_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                     upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff, 
                     lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R1_diff=r_lowerCI_trunc_R1_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff, 
                     upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R1_diff=r_upperCI_trunc_R1_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff, 
                     pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1,
                     pain_bq_reordered_R0=pain_bq_reordered_R0, 
                     containers$id_list=containers$id_list)
    }
  }
  
}



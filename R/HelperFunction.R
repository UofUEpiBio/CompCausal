gam.variables <- function(data){
  sapply(1:ncol(data), function(i){
    if (is.numeric(data[,i]) | is.integer(data[, i])){ paste("s", "(",names(data)[i], ")", sep="")}
    else{names(data[i]) }
  })
}

single.index.variables <- function(data){
  sapply(1:ncol(data), function(i){
    names(data[i])
  })
}

eq <- function(theta,z) {
  mean( pmin( (z)^2, theta^2 )/ theta^2) - log(length(z))/length(z)
}

## Induced estimates with single index model
counterfactual <- function(Y, M, R, X, t, trt, gamma, est, est_R1, est_R0){
  
  n <- length(t)
  Y[is.na(Y)] <- 0
  trt.ind <- as.numeric(t==trt) # create treatment indicator variable
  n_gamma <- length(gamma)
  
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model and 
  index.var.Y <- single.index.variables(X)
  X_adjust <- model.matrix(as.formula(paste("~", paste(index.var.Y, collapse = "+"))), data = X)[,-1]
  X_adjust_scale <- scale(X_adjust)
  
  X_adjust_scale_t_R0 <- X_adjust_scale[which(t==trt & R==0), ]
  X_adjust_scale_t_R1 <- X_adjust_scale[which(t==trt & R==1), ]
  X_R0 <- X[which(R==0), ]
  X_R1 <- X[which(R==1), ]
  Y_t_R0 <- Y[which(t==trt & R==0)]
  Y_t_R1 <- Y[which(t==trt & R==1)]
  M_t_R0 <- M[which(t==trt & R==0)]
  M_t_R1 <- M[which(t==trt & R==1)]
  trt.ind_R0 <- trt.ind[which(R==0)]
  trt.ind_R1 <- trt.ind[which(R==1)]
  
  ## fit outcome model
  requireNamespace('MAVE', quietly = TRUE)
  fit_t_R0_h <- fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R0[which(M_t_R0==1), ],
                                                         Y = Y_t_R0[which(M_t_R0==1)],
                                                         ids = 1:length(Y_t_R0[which(M_t_R0==1)]), 
                                                         kernel="dnorm", bw.selection="ise", bw.method="optim", use_mave=TRUE)
  
  fit_t_R1_h <- fit_SensIAT_single_index_norm1coef_model(X = X_adjust_scale_t_R1[which(M_t_R1==1), ],
                                                         Y = Y_t_R1[which(M_t_R1==1)],
                                                         ids = 1:length(Y_t_R1[which(M_t_R1==1)]), 
                                                         kernel="dnorm", bw.selection="ise", bw.method="optim", use_mave=TRUE)
  
  
  X_t_R0_beta_t_R0 <- as.vector(X_adjust_scale_t_R0[which(M_t_R0==1), ] %*% fit_t_R0_h$coef)
  X_beta_t_R0 <- as.vector(X_adjust_scale %*% fit_t_R0_h$coef)
  X_t_R1_beta_t_R1 <- as.vector(X_adjust_scale_t_R1[which(M_t_R1==1), ] %*% fit_t_R1_h$coef)
  X_beta_t_R1 <- as.vector(X_adjust_scale %*% fit_t_R1_h$coef)
  
  y_t_R0 = sort(unique(Y_t_R0[which(M_t_R0==1)]))    
  ny_t_R0 = length(y_t_R0) 
  F_X_t_R0 <- NW_new(Xb=X_t_R0_beta_t_R0, Y=Y_t_R0[which(M_t_R0==1)], 
                     xb=X_beta_t_R0, y=y_t_R0, h=fit_t_R0_h$bandwidth, 
                     kernel = "dnorm")
  y_t_R1 = sort(unique(Y_t_R1[which(M_t_R1==1)]))    
  ny_t_R1 = length(y_t_R1) 
  F_X_t_R1 <- NW_new(Xb=X_t_R1_beta_t_R1, Y=Y_t_R1[which(M_t_R1==1)], 
                     xb=X_beta_t_R1, y=y_t_R1, h=fit_t_R1_h$bandwidth, 
                     kernel = "dnorm")
  
  i1 = which(apply(F_X_t_R0==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X_beta_t_R0[i1], X_beta_t_R0[-i1], FUN = "-")), 1, which.min)
  F_X_t_R0[i1, ] <- F_X_t_R0[-i1, ][i1.closest, ]
  
  i1 = which(apply(F_X_t_R1==0,1,prod)==1)
  i1.closest <- apply(abs(outer(X_beta_t_R1[i1], X_beta_t_R1[-i1], FUN = "-")), 1, which.min)
  F_X_t_R1[i1, ] <- F_X_t_R1[-i1, ][i1.closest, ]
  
  dF_X_t_R0 <- F_X_t_R0[, -1, drop = FALSE]-F_X_t_R0[,-ny_t_R0, drop = FALSE]
  dF_X_t_R1 <- F_X_t_R1[, -1, drop = FALSE]-F_X_t_R1[,-ny_t_R1, drop = FALSE]
  
  ## conditional expectation of Y given R=0, T=t and X
  mu_X_t_R0 <- c(dF_X_t_R0 %*% (y_t_R0[-1]+y_t_R0[-ny_t_R0])/2)
  ## conditional expectation of Y given R=1, T=t and X
  mu_X_t_R1 <- c(dF_X_t_R1 %*% (y_t_R1[-1]+y_t_R1[-ny_t_R1])/2)
  
  ## fit g model
  g.fit <- mgcv::gam(as.formula(paste("R ~", gam.var)), data=cbind(R, X), family=binomial) 
  g1 <- stats::predict(g.fit, newdata=X, type="response") 
  g0 <- 1-g1
  
  ## fit treatment assignment model
  t_R0.fit <- mgcv::gam(as.formula(paste("trt.ind_R0 ~", gam.var)), data=X_R0, family=binomial) ## treatment model
  t_R1.fit <- mgcv::gam(as.formula(paste("trt.ind_R1 ~", gam.var)), data=X_R1, family=binomial) ## treatment model
  pi_R0 <- stats::predict(t_R0.fit, newdata=X, type="response") 
  pi_R1 <- stats::predict(t_R1.fit, newdata=X, type="response")
  
  ## empirical mean
  prop_t <- sum(trt.ind)/n
  prop_t_R0 <- mean(trt.ind[which(R==0)])
  prop_t_R1 <- mean(trt.ind[which(R==1)])
  prop_R1 <- sum(R)/n
  
  result_t0 <- (est-mean(mu_X_t_R0*pi_R0*g0+mu_X_t_R1*pi_R1*g1))/(1-prop_t)
  result_t <- mean(mu_X_t_R0*pi_R0*g0+mu_X_t_R1*pi_R1*g1)/prop_t

  result_R0_t0 <- (est_R0-mean(mu_X_t_R0*pi_R0*g0)/(1-prop_R1))/(1-prop_t_R0)
  result_R0_t <- mean(mu_X_t_R0*pi_R0*g0)/(prop_t_R0*(1-prop_R1))

  result_R1_t0 <- (est_R1-mean(mu_X_t_R1*pi_R1*g1)/prop_R1)/(1-prop_t_R1)
  result_R1_t <- mean(mu_X_t_R1*pi_R1*g1)/(prop_t_R1*prop_R1)

  return(data.frame(EY_t0=result_t0, EY_t=rep(result_t, n_gamma), 
                    EY_R0_t0=result_R0_t0, EY_R0_t=rep(result_R0_t, n_gamma), 
                    EY_R1_t0=result_R1_t0, EY_R1_t=rep(result_R1_t, n_gamma)))
  
}


## truth: beta regression
truth_beta <- function(Y, M, R, t, X, trt, gamma){
  
  truth_t <- c()
  truth_t_R1 <- c()
  truth_t_R0 <- c()
  
  n <- length(t)
  Y[is.na(Y)] <- 0
  trt.ind <- as.numeric(t==trt) # create treatment indicator variable
  
  gam.var <- paste(gam.variables(X), collapse = "+") ## gam variables for treatment assignment model and 
  index.var.Y <- single.index.variables(X)
  
  X_R0 <- X[which(R==0), ]
  X_R1 <- X[which(R==1), ]
  Y_t_R0 <- Y[which(t==trt & R==0)]
  Y_t_R1 <- Y[which(t==trt & R==1)]
  M_t_R0 <- M[which(t==trt & R==0)]
  M_t_R1 <- M[which(t==trt & R==1)]
  trt.ind_R0 <- trt.ind[which(R==0)]
  trt.ind_R1 <- trt.ind[which(R==1)]
  
  ## fit outcome model
  # rescale outcome to (0,1)
  eps <- 1e-6
  Y01 <- Y/100
  # nudge any 0/1 slightly inside (only needed if any exact 0 or 100)
  Y01[Y01 <= 0] <- eps
  Y01[Y01 >= 1] <- 1 - eps
  
  Y_t_R0.fit <- betareg::betareg(as.formula(paste("Y01 ~", paste(index.var.Y, collapse = "+"))), data=cbind(Y01, X)[which(t==trt & R==0 & M==1), ])
  Y_t_R1.fit <- betareg::betareg(as.formula(paste("Y01 ~", paste(index.var.Y, collapse = "+"))), data=cbind(Y01, X)[which(t==trt & R==1 & M==1), ])
  
  phi_t_R0 <- Y_t_R0.fit$coefficients$precision
  phi_t_R1 <- Y_t_R1.fit$coefficients$precision
  
  mu_X_t_R0 <- 100*betareg::predict(Y_t_R0.fit, newdata=X, type="response")
  mu_X_t_R1 <- 100*betareg::predict(Y_t_R1.fit, newdata=X, type="response")
  
  ## fit g model
  g.fit <- mgcv::gam(as.formula(paste("R ~", gam.var)), data=cbind(R, X), family=binomial) 
  g1 <- stats::predict(g.fit, newdata=X, type="response") 
  g0 <- 1-g1
  
  ## fit treatment assignment model
  t_R0.fit <- mgcv::gam(as.formula(paste("trt.ind_R0 ~", gam.var)), data=X_R0, family=binomial) ## treatment model
  pi_R0 <- stats::predict(t_R0.fit, newdata=X, type="response") 
  
  for(g in 1:length(gamma)){
    
    mu_Yexp_X_t_R0 <- mapply(function(a,b) {
      integrate(function(y) 100*y * exp(gamma[g]*pnorm((100*y-60)/25)) * dbeta(y, shape1 = a, shape2 = b), lower = 0, upper = 1)$value
    }, mu_X_t_R0/100*phi_t_R0, (1-mu_X_t_R0/100)*phi_t_R0)
    
    mu_exp_X_t_R0 <- mapply(function(a,b) {
      integrate(function(y) exp(gamma[g]*pnorm((100*y-60)/25)) * dbeta(y, shape1 = a, shape2 = b), lower = 0, upper = 1)$value
    }, mu_X_t_R0/100*phi_t_R0, (1-mu_X_t_R0/100)*phi_t_R0)
    
    
    truth_t <- c(truth_t, mean(g1*mu_X_t_R1+g0*(mu_X_t_R0*pi_R0+mu_Yexp_X_t_R0/mu_exp_X_t_R0*(1-pi_R0))))
    truth_t_R0 <- c(truth_t_R0, mean(g0*(mu_X_t_R0*pi_R0+mu_Yexp_X_t_R0/mu_exp_X_t_R0*(1-pi_R0)))/mean(g0))

  }
  
  truth_t_R1 <- mean(mu_X_t_R1*g1)/mean(g1)

  return(list(t=truth_t, t_R0=truth_t_R0, t_R1=truth_t_R1))
  
}






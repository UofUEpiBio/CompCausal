# Methods ----------------------------------------------------------------------

#' @rdname est_psi
#' @export
print.est_psi <- function(object, rounding=3, ...){
  
  if(object$simple_trunc){
    
    df <- data.frame(gamma=object$gamma, Estimates=round(object$est, rounding), Var=round(object$var, rounding), 
                     Lower_95CI=round(object$lowerCI, rounding), Upper_95CI=round(object$upperCI, rounding))
    
    df_R1 <- data.frame(Estimates=round(object$est_R1[1], rounding), Var=round(object$var_R1[1], rounding), 
                        Lower_95CI=round(object$lowerCI_R1[1], rounding), Upper_95CI=round(object$upperCI_R1[1], rounding))
    
    df_R0 <- data.frame(gamma=object$gamma, Estimates=round(object$est_R0, rounding), Var=round(object$var_R0, rounding), 
                        Lower_95CI=round(object$lowerCI_R0, rounding), Upper_95CI=round(object$upperCI_R0, rounding))
    
  }else{
    
    df <- data.frame(gamma=object$gamma, Estimates=round(object$est_trunc, rounding), Var=round(object$var_trunc, rounding), 
                     Lower_95CI=round(object$lowerCI_trunc, rounding), Upper_95CI=round(object$upperCI_trunc, rounding))
    
    df_R1 <- data.frame(Estimates=round(object$est_trunc_R1[1], rounding), Var=round(object$var_trunc_R1[1], rounding), 
                        Lower_95CI=round(object$lowerCI_trunc_R1[1], rounding), Upper_95CI=round(object$upperCI_trunc_R1[1], rounding))
    
    df_R0 <- data.frame(gamma=object$gamma, Estimates=round(object$est_trunc_R0, rounding), Var=round(object$var_trunc_R0, rounding), 
                        Lower_95CI=round(object$lowerCI_trunc_R0, rounding), Upper_95CI=round(object$upperCI_trunc_R0, rounding))
    
  }

  cat(sprintf("Estimation of E[Y(%s)]\n", object$trt))
  cat("=========================\n")
  print.data.frame(df, row.names = FALSE)
  cat("\n")
  
  cat(sprintf("Estimation of E[Y(%s)|R=1]\n", object$trt))
  cat("=========================\n")
  print.data.frame(df_R1, row.names = FALSE)
  cat("\n")
  
  cat(sprintf("Estimation of E[Y(%s)|R=0]\n", object$trt))
  cat("=========================\n")
  print.data.frame(df_R0, row.names = FALSE)
  cat("\n")
  
  invisible(object)
  
}




#' @rdname est_psi
#' @export
print_effects <- function(object_t1, object_t0, rounding=3, ...){
  
  n_gamma1 <- length(object_t1$gamma)
  n_gamma0 <- length(object_t0$gamma)
  
  ## data frame for results
  res_out_diff <- expand.grid(gamma1=object_t1$gamma, gamma0=object_t0$gamma, Type=c("CCCE", "PPCE"))
  res_out_diff <- rbind(res_out_diff, data.frame(gamma1=NA, gamma0=NA, Type="RTCE"))
  res_out_diff$Estimates <- 0
  res_out_diff$Var <- 0
  
  ## compute variance
  if((!is.null(object_t1$IF))&(!is.null(object_t0$IF))){
    scale_n <- 1 / length(object_t1$IF[[1]])
  }else{
    stop("Cannot estimate the variance of the treatment effect without influence functions for both treatment groups.")
  }
  
  fold_idx <- split(seq_along(object_t1$fold_index_l), object_t1$fold_index_l)
  
  if((!object_t1$simple_trunc)&(!object_t0$simple_trunc)){
    
    t1_IF_mat <- t(do.call(rbind, object_t1$IF_trunc))
    t1_IF_R0_mat <- t(do.call(rbind, object_t1$IF_trunc_R0))
    t1_IF_R1_mat <- t(do.call(rbind, object_t1$IF_trunc_R1))
    
    res_out_diff$Estimates[which(res_out_diff$type=="RTCE")] <- round(object_t1$est_trunc_R1[1]-object_t0$est_trunc_R1[1], rounding)
    
    IF_R1_diff <- t1_IF_R1_mat - object_t0$IF_trunc_R1[[1]]
    var_R1_temp <- vapply(fold_idx, function(id) colVars(IF_R1_diff[id, ]), numeric(n_gamma1))
    res_out_diff$Var[which(res_out_diff$Type=="RTCE")] <- round(rowMeans(var_R1_temp)[1]*scale_n, rounding)
    
    for (g_0 in seq_along(object_t0$gamma)) {
      
      IF_diff    <- t1_IF_mat - object_t0$IF_trunc[[g_0]]
      IF_R0_diff <- t1_IF_R0_mat - object_t0$IF_trunc_R0[[g_0]]
      
      var_temp    <- vapply(fold_idx, function(id) colVars(IF_diff[id, ]), numeric(n_gamma1))
      var_R0_temp <- vapply(fold_idx, function(id) colVars(IF_R0_diff[id, ]), numeric(n_gamma1))
      
      indx_CCCE <- which(res_out_diff$gamma0==object_t0$gamma[g_0]&res_out_diff$Type=="CCCE")
      res_out_diff$Estimates[indx_CCCE] <- round(object_t1$est_trunc-object_t0$est_trunc, rounding)
      res_out_diff$Var[indx_CCCE] <- round(rowMeans(var_temp)*scale_n, rounding)
      
      indx_PPCE <- which(res_out_diff$gamma0==object_t0$gamma[g_0]&res_out_diff$Type=="PPCE")
      res_out_diff$Estimates[indx_PPCE] <- round(object_t1$est_trunc_R0-object_t0$est_trunc_R0, rounding)
      res_out_diff$Var[indx_PPCE] <- round(rowMeans(var_R0_temp)*scale_n, rounding)
      
    }
    
    
  }else{
    
    t1_IF_mat <- t(do.call(rbind, object_t1$IF))
    t1_IF_R0_mat <- t(do.call(rbind, object_t1$IF_R0))
    t1_IF_R1_mat <- t(do.call(rbind, object_t1$IF_R1))
    
    res_out_diff$Estimates[which(res_out_diff$Type=="RTCE")] <- round(object_t1$est_R1[1]-object_t0$est_R1[1], rounding)
    
    IF_R1_diff <- t1_IF_R1_mat - object_t0$IF_R1[[1]]
    var_R1_temp <- vapply(fold_idx, function(id) colVars(IF_R1_diff[id, ]), numeric(n_gamma1))
    res_out_diff$Var[which(res_out_diff$Type=="RTCE")] <- round(rowMeans(var_R1_temp)[1]*scale_n, rounding)
    
    for (g_0 in seq_along(object_t0$gamma)) {
      
      IF_diff    <- t1_IF_mat - object_t0$IF[[g_0]]
      IF_R0_diff <- t1_IF_R0_mat - object_t0$IF_R0[[g_0]]
      
      var_temp    <- vapply(fold_idx, function(id) colVars(IF_diff[id, ]),    numeric(n_gamma1))
      var_R0_temp <- vapply(fold_idx, function(id) colVars(IF_R0_diff[id, ]), numeric(n_gamma1))
      
      indx_CCCE <- which(res_out_diff$gamma0==object_t0$gamma[g_0]&res_out_diff$Type=="CCCE")
      res_out_diff$Estimates[indx_CCCE] <- round(object_t1$est-object_t0$est, rounding)
      res_out_diff$Var[indx_CCCE] <- round(rowMeans(var_temp)*scale_n, rounding)
      
      indx_PPCE <- which(res_out_diff$gamma0==object_t0$gamma[g_0]&res_out_diff$Type=="PPCE")
      res_out_diff$Estimates[indx_PPCE] <- round(object_t1$est_R0-object_t0$est_R0, rounding)
      res_out_diff$Var[indx_PPCE] <- round(rowMeans(var_R0_temp)*scale_n, rounding)
      
    }
  
  }
  
  ## 95% CI
  res_out_diff$lowerCI <- round(res_out_diff$Estimates-qnorm(0.975)*sqrt(res_out_diff$Var), rounding)
  res_out_diff$upperCI <- round(res_out_diff$Estimates+qnorm(0.975)*sqrt(res_out_diff$Var), rounding)
  
  ## output table
  cat("Estimation of CCCE, PPCE, RTCE\n")
  cat("=========================\n")
  print.data.frame(res_out_diff, row.names = FALSE)
  cat("\n")
  
}





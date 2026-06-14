# Methods ----------------------------------------------------------------------

#' @rdname est_psi
#' @export
print.est_psi <- function(object, ...){
  
  if(object$simple_trunc){
    
    df <- data.frame(gamma=object$gamma, Estimates=object$est, Var=object$var, 
                     Lower_95CI=object$lowerCI, Upper_95CI=object$upperCI)
    
    df_R1 <- data.frame(Estimates=object$est_R1[1], Var=object$var_R1[1], 
                        Lower_95CI=object$lowerCI_R1[1], Upper_95CI=object$upperCI_R1[1])
    
    df_R0 <- data.frame(gamma=object$gamma, Estimates=object$est_R0, Var=object$var_R0, 
                        Lower_95CI=object$lowerCI_R0, Upper_95CI=object$upperCI_R0)
    
  }else{
    
    df <- data.frame(gamma=object$gamma, Estimates=object$est_trunc, Var=object$var_trunc, 
                     Lower_95CI=object$lowerCI_trunc, Upper_95CI=object$upperCI_trunc)
    
    df_R1 <- data.frame(Estimates=object$est_trunc_R1[1], Var=object$var_trunc_R1[1], 
                        Lower_95CI=object$lowerCI_trunc_R1[1], Upper_95CI=object$upperCI_trunc_R1[1])
    
    df_R0 <- data.frame(gamma=object$gamma, Estimates=object$est_trunc_R0, Var=object$var_trunc_R0, 
                        Lower_95CI=object$lowerCI_trunc_R0, Upper_95CI=object$upperCI_trunc_R0)
    
  }

  cat(sprintf("Estimation of E[Y(%s)]\n", object$trt))
  cat("================\n")
  print.data.frame(df, row.names = FALSE)
  cat("\n")
  
  cat(sprintf("Estimation of E[Y(%s)|R=1]\n", object$trt))
  cat("================\n")
  print.data.frame(df_R1, row.names = FALSE)
  cat("\n")
  
  cat(sprintf("Estimation of E[Y(%s)|R=0]\n", object$trt))
  cat("================\n")
  print.data.frame(df_R0, row.names = FALSE)
  cat("\n")
  
  invisible(object)
  
}











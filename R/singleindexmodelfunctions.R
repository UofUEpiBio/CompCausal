#' single index model's function
#' return the initial value for beta
#' @noRd
cumuSIR <- function(X, Y, eps = 1e-7){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  
  Y.CP <- matrix(
    Y[rep(1:number_n, times = number_n), ]<=
      Y[rep(1:number_n, each = number_n), ],
    nrow = number_n, ncol = number_n
  )
  
  # centralizing covariates
  X.cs <- t(t(X)-colMeans(X))
  
  # calculating m(y)=\E[X_i 1(Y_i\leq y)]
  m.y <- t(X.cs) %*% Y.CP/number_n
  # calculating K=\E[m(Y_i)m(Y_i)^T]
  Km <- m.y %*% t(m.y)/number_n
  
  Bhat <- eigen(solve(var(X) + eps*diag(number_p), Km))$vectors[, 1]
  
  return(Bhat)
}


#' find the optimal bandwidth, default set for kernel is 
#' @noRd
SIDRnew_h2 <- function(X, Y,
                       Y.CP = NULL,
                       beta = NULL,
                       kernel = "K2_Biweight",
                       method = "optim",
                       optim_method = "BFGS",
                       abs.tol = 1e-4,
                       initial_bandwidth = NULL,
                       wi.boot = NULL, penalty=0)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  number_n <- dim(X)[1]
  number_p <- dim(X)[2]
  
  if (is.null(beta))
  {
    beta <- c(1, rep(0, number_p-1))
  }else
  {
    beta <- as.vector(beta)
    beta <- beta/beta[1]
  }
  
  if (is.null(initial_bandwidth))
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h)
          b <- beta
          h <- exp(parameter)
          xb <- c(X%*%b) 
          y <- Y
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(xb-xb[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2+penalty*h)
        }
        
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
        #                                  h = h))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K2_Biweight kernel.")
      }
    }else if (kernel == "dnorm")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) dnorm(x/h, 0, 1)
          
          b <- beta
          h <- exp(parameter)
          
          x <- c(X%*%b) 
          y <- Y
          
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = h))^2)
        #   return(cv)
        # }
      }else
      {
        # wi.boot <- as.vector(wi.boot)
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = h, w = wi.boot))^2)
        #   return(cv)
        # }
        stop("There's no weighted version of the dnorm kernel.")
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) 
          
          b <- beta
          h <- exp(parameter)
          
          x <- c(X%*%b) 
          y <- Y
          
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.bh <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   h <- exp(parameter[number_p])
        #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
        #                                            h = h), 0), 1))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K4_Biweight kernel.")
      }
    }
    
    if(method == "nlminb")
    {
      esti <- nlminb(start = 0, 
                     objective = Eij3,
                     control = list(abs.tol = abs.tol))
      results <- list(coef = beta,
                      bandwidth = exp(esti$par),
                      details = esti)
    }else if (method == "optim")
    {
      # the new optimize function using optim, you can change the lower and upper
      esti <- optim(par = 0, 
                    fn = Eij3,
                    method = optim_method,
                    control = list(abstol = abs.tol))
      results <- list(coef = beta,
                      bandwidth = exp(esti$par),
                      details = esti)
    }
    else if (method == "optimize")
    {
      #   # esti <- nmk(par = 0, 
      #   #             fn = Eij3,
      #   #             control = list(tol = abs.tol))
      esti <- optimize(Eij3,
                       interval = c(-10, 10))
      results <- list(coef = beta,
                      bandwidth = exp(esti$minimum),
                      details = esti)
    }
    
  }else
  {
    if (kernel=="K2_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 15/16*(1-(x/h)^2)^2 * (abs(x) <= h)
          b <- beta
          h <- exp(parameter)
          xb <- c(X%*%b) 
          y <- Y
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(xb-xb[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2+penalty*h)
        }
        
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_K2B_rcpp(X = X %*% b, Y = Y.CP,
        #                                  h = bandwidth))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K2_Biweight kernel.")
      }
    }else if (kernel=="dnorm") 
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) dnorm(x/h,0,1) 
          
          b <- beta
          h <- exp(parameter)
          
          x <- c(X%*%b) 
          y <- Y
          
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_dnorm_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = bandwidth))^2)
        #   return(cv)
        # }
      }else
      {
        # wi.boot <- as.vector(wi.boot)
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-NWcv_K2B_w_rcpp(X = X %*% b, Y = Y.CP,
        #                                    h = bandwidth, w = wi.boot))^2)
        #   return(cv)
        # }
        stop("There's no weighted version of the dnorm kernel.")
      }
    }else if (kernel=="K4_Biweight")
    {
      if (is.null(wi.boot))
      {
        Eij3 <- function(parameter){
          K <- function(x, h) 105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) 
          
          b <- beta
          h <- exp(parameter)
          
          x <- c(X%*%b) 
          y <- Y
          
          n <- length(y)
          yo <- order(y)
          ys <- y[yo]
          uy <- rle(ys)[[1]]
          cols <- cumsum(uy)
          ei <- rep(0, n)
          for (i in 1:n){
            Kih <- K(x-x[i],h=h)
            Kih[i] <- 0
            denom <- sum(Kih)
            ei[i] <- sum(uy*(1*(y[i] <= ys)[cols] - (denom != 0)* cumsum(Kih[yo])[cols] / (denom + (denom == 0)))^2)
          }
          return(sum(ei)/n^2)
        }
        # cv.b <- function(parameter)
        # {
        #   b <- c(1, parameter[1:(number_p-1)])
        #   cv <- mean((Y.CP-pmin(pmax(NWcv_K4B_rcpp(X = X %*% b, Y = Y.CP,
        #                                            h = bandwidth), 0), 1))^2)
        #   return(cv)
        # }
      }else
      {
        stop("There's no weighted version of the K4_Biweight kernel.")
      }
    }
    
    if(method == "nlminb")
    {
      esti <- nlminb(start = initial_bandwidth, 
                     objective = Eij3,
                     control = list(abs.tol = abs.tol))
      results <- list(coef = beta,
                      bandwidth = exp(esti$par),
                      details = esti)
    }else if (method == "optim")
    {
      # the new optimize function using optim, you can change the lower and upper
      esti <- optim(par = initial_bandwidth, 
                    fn = Eij3,
                    method = optim_method,
                    control = list(abstol = abs.tol))
      results <- list(coef = beta,
                      bandwidth = exp(esti$par),
                      details = esti)
    }
    else if (method == "optimize")
    {
      #   # esti <- nmk(par = initial_bandwidth, 
      #   #             fn = Eij3,
      #   #             control = list(tol = abs.tol))
      esti <- optimize(Eij3,
                       interval = c(-10, 10))
      results <- list(coef = beta,
                      bandwidth = exp(esti$minimum),
                      details = esti)
      
    }
  }
  
  return(results)
}


# return the prediction using "K2_Biweight" kernel
NW_new <- function(Xb, Y, xb, y, h, kernel = "K2_Biweight"){
  
  if(kernel == "dnorm"){
    K <- function(x, h){dnorm(x/h, 0, 1)} # Gaussian 
  } else if(kernel == "K2_Biweight"){
    K <- function(x, h){15/16*(1-(x/h)^2)^2 * (abs(x) <= h)} # K2_Biweight
  } else if(kernel=="K4_Biweight"){
    K <- function(x, h){105/64*(1-3*((x/h)^2))*(1-(x/h)^2)^2 * (abs(x) <= h) }# K4_Biweight
  }
  
  Kxb <- sapply(xb, function(x, Xb) K(Xb-x, h), Xb=Xb)
  
  Ylty <- sapply(y, function(x, Y) 1*(Y <= x), Y=Y)
  
  denom <- colSums(Kxb)
  
  fyxb <- (denom!=0)*crossprod(Kxb, Ylty)/(denom + (denom==0))
  
  return(fyxb)
  
}


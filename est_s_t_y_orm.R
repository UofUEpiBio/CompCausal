
#' Create containers for cross-fitted ORM estimator
#' 
#' This helper function initializes an environment to hold influence-function
#' vectors, point estimates, variance estimates, and other intermediate 
#' quantities for the cross-fitted ordinal regression model (ORM) estimator. 
#' The structure is designed to facilitate storage and retrieval of these 
#' objects across multiple folds and sensitivity parameter values.
#' @noRd
est_s_t_y_orm_create_containers <- function(gamma, fold) {
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

#' Influence-function variance estimates and multiplier bootstrap quantiles
#'
#' Computes sandwich-type variance estimates from influence functions and,
#' optionally, runs a multiplier bootstrap for truncated influence functions to
#' produce bias estimates and pivotal quantiles.
#'
#' @param x For \code{var_estimator.default}: a list of overall influence-function
#'   vectors, one per element of \code{gamma} (i.e. the \code{IF} field of an
#'   \code{est_container}).  For \code{var_estimator.est_container}: the container
#'   environment itself — all IF fields are extracted automatically.
#' @param IF_R1,IF_R0 Lists of R1- and R0-stratum influence functions.
#' @param IF_diff,IF_R1_diff,IF_R0_diff Lists of difference-estimand influence
#'   functions (overall, R1-stratum, R0-stratum).
#' @param vector_R1,vector_R0 Numeric vectors of fold-level R1/R0 weights used
#'   to centre the stratum-specific influence functions.
#' @param n Total sample size.
#' @param IF_trunc,IF_R1_trunc,IF_R0_trunc,IF_trunc_diff,IF_R1_trunc_diff,IF_R0_trunc_diff
#'   Optional lists of truncated influence functions.  When supplied, truncated
#'   variance estimates and multiplier bootstrap quantiles are also returned.
#' @param r_est_trunc,r_est_trunc_R1,r_est_trunc_R0,r_est_trunc_diff,r_est_trunc_R1_diff,r_est_trunc_R0_diff
#'   Numeric vectors of truncated point estimates (length = number of gamma
#'   values), used as bootstrap centering constants.  Required when truncated IF
#'   arguments are supplied; computed automatically by the
#'   \code{est_container} method.
#' @param B Number of multiplier bootstrap replicates (default \code{2000}).
#' @param trunc For \code{var_estimator.est_container}: an \code{est_container}
#'   holding the truncated influence functions, or \code{NULL} (default) to skip
#'   the bootstrap.
#' @param ... Further arguments passed to \code{var_estimator.default}.
#'
#' @return A named list.  Always contains:
#'   \describe{
#'     \item{\code{var_2}, \code{var_R1_2}, \code{var_R0_2},
#'           \code{var_diff_2}, \code{var_R1_diff_2}, \code{var_R0_diff_2}}{%
#'       Sandwich variance estimates from the primary influence functions.}
#'   }
#'   When truncated influence functions are supplied, additionally contains:
#'   \describe{
#'     \item{\code{var_trunc_2}, \ldots, \code{var_trunc_R0_diff_2}}{%
#'       Sandwich variance estimates from the truncated influence functions.}
#'     \item{\code{var_trunc_star}, \ldots, \code{var_trunc_R0_diff_star}}{%
#'       Multiplier bootstrap variance estimates.}
#'     \item{\code{bias_star}, \ldots, \code{bias_R0_diff_star}}{%
#'       Multiplier bootstrap bias estimates.}
#'     \item{\code{q}, \code{q_R1}, \code{q_R0}, \code{q_diff},
#'           \code{q_R1_diff}, \code{q_R0_diff}}{%
#'       95th-percentile bootstrap quantiles of the pivotal statistic.}
#'   }
var_estimator <- function(x, ...) UseMethod("var_estimator")

#' @rdname var_estimator
var_estimator.default <- function(
  x,
  IF_R1, IF_R0,
  IF_diff, IF_R1_diff, IF_R0_diff,
  vector_R1, vector_R0,
  n,
  IF_trunc = NULL, IF_R1_trunc = NULL, IF_R0_trunc = NULL,
  IF_trunc_diff = NULL, IF_R1_trunc_diff = NULL, IF_R0_trunc_diff = NULL,
  r_est_trunc = NULL, r_est_trunc_R1 = NULL, r_est_trunc_R0 = NULL,
  r_est_trunc_diff = NULL, r_est_trunc_R1_diff = NULL, r_est_trunc_R0_diff = NULL,
  B = 2000,
  ...
) {
  IF       <- x
  ng       <- length(IF)
  do_trunc <- !is.null(IF_trunc)

  r_var_2 <- r_var_R1_2 <- r_var_R0_2 <-
    r_var_diff_2 <- r_var_R1_diff_2 <- r_var_R0_diff_2 <- numeric(ng)

  if (do_trunc) {
    r_var_trunc_2 <- r_var_trunc_R1_2 <- r_var_trunc_R0_2 <-
      r_var_trunc_diff_2 <- r_var_trunc_R1_diff_2 <- r_var_trunc_R0_diff_2 <- numeric(ng)
    r_var_trunc_star <- r_var_trunc_R1_star <- r_var_trunc_R0_star <-
      r_var_trunc_diff_star <- r_var_trunc_R1_diff_star <- r_var_trunc_R0_diff_star <- numeric(ng)
    r_bias_star <- r_bias_R1_star <- r_bias_R0_star <-
      r_bias_diff_star <- r_bias_R1_diff_star <- r_bias_R0_diff_star <- numeric(ng)
    q <- q_R1 <- q_R0 <- q_diff <- q_R1_diff <- q_R0_diff <- numeric(ng)
  }

  for (g in seq_len(ng)) {
    r_var_2[g]         <- sum((IF[[g]]          - mean(IF[[g]]))^2)                            / (n * (n - 1))
    r_var_R1_2[g]      <- sum((IF_R1[[g]]        - vector_R1 * mean(IF_R1[[g]]))^2)             / (n * (n - 1))
    r_var_R0_2[g]      <- sum((IF_R0[[g]]        - vector_R0 * mean(IF_R0[[g]]))^2)             / (n * (n - 1))
    r_var_diff_2[g]    <- sum((IF_diff[[g]]      - mean(IF_diff[[g]]))^2)                       / (n * (n - 1))
    r_var_R1_diff_2[g] <- sum((IF_R1_diff[[g]]   - vector_R1 * mean(IF_R1_diff[[g]]))^2)       / (n * (n - 1))
    r_var_R0_diff_2[g] <- sum((IF_R0_diff[[g]]   - vector_R0 * mean(IF_R0_diff[[g]]))^2)       / (n * (n - 1))

    if (do_trunc) {
      r_var_trunc_2[g]         <- sum((IF_trunc[[g]]         - mean(IF_trunc[[g]]))^2)                          / (n * (n - 1))
      r_var_trunc_R1_2[g]      <- sum((IF_R1_trunc[[g]]       - vector_R1 * mean(IF_R1_trunc[[g]]))^2)          / (n * (n - 1))
      r_var_trunc_R0_2[g]      <- sum((IF_R0_trunc[[g]]       - vector_R0 * mean(IF_R0_trunc[[g]]))^2)          / (n * (n - 1))
      r_var_trunc_diff_2[g]    <- sum((IF_trunc_diff[[g]]     - mean(IF_trunc_diff[[g]]))^2)                    / (n * (n - 1))
      r_var_trunc_R1_diff_2[g] <- sum((IF_R1_trunc_diff[[g]]  - vector_R1 * mean(IF_R1_trunc_diff[[g]]))^2)    / (n * (n - 1))
      r_var_trunc_R0_diff_2[g] <- sum((IF_R0_trunc_diff[[g]]  - vector_R0 * mean(IF_R0_trunc_diff[[g]]))^2)    / (n * (n - 1))

      wi <- matrix(sample(c(-1, 1), n * B, replace = TRUE), nrow = n, ncol = B)

      est_star         <- r_est_trunc[g]          + colMeans(wi * IF_trunc[[g]])
      est_R1_star      <- r_est_trunc_R1[g]       + colMeans(wi * (IF_R1_trunc[[g]]      - vector_R1 * mean(IF_R1_trunc[[g]])))
      est_R0_star      <- r_est_trunc_R0[g]       + colMeans(wi * (IF_R0_trunc[[g]]      - vector_R0 * mean(IF_R0_trunc[[g]])))
      est_diff_star    <- r_est_trunc_diff[g]     + colMeans(wi * (IF_trunc_diff[[g]]    - mean(IF_trunc_diff[[g]])))
      est_R1_diff_star <- r_est_trunc_R1_diff[g]  + colMeans(wi * (IF_R1_trunc_diff[[g]] - vector_R1 * mean(IF_R1_trunc_diff[[g]])))
      est_R0_diff_star <- r_est_trunc_R0_diff[g]  + colMeans(wi * (IF_R0_trunc_diff[[g]] - vector_R0 * mean(IF_R0_trunc_diff[[g]])))

      r_bias_star[g]         <- mean(est_star)         - r_est_trunc[g]
      r_bias_R1_star[g]      <- mean(est_R1_star)      - r_est_trunc_R1[g]
      r_bias_R0_star[g]      <- mean(est_R0_star)      - r_est_trunc_R0[g]
      r_bias_diff_star[g]    <- mean(est_diff_star)    - r_est_trunc_diff[g]
      r_bias_R1_diff_star[g] <- mean(est_R1_diff_star) - r_est_trunc_R1_diff[g]
      r_bias_R0_diff_star[g] <- mean(est_R0_diff_star) - r_est_trunc_R0_diff[g]

      var_star         <- rowSums((t(wi * IF_trunc[[g]])         - colMeans(wi * IF_trunc[[g]]))^2)         / (n * (n - 1))
      var_R1_star      <- rowSums((t(wi * IF_R1_trunc[[g]])      - colMeans(wi * IF_R1_trunc[[g]]))^2)      / (n * (n - 1))
      var_R0_star      <- rowSums((t(wi * IF_R0_trunc[[g]])      - colMeans(wi * IF_R0_trunc[[g]]))^2)      / (n * (n - 1))
      var_diff_star    <- rowSums((t(wi * IF_trunc_diff[[g]])    - colMeans(wi * IF_trunc_diff[[g]]))^2)    / (n * (n - 1))
      var_R1_diff_star <- rowSums((t(wi * IF_R1_trunc_diff[[g]]) - colMeans(wi * IF_R1_trunc_diff[[g]]))^2) / (n * (n - 1))
      var_R0_diff_star <- rowSums((t(wi * IF_R0_trunc_diff[[g]]) - colMeans(wi * IF_R0_trunc_diff[[g]]))^2) / (n * (n - 1))

      r_var_trunc_star[g]         <- var(est_star)
      r_var_trunc_R1_star[g]      <- var(est_R1_star)
      r_var_trunc_R0_star[g]      <- var(est_R0_star)
      r_var_trunc_diff_star[g]    <- var(est_diff_star)
      r_var_trunc_R1_diff_star[g] <- var(est_R1_diff_star)
      r_var_trunc_R0_diff_star[g] <- var(est_R0_diff_star)

      t_star         <- (est_star         - r_est_trunc[g])          / sqrt(var_star)
      t_R1_star      <- (est_R1_star      - r_est_trunc_R1[g])       / sqrt(var_R1_star)
      t_R0_star      <- (est_R0_star      - r_est_trunc_R0[g])       / sqrt(var_R0_star)
      t_diff_star    <- (est_diff_star    - r_est_trunc_diff[g])     / sqrt(var_diff_star)
      t_R1_diff_star <- (est_R1_diff_star - r_est_trunc_R1_diff[g])  / sqrt(var_R1_diff_star)
      t_R0_diff_star <- (est_R0_diff_star - r_est_trunc_R0_diff[g])  / sqrt(var_R0_diff_star)

      q[g]         <- quantile(abs(t_star),         1 - 0.05)
      q_R1[g]      <- quantile(abs(t_R1_star),      1 - 0.05)
      q_R0[g]      <- quantile(abs(t_R0_star),      1 - 0.05)
      q_diff[g]    <- quantile(abs(t_diff_star),    1 - 0.05)
      q_R1_diff[g] <- quantile(abs(t_R1_diff_star), 1 - 0.05)
      q_R0_diff[g] <- quantile(abs(t_R0_diff_star), 1 - 0.05)
    }
  }

  out <- list(
    var_2         = r_var_2,
    var_R1_2      = r_var_R1_2,
    var_R0_2      = r_var_R0_2,
    var_diff_2    = r_var_diff_2,
    var_R1_diff_2 = r_var_R1_diff_2,
    var_R0_diff_2 = r_var_R0_diff_2
  )
  if (do_trunc)
    out <- c(out, list(
      var_trunc_2          = r_var_trunc_2,
      var_trunc_R1_2       = r_var_trunc_R1_2,
      var_trunc_R0_2       = r_var_trunc_R0_2,
      var_trunc_diff_2     = r_var_trunc_diff_2,
      var_trunc_R1_diff_2  = r_var_trunc_R1_diff_2,
      var_trunc_R0_diff_2  = r_var_trunc_R0_diff_2,
      var_trunc_star       = r_var_trunc_star,
      var_trunc_R1_star    = r_var_trunc_R1_star,
      var_trunc_R0_star    = r_var_trunc_R0_star,
      var_trunc_diff_star  = r_var_trunc_diff_star,
      var_trunc_R1_diff_star = r_var_trunc_R1_diff_star,
      var_trunc_R0_diff_star = r_var_trunc_R0_diff_star,
      bias_star            = r_bias_star,
      bias_R1_star         = r_bias_R1_star,
      bias_R0_star         = r_bias_R0_star,
      bias_diff_star       = r_bias_diff_star,
      bias_R1_diff_star    = r_bias_R1_diff_star,
      bias_R0_diff_star    = r_bias_R0_diff_star,
      q         = q,
      q_R1      = q_R1,
      q_R0      = q_R0,
      q_diff    = q_diff,
      q_R1_diff = q_R1_diff,
      q_R0_diff = q_R0_diff
    ))
  out
}

#' @rdname var_estimator
var_estimator.est_container <- function(x, n, trunc = NULL, B = 2000, ...) {
  do_trunc <- !is.null(trunc)
  var_estimator.default(
    x          = x$IF,
    IF_R1      = x$IF_R1,
    IF_R0      = x$IF_R0,
    IF_diff    = x$IF_diff,
    IF_R1_diff = x$IF_R1_diff,
    IF_R0_diff = x$IF_R0_diff,
    vector_R1  = x$vector_R1,
    vector_R0  = x$vector_R0,
    n          = n,
    IF_trunc          = if (do_trunc) trunc$IF          else NULL,
    IF_R1_trunc       = if (do_trunc) trunc$IF_R1       else NULL,
    IF_R0_trunc       = if (do_trunc) trunc$IF_R0       else NULL,
    IF_trunc_diff     = if (do_trunc) trunc$IF_diff     else NULL,
    IF_R1_trunc_diff  = if (do_trunc) trunc$IF_R1_diff  else NULL,
    IF_R0_trunc_diff  = if (do_trunc) trunc$IF_R0_diff  else NULL,
    r_est_trunc          = if (do_trunc) colMeans(trunc$est_temp)         else NULL,
    r_est_trunc_R1       = if (do_trunc) colMeans(trunc$est_R1_temp)      else NULL,
    r_est_trunc_R0       = if (do_trunc) colMeans(trunc$est_R0_temp)      else NULL,
    r_est_trunc_diff     = if (do_trunc) colMeans(trunc$est_temp_diff)    else NULL,
    r_est_trunc_R1_diff  = if (do_trunc) colMeans(trunc$est_R1_temp_diff) else NULL,
    r_est_trunc_R0_diff  = if (do_trunc) colMeans(trunc$est_R0_temp_diff) else NULL,
    B = B,
    ...
  )
}

#' Cross-fitted ORM estimator for sensitivity-adjusted outcomes
#'
#' Estimates treatment-specific and overall sensitivity-adjusted outcome means
#' using cross-fitting with ordinal regression outcome models (`rms::orm`) and
#' nuisance models for treatment and missingness (`mgcv::gam`). The function
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
#'   probability weights. If `FALSE`, use fixed probability truncation and
#'   additional IF truncation diagnostics.
#' @param quant Numeric in `(0, 1)` used as the upper quantile for simple weight
#'   truncation when `simple_trunc = TRUE`.
#' @param coef_g.fit Optional starting values for a treatment model; currently
#'   retained for interface compatibility.
#' @param coef_t_R0.fit Optional starting coefficients for treatment model fit
#'   in `R = 0` stratum.
#' @param coef_t_R1.fit Optional starting coefficients for treatment model fit
#'   in `R = 1` stratum.
#' @param coef_M_R0.fit Optional starting coefficients for missingness model fit
#'   in `R = 0` stratum.
#' @param coef_M_R1.fit Optional starting coefficients for missingness model fit
#'   in `R = 1` stratum.
#'
#' @return A named list of estimates and uncertainty summaries for each value in
#'   `gamma`. Core elements include `est`, `est_R1`, `est_R0`, variance
#'   estimates (`var`, `var_R1`, `var_R0`, and `_2` analogs), and confidence
#'   interval bounds (`lowerCI*`, `upperCI*`). Additional components depend on
#'   `simple_trunc` and `IF_output`:
#'   \itemize{
#'   \item `simple_trunc = TRUE`: returns non-truncated summaries only.
#'   \item `simple_trunc = FALSE`: additionally returns truncated summaries,
#'   bootstrap-based diagnostics (`q*`, `var_trunc*_star`, `bias*_star`), and
#'   truncated IF objects when requested.
#'   \item `IF_output = TRUE`: includes influence-function lists (`IF*`) and,
#'   when relevant, truncated IF lists (`IF_trunc*`).
#'   }
#'
#' @details
#' The procedure uses sample-splitting/cross-fitting to reduce overfitting bias
#' in nuisance estimation. Outcome regressions are fit with ordinal logistic
#' models and then integrated over estimated conditional distributions to obtain
#' conditional means and sensitivity-adjusted moments.
#'
#' @examples
#' # out <- est_psi_orm(Y, M, R, X, t, trt = 1, gamma = c(0, 0.5),
#' #                   fold = 5, seed = 1, IF_output = FALSE,
#' #                   simple_trunc = TRUE, quant = 0.99)
est_psi_orm <- function(
  Y, M, R, X, t, trt, gamma, fold, seed, IF_output,
  simple_trunc, quant,
  coef_g.fit=NULL,
  coef_t_R0.fit=NULL,
  coef_t_R1.fit=NULL,
  coef_M_R0.fit=NULL,
  coef_M_R1.fit=NULL
  ){
  
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
  colnames(X_adjust) <- unlist(single.index.variables_colnames(X))
  
  ## empty containers
  containers <- est_s_t_y_orm_create_containers(gamma, fold)
  if (!simple_trunc)
    containers_trunc <- est_s_t_y_orm_create_containers(gamma, fold)
  
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
    X_adjust_out_fold <- X_adjust[out_fold_id_list, ]
    
    Y_out_fold_t_R0 <- Y_out_fold[which(t_out_fold==trt & R_out_fold==0)]
    Y_out_fold_t_R1 <- Y_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    X_adjust_out_fold_t_R0 <- X_adjust_out_fold[which(t_out_fold==trt & R_out_fold==0), ]
    X_adjust_out_fold_t_R1 <- X_adjust_out_fold[which(t_out_fold==trt & R_out_fold==1), ]
    M_out_fold_t_R0 <- M_out_fold[which(t_out_fold==trt & R_out_fold==0)]
    M_out_fold_t_R1 <- M_out_fold[which(t_out_fold==trt & R_out_fold==1)]
    
    mydata_out_fold_t_R0 <- cbind(Y_out_fold_t_R0, X_adjust_out_fold_t_R0)
    mydata_out_fold_t_R1 <- cbind(Y_out_fold_t_R1, X_adjust_out_fold_t_R1)
    
    ## in-fold data
    nk_in_fold <- length(fold_list[[k]])
    Y_in_fold <- Y[fold_list[[k]]]
    M_in_fold <- M[fold_list[[k]]]
    R_in_fold <- R[fold_list[[k]]]
    t_in_fold <- t[fold_list[[k]]]
    trt.ind_in_fold <- trt.ind[fold_list[[k]]]
    X_adjust_in_fold <- X_adjust[fold_list[[k]], ]

    X_adjust_in_fold_t_R0 <- X_adjust_in_fold[which(t_in_fold==trt & R_in_fold==0), ]
    X_adjust_in_fold_t0_R0 <- X_adjust_in_fold[which(t_in_fold!=trt & R_in_fold==0), ]
    X_adjust_in_fold_R1 <- X_adjust_in_fold[which(R_in_fold==1), ]
    
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
    fit_t_R0_h <- orm(as.formula(paste("Y_out_fold_t_R0~", paste(colnames(X_adjust), collapse = "+"))), 
                      data=data.frame(mydata_out_fold_t_R0[which(M_out_fold_t_R0==1), ]), family="logistic")
    
    fit_t_R1_h <- orm(as.formula(paste("Y_out_fold_t_R1~", paste(colnames(X_adjust), collapse = "+"))), 
                      data=data.frame(mydata_out_fold_t_R1[which(M_out_fold_t_R1==1), ]), family="logistic")
  
    ## compute CDF
    y_t_R0 = sort(unique(Y_out_fold_t_R0[which(M_out_fold_t_R0==1)]))    
    ny_t_R0 = length(y_t_R0) 
    F_t_R0_X_t_R0 <- 1-predict(fit_t_R0_h, newdata = X_adjust_in_fold_t_R0, type = "fitted")
    F_t_R0_X_t0_R0 <- 1-predict(fit_t_R0_h, newdata = X_adjust_in_fold_t0_R0, type = "fitted")

    F_t_R0_X_t_R0 <- cbind(0, F_t_R0_X_t_R0)
    F_t_R0_X_t0_R0 <- cbind(0, F_t_R0_X_t0_R0)
    
    y_t_R1 = sort(unique(Y_out_fold_t_R1[which(M_out_fold_t_R1==1)])) 
    ny_t_R1 = length(y_t_R1) 
    F_t_R1_X_R1 <- 1-predict(fit_t_R1_h, newdata = X_adjust_in_fold_R1, type = "fitted")

    F_t_R1_X_R1 <- cbind(0, F_t_R1_X_R1)
    
    #compute PDF
    dF_t_R0_X_t_R0 <- F_t_R0_X_t_R0[, -1, drop = FALSE]-F_t_R0_X_t_R0[, -ny_t_R0, drop = FALSE]
    dF_t_R0_X_t0_R0 <- F_t_R0_X_t0_R0[, -1, drop = FALSE]-F_t_R0_X_t0_R0[, -ny_t_R0, drop = FALSE]
    dF_t_R1_X_R1 <- F_t_R1_X_R1[,-1, drop = FALSE]-F_t_R1_X_R1[,-ny_t_R1, drop = FALSE]
    
    ## conditional expectation of Y given R=0, T=t and X, for I(R=0, T=t)
    mu_Y_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0[-1]+y_t_R0[-ny_t_R0])/2)
    
    ## conditional expectation of Y given R=1, T=t and X, for I(R=1)
    mu_Y_t_R1_X_R1 <- c(dF_t_R1_X_R1 %*% (y_t_R1[-1]+y_t_R1[-ny_t_R1])/2)
    
    ## P(R=1)
    prop.R1 <- mean(R_in_fold)
    
    ## containers$vector_R1 and containers$vector_R0
    containers$vector_R1 <- c(containers$vector_R1, c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                                rep(1/prop.R1, length(M_in_fold_R1))))
    containers$vector_R0 <- c(containers$vector_R0, c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                                rep(0, length(M_in_fold_R1))))
    
    ## Compute influence function for each gamma_t
    for (g in 1:length(gamma)){
      
      ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=t)
      mu_Yexp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (y_t_R0[-1]*exp(gamma[g]*pnorm((y_t_R0[-1]-60)/25))+
                                                     y_t_R0[-ny_t_R0]*exp(gamma[g]*pnorm((y_t_R0[-ny_t_R0]-60)/25)))/2)
      
      ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=t)
      mu_exp_t_R0_X_t_R0 <- c(dF_t_R0_X_t_R0 %*% (exp(gamma[g]*pnorm((y_t_R0[-1]-60)/25))+
                                                    exp(gamma[g]*pnorm((y_t_R0[-ny_t_R0]-60)/25)))/2)
      
      ## conditional expectation of Y * exp given R=0, T=t and X, for I(R=0, T=1-t)
      mu_Yexp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (y_t_R0[-1]*exp(gamma[g]*pnorm((y_t_R0[-1]-60)/25))+
                                                       y_t_R0[-ny_t_R0]*exp(gamma[g]*pnorm((y_t_R0[-ny_t_R0]-60)/25)))/2)
      
      ## conditional expectation of exp given R=0, T=t and X, for I(R=0, T=1-t)
      mu_exp_t_R0_X_t0_R0 <- c(dF_t_R0_X_t0_R0 %*% (exp(gamma[g]*pnorm((y_t_R0[-1]-60)/25))+
                                                      exp(gamma[g]*pnorm((y_t_R0[-ny_t_R0]-60)/25)))/2)
      
      if(!simple_trunc){
        mu_exp_t_R0_X_t_R0 <- pmax(mu_exp_t_R0_X_t_R0,0.01)
        mu_exp_t_R0_X_t0_R0 <- pmax(mu_exp_t_R0_X_t0_R0,0.01)
      }
      
    ## containers$IF+psi within each fold
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
    containers$var_R1_temp[k, g] <- sum((if_R1_temp - vector_mean_R1_temp)^2)/(nk_in_fold-1)
    containers$var_R0_temp[k, g] <- sum((if_R0_temp - vector_mean_R0_temp)^2)/(nk_in_fold-1)
    containers$est_temp_diff[k, g] <- mean(if_temp_diff)
    containers$est_R1_temp_diff[k, g] <- mean(if_R1_temp_diff)
    containers$est_R0_temp_diff[k, g] <- mean(if_R0_temp_diff)
    containers$var_temp_diff[k, g] <- sum((if_temp_diff - mean(if_temp_diff))^2)/(nk_in_fold-1)
    containers$var_R1_temp_diff[k, g] <- sum((if_R1_temp_diff - vector_mean_R1_temp_diff)^2)/(nk_in_fold-1)
    containers$var_R0_temp_diff[k, g] <- sum((if_R0_temp_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
    
    if(!simple_trunc){
      
      if (eq(max(abs(if_R1_temp)),z=if_R1_temp)>0) {
        trunc <- max(abs(if_R1_temp))
      } else {
        trunc <- uniroot(eq, z=if_R1_temp, interval=c(0.01, max(abs(if_R1_temp))))$root
      }
      if_R1_temp_trunc <- pmin(abs(if_R1_temp),trunc) *  sign(if_R1_temp)
      containers_trunc$IF_R1[[g]] <- c(containers_trunc$IF_R1[[g]], if_R1_temp_trunc)
      
      if (eq(max(abs(if_R1_temp_diff)),z=if_R1_temp_diff)>0) {
        trunc <- max(abs(if_R1_temp_diff))
      } else {
        trunc <- uniroot(eq, z=if_R1_temp_diff, interval=c(0.01, max(abs(if_R1_temp_diff))))$root
      }
      if_R1_temp_trunc_diff <- pmin(abs(if_R1_temp_diff),trunc) * sign(if_R1_temp_diff)
      containers_trunc$IF_R1_diff[[g]] <- c(containers_trunc$IF_R1_diff[[g]], if_R1_temp_trunc_diff)
      
      if (eq(max(abs(if_R0_temp)),z=if_R0_temp)>0) {
        trunc <- max(abs(if_R0_temp))
      } else {
        trunc <- uniroot(eq, z=if_R0_temp, interval=c(0.01, max(abs(if_R0_temp))))$root
      }
      if_R0_temp_trunc <- pmin(abs(if_R0_temp),trunc) * sign(if_R0_temp)
      containers_trunc$IF_R0[[g]] <- c(containers_trunc$IF_R0[[g]], if_R0_temp_trunc)
      
      if (eq(max(abs(if_R0_temp_diff)),z=if_R0_temp_diff)>0) {
        trunc <- max(abs(if_R0_temp_diff))
      } else {
        trunc <- uniroot(eq, z=if_R0_temp_diff, interval=c(0.01, max(abs(if_R0_temp_diff))))$root
      }
      if_R0_temp_trunc_diff <- pmin(abs(if_R0_temp_diff),trunc) * sign(if_R0_temp_diff)
      containers_trunc$IF_R0_diff[[g]] <- c(containers_trunc$IF_R0_diff[[g]], if_R0_temp_trunc_diff)
      
      if_temp_trunc <- if_R1_temp_trunc*prop.R1+if_R0_temp_trunc*(1-prop.R1)
      containers_trunc$IF[[g]] <- c(containers_trunc$IF[[g]], if_temp_trunc)
      if_temp_trunc_diff <- if_R1_temp_trunc_diff*prop.R1+if_R0_temp_trunc_diff*(1-prop.R1)
      containers_trunc$IF_diff[[g]] <- c(containers_trunc$IF_diff[[g]], if_temp_trunc_diff)
      
      vector_mean_R0_temp <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                               rep(0, length(M_in_fold_R1)))*mean(if_R0_temp_trunc)
      vector_mean_R0_temp_diff <- c(rep(1/(1-prop.R1), length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)), 
                                    rep(0, length(M_in_fold_R1)))*mean(if_R0_temp_trunc_diff)
      
      vector_mean_R1_temp <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                               rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp_trunc)
      vector_mean_R1_temp_diff <- c(rep(0, length(M_in_fold_t_R0)+length(M_in_fold_t0_R0)),
                                    rep(1/prop.R1, length(M_in_fold_R1)))*mean(if_R1_temp_trunc_diff)
      # 
      containers_trunc$est_temp[k, g] <- mean(if_temp_trunc)
      containers_trunc$est_R1_temp[k, g] <- mean(if_R1_temp_trunc)
      containers_trunc$est_R0_temp[k, g] <- mean(if_R0_temp_trunc)
      containers_trunc$var_temp[k, g] <- sum((if_temp_trunc - mean(if_temp_trunc))^2)/(nk_in_fold-1)
      containers_trunc$var_R1_temp[k, g] <- sum((if_R1_temp_trunc - vector_mean_R1_temp)^2)/(nk_in_fold-1)
      containers_trunc$var_R0_temp[k, g] <- sum((if_R0_temp_trunc - vector_mean_R0_temp)^2)/(nk_in_fold-1)
      
      containers_trunc$est_temp_diff[k, g] <- mean(if_temp_trunc_diff)
      containers_trunc$est_R1_temp_diff[k, g] <- mean(if_R1_temp_trunc_diff)
      containers_trunc$est_R0_temp_diff[k, g] <- mean(if_R0_temp_trunc_diff)
      containers_trunc$var_temp_diff[k, g] <- sum((if_temp_trunc_diff - mean(if_temp_trunc_diff))^2)/(nk_in_fold-1)
      containers_trunc$var_R1_temp_diff[k, g] <- sum((if_R1_temp_trunc_diff - vector_mean_R1_temp_diff)^2)/(nk_in_fold-1)
      containers_trunc$var_R0_temp_diff[k, g] <- sum((if_R0_temp_trunc_diff - vector_mean_R0_temp_diff)^2)/(nk_in_fold-1)
      
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
  
  ## another variance estimation & skewness & multiplier bootstrap
  ve <- var_estimator(
    containers, n = n,
    trunc = if (!simple_trunc) containers_trunc else NULL
  )
  r_var_2         <- ve$var_2
  r_var_R1_2      <- ve$var_R1_2
  r_var_R0_2      <- ve$var_R0_2
  r_var_diff_2    <- ve$var_diff_2
  r_var_R1_diff_2 <- ve$var_R1_diff_2
  r_var_R0_diff_2 <- ve$var_R0_diff_2
  r_bias_star         <- c()
  r_bias_R1_star      <- c()
  r_bias_R0_star      <- c()
  r_bias_diff_star    <- c()
  r_bias_R1_diff_star <- c()
  r_bias_R0_diff_star <- c()
  if (!simple_trunc) {
    r_var_trunc_2         <- ve$var_trunc_2
    r_var_trunc_R1_2      <- ve$var_trunc_R1_2
    r_var_trunc_R0_2      <- ve$var_trunc_R0_2
    r_var_trunc_diff_2    <- ve$var_trunc_diff_2
    r_var_trunc_R1_diff_2 <- ve$var_trunc_R1_diff_2
    r_var_trunc_R0_diff_2 <- ve$var_trunc_R0_diff_2
    r_var_trunc_star         <- ve$var_trunc_star
    r_var_trunc_R1_star      <- ve$var_trunc_R1_star
    r_var_trunc_R0_star      <- ve$var_trunc_R0_star
    r_var_trunc_diff_star    <- ve$var_trunc_diff_star
    r_var_trunc_R1_diff_star <- ve$var_trunc_R1_diff_star
    r_var_trunc_R0_diff_star <- ve$var_trunc_R0_diff_star
    r_bias_star         <- ve$bias_star
    r_bias_R1_star      <- ve$bias_R1_star
    r_bias_R0_star      <- ve$bias_R0_star
    r_bias_diff_star    <- ve$bias_diff_star
    r_bias_R1_diff_star <- ve$bias_R1_diff_star
    r_bias_R0_diff_star <- ve$bias_R0_diff_star
    q         <- ve$q
    q_R1      <- ve$q_R1
    q_R0      <- ve$q_R0
    q_diff    <- ve$q_diff
    q_R1_diff <- ve$q_R1_diff
    q_R0_diff <- ve$q_R0_diff
  }
  
  ## one-step jackknife
  # jack_est_trunc <- c()
  # jack_est_trunc_R0 <- c()
  # jack_est_trunc_R1 <- c()
  # for (g in 1:length(gamma)){
  #   jack_est_trunc <- cbind(jack_est_trunc, r_est_trunc[g]-containers_trunc$IF[[g]]/(n-1))
  #   jack_est_trunc_R0 <- cbind(jack_est_trunc_R0, r_est_trunc_R0[g]-containers_trunc$IF_R0[[g]]/(n-1))
  #   jack_est_trunc_R1 <- cbind(jack_est_trunc_R1, r_est_trunc_R1[g]-containers_trunc$IF_R1[[g]]/(n-1))
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
                     var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2,
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2,
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
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
                     var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2,
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                     var_trunc_2=r_var_trunc_2, var_trunc_R1_2=r_var_trunc_R1_2, var_trunc_R0_2=r_var_trunc_R0_2,
                     # skew_trunc=skew_trunc, skew_trunc_R1=skew_trunc_R1, skew_trunc_R0=skew_trunc_R0,
                     # jack_var_trunc=jack_var_trunc, jack_var_trunc_R0=jack_var_trunc_R0, jack_var_trunc_R1=jack_var_trunc_R1, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                     lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                     upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                     est_diff_mat=containers$est_temp_diff, est_R1_diff_mat=containers$est_R1_temp_diff, est_R0_diff_mat=containers$est_R0_temp_diff, 
                     est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                     est_trunc_diff_mat=containers_trunc$est_temp_diff, est_R1_trunc_diff_mat=containers_trunc$est_R1_temp_diff, est_R0_trunc_diff_mat=containers_trunc$est_R0_temp_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                     var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2,
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                     var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                     var_trunc_diff_2=r_var_trunc_diff_2, var_trunc_R1_diff_2=r_var_trunc_R1_diff_2, var_trunc_R0_diff_2=r_var_trunc_R0_diff_2,
                     # skew_trunc_diff=skew_trunc_diff, skew_trunc_R1_diff=skew_trunc_R1_diff, skew_trunc_R0_diff=skew_trunc_R0_diff, 
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
                     containers$id_list=containers$id_list, 
                     q=q, q_R1=q_R1, q_R0=q_R0, q_diff=q_diff, q_R1_diff=q_R1_diff, q_R0_diff=q_R0_diff, 
                     var_trunc_star=r_var_trunc_star, var_trunc_R1_star=r_var_trunc_R1_star, var_trunc_R0_star=r_var_trunc_R0_star, 
                     var_trunc_diff_star=r_var_trunc_diff_star, var_trunc_R1_diff_star=r_var_trunc_R1_diff_star, var_trunc_R0_diff_star=r_var_trunc_R0_diff_star, 
                     bias_star=r_bias_star, bias_R1_star=r_bias_R1_star, bias_R0_star=r_bias_R0_star, 
                     bias_diff_star=r_bias_diff_star, bias_R1_diff_star=r_bias_R1_diff_star, bias_R0_diff_star=r_bias_R0_diff_star)
    }
  }else{
    if(simple_trunc){
      result <- list(est=r_est, est_R1=r_est_R1, est_R0=r_est_R0, 
                     var=r_var, var_R1=r_var_R1, var_R0=r_var_R0, 
                     var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2,
                     # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                     lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0, 
                     upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0, 
                     est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff, 
                     var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff, 
                     var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2,
                     # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
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
                       var_2=r_var_2, var_R1_2=r_var_R1_2, var_R0_2=r_var_R0_2,
                       # skew=skew, skew_R1=skew_R1, skew_R0=skew_R0, 
                       var_trunc=r_var_trunc, var_trunc_R1=r_var_trunc_R1, var_trunc_R0=r_var_trunc_R0,
                       var_trunc_2=r_var_trunc_2, var_trunc_R1_2=r_var_trunc_R1_2, var_trunc_R0_2=r_var_trunc_R0_2,
                       # skew_trunc=skew_trunc, skew_trunc_R1=skew_trunc_R1, skew_trunc_R0=skew_trunc_R0,
                       # jack_var_trunc=jack_var_trunc, jack_var_trunc_R0=jack_var_trunc_R0, jack_var_trunc_R1=jack_var_trunc_R1, 
                       lowerCI=r_lowerCI, lowerCI_R1=r_lowerCI_R1, lowerCI_R0=r_lowerCI_R0,
                       upperCI=r_upperCI, upperCI_R1=r_upperCI_R1, upperCI_R0=r_upperCI_R0,
                       lowerCI_trunc=r_lowerCI_trunc, lowerCI_trunc_R1=r_lowerCI_trunc_R1, lowerCI_trunc_R0=r_lowerCI_trunc_R0,
                       upperCI_trunc=r_upperCI_trunc, upperCI_trunc_R1=r_upperCI_trunc_R1, upperCI_trunc_R0=r_upperCI_trunc_R0,
                       est_diff=r_est_diff, est_R1_diff=r_est_R1_diff, est_R0_diff=r_est_R0_diff,
                       est_diff_mat=containers$est_temp_diff, est_R1_diff_mat=containers$est_R1_temp_diff, est_R0_diff_mat=containers$est_R0_temp_diff, 
                       est_trunc_diff=r_est_trunc_diff, est_trunc_R1_diff=r_est_trunc_R1_diff, est_trunc_R0_diff=r_est_trunc_R0_diff,
                       est_trunc_diff_mat=containers_trunc$est_temp_diff, est_R1_trunc_diff_mat=containers_trunc$est_R1_temp_diff, est_R0_trunc_diff_mat=containers_trunc$est_R0_temp_diff, 
                       var_diff=r_var_diff, var_R1_diff=r_var_R1_diff, var_R0_diff=r_var_R0_diff,
                       var_diff_2=r_var_diff_2, var_R1_diff_2=r_var_R1_diff_2, var_R0_diff_2=r_var_R0_diff_2,
                       # skew_diff=skew_diff, skew_R1_diff=skew_R1_diff, skew_R0_diff=skew_R0_diff, 
                       var_trunc_diff=r_var_trunc_diff, var_trunc_R1_diff=r_var_trunc_R1_diff, var_trunc_R0_diff=r_var_trunc_R0_diff,
                       var_trunc_diff_2=r_var_trunc_diff_2, var_trunc_R1_diff_2=r_var_trunc_R1_diff_2, var_trunc_R0_diff_2=r_var_trunc_R0_diff_2,
                       # skew_trunc_diff=skew_trunc_diff, skew_trunc_R1_diff=skew_trunc_R1_diff, skew_trunc_R0_diff=skew_trunc_R0_diff, 
                       lowerCI_diff=r_lowerCI_diff, lowerCI_R1_diff=r_lowerCI_R1_diff, lowerCI_R0_diff=r_lowerCI_R0_diff,
                       upperCI_diff=r_upperCI_diff, upperCI_R1_diff=r_upperCI_R1_diff, upperCI_R0_diff=r_upperCI_R0_diff,
                       lowerCI_trunc_diff=r_lowerCI_trunc_diff, lowerCI_trunc_R1_diff=r_lowerCI_trunc_R1_diff, lowerCI_trunc_R0_diff=r_lowerCI_trunc_R0_diff,
                       upperCI_trunc_diff=r_upperCI_trunc_diff, upperCI_trunc_R1_diff=r_upperCI_trunc_R1_diff, upperCI_trunc_R0_diff=r_upperCI_trunc_R0_diff,
                       pain_bq_reordered=pain_bq_reordered, pain_bq_reordered_R1=pain_bq_reordered_R1,
                       pain_bq_reordered_R0=pain_bq_reordered_R0,
                       containers$id_list=containers$id_list, 
                       q=q, q_R1=q_R1, q_R0=q_R0, q_diff=q_diff, q_R1_diff=q_R1_diff, q_R0_diff=q_R0_diff, 
                       var_trunc_star=r_var_trunc_star, var_trunc_R1_star=r_var_trunc_R1_star, var_trunc_R0_star=r_var_trunc_R0_star, 
                       var_trunc_diff_star=r_var_trunc_diff_star, var_trunc_R1_diff_star=r_var_trunc_R1_diff_star, var_trunc_R0_diff_star=r_var_trunc_R0_diff_star, 
                       bias_star=r_bias_star, bias_R1_star=r_bias_R1_star, bias_R0_star=r_bias_R0_star, 
                       bias_diff_star=r_bias_diff_star, bias_R1_diff_star=r_bias_R1_diff_star, bias_R0_diff_star=r_bias_R0_diff_star)
    }
  }
  
}


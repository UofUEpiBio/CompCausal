#' comprehensivecohort: Inferring causal effects in comprehensive cohort studies
#'
#' This package provide R implementations of the estimation methods for 
#' inferring causal effects in comprehensive cohort studies, 
#' in the presence of unmeasured confounding and missing outcomes. The causal 
#' effects include comprehensive cohort causal effects, randomized trial causal effects, 
#' and patient preference causal effect. 
#' @import utils
#' @importFrom stats predict coef as.formula binomial dbeta dnorm integrate model.matrix nlminb optim optimize pnorm qnorm quantile sd uniroot var
#' @importFrom betareg betareg
#' @import splines
#' @import methods
#' @import dplyr
#' @import ggplot2
#' @import MAVE
#' @importFrom assertthat assert_that
#' @import dfoptim
#' @import ManifoldOptim
#' @import mgcv
#' @importFrom purrr map_dbl
#' @importFrom Rcpp Module
#' @import rlang
#' @keywords Comprehensive cohort, Unmeasured confounding, Missing at random, Causal effect
"_PACKAGE"

#' Simulated Data for the Comprehensive Cohort
#' @format 
#' A data frame wtih 563 observations and 8 variables:
#' \describe{
#'   \item{Y}{Continuous outcome ranging 0-100. Missing outcome as NA.}
#'   \item{M}{Binary missingness indicator for Y: M=1 if Y=Y_obs, M=0 if Y=NA.}
#'   \item{R}{Binary randomization consent indicator (1 for RCT, 0 for PPS).}
#'   \item{t}{Binary treatment assignment.}
#'   \item{age}{Continuous baseline covariate.}
#'   \item{womac_bq}{Baseline WOMAC pain score, ranging 0-100.}
#'   \item{expectationb}{Baseline covariate: Expected pain a year from now (much/a little worse, 
#'                       about the same, a little/much better/free of pain)}
#'   \item{ChronicPainb}{Baseline covariate: Chronic pain grade (grade 1-2, grade 3-4)}
#' }
#' 
#' @examples
#' data(ccohort)
"ccohort"

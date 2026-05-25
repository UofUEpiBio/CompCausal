

# comprehensivecohort: Inferring Comprehensive Cohort Causal Effects in the Presence of Unmeasured Confounding and Missing Outcomes

The `comprehensivecohort` package provides functions for estimating the
comprehensive cohort causal effects (CCCE) in comprehensive cohort
studies. We develop a semiparametric sensitivity analysis framework for
assessing the impact of unmeasured confounding in the observational arm,
and could also handle outcomes missing at random. Details about the
study design, assumptions and methodologies can be found in the
vignettes and paper.

## Installment

### Github version

``` r
remotes::install_github("UofUEpiBio/Comprehensive_cohort")
```

## Examples

Here is an example of how to use the main function `est_psi` with
quantile truncation of weights.

``` r
library(comprehensivecohort)

data(ccohort)

## simple truncation
out <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

We can also apply data adaptive truncation for the entire influence
function:

``` r
## data adaptive truncation
out <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

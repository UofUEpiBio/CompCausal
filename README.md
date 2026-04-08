

# Comprehensive Cohort R package

``` r
remotes::install_github("UofUEpiBio/Comprehensive_cohort")
```

## Examples

``` r
library(comprehensivecohort)
```

    Loading required package: splines

    Loading required package: betareg

``` r
data(ccohort)

out <- with(ccohort, {
  est_psi(Y, M, Y0 = NULL, R, X = data.frame(womac_bq), t, trt = 1, gamma = c(0, 0.5),
    fold = 5, seed = 1, IF_output = FALSE,
    simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
    single_index_method="norm1coef", method="optim")
})
```

    Error in `X_out_fold[which(R_out_fold == 0), ]`:
    ! incorrect number of dimensions

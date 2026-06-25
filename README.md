

<img src="man/figures/symbol4.png" align="right" width="220" alt="CompCausal logo" />

# CompCausal: <br> Inferring Comprehensive Cohort Causal Effects in the Presence of Unmeasured Confounding and Missing Outcomes

`CompCausal` (**Comp**rehensive Cohort **Causal** Effects) is an R
package for estimating comprehensive cohort causal effects (CCCEs) in
comprehensive cohort studies. The package implements a semi-parametric
sensitivity analysis framework for assessing the impact of unmeasured
confounding in the observational arm and accommodates outcomes that are
missing at random.

Detailed descriptions of the study design, identification assumptions,
estimation procedures and implementation are available in the package
vignettes and manuscript.

## Installation

Users can install `CompCausal` using the
<a href="https://cran.r-project.org/package=remotes"
target="_blank"><code>remotes</code></a> R package:

``` r
remotes::install_github("UofUEpiBio/CompCausal")
```

## Example

The package includes a simulated dataset based on the TOIB study, a
comprehensive cohort study that investigated whether older adults with
chronic knee pain should be advised to use topical or oral non-steroidal
anti-inflammatory drugs (NSAIDs) for pain management.

The dataset contains 563 observations. The primary outcome, $Y$, is the
Western Ontario and McMaster Universities Osteoarthritis Index (WOMAC)
pain score at 12 months, measured on a scale from 0 to 100. Some outcome
observations may be missing and are recorded as NA in $Y$. The variable
$M$ is the outcome missingness indicator, where $M=1$ indicates $Y$ is
observed and $M=0$ indicates $Y$ is missing.

Additional variables include the treatment indicator $t$, where $t=1$
corresponds to topical NSAIDs and $t=0$ corresponds to oral NSAIDs, and
the randomization consent indicator $R$, where $R=1$ denotes
participation in the randomized controlled trial (RCT) and $R=0$ denotes
participation in the parallel observational study (OBS). The remaining
variables are baseline covariates, including age, baseline WOMAC pain
score, expected pain level one year later, and chronic pain grade.

``` r
library(CompCausal)
## load in data
data(ccohort)
## data structure
str(ccohort)
```

    'data.frame':   563 obs. of  8 variables:
     $ Y           : num  26.7 72.1 55.3 NA 46 ...
     $ M           : num  1 1 1 0 1 1 1 0 1 1 ...
     $ R           : int  1 0 1 1 0 1 1 1 1 1 ...
     $ t           : num  1 1 1 0 0 0 0 1 0 0 ...
     $ age         : num  54 63 70 74 55 78 76 61 58 52 ...
     $ womac_bq    : num  22.2 100 25.8 30.4 27.4 31.6 45.2 73 39.4 39.4 ...
     $ expectationb: Factor w/ 3 levels "Much/A little worse",..: 3 2 1 2 3 1 1 3 1 1 ...
     $ ChronicPainb: Factor w/ 2 levels "1-2","3-4": 1 2 1 2 1 1 1 1 2 1 ...

The primary function, `est_psi()`, estimates $E[Y(t)]$, $E[Y(t)|R=0]$
and $E[Y(t)|R=1]$ for one or more user-specified values of $\gamma_t$.
In the following examples, we demonstrate how to use `est_psi()` and
summarize the resulting estimates.

The example below estimates $E[Y(1)]$, $E[Y(1)|R=0]$ and $E[Y(1)|R=1]$
under $\gamma_1=0$ and $\gamma_1=0.5$, using 5-fold sample splitting and
99th-percentile truncation of the estimated inverse probability weights.

``` r
## simple truncation
out_t1_simpleTrunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

As an alternative to quantile-based weight truncation, users can apply
the data-adaptive influence-function truncation procedure by setting
`simple_trunc = FALSE` and `quant = NULL`.

``` r
## data adaptive truncation
out_t1_ifTrunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

The `print()` method provides a concise summary of the estimation
results in three separate tables corresponding to $E[Y(t)]$,
$E[Y(t)|R=1]$ and $E[Y(t)|R=0]$. Each table includes the point
estimates, estimated variances, and 95% confidence intervals. For
$E[Y(t)]$ and $E[Y(t)|R=0]$, the associated values of $\gamma_t$ are
also reported.

``` r
print(out_t1_ifTrunc, rounding=2)
```

    Estimation of E[Y(1)]
    =========================
     gamma Estimates  Var Lower_95CI Upper_95CI
       0.0     40.12 2.16      37.24      43.00
       0.5     40.44 2.19      37.54      43.34

    Estimation of E[Y(1)|R=1]
    =========================
     Estimates  Var Lower_95CI Upper_95CI
          38.5 4.72      34.24      42.76

    Estimation of E[Y(1)|R=0]
    =========================
     gamma Estimates  Var Lower_95CI Upper_95CI
       0.0     41.70 3.95      37.80      45.60
       0.5     42.33 4.04      38.39      46.26

To estimate $E[Y(0)]$, $E[Y(0)|R=0]$ and $E[Y(0)|R=1]$, we repeat the
analysis with `trt = 0`. The example below considers $\gamma_0=0$ and
$\gamma_0=0.5$ and summarizes the resulting estimates.

``` r
## data adaptive truncation
out_t0_ifTrunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 0, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
print(out_t0_ifTrunc, rounding=2)
```

    Estimation of E[Y(0)]
    =========================
     gamma Estimates  Var Lower_95CI Upper_95CI
       0.0     40.81 2.67      37.61      44.02
       0.5     41.32 2.68      38.11      44.53

    Estimation of E[Y(0)|R=1]
    =========================
     Estimates  Var Lower_95CI Upper_95CI
         38.28 4.41      34.16       42.4

    Estimation of E[Y(0)|R=0]
    =========================
     gamma Estimates  Var Lower_95CI Upper_95CI
       0.0     43.15 6.14      38.29      48.01
       0.5     44.15 6.16      39.28      49.01

Three causal treatment effects can be estimated using the `est_psi()`
and `print_effects()` functions. First, obtain estimation results under
`trt=1` and `trt=0` using `est_psi()`. These results can then be passed
to `print_effects()` to compute treatment effects estimates, standard
errors, and confidence intervals. To ensure valid inference, users
should proceed as follows:

1.  Run `est_psi()` twice, once with `trt=1` and once with `trt=0`.
    Store the resulting objects separately.
2.  Use the same value of `seed` in both calls to ensure that
    sample-splitting assignments are aligned across treatment groups.
3.  Set `IF_output = TRUE` in both calls so that the estimated influence
    functions are returned.
4.  Use identical values for all other function arguments when fitting
    the models under `trt=1` and `trt=0`.

Aligning the sample splits and retaining the influence functions allows
`print_effects()` to properly account for the covariance between the two
estimators when computing standard errors and confidence intervals for
treatment effects.

``` r
out_t1_ifTrunc_IF <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = TRUE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})

out_t0_ifTrunc_IF <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 0, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = TRUE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
print_effects(out_t1_ifTrunc_IF, out_t0_ifTrunc_IF, rounding=2)
```

    Estimation of CCCE, PPCE, RTCE
    =========================
     gamma1 gamma0 Type Estimates  Var lowerCI upperCI
        0.0    0.0 CCCE     -0.69 4.14   -4.68    3.30
        0.5    0.0 CCCE     -0.88 4.16   -4.88    3.12
        0.0    0.5 CCCE     -0.69 4.15   -4.68    3.30
        0.5    0.5 CCCE     -0.88 4.17   -4.88    3.12
        0.0    0.0 PPCE     -1.45 8.40   -7.13    4.23
        0.5    0.0 PPCE     -1.82 8.45   -7.52    3.88
        0.0    0.5 PPCE     -1.45 8.43   -7.14    4.24
        0.5    0.5 PPCE     -1.82 8.49   -7.53    3.89
         NA     NA RTCE      0.22 8.13   -5.37    5.81

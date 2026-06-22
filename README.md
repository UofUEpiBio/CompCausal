

# comprehensivecohort: Inferring Comprehensive Cohort Causal Effects in the Presence of Unmeasured Confounding and Missing Outcomes

The `comprehensivecohort` package provides functions for estimating the
comprehensive cohort causal effects (CCCE) in comprehensive cohort
studies. We develop a semiparametric sensitivity analysis framework for
assessing the impact of unmeasured confounding in the observational arm.
Our methods can also handle outcomes missing at random. Details about
the study design, assumptions, methodology and implementation can be
found in the vignettes and paper.

## Installment

Users can install `comprehensivecohort` using the
<a href="https://cran.r-project.org/package=remotes"
target="_blank"><code>remotes</code></a> R package:

``` r
remotes::install_github("UofUEpiBio/Comprehensive_cohort", ref="cran")
```


    ── R CMD build ─────────────────────────────────────────────────────────────────
    * checking for file ‘/private/var/folders/kc/516dnvf974g17y535_20gss80000gn/T/RtmpeBVO7G/remotes1097750b9573c/UofUEpiBio-Comprehensive_cohort-ffcfab7/DESCRIPTION’ ... OK
    * preparing ‘comprehensivecohort’:
    * checking DESCRIPTION meta-information ... OK
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    Removed empty directory ‘comprehensivecohort/.devcontainer’
    Removed empty directory ‘comprehensivecohort/.github’
    Removed empty directory ‘comprehensivecohort/data-raw’
    Removed empty directory ‘comprehensivecohort/vignettes’
    * building ‘comprehensivecohort_0.0.9000.tar.gz’

Or install from CRAN:

## Example

The package includes a simulated dataset based on the TOIB study, a
comprehensive cohort study aiming to determine whether to advice older
adults with chronic knee pain to apply either topical or oral
non-steroidal anti-inflammatory drugs (NSAIDs) for knee pain management.
The dataset include 563 observations with outcome $Y$ as Western Ontario
and McMaster Universities Osteoarthritis Index (WOMAC) pain score
($0-100$) at 12 months. Some outcome observations might be missing
(denoted as NA in $Y$). Column $M$ is a outcome missingness indicator:
$1$ if $Y$ is observed, $0$ if $Y$ is missing. Other variables in the
dataset include $t$ the treatment indicator ($1$ for topical NSAIDs, $0$
for oral NSAIDs), and $R$ the randomization consent indicator ($1$ for
RCT, $0$ for OBS). Rest of the columns are baseline covariates (age,
baseline WOMAC pain score, expected pain one year later, chronic pain
grade).

``` r
library(comprehensivecohort)
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

The main function `est_psi` works to estimate $E[Y(t)]$, $E[Y(t)|R=0]$
and $E[Y(t)|R=1]$ under different $\gamma_t$ values. We will use the
following examples to demonstrate how to apply the function and present
results. Here is an example of how to estimate $E[Y(1)]$, $E[Y(1)|R=0]$
and $E[Y(1)|R=1]$ under $\gamma_1=0, 0.5$, using 5-fold sample splitting
and 99th quantile truncation of weights.

``` r
## simple truncation
out_t1_simpleTrunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

We can also apply data adaptive truncation to the influence functions by
setting `simple_trunc = FALSE` and `quant = NULL`.

``` r
## data adaptive truncation
out_t1_ifTrunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = c(0, 0.5), fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
```

Users can utilize the `print()` function to output the estimation
results in three separate tables, each containing $\gamma_t$ (except for
$E[Y(t)|R=1]$), estimates, variance, and 95% confidence intervals.

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

We can also estimate $E[Y(0)]$, $E[Y(0)|R=0]$ and $E[Y(0)|R=1]$ under
$\gamma_1=0, 0.5$ by setting `trt = 0`, and present the results.

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

Three types of causal treatment effects can be computed by inputting the
estimating results under `trt=1` and `trt=0` into the `print_effects()`
function. When the goal is to compute treatment effects, users need to
specify the parameters in the following ways:

1.  We need to run `est_psi()` twice (under `trt=1` and `trt=0`) and
    store the results separately.
2.  `seed` should be set to the same number to align observations.
3.  Set `IF_output = TRUE`.
4.  Other parameter specifications should be set to the same values when
    running `est_psi()` under `trt=1` and `trt=0`.

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
         NA     NA RTCE      0.00 8.13   -5.59    5.59

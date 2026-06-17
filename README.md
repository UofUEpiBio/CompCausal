

# comprehensivecohort: Inferring Comprehensive Cohort Causal Effects in the Presence of Unmeasured Confounding and Missing Outcomes

The `comprehensivecohort` package provides functions for estimating the
comprehensive cohort causal effects (CCCE) in comprehensive cohort
studies. We develop a semiparametric sensitivity analysis framework for
assessing the impact of unmeasured confounding in the observational arm,
and could also handle outcomes missing at random. Details about the
study design, assumptions, methodology and implementation can be found
in the vignettes and paper.

## Installment

Users can install `comprehensivecohort` using the
<a href="https://cran.r-project.org/package=remotes"
target="_blank"><code>remotes</code></a> R package:

``` r
remotes::install_github("UofUEpiBio/Comprehensive_cohort")
```

Or install from CRAN:

## Examples

The package includes a simulated dataset based on the TOIB study, a
comprehensive cohort study aiming to determine whether to advice older
adults with chronic knee pain to apply either topical or oral
non-steroidal anti-inflammatory drugs (NSAIDs) for knee pain management.
The dataset include 563 observations with outcome $Y$ as Western Ontario
and McMaster Universities Osteoarthritis Index (WOMAC) pain score
($0-100$) at 12 months. Some outcome observations might be missing
(denoted as NA in $Y$). Column $M$ is a outcome missingness indicator:
$1$ if $Y$ is observed, $0$ if $Y$ is missing. Other

We will demonstrate how to use the main functions in
`comprehensivecohort` by

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

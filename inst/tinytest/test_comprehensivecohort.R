library(comprehensivecohort)
data(ccohort)

## test estimation with est_psi, simple truncation
out_simple <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = 0.5, fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
expect_equal(out_simple$est, 40.60768, tolerance = 1e-4)
expect_equal(out_simple$est_R1, 40.02764, tolerance = 1e-4)
expect_equal(out_simple$est_R0, 40.90628, tolerance = 1e-4)
expect_equal(out_simple$var, 2.364001, tolerance = 1e-4)
expect_equal(out_simple$var_R1, 5.375859, tolerance = 1e-4)
expect_equal(out_simple$var_R0, 3.889298, tolerance = 1e-4)

## test estimation with est_psi, data adaptive truncation
out_if_trunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 0, gamma = 1, fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
expect_equal(out_if_trunc$est_trunc, 39.75124, tolerance = 1e-4)
expect_equal(out_if_trunc$est_trunc_R1, 34.18674, tolerance = 1e-4)
expect_equal(out_if_trunc$est_trunc_R0, 44.83376, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc, 2.6972, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc_R1, 3.54835, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc_R0, 7.046477, tolerance = 1e-4)

## test estimation with est_psi_exchange, simple truncation
out_exchange_simple <- with(ccohort, {
  est_psi_exchange(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
                   t, trt = 1, gamma = 0.5, fold = 5, seed = 1, IF_output = FALSE,
                   simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
                   single_index_method="norm1coef", method="optim")
})
expect_equal(out_exchange_simple$est, 41.90088, tolerance = 1e-4)
expect_equal(out_exchange_simple$est_R0, 43.60195, tolerance = 1e-4)
expect_equal(out_exchange_simple$var, 2.615983, tolerance = 1e-4)
expect_equal(out_exchange_simple$var_R0, 1.64063, tolerance = 1e-4)

## test estimation with est_psi_exchange, data adaptive truncation
out_exchange_if_trunc <- with(ccohort, {
  est_psi_exchange(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
                   t, trt = 0, gamma = -1, fold = 5, seed = 1, IF_output = FALSE,
                   simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
                   single_index_method="norm1coef", method="optim")
})
expect_equal(out_exchange_if_trunc$est_trunc, 33.53846, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$est_trunc_R0, 32.82498, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$var_trunc, 2.933737, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$var_trunc_R0, 1.403845, tolerance = 1e-4)



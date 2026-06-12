library(comprehensivecohort)
data(ccohort)

## test estimation with est_psi, simple truncation
out_simple <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 1, gamma = 0.5, fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
expect_equal(out_simple$est, 40.5326, tolerance = 1e-4)
expect_equal(out_simple$est_R1, 38.50266, tolerance = 1e-4)
expect_equal(out_simple$est_R0, 42.49174, tolerance = 1e-4)
expect_equal(out_simple$var, 2.270881, tolerance = 1e-4)
expect_equal(out_simple$var_R1, 4.721923, tolerance = 1e-4)
expect_equal(out_simple$var_R0, 4.330714, tolerance = 1e-4)

## test estimation with est_psi, data adaptive truncation
out_if_trunc <- with(ccohort, {
  est_psi(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
          t, trt = 0, gamma = 1, fold = 5, seed = 1, IF_output = FALSE,
          simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
          single_index_method="norm1coef", method="optim")
})
expect_equal(out_if_trunc$est_trunc, 41.74039, tolerance = 1e-4)
expect_equal(out_if_trunc$est_trunc_R1, 38.28075, tolerance = 1e-4)
expect_equal(out_if_trunc$est_trunc_R0, 44.96667, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc, 2.639573, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc_R1, 4.412838, tolerance = 1e-4)
expect_equal(out_if_trunc$var_trunc_R0, 5.974668, tolerance = 1e-4)

## test estimation with est_psi_exchange, simple truncation
out_exchange_simple <- with(ccohort, {
  est_psi_exchange(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
                   t, trt = 1, gamma = 0.5, fold = 5, seed = 1, IF_output = FALSE,
                   simple_trunc = TRUE, quant = 0.99, kernel="dnorm", 
                   single_index_method="norm1coef", method="optim")
})
expect_equal(out_exchange_simple$est, 39.25216, tolerance = 1e-4)
expect_equal(out_exchange_simple$est_R0, 39.97512, tolerance = 1e-4)
expect_equal(out_exchange_simple$var, 2.757351, tolerance = 1e-4)
expect_equal(out_exchange_simple$var_R0, 1.953512, tolerance = 1e-4)

## test estimation with est_psi_exchange, data adaptive truncation
out_exchange_if_trunc <- with(ccohort, {
  est_psi_exchange(Y, M, R, X = data.frame(age, womac_bq, expectationb, ChronicPainb), 
                   t, trt = 0, gamma = -1, fold = 5, seed = 1, IF_output = FALSE,
                   simple_trunc = FALSE, quant = NULL, kernel="dnorm", 
                   single_index_method="norm1coef", method="optim")
})
expect_equal(out_exchange_if_trunc$est_trunc, 37.38989, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$est_trunc_R0, 36.49553, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$var_trunc, 3.872997, tolerance = 1e-4)
expect_equal(out_exchange_if_trunc$var_trunc_R0, 2.217736, tolerance = 1e-4)



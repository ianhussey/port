test_that("non-square / non-numeric input errors", {
  expect_error(check_corr_psd(matrix(1, 2, 3)), "square")
  expect_error(check_corr_psd("not a matrix"), "numeric matrix")
})

test_that("p < 2 errors", {
  expect_error(check_corr_psd(matrix(1, 1, 1)), "at least 2")
})

test_that("NA entries error with a clear message", {
  R <- diag(3); R[1, 2] <- R[2, 1] <- NA
  expect_error(check_corr_psd(R), "NA or non-finite")
})

test_that("non-unit diagonal errors", {
  R <- diag(3); R[2, 2] <- 1.2
  expect_error(check_corr_psd(R), "Diagonal")
})

test_that("small asymmetry warns and is symmetrized; large asymmetry errors", {
  R <- matrix(c(1, 0.3, 0.2,
                0.30005, 1, 0.25,
                0.2, 0.25, 1), 3, 3)  # asymmetry 5e-5: in the warn band
  expect_warning(res <- check_corr_psd(R, decimals = 2), "not exactly symmetric")
  expect_equal(res$verdict, "possible")

  Rbig <- matrix(c(1, 0.3, 0.2,
                   0.9, 1, 0.25,
                   0.2, 0.25, 1), 3, 3)
  expect_error(check_corr_psd(Rbig, decimals = 2), "not symmetric")
})

test_that("data frames are accepted", {
  df <- as.data.frame(diag(3))
  expect_equal(check_corr_psd(df, decimals = 2)$verdict, "possible")
})

test_that("invalid decimals / delta / tau error", {
  expect_error(check_corr_psd(diag(3), decimals = -1), "non-negative integer")
  expect_error(check_corr_psd(diag(3), decimals = 1.5), "non-negative integer")
  expect_error(check_corr_psd(diag(3), delta = -0.1), "non-negative finite")
  expect_error(check_corr_psd(diag(3), tau = -1), "non-negative finite")
})

test_that("delta overrides decimals", {
  R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2, delta = 0.5)
  expect_equal(res$delta, 0.5)
})

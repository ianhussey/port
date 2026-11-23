test_that("check_witness() reproduces the bound of the reported witness", {
  triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  result <- check_corr_psd(triad, decimals = 2)
  checked <- check_witness(triad, result$witness, decimals = 2)
  expect_s3_class(checked, "tbl_df")
  expect_true(checked$inconsistent)
  expect_identical(checked$b_upper, result$b_upper)
  expect_identical(checked$margin, result$margin)

  # the general box (another rounding rule) uses the box bound, as
  # check_corr_psd() does
  result_box <- check_corr_psd(triad, decimals = 2, rounding = "truncate")
  checked_box <- check_witness(
    triad,
    result_box$witness,
    decimals = 2,
    rounding = "truncate"
  )
  expect_identical(checked_box$b_upper, result_box$b_upper)
})

test_that("check_witness() accepts a hand-made witness", {
  triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  # v'Rv = 3 - 5.4 = -2.4; rounding penalty 0.005 * (9 - 3) = 0.03
  checked <- check_witness(triad, c(1, -1, -1), decimals = 2)
  expect_true(checked$inconsistent)
  expect_equal(checked$margin, -2.37)
})

test_that("check_witness() never proves a genuine table inconsistent", {
  set.seed(1)
  genuine <- round(cor(matrix(rnorm(400), ncol = 4)), 2)
  checked <- check_witness(genuine, c(1, -1, 1, -1), decimals = 2)
  expect_false(checked$inconsistent)
  expect_gte(checked$b_upper, 0)
})

test_that("check_witness() frees missing cells to [-1, 1]", {
  triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  triad[2, 3] <- triad[3, 2] <- NA
  # with r23 free, r23 = 0.62 makes the table valid, so no witness can work
  expect_false(check_witness(triad, c(1, -1, -1), decimals = 2)$inconsistent)
})

test_that("check_witness() rejects malformed witnesses", {
  triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  expect_error(check_witness(triad, c(1, -1), decimals = 2), "length 3")
  expect_error(check_witness(triad, c(0, 0, 0), decimals = 2), "not all zero")
  expect_error(check_witness(triad, c(1, NA, 1), decimals = 2), "finite")
  expect_error(check_witness(triad, c(1, -1, -1) * 1e160, decimals = 2), "1e100")
  expect_error(check_witness(triad, c(1, -1, -1) * 1e-170, decimals = 2), "1e-100")
  expect_error(check_witness(triad, c(1, -1, -1), decimals = -1), "decimals")
})

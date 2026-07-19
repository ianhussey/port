test_that("an out-of-range entry is inconsistent via the precheck tier", {
  R <- matrix(c(1, 1.2, 0.2, 1.2, 1, 0.1, 0.2, 0.1, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_equal(res$tier, "precheck")
  expect_equal(res$detail$i, 1L)
  expect_equal(res$detail$j, 2L)
  expect_equal(res$detail$value, 1.2)
})

test_that("an entry just outside the box edge triggers the precheck", {
  # |r| - delta > 1  <=>  |r| > 1 + delta = 1.005 at 2dp
  R <- matrix(c(1, 1.01, 0, 1.01, 1, 0, 0, 0, 1), 3, 3)
  expect_equal(check_corr_psd(R, decimals = 2)$tier, "precheck")

  # |r| = 1.00 is in range (1.00 - 0.005 = 0.995 <= 1): NOT a precheck trigger.
  R2 <- matrix(c(1, 1.00, 0, 1.00, 1, 0, 0, 0, 1), 3, 3)
  res2 <- check_corr_psd(R2, decimals = 2)
  expect_false(identical(res2$tier, "precheck"))
})

test_that("precheck fires before eigdecomposition on a large out-of-range matrix", {
  R <- diag(20)
  R[3, 7] <- R[7, 3] <- 1.5
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_equal(res$tier, "precheck")
})

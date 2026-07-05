test_that("classic 3x3 (r12=r13=0.9, r23=-0.9) is inconsistent via witness", {
  R <- matrix(c(1, 0.9, 0.9,
                0.9, 1, -0.9,
                0.9, -0.9, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_equal(res$tier, "witness")
  expect_true(res$b_upper < 0)
  # margin is the plain B; b_upper adds nonnegative FP slack, so b_upper >= margin
  expect_true(res$b_upper >= res$margin)
})

test_that("the witness vector actually certifies inconsistency (independent recompute)", {
  R <- matrix(c(1, 0.9, 0.9,
                0.9, 1, -0.9,
                0.9, -0.9, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2)
  v <- res$witness
  delta <- res$delta
  # Recompute B = v'Rv + delta*(||v||_1^2 - ||v||_2^2) from scratch.
  B <- as.numeric(t(v) %*% R %*% v) + delta * (sum(abs(v))^2 - sum(v^2))
  expect_true(B < 0)
  expect_equal(B, res$margin, tolerance = 1e-10)
})

test_that("inconsistency is robust to the 2dp rounding radius", {
  # The 3x3 determinant is strongly negative; it should stay inconsistent.
  R <- matrix(c(1, 0.9, 0.9,
                0.9, 1, -0.9,
                0.9, -0.9, 1), 3, 3)
  # det = 1 + 2*.9*.9*(-.9) - .9^2*3 = 1 - 1.458 - 2.43 < 0
  expect_lt(det(R), 0)
  expect_equal(check_corr_psd(R, decimals = 2)$verdict, "inconsistent")
  # Even with a generous absolute delta the classic case remains inconsistent.
  expect_equal(check_corr_psd(R, delta = 0.02)$verdict, "inconsistent")
})

test_that("a larger clearly-indefinite matrix is inconsistent", {
  # All-pairs high positive except one strong negative triangle inconsistency.
  R <- matrix(0.8, 5, 5)
  diag(R) <- 1
  R[1, 2] <- R[2, 1] <- -0.8
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_true(res$b_upper < 0)
})

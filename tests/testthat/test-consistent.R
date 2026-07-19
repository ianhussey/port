test_that("identity matrices are consistent with a certified PSD witness matrix", {
  for (p in 2:6) {
    res <- check_corr_psd(diag(p), decimals = 2)
    expect_equal(res$verdict, "consistent", info = paste("p =", p))
    expect_false(is.null(res$certified_matrix))
  }
})

test_that("a mild well-conditioned matrix is consistent", {
  R <- matrix(c(1, 0.3, 0.2, 0.3, 1, 0.25, 0.2, 0.25, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "consistent")
  expect_false(is.null(res$certified_matrix))
  # the certified matrix must be in-box and pass the Rump verification
  expect_true(port:::.in_box(res$certified_matrix, R, res$delta))
  expect_true(port:::.verify_psd(res$certified_matrix))
})

test_that("consistent verdicts never carry a negative b_upper", {
  R <- matrix(c(1, 0.5, 0.4, 0.5, 1, 0.45, 0.4, 0.45, 1), 3, 3)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "consistent")
  expect_true(is.na(res$b_upper) || res$b_upper >= 0)
})

test_that("certificate() accessor returns the expected shape", {
  R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  cert <- certificate(check_corr_psd(R, decimals = 2))
  expect_named(cert, c("witness", "margin", "b_upper", "verdict", "tier"))
  expect_equal(cert$verdict, "inconsistent")
  expect_length(cert$witness, 3)
  # on a consistent case the witness is NULL
  cert2 <- certificate(check_corr_psd(diag(3), decimals = 2))
  expect_null(cert2$witness)
  expect_error(certificate(42), "corr_psd_check")
})

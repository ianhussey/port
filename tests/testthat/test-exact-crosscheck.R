# Cross-check the exact small-p path (principal-minor maximization) against the
# witness verdict. The witness is the source of truth; the exact check is an
# independent oracle. Soundness requires: whenever the witness says
# "impossible", the exact check also says "impossible".

test_that("exact p=2 check matches the range condition", {
  # |r| <= 1 + delta  is possible; beyond is impossible
  expect_equal(psdness:::.exact_possible(matrix(c(1, 1.00, 1.00, 1), 2), 0.005), "possible")
  expect_equal(psdness:::.exact_possible(matrix(c(1, 1.01, 1.01, 1), 2), 0.005), "impossible")
})

test_that("exact p=3 check matches the classic impossible case", {
  R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  expect_equal(psdness:::.exact_possible(R, 0.005), "impossible")
})

test_that("witness never contradicts the exact verdict on random p in {2,3}", {
  set.seed(11)
  n_undecided <- 0L
  for (t in seq_len(3000)) {
    p <- sample(2:3, 1)
    M <- matrix(runif(p * p, -1, 1), p, p)
    M <- (M + t(M)) / 2
    diag(M) <- 1
    Mr <- round(M, 2)
    ex <- psdness:::.exact_possible(Mr, 0.005)
    res <- check_corr_psd(Mr, decimals = 2)
    if (res$verdict == "undecided") {
      n_undecided <- n_undecided + 1L
      next
    }
    # Soundness of impossibility: witness-impossible => exact-impossible.
    if (res$verdict == "impossible") {
      expect_equal(ex, "impossible",
                   info = paste("witness impossible but exact possible, trial", t))
    }
    # And a certified-possible (has an in-box PSD matrix) must be exact-possible.
    if (res$verdict == "possible" && !is.null(res$certified_matrix)) {
      expect_equal(ex, "possible",
                   info = paste("certified possible but exact impossible, trial", t))
    }
  }
  # There should be at most a handful of genuinely ambiguous (undecided) cases.
  expect_lt(n_undecided, 30L)
})

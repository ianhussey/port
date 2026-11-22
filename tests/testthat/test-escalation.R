# The POCS escalation tier fires only in the precision-limited zone [0, tau]
# when the cheap constructions did not settle consistency. It is self-contained
# (no external solver) and every matrix it returns is independently Rump-verified.

# A near-boundary 4x4 whose witness margin lands in [0, tau] and for which the
# cheap constructions fail: alternating projections still find an in-box PSD
# matrix, so it resolves to a *certified* consistent.
ambiguous_R <- matrix(
  c(
    1.00,
    0.96,
    0.23,
    -0.21,
    0.96,
    1.00,
    0.03,
    -0.04,
    0.23,
    0.03,
    1.00,
    0.26,
    -0.21,
    -0.04,
    0.26,
    1.00
  ),
  4,
  4,
  byrow = TRUE
)

test_that("near-boundary matrix reaches the escalation branch", {
  res <- check_corr_psd(ambiguous_R, decimals = 2)
  expect_true(res$b_upper >= 0)
  expect_true(res$b_upper <= 10 * res$delta) # inside tau default
  expect_false(identical(res$verdict, "inconsistent")) # witness did not fire
})

test_that("POCS resolves the ambiguous case to a certified consistent", {
  res <- check_corr_psd(ambiguous_R, decimals = 2)
  expect_equal(res$verdict, "consistent")
  expect_equal(res$tier, "pocs")
  expect_false(is.null(res$certified_matrix))
  # the returned matrix is genuinely in-box and rigorously PSD
  expect_true(port:::.in_box(res$certified_matrix, ambiguous_R, res$delta))
  expect_true(port:::.verify_psd(res$certified_matrix))
})

test_that("a near-boundary case the dense witness search missed is certified", {
  # The dense directions leave b_upper ~ 6e-5 > 0 and no strictly-PD in-box
  # point exists, so this was once "undecided". The greedy sparse search finds a
  # witness on variables 1-4 (exact integer check: max_box v'Xv < 0).
  R <- matrix(
    c(
      1,
      -0.24,
      -0.27,
      0.80,
      -0.13,
      -0.24,
      1,
      -0.16,
      -0.58,
      0.18,
      -0.27,
      -0.16,
      1,
      0.32,
      0.11,
      0.80,
      -0.58,
      0.32,
      1,
      -0.10,
      -0.13,
      0.18,
      0.11,
      -0.10,
      1
    ),
    5,
    5,
    byrow = TRUE
  )
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_equal(res$tier, "witness")
  expect_equal(res$witness[5], 0)
})

test_that("an unresolved near-boundary case is honestly undecided, not guessed", {
  # Neither a witness (b_upper ~ 2e-4) nor an in-box PSD matrix is found.
  R <- matrix(c(
    1, 0.57, -0.51, -0.95, 0.98,
    0.57, 1, 0.34, -0.3, 0.53,
    -0.51, 0.34, 1, 0.72, -0.6,
    -0.95, -0.3, 0.72, 1, -0.97,
    0.98, 0.53, -0.6, -0.97, 1
  ), 5, 5)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "undecided")
  expect_equal(res$tier, "pocs")
  expect_false(res$certified)
})

test_that(".pocs_consistent returns a verified in-box matrix when one exists", {
  # exactly-singular (rank-2) PSD matrix: box contains a PSD matrix
  L <- matrix(c(1, 0, 0.6, 0.8, 0.6, -0.8), 3, 2, byrow = TRUE)
  Rs <- L %*% t(L)
  diag(Rs) <- 1
  hit <- port:::.pocs_consistent(Rs, 0.005)
  expect_false(is.null(hit))
  expect_true(port:::.in_box(hit$X, Rs, 0.005))
  expect_true(port:::.verify_psd(hit$X))
  # and the top-level verdict is never inconsistent
  expect_false(identical(
    check_corr_psd(Rs, decimals = 2)$verdict,
    "inconsistent"
  ))
})

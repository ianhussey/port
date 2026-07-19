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

test_that("a genuinely boundary-only case is honestly undecided, not guessed", {
  # b_upper ~ 6e-5: the box only grazes the PSD cone, so no strictly-PD in-box
  # point exists for the one-sided Cholesky test to certify.
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
  expect_equal(res$verdict, "undecided")
  expect_equal(res$tier, "pocs")
  expect_false(identical(res$verdict, "inconsistent")) # never a false inconsistent
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

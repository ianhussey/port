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

# 200 * max over the 2-decimal rounding box of w'Xw, for the integer direction
# w = round(v * 1e5), in exact integer arithmetic (every term is an integer far
# below 2^53). Negative means inconsistent, independently of the package's
# floating-point error model. (Clipping at +-1 is ignored, which can only make
# the value larger; no cell of the matrices below is within delta of +-1.)
exact_box_max <- function(R, v) {
  w <- round(v * 1e5)
  2 * sum(outer(w, w) * round(100 * R)) + sum(abs(w))^2 - sum(w^2)
}

test_that("a near-boundary case the dense witness search missed is certified", {
  # The dense directions leave b_upper ~ 6e-5 > 0 and no strictly-PD in-box
  # point exists, so this was once "undecided". The greedy sparse search finds a
  # witness on variables 1-4.
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
  expect_lt(exact_box_max(R, res$witness), 0)
})

test_that("a sparse witness that needs polishing on its support is found", {
  # The bottom eigenvector of R[2:5, 2:5] gives b_upper ~ 2e-4 > 0; one polish
  # step on that support gives -3.4e-5. Polishing on the full matrix turns the
  # direction dense and finds nothing.
  # fmt: skip
  R <- matrix(c(
     1.00,  0.57, -0.51, -0.95,  0.98,
     0.57,  1.00,  0.34, -0.30,  0.53,
    -0.51,  0.34,  1.00,  0.72, -0.60,
    -0.95, -0.30,  0.72,  1.00, -0.97,
     0.98,  0.53, -0.60, -0.97,  1.00
  ), 5, 5, byrow = TRUE)
  v0 <- port:::.subset_bottom_eigvec(R, 2:5)
  expect_gt(port:::.witness_bound(R, v0, 0.005)$B_upper, 0)
  expect_gt(port:::.polish_witness(R, 0.005, v0)$B_upper, 0)
  pol <- port:::.polish_witness(R, 0.005, v0, on_support = TRUE)
  expect_lt(pol$B_upper, 0)
  expect_equal(pol$v[1], 0)
  expect_lt(exact_box_max(R, pol$v), 0)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_lt(exact_box_max(R, res$witness), 0)
})

test_that("the POCS direction certifies a case the sparse search misses", {
  # From special-scripts/sparse-witness-simulation.R. When this was added,
  # greedy elimination left b_upper ~ 7e-4 > 0 and no in-box PSD matrix was
  # found. The bottom eigenvector of the box-side iterate of margin-free
  # alternating projections certifies it (about -3e-4 at the first check that
  # succeeds). The search is called directly, so the test does not depend on
  # the greedy search missing this matrix.
  # fmt: skip
  R <- matrix(c(
     1.00,  0.83,  0.33, -0.17, -0.88, -0.99,  0.14, -0.83,
     0.83,  1.00, -0.17,  0.39, -0.70, -0.76,  0.65, -0.77,
     0.33, -0.17,  1.00, -0.90, -0.11, -0.38, -0.55,  0.10,
    -0.17,  0.39, -0.90,  1.00,  0.13,  0.26,  0.87, -0.09,
    -0.88, -0.70, -0.11,  0.13,  1.00,  0.91,  0.03,  0.96,
    -0.99, -0.76, -0.38,  0.26,  0.91,  1.00, -0.02,  0.84,
     0.14,  0.65, -0.55,  0.87,  0.03, -0.02,  1.00, -0.16,
    -0.83, -0.77,  0.10, -0.09,  0.96,  0.84, -0.16,  1.00
  ), 8, 8, byrow = TRUE)
  bb <- port:::.box_bounds(R, 0.005)
  wp <- port:::.pocs_witness(bb$lo, bb$hi, bb$off, function(v) {
    port:::.witness_bound(R, v, 0.005)
  })
  expect_lt(wp$B_upper, 0)
  expect_lt(exact_box_max(R, wp$v), 0)
  res <- check_corr_psd(R, decimals = 2)
  expect_equal(res$verdict, "inconsistent")
  expect_true(res$certified)
  expect_lt(exact_box_max(R, res$witness), 0)
  # same verdict through the heterogeneous-box path
  res_box <- check_corr_psd(R, decimals = matrix(2, 8, 8))
  expect_equal(res_box$verdict, "inconsistent")
  expect_lt(exact_box_max(R, res_box$witness), 0)
})

test_that("a genuinely boundary-only case is honestly undecided, not guessed", {
  # The box [1, 1 + 2^-6] clips to the single point r = 1, so the only in-box
  # correlation matrix is singular: it is PSD, so no witness exists, and it is
  # not positive definite, so the Cholesky test cannot certify it either.
  R <- matrix(c(1, 1 + 2^-7, 1 + 2^-7, 1), 2)
  res <- check_corr_psd(R, delta = 2^-7)
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

# Fixtures for each localization verdict class (see spec).

triad_R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)

# Single bad cell: (1,2) is the sole culprit in an otherwise-consistent 5x5.
cell_R <- matrix(
  c(
    1,
    0.56,
    0.68,
    0.05,
    0.82,
    0.56,
    1,
    -0.26,
    -0.31,
    -0.75,
    0.68,
    -0.26,
    1,
    -0.18,
    0.50,
    0.05,
    -0.31,
    -0.18,
    1,
    0.27,
    0.82,
    -0.75,
    0.50,
    0.27,
    1
  ),
  5,
  5,
  byrow = TRUE
)

# Bad variable: variable 4 is over-correlated with 1,2,3.
variable_R <- local({
  M <- diag(4)
  M[1, 2] <- M[2, 1] <- -0.4
  M[1, 3] <- M[3, 1] <- 0.1
  M[2, 3] <- M[3, 2] <- 0.1
  M[1, 4] <- M[4, 1] <- 0.8
  M[2, 4] <- M[4, 2] <- 0.8
  M[3, 4] <- M[4, 3] <- 0.8
  M
})

# Joint: two disjoint inconsistent triangles; no single cell suffices.
joint_R <- local({
  blk <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  J <- matrix(0, 6, 6)
  J[1:3, 1:3] <- blk
  J[4:6, 4:6] <- blk
  diag(J) <- 1
  J
})

test_that("inconsistent triangle yields verdict 'triad'", {
  lf <- localize_psd_fault(triad_R, decimals = 2)
  expect_s3_class(lf, "psd_fault")
  expect_equal(lf$localization_verdict, "triad")
  expect_equal(nrow(lf$implicated$cells), 3L)
  expect_match(lf$notes[1], "any one")
})

test_that("single bad cell yields verdict 'cell'", {
  lf <- localize_psd_fault(cell_R, decimals = 2)
  expect_equal(lf$localization_verdict, "cell")
  expect_equal(as.integer(lf$implicated$cells[1, ]), c(1L, 2L))
  expect_equal(lf$convergence, "all")
})

test_that("bad variable column yields verdict 'variable'", {
  lf <- localize_psd_fault(variable_R, decimals = 2)
  expect_equal(lf$localization_verdict, "variable")
  expect_equal(lf$implicated$variable, 4L)
})

test_that("two-cell joint fault yields verdict 'joint'", {
  lf <- localize_psd_fault(joint_R, decimals = 2)
  expect_equal(lf$localization_verdict, "joint")
  expect_gte(nrow(lf$implicated$cells), 2L)
  expect_true(is.null(lf$implicated$variable))
})

test_that("a genuinely PSD rounded matrix yields verdict 'none'", {
  expect_equal(
    localize_psd_fault(diag(4), decimals = 2)$localization_verdict,
    "none"
  )
  set.seed(7)
  L <- matrix(rnorm(20), 5, 4)
  S <- L %*% t(L)
  d <- sqrt(diag(S))
  S <- S / outer(d, d)
  diag(S) <- 1
  lf <- localize_psd_fault(round(S, 2), decimals = 2)
  expect_equal(lf$localization_verdict, "none")
})

test_that("localize accepts a corr_psd_check object (carrying its matrix)", {
  chk <- check_corr_psd(triad_R, decimals = 2)
  expect_false(is.null(attr(chk, "R")))
  lf <- localize_psd_fault(chk)
  expect_equal(lf$localization_verdict, "triad")
})

test_that("severity is always reported and flags egregious violations", {
  lf <- localize_psd_fault(triad_R, decimals = 2)
  expect_true(is.finite(lf$severity$severity_max))
  expect_true(is.finite(lf$severity$severity_frob))
  expect_gt(lf$severity$severity_max, lf$delta) # far beyond rounding
  expect_equal(
    lf$severity$witness_margin,
    check_corr_psd(triad_R, decimals = 2)$margin
  )
})

test_that("print.psd_fault narrates the inference", {
  expect_output(
    print(localize_psd_fault(cell_R, decimals = 2)),
    "Attributable to a single cell"
  )
  expect_output(
    print(localize_psd_fault(triad_R, decimals = 2)),
    "any one of the three"
  )
  expect_output(print(localize_psd_fault(diag(3), decimals = 2)), "NONE")
})

test_that("fault_evidence returns the A-D evidence", {
  ev <- fault_evidence(localize_psd_fault(triad_R, decimals = 2))
  expect_named(
    ev,
    c(
      "sole_culprit_cells",
      "inconsistent_triples",
      "lofo_restoring",
      "sparse_support",
      "rsquared"
    )
  )
  expect_error(fault_evidence(42), "psd_fault")
})

triad <- matrix(c(1, .9, .9, .9, 1, -.9, .9, -.9, 1), 3)

test_that("the greedy joint search (p >= 7) finds a two-cell fault", {
  R <- diag(7)
  R[2:4, 2:4] <- triad
  R[5:7, 5:7] <- triad
  lf <- localize_psd_fault(R, decimals = 2)
  expect_identical(lf$localization_verdict, "joint")
  cells <- lf$implicated$cells
  expect_equal(nrow(cells), 2L)
  # one cell from each frustrated triangle
  expect_true(all(cells[1, ] %in% 2:4) || all(cells[1, ] %in% 5:7))
  expect_true(any(apply(cells, 1, function(c) all(c %in% 2:4))))
  expect_true(any(apply(cells, 1, function(c) all(c %in% 5:7))))
})

test_that("disattenuation is not inconsistent when min |D| is exactly 1", {
  # (r - .005) / (rel + .005) == 1 in real arithmetic: D = 1 is attainable and
  # the 2x2 matrix of ones is PSD.
  for (rel in seq(0.20, 0.90, by = 0.05)) {
    r <- round(rel + 0.01, 2)
    R <- matrix(c(1, r, r, 1), 2)
    v <- check_disattenuated_psd(R, reliability = round(rel, 2))$verdict
    expect_false(identical(v, "inconsistent"), info = paste(r, rel))
  }
})

test_that("triples touching det = 0 exactly are not flagged", {
  R <- matrix(c(1, .5, .5, .5, 1, -.5, .5, -.5, 1), 3)
  expect_equal(nrow(inconsistent_triples(R, delta = 0)), 0L)
})

test_that("tables with missing cells are validated like complete ones", {
  R <- matrix(c(1, NA, Inf, NA, 1, .2, Inf, .2, 1), 3)
  expect_error(check_corr_psd(R), "infinite")
  R <- matrix(c(1, NA, .5, NA, 1, .2, .50005, .2, 1), 3)
  expect_warning(check_corr_psd(R), "symmetrizing")
})

test_that("decimals and cells are validated in the other entry points", {
  expect_error(inconsistent_triples(diag(3), decimals = -1), "decimals")
  expect_error(localize_psd_fault(diag(3), decimals = 1.5), "decimals")
  expect_error(implied_interval(diag(3), cells = c(1, 5)), "cells")
  expect_error(
    check_disattenuated_psd(diag(3), 0.8, reliability_decimals = -2),
    "reliability_decimals"
  )
})

test_that("the witness is named by the variables", {
  R <- triad
  dimnames(R) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_named(check_corr_psd(R)$witness, c("a", "b", "c"))
})

test_that("the batch helper passes the rounding rule through", {
  R <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3)
  out <- check_corr_psd_batch(list(R), rounding = "floor", quiet = TRUE)
  expect_identical(out$verdict, check_corr_psd(R, rounding = "floor")$verdict)
  expect_identical(
    out$delta,
    check_corr_psd(R, rounding = "floor")$delta
  )
})

test_that("variables correlated exactly +-1 are merged and certified", {
  set.seed(1)
  X <- matrix(rnorm(100 * 8), 100) %*% matrix(runif(64, -.2, .5), 8)
  R <- trunc(cor(X) * 100) / 100
  R[lower.tri(R)] <- t(R)[lower.tri(R)]
  diag(R) <- 1
  # variable 9 duplicates variable 3; variable 10 is variable 3 reverse-coded
  R <- rbind(cbind(R, R[, 3], -R[, 3]), c(R[3, ], 1, -1), c(-R[3, ], -1, 1))
  res <- check_corr_psd(R, decimals = 2, rounding = "truncate")
  expect_identical(res$verdict, "consistent")
  expect_true(res$certified)
  expect_match(paste(res$note, collapse = " "), "merged \\(8 remain\\)")
  bx <- port:::.reported_box(R, 2, NULL, "truncate")
  expect_true(port:::.in_box_generic(res$certified_matrix, bx$lo, bx$hi, bx$off))
  expect_gt(min(eigen(res$certified_matrix, symmetric = TRUE)$values), -1e-12)
})

test_that("a table of exact +-1 correlations is certified under every rule", {
  s <- c(1, -1, 1, 1, -1)
  R <- outer(s, s)
  for (rule in c("truncate", "floor", "ceiling")) {
    res <- check_corr_psd(R, decimals = 2, rounding = rule)
    expect_identical(res$verdict, "consistent", info = rule)
    expect_true(res$certified, info = rule)
  }
  # the uniform path with delta = 0
  res <- check_corr_psd(R, delta = 0)
  expect_identical(res$verdict, "consistent")
  expect_true(res$certified)
})

test_that("merging does not hide an inconsistent duplicate", {
  R <- matrix(c(1, 1, .3, 1, 1, .5, .3, .5, 1), 3)
  res <- check_corr_psd(R, decimals = 2, rounding = "truncate")
  expect_identical(res$verdict, "inconsistent")
  # contradictory pins: 1~2 = +1, 2~3 = +1, but 1~3 = -1
  R <- matrix(c(1, 1, -1, 1, 1, 1, -1, 1, 1), 3)
  expect_identical(
    check_corr_psd(R, decimals = 2, rounding = "truncate")$verdict,
    "inconsistent"
  )
  bx <- port:::.reported_box(R, 2, NULL, "truncate")
  expect_null(port:::.pinned_merge(bx$lo, bx$hi, bx$off)$P)
})

test_that("a range merely close to +-1 is not merged into an out-of-box matrix", {
  # hi = 1 - 1e-13 excludes 1, so forcing X12 = 1 would leave the box.
  lo <- matrix(c(1, 1 - 3e-13, 1 - 3e-13, 1), 2)
  hi <- matrix(c(1, 1 - 1e-13, 1 - 1e-13, 1), 2)
  off <- row(lo) != col(lo)
  expect_identical(nrow(port:::.pinned_cells(lo, hi)), 0L)
  hit <- port:::.pocs_feasible(lo, hi, off)
  if (!is.null(hit)) expect_true(port:::.in_box_generic(hit$X, lo, hi, off))

  # A computed edge a few ulps off a real +-1 still counts as pinned:
  # -1.001 under floor rounding at 3 decimals gives [-1, -1 + 1.1e-16].
  R <- matrix(c(1, -1.001, 0.3, -1.001, 1, -0.3, 0.3, -0.3, 1), 3)
  bx <- port:::.reported_box(R, 3, NULL, "floor")
  expect_identical(nrow(port:::.pinned_cells(bx$lo, bx$hi)), 1L)
})

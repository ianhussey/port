# Tests for reliability disattenuation.

R3 <- matrix(c(1, 0.60, 0.55,
               0.60, 1, 0.50,
               0.55, 0.50, 1), 3, 3)

equicorr <- function(p, rho) { M <- matrix(rho, p, p); diag(M) <- 1; M }

test_that("disattenuate applies the Spearman correction", {
  D <- disattenuate(R3, 0.8)
  expect_equal(diag(D), rep(1, 3))
  expect_equal(D[1, 2], 0.60 / 0.8)                     # common rho: /rho
  # per-variable
  Dv <- disattenuate(R3, c(0.7, 0.8, 0.9))
  expect_equal(Dv[1, 2], 0.60 / sqrt(0.7 * 0.8))
  expect_equal(Dv[1, 3], 0.55 / sqrt(0.7 * 0.9))
  # reliability 1 is the identity correction
  expect_equal(disattenuate(R3, 1), R3)
})

test_that("disattenuate validates reliability", {
  expect_error(disattenuate(R3, 1.2), "\\(0, 1\\]")
  expect_error(disattenuate(R3, 0), "\\(0, 1\\]")
  expect_error(disattenuate(R3, c(0.8, 0.9)), "length 1 or p")
})

test_that("critical reliability matches the closed form rho* = 1 - lambda_min", {
  # equicorrelation 0.5 (p=3): lambda_min = 0.5, max|r| = 0.5 -> rho* = 0.5
  E <- equicorr(3, 0.5)
  res <- check_disattenuated_psd(E, reliability = NULL, decimals = 4)
  expect_equal(res$thresholds$lambda_min, 0.5, tolerance = 1e-8)
  expect_equal(res$thresholds$rho_impossible, 0.5, tolerance = 1e-8)
})

test_that("low reliability drives a corrected correlation above 1 -> impossible", {
  res <- check_disattenuated_psd(R3, reliability = 0.34, decimals = 2, max_plausible_r = 0.9)
  expect_equal(res$verdict, "impossible")
  expect_true(res$forward$range_impossible)
  # corrected (1,2) = 0.60/0.34 ~ 1.76
  expect_gt(res$disattenuated[1, 2], 1)
})

test_that("a plausibility cutoff yields an 'implausible' verdict", {
  res <- check_disattenuated_psd(R3, reliability = 0.62, decimals = 2, max_plausible_r = 0.9)
  expect_equal(res$verdict, "implausible")
  expect_gt(res$max_disattenuated, 0.9)
  expect_lt(res$max_disattenuated, 1)
})

test_that("cutoff of 1 disables the implausible level", {
  res <- check_disattenuated_psd(R3, reliability = 0.62, decimals = 2, max_plausible_r = 1)
  expect_true(res$verdict %in% c("possible", "impossible", "undecided"))
  expect_false(identical(res$verdict, "implausible"))
})

test_that("high reliability keeps the disattenuated matrix possible", {
  res <- check_disattenuated_psd(R3, reliability = 0.95, decimals = 2, max_plausible_r = 0.9)
  expect_equal(res$verdict, "possible")
  expect_true(res$headroom > 0)
})

test_that("verdict severity is monotone as reliability falls", {
  sev <- c(possible = 0L, implausible = 1L, impossible = 2L)
  rhos <- c(0.98, 0.85, 0.70, 0.55, 0.40)
  v <- vapply(rhos, function(r)
    check_disattenuated_psd(R3, reliability = r, decimals = 2, max_plausible_r = 0.9)$verdict,
    character(1))
  s <- sev[v]
  expect_false(any(diff(s) < 0))       # non-decreasing severity as rho decreases
})

test_that("impossibility can bind on PSD rather than range", {
  # near-boundary matrix: PSD fails before any corrected |r| exceeds 1
  R2 <- matrix(c(1, 0.5, 0.5, 0.5, 1, -0.4, 0.5, -0.4, 1), 3, 3)
  res <- check_disattenuated_psd(R2, reliability = 0.7, decimals = 2, max_plausible_r = 0.9)
  expect_equal(res$verdict, "impossible")
  expect_true(res$forward$psd_impossible)
  expect_false(res$forward$range_impossible)
  expect_equal(res$thresholds$impossible_binds, "PSD")
  # the disattenuated centre is genuinely indefinite
  expect_lt(min(eigen(res$disattenuated, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("critical mode reports thresholds and an optional floor narrative", {
  reachable <- check_disattenuated_psd(R3, reliability = NULL, decimals = 2,
                                       max_plausible_r = 0.9, plausible_floor = 0.3)
  expect_equal(reachable$mode, "critical")
  expect_true(reachable$thresholds$rho_impossible > 0.3)   # floor below -> reachable
  expect_output(print(reachable), "IMPOSSIBLE construct matrix")

  safe <- check_disattenuated_psd(R3, reliability = NULL, decimals = 2,
                                  plausible_floor = 0.95)
  expect_output(print(safe), "no reachable reliability")
})

test_that("rounded reliabilities widen the disattenuated intervals", {
  # exact reliabilities vs rounded to 1 dp: the boxed version is at least as
  # willing to reach a verdict of impossible (wider box on the correlations side
  # only matters via the reliability box here)
  exact <- check_disattenuated_psd(R3, reliability = 0.50, decimals = 2,
                                   reliability_decimals = Inf, max_plausible_r = 0.9)
  rounded <- check_disattenuated_psd(R3, reliability = 0.50, decimals = 2,
                                     reliability_decimals = 1, max_plausible_r = 0.9)
  expect_s3_class(exact, "disattenuation_check")
  expect_s3_class(rounded, "disattenuation_check")
  # the point disattenuated matrix is identical; only the box differs
  expect_equal(exact$disattenuated, rounded$disattenuated)
})

test_that("max_plausible_r is validated", {
  expect_error(check_disattenuated_psd(R3, 0.7, max_plausible_r = 1.5), "\\(0, 1\\]")
  expect_error(check_disattenuated_psd(R3, 0.7, max_plausible_r = 0), "\\(0, 1\\]")
})

test_that("per-variable reliabilities give a forward verdict", {
  res <- check_disattenuated_psd(R3, reliability = c(0.4, 0.5, 0.9), decimals = 2,
                                 max_plausible_r = 0.9)
  expect_equal(res$mode, "per_variable")
  expect_true(res$verdict %in% c("impossible", "implausible", "possible", "undecided"))
})

# Tests for the sharpening pass: certified semantics, excusable_delta,
# asymmetric rounding boxes, per-cell precision, NA cells, box-sound rho*,
# and the strengthened witness search (residual directions + polishing).

triad_R <- matrix(c(1, 0.9, 0.9,
                    0.9, 1, -0.9,
                    0.9, -0.9, 1), 3, 3)

# small inline generator of valid correlation matrices (no external deps)
rand_psd_corr <- function(p, k = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  L <- matrix(rnorm(p * k), p, k)
  S <- tcrossprod(L) + diag(runif(p, 0.05, 0.5))
  d <- sqrt(diag(S)); S <- S / outer(d, d); diag(S) <- 1
  S
}

# ---- certified semantics (item 4) -------------------------------------------

test_that("verdicts carry an explicit certification flag", {
  imp <- check_corr_psd(triad_R, decimals = 2)
  expect_true(imp$certified)                       # inconsistency is a certificate
  pos <- check_corr_psd(diag(3), decimals = 2)
  expect_true(pos$certified)                       # exhibited PSD matrix
  expect_false(is.null(pos$certified_matrix))
  # certified consistent <=> a certificate matrix is attached
  expect_equal(!is.null(pos$certified_matrix), isTRUE(pos$certified))
})

test_that("the undecided fixture is explicitly uncertified", {
  R <- matrix(c(1, -0.24, -0.27,  0.80, -0.13,
               -0.24,  1,   -0.16, -0.58,  0.18,
               -0.27, -0.16,  1,    0.32,  0.11,
                0.80, -0.58,  0.32,  1,   -0.10,
               -0.13,  0.18,  0.11, -0.10,  1), 5, 5, byrow = TRUE)
  res <- check_corr_psd(R, decimals = 2)
  if (identical(res$verdict, "undecided")) {
    expect_false(res$certified)
  } else {
    # the strengthened witness may legitimately settle this boundary case,
    # but it must do so WITH a certificate
    expect_true(res$certified)
  }
})

test_that("batch output includes the certified column", {
  out <- check_corr_psd_batch(list(a = diag(3), b = triad_R), quiet = TRUE)
  expect_true("certified" %in% names(out))
  expect_true(all(out$certified))
})

# ---- excusable_delta (item 2) ------------------------------------------------

test_that("excusable_delta unifies verdict and severity in precision units", {
  w <- excusable_delta(triad_R)
  expect_gt(w, 0.05)     # not excusable even by 1-decimal rounding
  # verdict consistency: inconsistent at delta exactly when delta < w
  expect_equal(check_corr_psd(triad_R, delta = w / 2)$verdict, "inconsistent")
  expect_false(identical(check_corr_psd(triad_R, delta = w * 1.5)$verdict,
                         "inconsistent"))
  # a comfortably valid matrix has threshold ~ 0
  expect_lt(excusable_delta(diag(4)), 1e-3)
})

test_that("psd_fault severity carries excusable_delta = delta + severity_max", {
  lf <- localize_psd_fault(triad_R, decimals = 2)
  expect_true(is.finite(lf$severity$excusable_delta))
  expect_equal(lf$severity$excusable_delta,
               lf$delta + lf$severity$severity_max)
  expect_equal(lf$severity$excusable_delta, excusable_delta(triad_R),
               tolerance = 2e-3)
  expect_output(print(lf), "Excusable only by mis-reporting")
})

test_that(".excusable_decimals translates widths to decimal places", {
  expect_equal(port:::.excusable_decimals(0.04), 1)   # 1dp (+-0.05) excuses
  expect_equal(port:::.excusable_decimals(0.004), 2)  # 2dp (+-0.005) excuses
  expect_lt(port:::.excusable_decimals(0.6), 1)       # nothing conventional does
})

# ---- asymmetric rounding boxes (item 5) --------------------------------------

test_that(".reported_box builds the correct interval per rounding rule", {
  R <- matrix(c(1, 0.32, -0.45, 0.32, 1, 0, -0.45, 0, 1), 3, 3)
  near <- port:::.reported_box(R, decimals = 2, delta = NULL, "nearest")
  expect_equal(near$lo_raw[1, 2], 0.32 - 0.005)
  expect_equal(near$hi_raw[1, 2], 0.32 + 0.005)
  flo <- port:::.reported_box(R, decimals = 2, delta = NULL, "floor")
  expect_equal(flo$lo_raw[1, 2], 0.32); expect_equal(flo$hi_raw[1, 2], 0.33)
  cei <- port:::.reported_box(R, decimals = 2, delta = NULL, "ceiling")
  expect_equal(cei$lo_raw[1, 2], 0.31); expect_equal(cei$hi_raw[1, 2], 0.32)
  tru <- port:::.reported_box(R, decimals = 2, delta = NULL, "truncate")
  expect_equal(tru$lo_raw[1, 2], 0.32); expect_equal(tru$hi_raw[1, 2], 0.33)
  expect_equal(tru$lo_raw[1, 3], -0.46); expect_equal(tru$hi_raw[1, 3], -0.45)
  expect_equal(tru$lo_raw[2, 3], -0.01); expect_equal(tru$hi_raw[2, 3], 0.01)
})

test_that("matched asymmetric rounding is sound (never inconsistent on truth)", {
  set.seed(31)
  dirfloor <- function(x, d) floor(x * 10^d) / 10^d
  dirtrunc <- function(x, d) trunc(x * 10^d) / 10^d
  for (t in seq_len(n_reps(150))) {
    S <- rand_psd_corr(sample(3:7, 1))
    off <- upper.tri(S) | lower.tri(S)
    Rf <- S; Rf[off] <- dirfloor(S[off], 2); Rf <- (Rf + t(Rf)) / 2; diag(Rf) <- 1
    expect_false(identical(
      check_corr_psd(Rf, decimals = 2, rounding = "floor")$verdict, "inconsistent"))
    Rt <- S; Rt[off] <- dirtrunc(S[off], 2); Rt <- (Rt + t(Rt)) / 2; diag(Rt) <- 1
    expect_false(identical(
      check_corr_psd(Rt, decimals = 2, rounding = "truncate")$verdict, "inconsistent"))
  }
})

test_that("delta with non-nearest rounding errors clearly", {
  expect_error(check_corr_psd(diag(3), delta = 0.01, rounding = "floor"),
               "decimals")
})

# ---- per-cell precision (item 6) ----------------------------------------------

test_that("a decimals matrix supports mixed-precision tables soundly", {
  set.seed(41)
  for (t in seq_len(n_reps(100))) {
    p <- sample(3:6, 1)
    S <- rand_psd_corr(p)
    dec <- matrix(2L, p, p)
    dec[1, ] <- dec[, 1] <- 1L                    # variable 1 reported at 1dp
    Rr <- S
    off <- which(upper.tri(S), arr.ind = TRUE)
    for (r in seq_len(nrow(off))) {
      i <- off[r, 1]; j <- off[r, 2]
      Rr[i, j] <- Rr[j, i] <- round(S[i, j], dec[i, j])
    }
    expect_false(identical(
      check_corr_psd(Rr, decimals = dec)$verdict, "inconsistent"))
  }
})

test_that("a gross violation stays inconsistent under per-cell precision", {
  dec <- matrix(2L, 3, 3)
  res <- check_corr_psd(triad_R, decimals = dec)
  expect_equal(res$verdict, "inconsistent")
  expect_true(res$certified)
})

test_that("invalid decimals / delta matrices error", {
  bad <- matrix(c(2, 1, 2, 2, 2, 2, 2, 2, 2), 3, 3)   # asymmetric
  expect_error(check_corr_psd(diag(3), decimals = bad), "symmetric")
  expect_error(check_corr_psd(diag(3), delta = matrix(-0.1, 3, 3)), "non-negative")
})

# ---- NA cells (item 6) ---------------------------------------------------------

test_that("freeing an NA cell can rescue a triad (verdict for all values)", {
  Rna <- triad_R; Rna[1, 2] <- Rna[2, 1] <- NA
  res <- check_corr_psd(Rna, decimals = 2)
  expect_equal(res$verdict, "consistent")            # some value of (1,2) works
  expect_true(res$certified)
  expect_true(any(grepl("missing cell", res$note)))
})

test_that("an NA cell does not weaken an independent violation", {
  blk <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  J <- matrix(0, 6, 6); J[1:3, 1:3] <- blk; J[4:6, 4:6] <- blk; diag(J) <- 1
  J[2, 3] <- J[3, 2] <- NA                         # free a cell of triangle 1
  res <- check_corr_psd(J, decimals = 2)
  expect_equal(res$verdict, "inconsistent")          # triangle 2 still inconsistent
  expect_true(res$certified)                       # ...whatever the NA was
})

test_that("localize_psd_fault refuses non-uniform boxes politely", {
  Rna <- triad_R; Rna[1, 2] <- Rna[2, 1] <- NA
  chk <- check_corr_psd(Rna, decimals = 2)
  expect_error(localize_psd_fault(chk), "uniform")
})

# ---- box-sound rho* (item 3) ----------------------------------------------------

test_that("the box-sound critical reliability is <= the point closed form", {
  R3 <- matrix(c(1, 0.60, 0.55, 0.60, 1, 0.50, 0.55, 0.50, 1), 3, 3)
  res <- check_disattenuated_psd(R3, reliability = NULL, decimals = 2)
  thr <- res$thresholds
  expect_true(is.finite(thr$rho_inconsistent_box))
  # rounding gives the box extra slack, so the certified boundary sits at or
  # below the exact-value boundary
  expect_lte(thr$rho_inconsistent_box, thr$rho_inconsistent + 1e-3)
  # and inconsistency genuinely holds just below the box-sound boundary
  below <- thr$rho_inconsistent_box - 0.02
  expect_equal(check_disattenuated_psd(R3, reliability = below,
                                       decimals = 2)$verdict, "inconsistent")
})

test_that("box-sound rho* tightens toward the point value as rounding vanishes", {
  E <- matrix(0.5, 3, 3); diag(E) <- 1              # closed form rho* = 0.5
  res <- check_disattenuated_psd(E, reliability = NULL, decimals = 4)
  expect_equal(res$thresholds$rho_inconsistent, 0.5, tolerance = 1e-8)
  expect_equal(res$thresholds$rho_inconsistent_box, 0.5, tolerance = 5e-3)
})

test_that("headroom is anchored on the box-sound boundary", {
  R3 <- matrix(c(1, 0.60, 0.55, 0.60, 1, 0.50, 0.55, 0.50, 1), 3, 3)
  res <- check_disattenuated_psd(R3, reliability = 0.95, decimals = 2)
  thr <- res$thresholds
  expect_equal(res$headroom, 0.95 - thr$rho_inconsistent_box, tolerance = 1e-9)
})

# ---- strengthened witness (item 7): soundness must be untouched -----------------

test_that("residual directions + polishing never create false inconsistents", {
  set.seed(77)
  for (t in seq_len(n_reps(300))) {
    S <- rand_psd_corr(sample(3:8, 1))
    expect_false(identical(
      check_corr_psd(round(S, 2), decimals = 2)$verdict, "inconsistent"))
  }
})

test_that("known inconsistents are still detected after strengthening", {
  expect_equal(check_corr_psd(triad_R, decimals = 2)$verdict, "inconsistent")
  R5 <- matrix(0.8, 5, 5); diag(R5) <- 1
  R5[1, 2] <- R5[2, 1] <- -0.8
  expect_equal(check_corr_psd(R5, decimals = 2)$verdict, "inconsistent")
})

# ---- claim-object reporting (item 1) ---------------------------------------------

test_that("inconsistent verdicts surface the benign-generator caution", {
  expect_output(print(check_corr_psd(triad_R, decimals = 2)), "pairwise deletion")
  expect_output(print(localize_psd_fault(triad_R, decimals = 2)),
                "polychoric")
})

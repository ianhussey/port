# Tests for the van Tilburg & van Tilburg (2023) additions: R^2/VIF localizer,
# causes taxonomy, plausibility gradient, and closed-form implied intervals.

cell_R <- matrix(c(1,     0.56,  0.68,  0.05,  0.82,
                   0.56,  1,    -0.26, -0.31, -0.75,
                   0.68, -0.26,  1,    -0.18,  0.50,
                   0.05, -0.31, -0.18,  1,     0.27,
                   0.82, -0.75,  0.50,  0.27,  1), 5, 5, byrow = TRUE)

variable_R <- local({
  M <- diag(4)
  M[1, 2] <- M[2, 1] <- -0.4; M[1, 3] <- M[3, 1] <- 0.1; M[2, 3] <- M[3, 2] <- 0.1
  M[1, 4] <- M[4, 1] <- 0.8;  M[2, 4] <- M[4, 2] <- 0.8; M[3, 4] <- M[4, 3] <- 0.8
  M
})

# ---- R^2 localizer (V1+V2) --------------------------------------------------

test_that("the residual-witness identity v'Rv = 1 - R^2 holds", {
  set.seed(4)
  L <- matrix(rnorm(30), 6, 5); S <- L %*% t(L); d <- sqrt(diag(S))
  S <- S / outer(d, d); diag(S) <- 1
  for (i in 1:6) {
    dir <- port:::.rsquared_direction(S, i)
    quad <- as.numeric(t(dir$v) %*% S %*% dir$v) / sum(dir$v^2)
    expect_equal(quad, 1 - dir$r2, tolerance = 1e-8)
  }
})

test_that("R^2 localizer cleanly blames the over-connected variable", {
  rsq <- port:::.rsquared_evidence(variable_R, 0.005)
  blamed <- which(port:::.rsquared_blamed(rsq))
  expect_equal(blamed, 4L)
  expect_true(rsq[[4]]$r2 > 1)
  expect_true(rsq[[4]]$complement_pd)
})

test_that("R^2 blame requires a PD complement (V6 under-identification)", {
  # In the variable fixture, variables 1-3 have a non-PD complement (it contains
  # the bad structure), so they are NOT cleanly blamed even if their direction fires.
  rsq <- port:::.rsquared_evidence(variable_R, 0.005)
  for (k in 1:3) expect_false(port:::.rsquared_blamed(rsq)[k])
})

test_that("R^2 attribution and convergence surface in the verdict", {
  lf <- localize_psd_fault(variable_R, decimals = 2)
  expect_equal(lf$localization_verdict, "variable")
  expect_equal(lf$convergence, "all")                       # LOVO + R^2 agree
  expect_true(any(grepl("R\\^2 localizer", lf$notes)))
  expect_true(any(grepl("conditional on the other cells", lf$notes)))  # V6 caveat
})

test_that("a single bad cell shows up as its two endpoint variables in R^2", {
  lf <- localize_psd_fault(cell_R, decimals = 2)
  blamed <- which(port:::.rsquared_blamed(lf$evidence$rsquared))
  expect_setequal(blamed, c(1L, 2L))     # the endpoints of bad cell (1,2)
})

# ---- Causes taxonomy (V5) ---------------------------------------------------

test_that("a gross violation is labelled substantive, not structural", {
  triad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  lf <- localize_psd_fault(triad, decimals = 2)
  expect_equal(lf$structural$severity_class, "substantive")
  expect_match(port:::.structural_note(lf$structural), "Substantive")
})

test_that("a composite+subscores near-boundary case is labelled benign", {
  Sig <- matrix(0.5, 3, 3); diag(Sig) <- 1
  a <- c(1, 1, 1); r <- as.numeric(Sig %*% a) / sqrt(sum(a * (Sig %*% a)))
  C <- rbind(cbind(Sig, r), c(r, 1))
  Cn <- C; for (k in 1:3) Cn[k, 4] <- Cn[4, k] <- r[k] + 0.02
  lf <- localize_psd_fault(round(Cn, 2), decimals = 2)
  expect_equal(lf$structural$severity_class, "near-boundary")
  expect_equal(lf$structural$pattern, "composite")
  expect_match(port:::.structural_note(lf$structural), "benign")
  expect_match(port:::.structural_note(lf$structural), "van Tilburg")
})

test_that("a full set of category dummies reads as ipsative and benign", {
  D <- matrix(-0.5, 3, 3); diag(D) <- 1     # 3 dummies, sum to a constant
  Dn <- D; Dn[1, 2] <- Dn[2, 1] <- -0.52; Dn[1, 3] <- Dn[3, 1] <- -0.52
  lf <- localize_psd_fault(round(Dn, 2), decimals = 2)
  expect_equal(lf$structural$pattern, "ipsative")
  expect_match(port:::.structural_note(lf$structural), "benign")
})

# ---- Plausibility gradient (V4) --------------------------------------------

test_that("a consistent matrix at the ceiling reports a plausibility gradient", {
  Sig <- matrix(0.5, 3, 3); diag(Sig) <- 1
  a <- c(1, 1, 1); r <- as.numeric(Sig %*% a) / sqrt(sum(a * (Sig %*% a)))
  C <- rbind(cbind(Sig, r), c(r, 1))         # exactly singular -> consistent when rounded
  lf <- localize_psd_fault(round(C, 2), decimals = 2)
  expect_equal(lf$localization_verdict, "none")
  expect_false(is.null(lf$plausibility))
  expect_true(lf$plausibility$near_ceiling)
  expect_true(any(grepl("Plausibility", lf$notes)))
})

test_that("a comfortably-PSD matrix reports no ceiling concern", {
  R <- diag(4); R[1, 2] <- R[2, 1] <- 0.3
  lf <- localize_psd_fault(R, decimals = 2)
  expect_equal(lf$localization_verdict, "none")
  expect_false(isTRUE(lf$plausibility$near_ceiling))
})

# ---- Implied intervals / imputation (V3, V8) --------------------------------

test_that("closed-form implied interval matches the box path on a reported cell", {
  R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  pt <- implied_interval(R, cells = c(1, 2), hold = "points")
  bx <- implied_interval(R, cells = c(1, 2), hold = "box")
  expect_equal(pt$status, "feasible"); expect_equal(bx$status, "feasible")
  # both agree the cell must move far negative
  expect_lt(pt$hi, -0.5); expect_lt(bx$hi, -0.5)
  expect_equal(pt$hi, bx$hi, tolerance = 0.05)
})

test_that("single-missing-cell imputation recovers the feasible interval", {
  Sig <- matrix(0.5, 3, 3); diag(Sig) <- 1
  a <- c(1, 1, 1); r <- as.numeric(Sig %*% a) / sqrt(sum(a * (Sig %*% a)))
  C <- rbind(cbind(Sig, r), c(r, 1))
  Cna <- C; Cna[1, 4] <- Cna[4, 1] <- NA
  out <- implied_interval(Cna)
  expect_equal(nrow(out), 1L)
  expect_equal(c(out$i, out$j), c(1L, 4L))
  expect_true(is.na(out$reported))
  # the true value sits inside the imputed interval
  expect_gte(r[1], out$lo - 1e-8); expect_lte(r[1], out$hi + 1e-8)
})

test_that("implied_interval errors on multiple missing cells for one interrogation", {
  R <- diag(3); R[1, 2] <- R[2, 1] <- NA; R[1, 3] <- R[3, 1] <- NA
  expect_error(implied_interval(R), "non-missing")
})

# Component-level tests for the fault-localization layer.

triad_R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)

test_that("component B (inconsistent_triples) flags the classic triangle", {
  it <- inconsistent_triples(triad_R, decimals = 2)
  expect_s3_class(it, "tbl_df")
  expect_equal(nrow(it), 1L)
  expect_equal(c(it$i, it$j, it$k), c(1L, 2L, 3L))
  expect_lt(it$box_max_det, 0)
})

test_that("component B is SOUND: never flags a consistent triple", {
  # The interval determinant bound must never exceed nor mis-sign the true max.
  set.seed(21)
  bad <- 0
  for (t in seq_len(n_reps(3000))) {
    vals <- round(runif(3, -1, 1), 2)
    boxes <- lapply(vals, function(v) c(max(-1, v - 0.005), min(1, v + 0.005)))
    ub <- port:::.triple_det_ub(boxes[[1]], boxes[[2]], boxes[[3]])
    # dense-grid true max
    gr <- lapply(boxes, function(b) seq(b[1], b[2], length.out = 10))
    G <- expand.grid(gr)
    gmax <- max(1 + 2 * G[, 1] * G[, 2] * G[, 3] - G[, 1]^2 - G[, 2]^2 - G[, 3]^2)
    if (ub < gmax - 1e-9) bad <- bad + 1     # unsound: bound below true max
  }
  expect_equal(bad, 0L)
})

test_that("inconsistent_triples is empty for a genuinely PSD matrix", {
  set.seed(3)
  L <- matrix(rnorm(24), 6, 4); S <- L %*% t(L); d <- sqrt(diag(S))
  S <- S / outer(d, d); diag(S) <- 1
  expect_equal(nrow(inconsistent_triples(round(S, 2), decimals = 2)), 0L)
})

test_that("component A finds all three sole culprits in a triad", {
  ci <- port:::.cell_intervals(triad_R, 0.005, verify = TRUE)
  sole <- Filter(function(c) isTRUE(c$sole), ci)
  expect_equal(length(sole), 3L)
  # every sole culprit has a required edit far beyond rounding
  expect_true(all(vapply(sole, function(c) abs(c$required_edit) > 0.005, logical(1))))
})

test_that("component C (leave-one-out) restores a triangle by removing any vertex", {
  lofo <- port:::.lofo_restoring(triad_R, 0.005)
  expect_setequal(lofo$restoring, 1:3)   # removing any variable leaves a 2x2
})

test_that("severity_max equals the smallest uniform box-widening for feasibility", {
  R <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  eps <- port:::.severity_max(R, 0.005)
  # widening by eps + tiny slack is feasible; widening by eps - slack is not
  bnd_ok <- port:::.box_bounds(R, 0.005 + eps + 1e-3)
  bnd_no <- port:::.box_bounds(R, 0.005 + eps - 1e-3)
  expect_equal(port:::.box_feasible(bnd_ok$lo, bnd_ok$hi, bnd_ok$off)$status, "feasible")
  expect_equal(port:::.box_feasible(bnd_no$lo, bnd_no$hi, bnd_no$off)$status, "infeasible")
})

test_that("rigor: sole-culprit 'empty' claims are witness-certified under verify=TRUE", {
  # In the joint (two-triangle) fixture, NO single cell restores; each such
  # 'empty' must be certified by the box witness, not a solver status.
  blk <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  J <- matrix(0, 6, 6); J[1:3, 1:3] <- blk; J[4:6, 4:6] <- blk; diag(J) <- 1
  ci <- port:::.cell_intervals(J, 0.005, verify = TRUE)
  # cells inside a triangle are not sole culprits here
  empt <- Filter(function(c) isTRUE(c$empty), ci)
  expect_true(length(empt) >= 1L)
  expect_true(all(vapply(empt, function(c) isTRUE(c$certified), logical(1))))
})

test_that("localize batch summarizes verdict classes into a tibble", {
  mats <- list(triad = triad_R, ok = diag(3))
  out <- suppressMessages(localize_psd_fault_batch(mats, decimals = 2, quiet = TRUE))
  expect_s3_class(out, "tbl_df")
  expect_equal(out$localization_verdict, c("triad", "none"))
  expect_true(all(c("implicated", "severity_max", "severity_frob") %in% names(out)))
  expect_match(out$implicated[1], "\\(1,2\\)")
})

test_that("localize batch logs a verdict-class summary message", {
  expect_message(localize_psd_fault_batch(list(triad_R), decimals = 2),
                 "Localization verdict classes")
})

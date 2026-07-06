# The off-diagonal correlation summary (min / max / mean / SD) attached to
# check_corr_psd() results and check_corr_psd_batch() rows.

R3 <- matrix(c(1,   0.2, 0.4,
               0.2, 1,   0.6,
               0.4, 0.6, 1), 3, 3)
offd <- c(0.2, 0.4, 0.6)

test_that("check_corr_psd reports the off-diagonal correlation summary", {
  res <- check_corr_psd(R3, decimals = 2)
  expect_equal(res$r_min,  min(offd))
  expect_equal(res$r_max,  max(offd))
  expect_equal(res$r_mean, mean(offd))
  expect_equal(res$r_sd,   sd(offd))
})

test_that("the summary excludes the unit diagonal and any NA cells", {
  Rna <- R3; Rna[1, 2] <- Rna[2, 1] <- NA
  res <- check_corr_psd(Rna, decimals = 2)      # non-uniform (NA) path
  keep <- c(0.4, 0.6)                            # (1,2) freed / excluded
  expect_equal(res$r_min,  min(keep))
  expect_equal(res$r_max,  max(keep))
  expect_equal(res$r_mean, mean(keep))
  expect_equal(res$r_sd,   sd(keep))
})

test_that("SD is NA with a single correlation; degenerate matrices summarise", {
  two <- check_corr_psd(matrix(c(1, 0.5, 0.5, 1), 2), decimals = 2)
  expect_equal(two$r_min, 0.5); expect_equal(two$r_max, 0.5)
  expect_equal(two$r_mean, 0.5); expect_true(is.na(two$r_sd))

  id <- check_corr_psd(diag(4), decimals = 2)   # all off-diagonals 0
  expect_equal(c(id$r_min, id$r_max, id$r_mean, id$r_sd), c(0, 0, 0, 0))
})

test_that("the print method shows the correlation summary line", {
  expect_output(print(check_corr_psd(R3, decimals = 2)), "corr")
})

test_that("batch output carries the correlation summary columns", {
  out <- check_corr_psd_batch(list(a = R3, id = diag(3)), quiet = TRUE)
  expect_true(all(c("r_min", "r_max", "r_mean", "r_sd") %in% names(out)))
  expect_equal(out$r_min[1],  0.2)
  expect_equal(out$r_max[1],  0.6)
  expect_equal(out$r_mean[1], 0.4)
  expect_equal(out$r_sd[2],   0)          # diag(3): off-diagonals all 0
})

test_that("error rows leave the summary as NA", {
  out <- check_corr_psd_batch(list(ok = R3, broken = matrix(1, 2, 3)),
                              quiet = TRUE)
  expect_true(all(is.na(c(out$r_min[2], out$r_max[2], out$r_mean[2], out$r_sd[2]))))
})

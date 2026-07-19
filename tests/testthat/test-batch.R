test_that("batch helper returns a tidy tibble with one row per matrix", {
  good <- diag(3)
  bad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3)
  oob <- matrix(c(1, 1.3, 0, 1.3, 1, 0, 0, 0, 1), 3, 3)
  out <- check_corr_psd_batch(
    list(good = good, bad = bad, oob = oob),
    decimals = 2,
    quiet = TRUE
  )
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 3L)
  expect_equal(out$id, c("good", "bad", "oob"))
  expect_equal(out$verdict, c("consistent", "inconsistent", "inconsistent"))
  expect_equal(out$tier, c("witness", "witness", "precheck"))
  expect_true(all(
    c("margin", "b_upper", "delta", "p", "message") %in% names(out)
  ))
})

test_that("batch helper indexes unnamed lists positionally", {
  out <- check_corr_psd_batch(list(diag(2), diag(3)), quiet = TRUE)
  expect_equal(out$id, c("1", "2"))
  expect_equal(out$p, c(2L, 3L))
})

test_that("batch helper records validation errors as rows by default", {
  out <- check_corr_psd_batch(
    list(ok = diag(3), broken = matrix(1, 2, 3)),
    quiet = TRUE
  )
  expect_equal(out$verdict, c("consistent", "error"))
  expect_true(!is.na(out$message[2]))
})

test_that("batch helper can re-raise errors when asked", {
  expect_error(
    check_corr_psd_batch(
      list(matrix(1, 2, 3)),
      on_error = "stop",
      quiet = TRUE
    ),
    "square"
  )
})

test_that("batch helper logs escalation-tier frequency", {
  expect_message(
    check_corr_psd_batch(list(diag(3)), decimals = 2),
    "POCS escalation tier fired"
  )
})

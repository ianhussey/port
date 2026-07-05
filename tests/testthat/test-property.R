# Anti-false-positive property test.
#
# A genuinely PSD correlation matrix rounded to d decimals always leaves the
# ORIGINAL matrix inside the rounding box (each entry moves by at most
# delta = 0.5*10^(-d)). Hence a PSD matrix exists in the box and the tool must
# NEVER return "impossible". Because the witness bound is a sound over-estimate,
# this is guaranteed in exact arithmetic; the test guards against an
# over-firing floating-point error model.

test_that("genuine PSD matrices rounded to 2dp are never impossible", {
  set.seed(2024)
  n_trials <- 1500
  seen <- c(impossible = 0L, possible = 0L, undecided = 0L)
  for (t in seq_len(n_trials)) {
    p <- sample(2:8, 1)
    k <- sample(1:p, 1)
    L <- matrix(rnorm(p * k), p, k)
    S <- L %*% t(L) + diag(runif(p, 0.01, 0.5))
    d <- sqrt(diag(S))
    S <- S / outer(d, d)
    diag(S) <- 1
    stopifnot(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
    Rr <- round(S, 2)
    res <- check_corr_psd(Rr, decimals = 2)
    seen[res$verdict] <- seen[res$verdict] + 1L
    expect_false(identical(res$verdict, "impossible"),
                 info = paste("false positive at trial", t))
  }
  # sanity: the overwhelming majority resolve to a definite "possible"
  expect_gt(seen["possible"], 0.9 * n_trials)
})

test_that("rounding-induced near-boundary cases are handled without false impossibility", {
  # Deliberately construct matrices whose smallest eigenvalue is ~0 before
  # rounding (the hardest anti-false-positive regime).
  set.seed(99)
  for (t in seq_len(400)) {
    p <- sample(3:6, 1)
    A <- matrix(rnorm(p * p), p, p)
    A <- crossprod(A)                       # PSD
    eg <- eigen(A, symmetric = TRUE)
    eg$values[p] <- 0                       # force exact singularity
    S <- eg$vectors %*% (eg$values * t(eg$vectors))
    dd <- diag(S)
    if (any(dd <= 1e-8)) next
    S <- S / outer(sqrt(dd), sqrt(dd)); diag(S) <- 1
    res <- check_corr_psd(round(S, 2), decimals = 2)
    expect_false(identical(res$verdict, "impossible"))
  }
})

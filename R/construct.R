# -----------------------------------------------------------------------------
# Constructive consistency certificates.
#
# We try to exhibit an explicit in-box matrix and *verify* it PSD with
# .verify_psd(). Any success is a rigorous proof of consistency (an in-box PSD
# matrix exists). Failure of all candidates is inconclusive, not a proof of
# inconsistency.
# -----------------------------------------------------------------------------

# Is X inside the rounding box around R (off-diagonals within delta, clipped to
# [-1, 1]; unit diagonal)? Used to double-check solver / constructed solutions.
.in_box <- function(X, R, delta, tol = 1e-9) {
  p <- nrow(R)
  if (any(abs(diag(X) - 1) > tol)) return(FALSE)
  off <- upper.tri(R) | lower.tri(R)
  lo <- pmax(R[off] - delta, -1) - tol
  hi <- pmin(R[off] + delta, 1) + tol
  x <- X[off]
  all(x >= lo & x <= hi)
}

# Shrink every off-diagonal toward 0 by up to `amount` (staying in the box):
# the most diagonally-dominant in-box matrix for amount = delta.
.shrink_offdiag <- function(R, amount) {
  X <- R
  off <- upper.tri(R) | lower.tri(R)
  X[off] <- sign(R[off]) * pmax(0, abs(R[off]) - amount)
  X
}

# Project R onto the PSD cone (clip negative eigenvalues), restore the unit
# diagonal, then clamp off-diagonals back into the box. The result is in-box by
# construction; whether it is PSD is decided by .verify_psd().
.project_clip <- function(R, delta) {
  eg <- eigen(R, symmetric = TRUE)
  lam <- pmax(eg$values, 0)
  X <- eg$vectors %*% (lam * t(eg$vectors))
  X <- (X + t(X)) / 2
  diag(X) <- 1
  off <- upper.tri(R) | lower.tri(R)
  lo <- pmax(R[off] - delta, -1)
  hi <- pmin(R[off] + delta, 1)
  X[off] <- pmin(pmax(X[off], lo), hi)
  X
}

# Try a small family of in-box candidates and return the first verified-PSD one
# as a consistency certificate, or NULL.
.construct_consistent <- function(R, delta) {
  candidates <- list(
    list(X = R,                          how = "the reported matrix itself is PSD"),
    list(X = .shrink_offdiag(R, delta),  how = "an in-box matrix shrunk toward the diagonal is PSD"),
    list(X = .shrink_offdiag(R, delta / 2), how = "a partially shrunk in-box matrix is PSD"),
    list(X = .project_clip(R, delta),    how = "the box-clipped PSD projection of R is PSD")
  )
  for (cand in candidates) {
    if (.in_box(cand$X, R, delta) && .verify_psd(cand$X)) return(cand)
  }
  NULL
}

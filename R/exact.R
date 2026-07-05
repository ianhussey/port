# -----------------------------------------------------------------------------
# Exact small-p check (cross-validation only).
#
# For p in {2, 3} we can decide consistency exactly by maximizing the relevant
# principal minors over the rounding box. This is used in the test suite to
# cross-check the witness verdict; the witness path remains the source of truth
# in production. Returns "consistent" or "inconsistent".
# -----------------------------------------------------------------------------

# Per-variable box (off-diagonal), clipped to [-1, 1].
.offdiag_box <- function(r, delta) {
  c(lo = max(-1, r - delta), hi = min(1, r + delta))
}

# Maximize the 3x3 correlation determinant
#   f(x, y, z) = 1 + 2xyz - x^2 - y^2 - z^2
# over the product box. f is a downward parabola in each variable (vertex at the
# product of the other two), so coordinate ascent from every corner plus the
# centre reliably finds the global maximum over these tiny boxes.
.max_det_p3 <- function(lo, hi) {
  fval <- function(v) 1 + 2 * v[1] * v[2] * v[3] - sum(v^2)
  clamp1 <- function(t, k) min(max(t, lo[k]), hi[k])
  corners <- as.matrix(expand.grid(c(lo[1], hi[1]), c(lo[2], hi[2]), c(lo[3], hi[3])))
  starts <- rbind(corners, (lo + hi) / 2)
  best <- -Inf
  for (r in seq_len(nrow(starts))) {
    v <- as.numeric(starts[r, ])
    for (iter in seq_len(200L)) {
      vprev <- v
      v[1] <- clamp1(v[2] * v[3], 1L)
      v[2] <- clamp1(v[1] * v[3], 2L)
      v[3] <- clamp1(v[1] * v[2], 3L)
      if (max(abs(v - vprev)) < 1e-15) break
    }
    best <- max(best, fval(v))
  }
  for (r in seq_len(nrow(corners))) best <- max(best, fval(as.numeric(corners[r, ])))
  best
}

# Exact consistency verdict for p in {2, 3}. Errors for other p.
.exact_consistent <- function(R, delta, tol = 1e-12) {
  p <- nrow(R)
  if (p == 2L) {
    b <- .offdiag_box(R[1, 2], delta)
    return(if (b["lo"] <= b["hi"] + tol) "consistent" else "inconsistent")
  }
  if (p == 3L) {
    bx <- .offdiag_box(R[1, 2], delta)
    by <- .offdiag_box(R[1, 3], delta)
    bz <- .offdiag_box(R[2, 3], delta)
    # Any empty per-variable box already forces inconsistency.
    if (bx["lo"] > bx["hi"] || by["lo"] > by["hi"] || bz["lo"] > bz["hi"]) {
      return("inconsistent")
    }
    lo <- c(bx["lo"], by["lo"], bz["lo"])
    hi <- c(bx["hi"], by["hi"], bz["hi"])
    best <- .max_det_p3(lo, hi)
    return(if (best >= -tol) "consistent" else "inconsistent")
  }
  stop("The exact cross-check supports only p in {2, 3}.", call. = FALSE)
}

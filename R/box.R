# -----------------------------------------------------------------------------
# Generalized rounding-box infrastructure for fault localization.
#
# The base package works with a uniform box (centre R, half-width delta). The
# localization layer needs *heterogeneous* boxes: some cells freed to [-1, 1],
# some pinned to a value, others held to their rounding interval. Everything
# below accepts explicit per-cell bounds `lo`, `hi` (p x p matrices, diagonal
# ignored) plus the off-diagonal mask `off`.
#
# Two kinds of claim, two kinds of evidence (mirroring the base package):
#   * feasibility ("an in-box PSD matrix exists") is an EXISTENCE claim, proven
#     by exhibiting a Rump-verified matrix (found by POCS);
#   * infeasibility ("no in-box matrix is PSD") is proven by a sound witness
#     vector via .witness_box_bound(), never by a solver status.
# POCS is only ever a search; the certificates are witness + Rump.
# -----------------------------------------------------------------------------

# Build the reported-value box for the general front door: per-cell decimals or
# delta (scalar or p x p matrix), a rounding rule, and NA cells freed to [-1,1].
# Interval per reported value r at width w = 10^(-d):
#   nearest : [r - w/2, r + w/2]        (error of round-half at most w/2)
#   floor   : [r, r + w]
#   ceiling : [r - w, r]
#   truncate: toward zero -- r >= 0 gives [r, r + w]; r < 0 gives [r - w, r];
#             r == 0 gives [-w, w] (either sign truncates to zero).
# Open interval ends are closed (a sound superset). Returns raw (unclipped) and
# clipped bounds, the NA mask, per-cell emptiness (raw interval entirely outside
# [-1, 1]), and the maximum half-width (used for tau and reporting).
.reported_box <- function(R, decimals, delta, rounding) {
  p <- nrow(R)
  off <- upper.tri(R) | lower.tri(R)
  na_mask <- is.na(R) & off

  if (!is.null(delta)) {
    dmat <- if (is.matrix(delta)) delta else matrix(delta, p, p)
    lo_raw <- R - dmat
    hi_raw <- R + dmat
  } else {
    dec <- if (is.matrix(decimals)) decimals else matrix(decimals, p, p)
    w <- 10^(-dec)
    if (identical(rounding, "nearest")) {
      lo_raw <- R - w / 2
      hi_raw <- R + w / 2
    } else if (identical(rounding, "floor")) {
      lo_raw <- R
      hi_raw <- R + w
    } else if (identical(rounding, "ceiling")) {
      lo_raw <- R - w
      hi_raw <- R
    } else {
      # truncate (toward zero)
      lo_raw <- ifelse(R > 0, R, R - w)
      hi_raw <- ifelse(R < 0, R, R + w)
    }
  }
  lo_raw[na_mask] <- -1
  hi_raw[na_mask] <- 1

  empty <- off & !na_mask & (lo_raw > 1 | hi_raw < -1)
  lo <- pmax(lo_raw, -1)
  hi <- pmin(hi_raw, 1)
  diag(lo) <- 1
  diag(hi) <- 1
  # reporting half-width: over the REPORTED cells (freed NA cells excluded)
  reported <- off & !na_mask
  half_max <- if (any(reported)) max((hi_raw - lo_raw)[reported]) / 2 else 1

  list(
    lo = lo,
    hi = hi,
    lo_raw = lo_raw,
    hi_raw = hi_raw,
    off = off,
    na_mask = na_mask,
    empty = empty,
    half_max = half_max
  )
}

# Off-diagonal bounds for the uniform rounding box around R.
.box_bounds <- function(R, delta) {
  p <- nrow(R)
  off <- upper.tri(R) | lower.tri(R)
  lo <- matrix(0, p, p)
  hi <- matrix(0, p, p)
  lo[off] <- pmax(R[off] - delta, -1)
  hi[off] <- pmin(R[off] + delta, 1)
  list(lo = lo, hi = hi, off = off)
}

# Is X inside the (heterogeneous) box? Unit diagonal, off-diagonals in [lo, hi].
.in_box_generic <- function(X, lo, hi, off, tol = 1e-9) {
  if (any(abs(diag(X) - 1) > tol)) {
    return(FALSE)
  }
  x <- X[off]
  all(x >= lo[off] - tol & x <= hi[off] + tol)
}

# POCS search for an in-box PSD matrix in a heterogeneous box. Returns a
# Rump-verified in-box matrix (list with X, mu) or NULL. NULL means "not found",
# NOT "infeasible"; infeasibility must be established by .box_inconsistent().
.pocs_feasible <- function(
  lo,
  hi,
  off,
  mus = c(1e-2, 1e-3, 1e-4, 1e-6),
  max_iter = 1000L,
  tol = 1e-13
) {
  p <- nrow(lo)
  start <- matrix(0, p, p)
  start[off] <- (lo[off] + hi[off]) / 2
  diag(start) <- 1
  for (mu in mus) {
    X <- .proj_box(start, lo, hi, off)
    for (k in seq_len(max_iter)) {
      Xn <- .proj_box(.proj_psd_margin(X, mu), lo, hi, off)
      change <- max(abs(Xn - X))
      X <- Xn
      if (change < tol) break
    }
    if (.in_box_generic(X, lo, hi, off) && .verify_psd(X)) {
      return(list(X = X, mu = mu))
    }
  }
  NULL
}

# Sound upper bound on max_{X in box} v'Xv for a heterogeneous box (diag = 1).
# For each off-diagonal pair the box-optimal entry is hi_st when v_s v_t >= 0
# and lo_st otherwise. If B_upper < 0 then no in-box matrix is PSD along v.
# A rigorous a-priori forward-error slack is added (Higham gamma bounds), so the
# verdict is floating-point-safe (same model as the base witness bound).
.witness_box_bound <- function(lo, hi, off, v) {
  p <- length(v)
  u <- .unit_roundoff()
  vv <- outer(v, v) # v_s v_t
  best <- ifelse(vv >= 0, hi, lo) # box-optimal entry per pair
  best[!off] <- 0
  diag_term <- sum(v * v) # sum_i v_i^2 (diagonal = 1)
  cross <- sum(vv[off] * best[off]) # sum over ALL ordered off pairs = 2 * sum_{s<t}
  M_hat <- diag_term + cross
  # Magnitude scaffold for the error bound.
  P <- sum(v * v) + sum(abs(vv[off]) * abs(best[off]))
  E <- (.gamma(2L * p * p + 2L) + u) * P
  E <- E * 1.000001
  list(v = v, M_hat = M_hat, E = E, B_upper = M_hat + E)
}

# Fixed-point polish for the heterogeneous box (see .polish_witness in
# witness.R): X(v) picks the box-optimal entry per pair, the next direction is
# its bottom eigenvector. Sound by construction (every iterate is evaluated
# through the rigorous box bound).
.polish_witness_box <- function(lo, hi, off, v, max_iter = 8L) {
  best <- .witness_box_bound(lo, hi, off, v)
  p <- nrow(lo)
  for (k in seq_len(max_iter)) {
    vv <- outer(v, v)
    X <- ifelse(vv >= 0, hi, lo)
    X[!off] <- 0
    diag(X) <- 1
    v_new <- eigen(X, symmetric = TRUE)$vectors[, p]
    wb <- .witness_box_bound(lo, hi, off, v_new)
    if (wb$B_upper < best$B_upper - 1e-15) {
      best <- wb
      v <- v_new
    } else {
      break
    }
  }
  best
}

# Decide whether a heterogeneous box is inconsistent (no in-box PSD matrix), using
# a search over sound witness directions. Returns list(inconsistent, witness,
# b_upper, margin). Directions: bottom eigenvectors of the box midpoint, all
# coordinate pairs, (small p) 3x3 submatrix eigenvectors, and the p regression-
# residual directions of the midpoint -- the same family that powers the base
# checker, all evaluated with the rigorous box bound, then polished.
.box_inconsistent <- function(lo, hi, off, triples_max_p = 12L) {
  p <- nrow(lo)
  mid <- matrix(0, p, p)
  mid[off] <- (lo[off] + hi[off]) / 2
  diag(mid) <- 1

  cands <- list()
  eg <- eigen(mid, symmetric = TRUE)
  for (j in seq_len(min(p, 4L))) {
    cands[[length(cands) + 1L]] <- eg$vectors[, p - j + 1L]
  }
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      s <- if (mid[i, j] >= 0) 1 else -1
      vv <- numeric(p)
      vv[i] <- 1 / sqrt(2)
      vv[j] <- -s / sqrt(2)
      cands[[length(cands) + 1L]] <- vv
    }
  }
  if (p >= 3L && p <= triples_max_p) {
    combs <- utils::combn(p, 3L)
    for (c in seq_len(ncol(combs))) {
      S <- combs[, c]
      es <- eigen(mid[S, S, drop = FALSE], symmetric = TRUE)
      vv <- numeric(p)
      vv[S] <- es$vectors[, 3L]
      cands[[length(cands) + 1L]] <- vv
    }
  }
  for (i in seq_len(p)) {
    d <- .rsquared_direction(mid, i)
    if (!is.null(d)) cands[[length(cands) + 1L]] <- d$v
  }

  best <- NULL
  for (v in cands) {
    if (all(v == 0)) {
      next
    }
    wb <- .witness_box_bound(lo, hi, off, v)
    if (is.null(best) || wb$B_upper < best$B_upper) best <- wb
  }
  if (!is.null(best)) {
    pol <- .polish_witness_box(lo, hi, off, best$v)
    if (pol$B_upper < best$B_upper) best <- pol
  }
  list(
    inconsistent = !is.null(best) && best$B_upper < 0,
    witness = if (is.null(best)) NULL else best$v,
    b_upper = if (is.null(best)) NA_real_ else best$B_upper,
    margin = if (is.null(best)) NA_real_ else best$M_hat
  )
}

# Feasibility oracle used by the bisection localizers. Combines the two sound
# certificates: try to EXHIBIT a verified in-box PSD matrix (POCS), else try to
# REFUTE via the box witness. Returns:
#   status "feasible" (with $X), "infeasible" (with $witness), or "undecided".
# When verify = FALSE, a POCS point that fails Rump is still reported feasible
# (labelled solver/search-based) so callers can run cheaply.
.box_feasible <- function(lo, hi, off, verify = TRUE) {
  # Witness first: a sound refutation is cheap (O(p^2) directions) and lets the
  # bisection localizers reject infeasible pins without a full POCS search.
  imp <- .box_inconsistent(lo, hi, off)
  if (isTRUE(imp$inconsistent)) {
    return(list(
      status = "infeasible",
      witness = imp$witness,
      b_upper = imp$b_upper,
      certified = TRUE
    ))
  }
  hit <- .pocs_feasible(lo, hi, off)
  if (!is.null(hit)) {
    return(list(status = "feasible", X = hit$X, certified = TRUE))
  }
  if (!verify) {
    # Last-resort unverified search: a single low-margin POCS pass.
    p <- nrow(lo)
    start <- matrix(0, p, p)
    start[off] <- (lo[off] + hi[off]) / 2
    diag(start) <- 1
    X <- start
    for (k in seq_len(2000L)) {
      Xn <- .proj_box(.proj_psd_margin(X, 0), lo, hi, off)
      if (max(abs(Xn - X)) < 1e-13) {
        X <- Xn
        break
      }
      X <- Xn
    }
    gap <- max(abs(.proj_psd_margin(X, 0) - X))
    if (gap < 1e-7) return(list(status = "feasible", X = X, certified = FALSE))
  }
  list(status = "undecided", X = NULL)
}

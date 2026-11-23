# -----------------------------------------------------------------------------
# Self-contained consistency search by projections onto convex sets (POCS).
#
# Escalation tier used only in the ambiguous zone, when neither the witness
# search nor the cheap constructions in construct.R settled the verdict. It is
# fully self-contained: no external solver or optional dependency.
#
# The feasible set is the intersection of two closed convex sets:
#   * the box  B = { X symmetric : X_ii = 1, X_ij in [lo_ij, hi_ij] };
#   * the tightened PSD cone C_mu = { X : X - mu*I is PSD } for a margin mu >= 0.
# Alternating projections X <- P_B(P_{C_mu}(X)) converge to a point in B intersect C_mu
# whenever that intersection is nonempty (von Neumann / Cheney-Goldstein for two
# closed convex sets). A strictly-positive margin mu makes the recovered point
# strictly positive definite, so the one-sided Cholesky test of verify.R can
# certify it. The witness path remains the source of truth: POCS is used only to
# *find* a candidate, which is then independently Rump-verified. A failure to
# find a verifiable point is inconclusive, never an inconsistency certificate
# by itself. .pocs_witness() (below) can still certify inconsistency, but only
# through the rigorous witness bound, like any other test direction.
# -----------------------------------------------------------------------------

# Euclidean projection onto the box (unit diagonal; off-diagonals clipped).
.proj_box <- function(X, lo, hi, off) {
  X[off] <- pmin(pmax(X[off], lo[off]), hi[off])
  diag(X) <- 1
  X
}

# Euclidean projection onto { Y : Y - mu*I is PSD } (eigenvalue thresholding).
.proj_psd_margin <- function(X, mu) {
  X <- (X + t(X)) / 2
  eg <- eigen(X, symmetric = TRUE)
  Y <- eg$vectors %*% (pmax(eg$values, mu) * t(eg$vectors))
  (Y + t(Y)) / 2
}

# Try to certify consistency via POCS. Returns a list(X, mu, how) with a
# verified in-box PSD matrix, or NULL. `mus` is a decreasing sequence of PD
# margins: a larger feasible margin yields a comfortably-verifiable point, and
# smaller margins are tried when the box only just reaches the PSD cone.
.pocs_consistent <- function(
  R,
  delta,
  mus = c(1e-2, 1e-3, 1e-4, 1e-6),
  max_iter = 1000L,
  tol = 1e-12
) {
  p <- nrow(R)
  off <- upper.tri(R) | lower.tri(R)
  lo <- matrix(0, p, p)
  hi <- matrix(0, p, p)
  lo[off] <- pmax(R[off] - delta, -1)
  hi[off] <- pmin(R[off] + delta, 1)

  for (mu in mus) {
    X <- .proj_box(R, lo, hi, off)
    for (k in seq_len(max_iter)) {
      Xn <- .proj_box(.proj_psd_margin(X, mu), lo, hi, off)
      change <- max(abs(Xn - X))
      X <- Xn
      if (change < tol) break
    }
    # verify.R's one-sided Cholesky test is the rigorous gate: only an in-box
    # matrix it certifies PSD is accepted as a consistency certificate.
    if (.in_box(X, R, delta) && .verify_psd(X)) {
      return(list(
        X = X,
        mu = mu,
        how = sprintf(
          paste(
            "alternating projections found an in-box matrix",
            "(PD margin mu = %g), independently verified PSD"
          ),
          mu
        )
      ))
    }
  }
  NULL
}

# Witness from alternating projections between the box and the PSD cone, run
# without a margin. If the box is inconsistent, the iterates approach a closest
# pair (X in the box, its projection Y on the cone), and Y - X, the negative
# eigen-part of X with the sign flipped, separates the box from the cone. The
# bottom eigenvector of X is its leading rank-one term. It is only a test
# direction: `bound(v)` must return a rigorous bound record with element
# `B_upper`, and every direction goes through it, so failure to converge costs
# power, never soundness. The direction is tested every `every` iterations and
# the search stops at the first certificate. Returns the last record tested.
# The iteration converges slowly near the boundary: on held-out matrices a cap
# of 1000 iterations would leave 9 of 546 open cases undecided, and 10000
# leaves 6 (special-scripts/sparse-witness-simulation.R).
# ponytail: up to max_iter eigendecompositions (about 17 s at p = 100 if no
# certificate turns up), so check_corr_psd() uses it only
# where the verdict would otherwise be "undecided", and the localizers' oracle
# does not use it at all. Add it there if undecided oracle calls turn out to
# matter. The full matrix Y - X is a stronger certificate than its leading
# term, but the result object holds a single witness vector.
.pocs_witness <- function(
  lo,
  hi,
  off,
  bound,
  max_iter = 10000L,
  every = 100L,
  tol = 1e-13
) {
  p <- nrow(lo)
  X <- matrix(0, p, p)
  X[off] <- (lo[off] + hi[off]) / 2
  diag(X) <- 1
  for (k in seq_len(max_iter)) {
    Xn <- .proj_box(.proj_psd_margin(X, 0), lo, hi, off)
    done <- k == max_iter || max(abs(Xn - X)) < tol
    X <- Xn
    if (done || k %% every == 0L) {
      wb <- bound(eigen(X, symmetric = TRUE)$vectors[, p])
      if (done || wb$B_upper < 0) break
    }
  }
  wb
}

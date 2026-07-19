# -----------------------------------------------------------------------------
# Component D: sparse minimal-cardinality correction + severity (all POCS-based).
#
# The spec's D1/D2/D3 are convex programs over box-excess; we realize their
# INTENT with the feasibility oracle instead of an SDP solver:
#   * sparse_support (D1, "minimal-cardinality correction"): the smallest set of
#     cells whose freeing beyond their boxes restores feasibility. A greedy /
#     bounded-exhaustive search over the POCS+witness oracle. Its cardinality-1
#     pass is exactly component A's sole-culprit test, so A and D1 stay
#     consistent. Minimality is heuristic; every feasibility fact in it is
#     rigorously certified (Rump-verified restoring matrix).
#   * severity_max (D2, min-max box-excess): the smallest UNIFORM box-widening
#     eps that makes the box feasible -- located exactly by bisection. This is
#     "the largest single edit required, minimized over corrections".
#   * severity_frob (D3, total mass): an achievable total beyond-rounding
#     movement, from the nearest PSD correlation matrix; an upper bound on the
#     true minimum, reported as such.
# -----------------------------------------------------------------------------

# Nearest-correlation-matrix by alternating PSD / unit-diagonal projection
# (Higham-style, simple variant). Used only to report an achievable Frobenius
# severity; not on any rigorous path.
.nearcorr <- function(R, max_iter = 500L, tol = 1e-10) {
  X <- R
  off <- upper.tri(R) | lower.tri(R)
  for (k in seq_len(max_iter)) {
    Y <- .proj_psd_margin(X, 0)
    Y[off] <- pmin(pmax(Y[off], -1), 1)
    diag(Y) <- 1
    if (max(abs(Y - X)) < tol) {
      X <- Y
      break
    }
    X <- Y
  }
  X
}

# excess_ij(X) = distance of X_ij from the rounding interval [R_ij +/- delta].
.excess_offdiag <- function(X, R, delta) {
  p <- nrow(R)
  idx <- which(upper.tri(R), arr.ind = TRUE)
  ex <- numeric(nrow(idx))
  for (r in seq_len(nrow(idx))) {
    i <- idx[r, 1]
    j <- idx[r, 2]
    ex[r] <- max(0, X[i, j] - (R[i, j] + delta), (R[i, j] - delta) - X[i, j])
  }
  list(idx = idx, excess = ex)
}

# severity_frob (D3): achievable total beyond-rounding mass (Frobenius).
.severity_frob <- function(R, delta) {
  Xnc <- .nearcorr(R)
  ex <- .excess_offdiag(Xnc, R, delta)$excess
  sqrt(sum(ex^2))
}

# severity_max (D2): smallest uniform box-widening eps giving feasibility.
.severity_max <- function(R, delta, verify = TRUE, tol = 1e-6, max_it = 40L) {
  p <- nrow(R)
  off <- upper.tri(R) | lower.tri(R)
  feas <- function(eps) {
    lo <- matrix(0, p, p)
    hi <- matrix(0, p, p)
    lo[off] <- pmax(R[off] - delta - eps, -1)
    hi[off] <- pmin(R[off] + delta + eps, 1)
    identical(.box_feasible(lo, hi, off, verify = verify)$status, "feasible")
  }
  lo_e <- 0
  hi_e <- 2 # eps = 2 makes every box [-1,1]
  for (it in seq_len(max_it)) {
    if (hi_e - lo_e < tol) {
      break
    }
    m <- (lo_e + hi_e) / 2
    if (feas(m)) hi_e <- m else lo_e <- m
  }
  hi_e
}

# best achievable lambda_min over the box (search-based, POCS): the largest
# margin mu for which the box intersects { X : X - mu I is PSD }. Reported as a
# secondary, less-intuitive severity.
.best_lambda_min <- function(R, delta, tol = 1e-6, max_it = 40L) {
  p <- nrow(R)
  bnd <- .box_bounds(R, delta)
  lo <- bnd$lo
  hi <- bnd$hi
  off <- bnd$off
  feas_margin <- function(mu) {
    start <- matrix(0, p, p)
    start[off] <- (lo[off] + hi[off]) / 2
    diag(start) <- 1
    X <- .proj_box(start, lo, hi, off)
    for (k in seq_len(1000L)) {
      Xn <- .proj_box(.proj_psd_margin(X, mu), lo, hi, off)
      if (max(abs(Xn - X)) < 1e-12) {
        X <- Xn
        break
      }
      X <- Xn
    }
    max(abs(.proj_psd_margin(X, mu) - X)) < 1e-7
  }
  lo_m <- -1
  hi_m <- 0 # feasible for very negative mu
  for (it in seq_len(max_it)) {
    if (hi_m - lo_m < tol) {
      break
    }
    m <- (lo_m + hi_m) / 2
    if (feas_margin(m)) lo_m <- m else hi_m <- m
  }
  lo_m
}

# Does freeing the given set of cells (to [-1,1]) restore feasibility? Returns
# the Rump-verified restoring matrix or NULL.
.free_set_box <- function(bnd, cells) {
  lo <- bnd$lo
  hi <- bnd$hi
  for (r in seq_len(nrow(cells))) {
    i <- cells[r, 1]
    j <- cells[r, 2]
    lo[i, j] <- lo[j, i] <- -1
    hi[i, j] <- hi[j, i] <- 1
  }
  list(lo = lo, hi = hi, off = bnd$off)
}
.set_restores <- function(bnd, cells, verify) {
  b <- .free_set_box(bnd, cells)
  fo <- .box_feasible(b$lo, b$hi, b$off, verify = verify)
  if (identical(fo$status, "feasible")) fo$X else NULL
}

# Sparse minimal-cardinality support. `sole_cells` is a matrix of the cells whose
# singleton freeing already restores (from component A). Returns list(cells,
# magnitudes, cardinality) or NULL if no set of size <= sparse_k restores.
.sparse_support <- function(
  R,
  delta,
  sole_cells,
  sole_edits,
  verify = TRUE,
  sparse_k = 3L,
  max_cells_exhaustive = 15L
) {
  bnd <- .box_bounds(R, delta)
  allcells <- which(upper.tri(R), arr.ind = TRUE)
  m <- nrow(allcells)

  score_set <- function(cells, X) {
    # total beyond-box movement of the restoring matrix on the freed cells
    tot <- 0
    for (r in seq_len(nrow(cells))) {
      i <- cells[r, 1]
      j <- cells[r, 2]
      tot <- tot +
        max(0, X[i, j] - (R[i, j] + delta), (R[i, j] - delta) - X[i, j])
    }
    tot
  }

  # Cardinality 1: use the sole culprits found in A (avoid recomputation).
  if (nrow(sole_cells) > 0L) {
    best_r <- which.min(abs(sole_edits))
    c1 <- sole_cells[best_r, , drop = FALSE]
    return(list(
      cells = c1,
      magnitudes = abs(sole_edits[best_r]),
      cardinality = 1L
    ))
  }

  # Cardinality 2..sparse_k: bounded-exhaustive when the cell count is small,
  # greedy otherwise.
  if (m <= max_cells_exhaustive) {
    for (kk in 2:max(2L, sparse_k)) {
      if (kk > m) {
        break
      }
      combs <- utils::combn(m, kk)
      best <- NULL
      best_score <- Inf
      for (cc in seq_len(ncol(combs))) {
        cells <- allcells[combs[, cc], , drop = FALSE]
        X <- .set_restores(bnd, cells, verify)
        if (!is.null(X)) {
          sc <- score_set(cells, X)
          if (sc < best_score) {
            best_score <- sc
            best <- list(cells = cells, X = X)
          }
        }
      }
      if (!is.null(best)) {
        mags <- vapply(
          seq_len(nrow(best$cells)),
          function(r) {
            i <- best$cells[r, 1]
            j <- best$cells[r, 2]
            max(
              0,
              best$X[i, j] - (R[i, j] + delta),
              (R[i, j] - delta) - best$X[i, j]
            )
          },
          numeric(1)
        )
        return(list(cells = best$cells, magnitudes = mags, cardinality = kk))
      }
    }
    return(NULL)
  }

  # Greedy for larger matrices: add the cell that most reduces the best achievable
  # margin deficit until feasible or sparse_k reached.
  chosen <- matrix(integer(0), ncol = 2)
  remaining <- allcells
  for (kk in seq_len(sparse_k)) {
    best_cell <- NULL
    best_X <- NULL
    best_score <- Inf
    for (r in seq_len(nrow(remaining))) {
      cand <- rbind(chosen, remaining[r, , drop = FALSE])
      X <- .set_restores(bnd, cand, verify)
      if (!is.null(X)) {
        best_cell <- remaining[r, , drop = FALSE]
        best_X <- X
        best_score <- 0
        break
      }
    }
    if (!is.null(best_X)) {
      cells <- rbind(chosen, best_cell)
      mags <- vapply(
        seq_len(nrow(cells)),
        function(r) {
          i <- cells[r, 1]
          j <- cells[r, 2]
          max(
            0,
            best_X[i, j] - (R[i, j] + delta),
            (R[i, j] - delta) - best_X[i, j]
          )
        },
        numeric(1)
      )
      return(list(cells = cells, magnitudes = mags, cardinality = nrow(cells)))
    }
    # otherwise greedily commit the cell whose addition most lowers the witness
    # violation (approximate); pick the first remaining as a simple heuristic
    chosen <- rbind(chosen, remaining[1, , drop = FALSE])
    remaining <- remaining[-1, , drop = FALSE]
    if (nrow(remaining) == 0L) break
  }
  NULL
}

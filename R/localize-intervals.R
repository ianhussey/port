# -----------------------------------------------------------------------------
# Component A: per-cell implied interval + single-cell sufficiency (box-aware).
#
# For each off-diagonal (i,j): FREE that cell to [-1, 1], hold all OTHER
# off-diagonals to their rounding boxes, diag = 1, and compute the interval of
# X_ij for which some in-box setting of the others makes X PSD:
#     [lo_ij, hi_ij] = { x : exists in-box others with X_ij = x and X PSD }.
# The box+PSD feasible set is convex, so its projection onto axis (i,j) is this
# interval; we locate its endpoints by bisection on the feasibility oracle.
#   * Empty   -> freeing (i,j) alone cannot restore PSD; (i,j) is not a sole culprit.
#   * Nonempty-> (i,j) is a sole-culprit candidate. required_edit is the signed
#     distance from the cell's rounding box to [lo_ij, hi_ij].
#
# Rigor: "empty" is certified by the sound box witness; "nonempty" exhibits a
# Rump-verified matrix. Endpoints are converged from the certified-feasible side,
# so [lo_ij, hi_ij] is a (slightly conservative) subset of the true interval and
# required_edit is never understated.
#
# Note: for an impossible-given-rounding matrix, no cell's rounding box can
# overlap its own feasible interval (an overlap would exhibit a fully in-box PSD
# matrix), so every sole-culprit cell has a strictly non-zero required_edit.
# -----------------------------------------------------------------------------

.cell_free_box <- function(bnd, i, j) {
  lo <- bnd$lo; hi <- bnd$hi
  lo[i, j] <- lo[j, i] <- -1
  hi[i, j] <- hi[j, i] <-  1
  list(lo = lo, hi = hi, off = bnd$off)
}

.cell_pin_box <- function(bnd, i, j, x) {
  lo <- bnd$lo; hi <- bnd$hi
  lo[i, j] <- lo[j, i] <- x
  hi[i, j] <- hi[j, i] <- x
  list(lo = lo, hi = hi, off = bnd$off)
}

.feasible_status_at <- function(bnd, i, j, x, verify) {
  b <- .cell_pin_box(bnd, i, j, x)
  .box_feasible(b$lo, b$hi, b$off, verify = verify)$status
}

# Bisect for the feasibility boundary between a known-infeasible x and a
# known-feasible x; returns the certified-feasible endpoint (undecided treated
# as the infeasible side, so the interval stays a subset of the true one).
.bisect_boundary <- function(bnd, i, j, x_infeas, x_feas, verify,
                             tol = 1e-6, max_it = 40L) {
  for (it in seq_len(max_it)) {
    if (abs(x_feas - x_infeas) < tol) break
    m <- (x_infeas + x_feas) / 2
    if (identical(.feasible_status_at(bnd, i, j, m, verify), "feasible")) {
      x_feas <- m
    } else {
      x_infeas <- m
    }
  }
  x_feas
}

# Signed distance from box [box_lo, box_hi] to interval [lo, hi] (0 if overlap;
# positive if the interval lies above the box, negative if below).
.signed_gap <- function(box_lo, box_hi, lo, hi) {
  if (lo > box_hi) return(lo - box_hi)
  if (hi < box_lo) return(hi - box_lo)
  0
}

# Compute per-cell intervals for every off-diagonal. Returns a list of per-cell
# records; sole-culprit cells have empty = FALSE and an interval + required_edit.
.cell_intervals <- function(R, delta, verify = TRUE) {
  p <- nrow(R)
  bnd <- .box_bounds(R, delta)
  cells <- which(upper.tri(R), arr.ind = TRUE)
  out <- vector("list", nrow(cells))
  for (r in seq_len(nrow(cells))) {
    i <- cells[r, 1]; j <- cells[r, 2]
    fb <- .cell_free_box(bnd, i, j)
    fo <- .box_feasible(fb$lo, fb$hi, fb$off, verify = verify)

    if (identical(fo$status, "infeasible")) {
      out[[r]] <- list(i = i, j = j, empty = TRUE, sole = FALSE,
                       certified = fo$certified)
      next
    }
    if (identical(fo$status, "undecided")) {
      out[[r]] <- list(i = i, j = j, empty = NA, sole = NA, undecided = TRUE)
      next
    }

    # Feasible: the interval is nonempty. Locate its endpoints.
    x_feas <- fo$X[i, j]
    lo_ij <- if (identical(.feasible_status_at(bnd, i, j, -1, verify), "feasible")) {
      -1
    } else {
      .bisect_boundary(bnd, i, j, x_infeas = -1, x_feas = x_feas, verify = verify)
    }
    hi_ij <- if (identical(.feasible_status_at(bnd, i, j, 1, verify), "feasible")) {
      1
    } else {
      .bisect_boundary(bnd, i, j, x_infeas = 1, x_feas = x_feas, verify = verify)
    }

    box_lo <- R[i, j] - delta
    box_hi <- R[i, j] + delta
    edit <- .signed_gap(box_lo, box_hi, lo_ij, hi_ij)
    implied <- min(max(R[i, j], lo_ij), hi_ij)
    out[[r]] <- list(i = i, j = j, empty = FALSE, sole = TRUE,
                     lo = lo_ij, hi = hi_ij, required_edit = edit,
                     implied = implied, reported = R[i, j],
                     witness_matrix = fo$X, certified = fo$certified)
  }
  out
}
